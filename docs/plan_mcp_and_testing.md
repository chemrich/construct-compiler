# Plan: MCP Server + Automated Testing Harness

**Date:** 2026-03-24
**Status:** Draft

---

## Vision

Enable a conversational design workflow where a user describes a construct in natural language and Claude iteratively builds, validates, and refines a YAML spec — using construct_compiler as an MCP tool server. Separately, establish CI-level automated testing so that pipeline changes are gated by the validation harness.

These two goals reinforce each other: the MCP tools call the same harness that CI uses, so an agent-driven design session and a CI run exercise the same validation logic.

---

## Part 1: MCP Server (`construct-compiler-mcp`)

### What is an MCP server?

The Model Context Protocol (MCP) is a JSON-RPC-over-stdio protocol that lets AI assistants (like Claude Code) discover and call tools exposed by external programs. The server runs as a subprocess, reads JSON-RPC from stdin, writes responses to stdout. No HTTP, no ports — just piped I/O.

For construct_compiler this means: Claude Code spawns `construct-compiler-mcp`, discovers three tools (compile, check, evaluate_variants), and can call them mid-conversation.

### Tools to expose

#### 1. `compile_spec`

Compiles a YAML spec dict through the full pipeline. Returns the construct summary, assembly strategies with costs, any errors/warnings, and optionally the GenBank output.

**Input schema:**
```json
{
  "spec": { "construct": { ... } },
  "skip_constraints": false,
  "include_genbank": false
}
```

**Output:** construct name, host, insert length, parts list, assembly strategies (cost breakdown, days, fragments), errors, warnings, and optionally a GenBank string.

**Why this tool:** This is the core compilation step. The agent drafts a YAML spec from the user's natural language description and compiles it to see if it works and what it costs.

#### 2. `check_spec`

Runs the four validity checks (reading frames, start codons, translation fidelity, internal stops) and returns a structured result with pass/fail, score, and individual check details.

**Input schema:**
```json
{
  "spec": { "construct": { ... } },
  "skip_constraints": false,
  "validate_intermediate": false
}
```

**Output:** passed (bool), score (0.0–1.0), error_count, warning_count, individual check results with severity and messages, compile_time_s.

**Why this tool:** After compiling, the agent calls check to verify biological correctness. If checks fail, the agent reads the error messages and modifies the spec accordingly. This is the core of the turn-based iteration loop.

#### 3. `evaluate_variants`

Generates a combinatorial design space from a base spec and a set of axes, evaluates all variants, and returns results ranked by score.

**Input schema:**
```json
{
  "spec": { "construct": { ... } },
  "axes": [
    { "name": "expression", "field_path": "cassette.1.cistron.expression", "values": ["high", "medium", "low"] },
    { "name": "spacer", "field_path": "cassette.2.spacer", "values": [20, 30, 50] }
  ],
  "skip_constraints": true
}
```

**Output:** list of `{ variant_name, parameters, passed, score, error_count, warning_count }` sorted by score descending, plus a batch summary.

**Why this tool:** Once a base design passes checks, the agent can explore the parameter space to find optimal configurations. This is the "design exploration" phase of the conversation.

### Implementation plan

#### File: `src/construct_compiler/mcp_server.py`

The MCP server is a single Python module (~200–300 lines). It uses the `mcp` Python SDK (the official reference implementation from Anthropic) which handles the JSON-RPC/stdio transport. Our code just registers tool handlers.

```
pip install mcp
```

**Structure:**

```python
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

server = Server("construct-compiler")

@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(name="compile_spec", description="...", inputSchema={...}),
        Tool(name="check_spec",   description="...", inputSchema={...}),
        Tool(name="evaluate_variants", description="...", inputSchema={...}),
    ]

@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    if name == "compile_spec":
        return _handle_compile(arguments)
    elif name == "check_spec":
        return _handle_check(arguments)
    elif name == "evaluate_variants":
        return _handle_variants(arguments)

async def main():
    async with stdio_server() as (read, write):
        await server.run(read, write, server.create_initialization_options())
```

Each handler is a thin wrapper around the existing Python API:

- `_handle_compile` → calls `parse_spec()` + `run_pipeline()` + `plan_assembly()` (same logic as `server.py`'s `/api/compile`)
- `_handle_check` → calls `evaluate_spec()` from `validation.harness`
- `_handle_variants` → calls `vary_spec()` + `evaluate_batch()`

**Entry point** (in pyproject.toml):
```toml
[project.scripts]
construct-compiler-mcp = "construct_compiler.mcp_server:main"
```

#### Claude Code configuration

Users add this to their `.claude/settings.json` or project-level `.mcp.json`:

```json
{
  "mcpServers": {
    "construct-compiler": {
      "command": "construct-compiler-mcp",
      "args": []
    }
  }
}
```

Or if running from the repo checkout:

```json
{
  "mcpServers": {
    "construct-compiler": {
      "command": "python",
      "args": ["-m", "construct_compiler.mcp_server"]
    }
  }
}
```

### The conversational workflow

Here's what a typical session looks like once the MCP server is wired up:

```
User:  "I need an E. coli expression vector for producing a GFP-tagged
        version of human lysozyme, with a His tag for purification and
        a TEV cleavage site. Use Golden Gate assembly."

Claude: [drafts a YAML spec based on the description]
        [calls compile_spec tool]
        → "Compiled successfully. 2 assembly strategies found.
           Cheapest: 2-fragment Golden Gate, $127, 14 days."

Claude: [calls check_spec tool]
        → "PASSED (score 1.0). All 4 checks clean."

Claude: "Your construct compiled and passed all checks. Here's the
         spec I generated: [shows YAML]. The cheapest assembly
         strategy is a 2-fragment Golden Gate at $127. Want me to
         explore different expression levels or spacer lengths?"

User:  "Yes, try high/medium/low expression and spacers 20–100."

Claude: [calls evaluate_variants tool with 3×4 = 12 variants]
        → "All 12 passed. Best: high expression + 30bp spacer (score 1.0).
           Cheapest: low expression + 20bp spacer ($89)."
```

### Implementation phases

**Phase 1 — Minimal viable server (1–2 days)**
- Implement `mcp_server.py` with `compile_spec` and `check_spec`
- Add `mcp` to optional dependencies: `[project.optional-dependencies] mcp = ["mcp>=1.0"]`
- Add entry point to pyproject.toml
- Write 3–5 integration tests that spawn the server and exercise the tools
- Test with Claude Code manually

**Phase 2 — Variant exploration (1 day)**
- Add `evaluate_variants` tool
- Add input validation (cap max variants at ~100 to avoid long waits)
- Test variant workflow end-to-end

**Phase 3 — Polish (1 day)**
- Tool descriptions tuned for Claude (clear parameter docs, examples in descriptions)
- Error messages designed for agent consumption (actionable, not stack traces)
- Add a `resources` endpoint exposing the example YAML specs as MCP resources (so Claude can read them as templates)
- Add a `prompts` endpoint with a "design a construct" prompt template

---

## Part 2: Automated Testing with the Harness

### Goal

Two levels of automated testing:

1. **CI regression suite** — pytest tests that compile every example spec, run all checks, and fail the build if anything regresses.
2. **Agent-driven design loop** — Claude uses the MCP tools to iterate on designs, with the harness as the source of truth for correctness.

The CI suite ensures pipeline changes don't break existing specs. The agent loop enables design exploration that feeds back into the example/test corpus.

### CI regression suite

#### File: `tests/test_harness_regression.py`

```python
"""
Regression tests: every example spec must compile and pass all checks.

These tests ensure that pipeline changes don't silently break existing
construct designs. Run as part of the normal pytest suite.
"""
import pytest
from pathlib import Path
from construct_compiler.validation.harness import evaluate_spec

EXAMPLES_DIR = Path(__file__).parent.parent / "examples"
EXAMPLE_SPECS = sorted(EXAMPLES_DIR.glob("*.yaml"))

@pytest.fixture(params=EXAMPLE_SPECS, ids=lambda p: p.stem)
def example_spec(request):
    return request.param

def test_example_compiles_and_passes(example_spec):
    """Every example spec must compile cleanly and pass all 4 checks."""
    result = evaluate_spec(example_spec, skip_constraints=True)
    assert result.passed, (
        f"{example_spec.name} failed with score {result.score}:\n"
        f"{result.summary()}"
    )
    assert result.error_count == 0

def test_example_with_constraints(example_spec):
    """Full pipeline including codon optimization — may be slow."""
    result = evaluate_spec(example_spec, skip_constraints=False)
    assert result.passed, result.summary()
```

#### File: `tests/test_harness_variants.py`

```python
"""
Variant evaluation smoke tests.

Ensures the variant generator + batch evaluator work together,
and that known-good parameter combinations still pass.
"""
from construct_compiler.validation.harness import evaluate_batch
from construct_compiler.validation.variants import DesignAxis, vary_spec_dicts

def test_expression_variants():
    """All expression levels should produce valid constructs."""
    axes = [
        DesignAxis("expression", "cassette.1.cistron.expression",
                   ["high", "medium", "low"]),
    ]
    specs = vary_spec_dicts("examples/his_tev_mbp_egfp.yaml", axes)
    results = evaluate_batch(specs, skip_constraints=True)
    for r in results:
        assert r.passed, r.summary()

def test_spacer_variants():
    """All reasonable spacer lengths should produce valid constructs."""
    axes = [
        DesignAxis("spacer", "cassette.2.spacer", [20, 30, 50, 100]),
    ]
    specs = vary_spec_dicts("examples/his_tev_mbp_egfp.yaml", axes)
    results = evaluate_batch(specs, skip_constraints=True)
    for r in results:
        assert r.passed, r.summary()
```

#### pytest markers for speed control

```ini
# pyproject.toml additions
[tool.pytest.ini_options]
markers = [
    "slow: tests that run codon optimization (deselect with -m 'not slow')",
]
```

The `test_example_with_constraints` tests get the `@pytest.mark.slow` marker since codon optimization takes seconds per spec. CI can run fast tests on every push and slow tests nightly.

### Agent-driven testing loop

This is where the MCP server and the harness connect. The workflow:

1. Claude generates or modifies a YAML spec
2. Claude calls `check_spec` → harness runs 4 validators
3. If checks fail, Claude reads the errors and edits the spec
4. Repeat until all checks pass
5. Optionally, Claude calls `evaluate_variants` to explore the design space
6. Best-scoring variants can be saved as new example specs (expanding the regression corpus)

The key insight: **every design the agent validates in conversation also becomes a potential test case.** We can add a convention where specs that Claude generates and validates get saved to `examples/agent_generated/` and automatically picked up by the regression suite.

### Implementation phases

**Phase 1 — Regression suite (half day)**
- Write `test_harness_regression.py` and `test_harness_variants.py`
- Add pytest markers for slow tests
- Verify all existing example specs pass
- Add to CI (GitHub Actions or similar)

**Phase 2 — Golden spec corpus (half day)**
- Create `examples/agent_generated/` directory
- Add a `save_spec` helper that the MCP server can call to persist validated specs
- Regression suite auto-discovers new specs

**Phase 3 — CI pipeline (half day)**
- GitHub Actions workflow: fast tests on push, full suite (with constraints) nightly
- Coverage reporting on the harness and variant code
- Badge in README

---

## Dependency additions

```toml
[project.optional-dependencies]
mcp = [
    "mcp>=1.0",
]
```

The `mcp` SDK is the only new dependency. Everything else (validation harness, variant engine, compilation pipeline) already exists.

---

## Summary of deliverables

| Deliverable | New files | Effort |
|---|---|---|
| MCP server (compile + check) | `mcp_server.py`, entry point | 1–2 days |
| MCP variant tool | addition to `mcp_server.py` | 1 day |
| MCP resources + prompts | addition to `mcp_server.py` | half day |
| Regression test suite | `tests/test_harness_regression.py` | half day |
| Variant test suite | `tests/test_harness_variants.py` | half day |
| CI workflow | `.github/workflows/test.yml` | half day |
| **Total** | | **~4–5 days** |

---

## Decisions (2026-03-24)

1. **GenBank output: inline.** Max construct size ~7,500 bp (~15–20 KB as GenBank text). Fine for inline return.
2. **Variant cap: 50.** `evaluate_variants` rejects requests that would produce >50 combinations.
3. **Spec persistence: auto-save.** Validated specs from agent sessions auto-save to `examples/agent_generated/`. The regression suite discovers them automatically.
