# construct-compiler

Genetic construct design compiler for polycistronic expression vectors. YAML spec in, annotated DNA + assembly plan + cost estimate out.

## Project structure

```
src/construct_compiler/
  core/           # IR: types.py (PortKind, Constraint), parts.py (12 part types), graph.py (ConstructGraph)
  frontend/       # parser.py — YAML → IR
  passes/         # 4-pass pipeline: part_resolution → reverse_translation → constraint_resolution → assembly_planning
  validation/     # construct_checks.py (4 validators), harness.py (eval engine), variants.py (design space)
  backends/       # genbank.py — GenBank export
  vendors/        # twist.py, idt.py — synthesis vendor APIs
  data/           # parts_db.py — 23 vectors, promoters, RBS, tags, codons, overhangs
  server.py       # FastAPI web server
  __main__.py     # CLI (compile, validate, check, parts)
tests/            # pytest suite (39 tests)
examples/         # YAML spec examples
```

## Key commands

```bash
# Compile a spec to GenBank + assembly plan
construct-compiler compile examples/his_tev_mbp_egfp.yaml -o output/

# Run validity checks (reading frame, start codons, translation fidelity, internal stops)
construct-compiler check examples/his_tev_mbp_egfp.yaml

# JSON output for machine consumption
construct-compiler check spec.yaml --json

# Run tests
pytest tests/ -v
```

## Validation harness

The harness is the primary tool for automated design evaluation. Use it whenever you create or modify a construct spec.

### After editing any spec or pipeline code, always run:

```bash
construct-compiler check examples/his_tev_mbp_egfp.yaml
```

Or from Python:

```python
from construct_compiler.validation import evaluate_spec
result = evaluate_spec("path/to/spec.yaml")
assert result.passed, result.summary()
```

### To evaluate design variants:

```python
from construct_compiler.validation import evaluate_batch
from construct_compiler.validation.variants import DesignAxis, vary_spec_dicts

axes = [
    DesignAxis("expression", "cassette.1.cistron.expression", ["high", "medium", "low"]),
    DesignAxis("spacer", "cassette.2.spacer", [20, 30, 50, 100]),
]
specs = vary_spec_dicts("examples/his_tev_mbp_egfp.yaml", axes)
results = evaluate_batch(specs, skip_constraints=True)
best = results[0]  # sorted by score descending
```

### There is also a standalone runner script:

```bash
python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml
python scripts/design_evaluate.py examples/his_tev_mbp_egfp.yaml --axes expression=high,medium,low spacer=20,30,50
```

### DesignAxis field paths

Field paths point into the YAML `construct` block using dot-separated keys and integer list indices:

- `cassette.0` — first cassette element (usually the promoter)
- `cassette.1.cistron.expression` — first cistron's expression level
- `cassette.2.spacer` — spacer length (if spacer is at index 2)
- `cassette.3.cistron.expression` — second cistron's expression level
- `cassette.4.terminator` — terminator choice

## What the validators check

1. **Reading frame continuity** — every coding part's DNA length is divisible by 3, no frame drift across fusion chains
2. **Start codon placement** — first coding element in each cistron starts with ATG
3. **Translation fidelity** — translating the DNA back gives the expected protein sequence
4. **Internal stop codons** — no premature stops in coding regions or at part junctions

## Design rules to follow

- Every cistron needs an RBS before the first coding element
- The first coding element after the RBS gets the start codon (ATG) — this is handled automatically by the compiler
- Fusion chains (tag → cleavage → solubility_tag → gene) are in-frame with no stop codons between them
- The last CDS in a cistron should have `has_stop: true` (set automatically by the parser for the last element in a chain)
- Spacers between cistrons should be at least 20 bp
- Terminators go at the end of the cassette (after all cistrons)

## Gene resolution data gap

**Current state:** `passes/part_resolution.py::_resolve_cds` can resolve CDS parts from:
1. A protein sequence already in the spec
2. `source: fpbase` — fetches by FPbase slug (e.g. `mEGFP`, `mCherry`) via REST API, with a 6-protein built-in fallback
3. `source: uniprot` — fetches by accession (e.g. `P42212`) via UniProt REST API, with a 1-protein built-in fallback

**The gap:** LLM-generated specs use gene names without accessions (`lacI`, `tetR`, `araC`, `luxR`, `luxI`, `dCas9`, `rtTA3`, `cI`, `araC`) and common aliases (`GFP`, `EGFP`, `mVenus`, `YFP`). These all silently fail with "CDS could not be resolved." In the 100-spec panel test this accounts for ~93% of failures.

**Planned approach:**

1. **Expand built-in FP sequences** in `_BUILTIN_FP_SEQUENCES` — add `mVenus`, `mTurquoise2`, `mTagBFP2`, `mCerulean3`, `tdTomato`, `mCitrine`, `iRFP713`, and common aliases (`GFP`/`EGFP` → mEGFP). Sources: FPbase slugs.

2. **Add built-in regulatory/enzyme sequences** in a new `_BUILTIN_COMMON_SEQUENCES` dict. Priority targets (E. coli expression biology):
   - Repressors/activators: `lacI` (P03023), `tetR` (P0ACT4), `araC` (P0A9E0), `cI` (P03034), `luxR` (P12746), `luxI` (P12747), `rtTA3` (synthetic), `KRAB` (human KRAB domain)
   - Enzymes: `dCas9` (synthetic S. pyogenes), common metabolic enzymes from the eval corpus

3. **UniProt gene-name search fallback** — when `source: uniprot` but the ID looks like a gene name (not an accession), call `https://rest.uniprot.org/uniprotkb/search?query=gene:{name}+AND+organism_id:83333&format=fasta` and take the first reviewed hit. Add a small in-memory cache to avoid duplicate requests.

4. **Case-insensitive lookup** — normalize gene names before dict lookup (e.g. `TetR` → `tetR`, `GFP` → `EGFP`).

**Files to change:**
- `src/construct_compiler/passes/part_resolution.py` — main implementation
- `tests/test_part_resolution.py` — add tests for each new resolution path
- `scripts/omega_panel_test.py` — re-run with `--n 100` to measure improvement in pass rate

**Success metric:** panel pass rate ≥50% without requiring any spec changes.

## Dependencies

Python 3.10+. Key libraries: biopython, dnachisel, pyyaml, click. Install with `pip install -e ".[dev]"`.
