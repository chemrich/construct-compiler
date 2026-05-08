# Twist Eval Pipeline

End-to-end harness for evaluating spec corpora through Twist's codon-optimization
and pricing APIs. Three scripts in `scripts/`:

| Script | Role |
|---|---|
| `strip_yaml_fences.py`     | One-shot corpus cleanup — strips markdown code fences and other LLM artifacts from generated specs so they parse. |
| `twist_codon_eval.py`      | Per-spec sidecar — compiles each spec, swaps Twist-optimized codons in for DNA Chisel, runs the validators, writes JSONL. |
| `twist_pricing_annotate.py`| Annotates a JSONL output of the codon eval with heuristic + live Twist quote pricing. |

All three are rate-limit-aware: a `RateLimitedTwist` subclass holds a global lock + min-interval gate on every `TwistVendor._request`, so chained async-job polling stays under the configured rate.

## Prereqs

Set in `~/.zshrc` or `.env`:

```bash
TWIST_JWT_TOKEN=...
TWIST_END_USER_TOKEN=...
TWIST_USER_EMAIL=you@example.com
TWIST_USER_PHONE=5551234567        # only for twist_pricing_annotate.py (Twist
                                   # quotes require a phone field; reads from env
                                   # so it doesn't need to be in code)
```

`.env` is gitignored — keep PII (phone, address IDs) there or pass via flags.

For pricing, the user's Twist account also needs:
- a verified Shipping address (`pick_shipping_address` finds the first VERIFIED/APPROVED one; pass `--allow-pending-address` to also accept `PENDING_REVIEW`)
- a phone (from profile, env var, or `--phone` flag)

## 1. Corpus cleanup — `strip_yaml_fences.py`

LLM-generated specs in `evals/generated_specs/` were captured with surrounding markdown fences and stray prose, breaking `yaml.safe_load`. The script:

- keeps only the content between the first ` ```yaml ` opener and the next ` ``` ` closer (drops anything before/after — e.g. `**Note:** ...` paragraphs)
- when a file contains multiple `---`-separated documents, falls back to `safe_load_all` and keeps the first non-null doc (the harness expects one construct per file)
- verifies each cleaned file parses before writing
- is idempotent (clean files are untouched) and supports `--dry-run`

```bash
uv run python scripts/strip_yaml_fences.py --dry-run    # preview
uv run python scripts/strip_yaml_fences.py              # apply
```

Result on the 6745-spec corpus:

| Bucket | Count |
|---|---|
| Already clean | 1,289 |
| Fenced & fixed | 5,381 |
| Multi-doc collapsed to first doc | 47 |
| Still failing parse | 28 |

The remaining 28 are LLM refusals or clarifying-question responses with no YAML content — not auto-fixable.

## 2. Codon-opt evaluation — `twist_codon_eval.py`

For each spec:

1. parse → resolve_parts → reverse_translate (graph with local-codon DNA in every coding part)
2. for each coding part with a known protein, calls `TwistVendor.optimize_codons()` (cached by `(protein, organism)`) and replaces the part's sequence; preserves leading ATG and trailing stop codons that revtrans had added
3. runs `run_all_checks` on the assembled graph
4. appends a JSONL record (resumable: re-running skips specs already in the output)

```bash
# 100 random specs, seed=0, default 5s rate limit
uv run python scripts/twist_codon_eval.py

# 400 specs at 2 req/sec
uv run python scripts/twist_codon_eval.py -n 400 --seed 1 --interval 0.5 \
    --output evals/results/twist_codon_eval_v3.jsonl

# Specific specs (e.g. for a sanity test)
uv run python scripts/twist_codon_eval.py --specs evals/generated_specs/basic_gfp.yaml \
    --output /tmp/sanity.jsonl
```

Each JSONL record:

```json
{
  "spec_path": "...", "spec_name": "...", "construct_name": "...",
  "passed": true, "compile_error": null,
  "twist_calls": 9, "proteins_optimized": 2, "proteins_cached_hits": 0,
  "errors": [], "warnings": [], "error_count": 0, "warning_count": 0,
  "insert_length_bp": 960, "cistron_count": 1, "elapsed_s": 41.77
}
```

### Findings (500 specs across two batches)

| | seed=0 (n=100, 5s gate) | seed=1 (n=400, 0.5s gate) |
|---|---|---|
| Compiled (post fence-fix) | 100/100 | 398/400 |
| Validation pass rate among compiled | 100% | 100% |
| Twist HTTP calls | 229 | 493 |
| Unique proteins optimized | 42 | 85 |
| Cache hits | 131 | 635 |
| Cache hit ratio | 76% | 88% |
| Twist API errors | 0 | 0 |
| Validator errors after Twist swap | 0 | 0 |
| Wall time | 7.4 min | 17.6 min |

The cache pays off heavily — common parts (His-tag, MBP, eGFP, etc.) appear across many specs.

## 3. Pricing annotation — `twist_pricing_annotate.py`

Reads a codon-eval JSONL and writes a new file with two new fields per record: a heuristic price ($/bp tiered by Twist's complexity rating) and a real Twist quote price. The pipeline is two phases:

1. **Per-record screening.** For each compilable record:
   - recompile (parse → resolve → revtrans), assemble the full insert (concat of every non-Backbone part, matching the GenBank backend)
   - if the insert is in Twist's NON_CLONED_GENE window (300–5000 bp), submit it via `_create_construct` and poll `_bulk_retrieve_construct` for scoring → captures `twist_difficulty` + `twist_score`
   - heuristic price = `length × $0.07 × difficulty_multiplier` (STANDARD 1.0, MODERATE 1.2, DIFFICULT 1.4, COMPLEX 1.6, VERY_COMPLEX 1.85)

2. **Bundled quote.** Files **one** `create_quote` covering every BUILDABLE construct (saves account clutter and gets bulk pricing). Polls until `status_info.status == "SUCCESS"`. Twist requires `order_settings`; the script defaults to `[{"name": "Delivery Format", "product_code": "SER_PKG_TUBE"}]`.

3. **Tier extraction.** Twist prices by length tier, not per-line-item. The script parses tiers out of `quote.quote_lines` (each tier has a setup line and a sibling base-pairs line) and computes per-construct prices as `tier_setup + (per_bp × length)`.

```bash
# Full run
uv run python scripts/twist_pricing_annotate.py \
    --input  evals/results/twist_codon_eval_v3.jsonl \
    --output evals/results/twist_codon_eval_v3_priced.jsonl \
    --interval 0.5

# Heuristic only — no API calls, no address required
uv run python scripts/twist_pricing_annotate.py \
    --input  in.jsonl --output out.jsonl --dry-run-heuristic-only

# Skip Phase 1 — reuse twist_construct_id values already in the input file
uv run python scripts/twist_pricing_annotate.py \
    --input priced.jsonl --output priced.jsonl --quote-only

# Skip create_quote — re-parse an existing quote (e.g., after fixing the parser)
uv run python scripts/twist_pricing_annotate.py \
    --input priced.jsonl --output priced.jsonl \
    --quote-only --existing-quote-id 74d638d1-d445-4c0c-80b5-3a18be302cad
```

### Twist's pricing tiers (extracted from a live quote, May 2026)

| Length window | Setup | Per bp |
|---|---|---|
| 300 – 500 bp   | $35.00 | (flat) |
| 501 – 1,800 bp | $0     | $0.07 |
| 1,801 – 3,200 bp | $0   | $0.08 |
| 3,201 – 5,000 bp | $0   | $0.09 |

Inserts outside 300–5000 bp can't go through `NON_CLONED_GENE` and are skipped (heuristic-only). Constructs that Twist scores `UNBUILDABLE` are likewise skipped.

### Findings (the v3 batch, 400 specs)

| Bucket | n |
|---|---|
| Priced (Twist + heuristic side-by-side) | 108 |
| UNBUILDABLE per Twist screening | 107 |
| Length out of NON_CLONED_GENE window | 176 |
| Compile/recompile failures | 9 |

For the 108 priced specs:

- Heuristic total: **$9,046.32** (mean $83.76)
- Twist total:    **$8,006.57** (mean $74.13)
- **Twist / Heuristic ratio: 0.885x** — heuristic over-estimates by ~13%
- STANDARD difficulty: avg $63.68 (n=88); COMPLEX: avg $120.15 (n=20)

The heuristic over-estimates because it applies a flat 1.0× multiplier to STANDARD constructs even though Twist's smaller-tier per-bp rate ($0.07) is at the heuristic baseline — but most of our constructs land in the 501–1800 bp tier where there's no setup fee, so the heuristic includes phantom complexity for constructs Twist treats as cheap.

## Rate-limit guidance

The shared `RateLimitedTwist` gate counts every HTTP call (POSTs and the polling GETs that come with each async job). Numbers from the live runs:

- 1 req / 5 sec (`--interval 5.0`): safe default, never trips IP/API rate limits even on long batches
- 1 req / 1 sec: comfortable for batches up to a few hundred constructs
- 1 req / 0.5 sec: clean against 400-spec batches; tested production, zero API errors

If you scale past a few hundred constructs in a single run, consider 1 req / 0.5 sec with regular checkpoints (the codon-eval script writes JSONL incrementally and is resumable; the pricing script writes once at the end so use `--quote-only` / `--existing-quote-id` for retries).

## File layout reference

```
scripts/
  strip_yaml_fences.py         # corpus cleanup
  twist_codon_eval.py          # per-spec codon-opt + harness validation
  twist_pricing_annotate.py    # heuristic + Twist quote pricing
evals/
  generated_specs/             # 6745 LLM-generated specs (cleaned)
  results/                     # JSONL outputs (gitignored)
docs/
  twist_integration.md         # vendor / TAPI reference
  twist_eval_pipeline.md       # this file
```
