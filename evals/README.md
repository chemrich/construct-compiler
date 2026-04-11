# Construct Compiler Evaluations

This directory contains the evaluation harness for the `construct-compiler`. It is used to test the full pipeline from natural language prompt to compiled and validated construct.

## Overview

The evaluation harness sends prompts to an LLM, parses the generated YAML specs, compiles them, and checks both harness validity and specific expectation properties (host, cistron count, required parts).

## Files

*   `run_eval.py`: The main execution script for running evaluations.
*   `generate_persona_prompts.py`: Script to generate prompts for specific personas.
*   `thesaurus.yaml`: Normalization mapping for biological terms to reduce false negatives in expectation matching.
*   `spec_generation_prompt.txt`: The system prompt used for the generator model.

## Corpora

*   `prompt_corpus_baseline_1000.yaml`: A combined corpus of 1000 prompts, taking the first 100 prompts from each of the 10 personas. Used for baseline calibration.
*   Persona-specific corpora (e.g., `prompt_corpus_undergrad.yaml`, `prompt_corpus_protein_engineer_bp.yaml`, etc.).

## Baseline Evaluation Results (2026-04-10)

A comprehensive baseline evaluation was run on `prompt_corpus_baseline_1000.yaml` using `gemini-3-flash-preview` as the generator and judge.

### Summary

*   **Total Prompts**: 1000
*   **Specs Generated**: 999 (99.9%)
*   **Harness Passed**: 995 (99.5%)
*   **All Expectations Met**: **68.3%**

### Persona Breakdown

| Persona | Total | Expectations Met | Rate |
| :--- | :--- | :--- | :--- |
| `minimalist` | 100 | 95 | 95.0% |
| `protein_engineer_bp` | 100 | 86 | 86.0% |
| `old_school` | 100 | 80 | 80.0% |
| `industry_vet` | 100 | 76 | 76.0% |
| `undergrad` | 100 | 75 | 75.0% |
| `metabolic_engineer` | 100 | 68 | 68.0% |
| `postdoc` | 100 | 64 | 64.0% |
| `therapeutic` | 100 | 58 | 58.0% |
| `repetition_tester` | 100 | 51 | 51.0% |
| `pure_logic` | 100 | 30 | 30.0% |

### Key Observations

*   **`minimalist`** and **`protein_engineer_bp`** personas performed best.
*   **`pure_logic`** performed worst, indicating difficulty with abstract logical constraints.
*   Common failures included `logic_verification` in `therapeutic` and `host_compatibility` in `protein_engineer_bp` (some of which may be judge strictness issues).

## Haiku Full-Persona Evaluation (2026-04-11)

A cost-optimised run across all 10 persona corpora (500 prompts each, 5000 total) using `claude-haiku-4-5` as the generator and LLM judge.

### Summary

| Metric | Value |
| :--- | :--- |
| Total prompts | 5,000 |
| Harness pass rate | **97.1%** (4,856 / 5,000) |
| All expectations met | **31.8%** (1,590 / 5,000) |
| Individual expectation pass rate | **79.1%** |

### Persona Breakdown

| Persona | Pass% | All-exp-met% |
| :--- | :--- | :--- |
| `metabolic_engineer` | 100.0% | 16.2% |
| `protein_engineer_bp` | 99.4% | 43.3% |
| `old_school` | 99.2% | 43.8% |
| `undergrad` | 98.8% | 26.9% |
| `minimalist` | 98.6% | **71.6%** |
| `industry_vet` | 97.0% | 45.4% |
| `therapeutic` | 96.4% | 29.7% |
| `pure_logic` | 96.2% | 5.6% |
| `repetition_tester` | 93.6% | 25.0% |
| `postdoc` | 92.0% | 18.3% |

### Top Expectation Failures

| Expectation | Fail rate |
| :--- | :--- |
| `has_KozakSequence` | 92.6% |
| `has_IRES` | 85.4% |
| `has_RBS` | 82.7% |
| `logic_verification` | 64.9% |
| `has_PolyA_Signal` | 65.4% |
| `host_compatibility` | 35.0% |

---

## Model Selection — Impact on Results

**The model you choose for generation and evaluation has significant, sometimes serious implications on your results.** This is not just a speed/cost trade-off.

### Generator model

The generator model is responsible for turning natural language prompts into valid YAML specs. Differences between model tiers show up in two ways:

**Structural quality (harness pass rate):** Smaller/cheaper models (e.g. Haiku) are more likely to emit malformed output — most commonly, Markdown formatting syntax (e.g. `**bold text**`) embedded in or after the YAML block, which causes YAML parse failures. They also more frequently output `null` for required fields like `gene.id` or `promoter`, which can crash the compiler pipeline. Frontier models (Sonnet, Gemini Flash) are considerably more reliable here.

**Semantic quality (expectation rate):** This is where the gap is most pronounced. Haiku achieves only **31.8% all-expectations-met** vs. **68.3% for Gemini Flash** on the same baseline corpus. Haiku systematically omits biologically important elements that the prompt implies but does not explicitly name — RBS for prokaryotic constructs, Kozak sequences for mammalian constructs, IRES for bicistronic designs, correct host-terminator pairing, and so on. It also cannot reliably implement abstract logic gate circuits. Expect Haiku results to understate the system's true capability.

**Cost guidance:** Use Haiku for rapid iteration on corpus structure, prompt wording, and pipeline mechanics. Switch to Sonnet (or equivalent frontier model) for any results you intend to treat as ground truth or share externally.

### Judge model

The LLM judge evaluates whether a generated spec meets the semantic expectations defined in the corpus. Judge model choice affects:

- **False negatives:** A weaker judge may miss a valid design that satisfies an expectation through equivalent but non-obvious means (e.g. a different RBS name that it doesn't recognise).
- **False positives:** A weaker judge may pass a spec that doesn't actually meet a biological requirement, inflating scores.
- **`logic_verification` reliability:** Abstract circuit expectations (toggle switches, repressilators, Boolean gates) require strong reasoning. Haiku-as-judge gives unreliable verdicts on these; the ~65% fail rate for `logic_verification` should be interpreted cautiously when Haiku is the judge.

**Recommendation:** Always use a frontier model (Sonnet, Gemini Flash, GPT-4o) as the judge when expectation accuracy matters. The `--judge-model` flag controls this independently of `--model`.

---

## Running Evaluations

To run the baseline evaluation again:

```bash
uv run python run_eval.py --corpus prompt_corpus_baseline_1000.yaml --model gemini-3-flash-preview --concurrency 20 --rpm 300 --llm-judge
```
