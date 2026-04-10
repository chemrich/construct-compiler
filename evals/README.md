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

## Running Evaluations

To run the baseline evaluation again:

```bash
uv run python run_eval.py --corpus prompt_corpus_baseline_1000.yaml --model gemini-3-flash-preview --concurrency 20 --rpm 300 --llm-judge
```
