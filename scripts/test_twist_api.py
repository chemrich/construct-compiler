#!/usr/bin/env python3
"""
Quick standalone test for Twist TAPI integration.

Confirms credentials are visible, then exercises the lookup and async
job endpoints against the live API. Writes results to a JSON file.

Note: requires TWIST_USER_EMAIL in addition to the two token env vars.

Usage:
    uv run python scripts/test_twist_api.py
    uv run python scripts/test_twist_api.py --output /tmp/twist_results.json
"""

import json
import os
import sys
import time
import argparse
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from dotenv import load_dotenv
load_dotenv(Path(__file__).resolve().parent.parent / ".env")

from construct_compiler.vendors.twist import TwistVendor


def run_tests() -> dict:
    results = {"timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"), "tests": {}}

    user_email = os.environ.get("TWIST_USER_EMAIL")
    if not user_email:
        raise RuntimeError("TWIST_USER_EMAIL environment variable is required")
    vendor = TwistVendor(user_email=user_email)

    results["credentials"] = {
        "jwt_token_chars": len(vendor.jwt_token),
        "end_user_token_chars": len(vendor.end_user_token),
        "user_email": vendor.user_email,
        "authenticated": vendor.authenticated,
    }

    if not vendor.authenticated:
        results["summary"] = "FAILED — credentials not in env. Need TWIST_JWT_TOKEN, TWIST_END_USER_TOKEN, TWIST_USER_EMAIL."
        return results

    # --- Test 0: Auth Probe (GET /users/{email}/) ---
    import requests
    test = {"name": f"Auth Probe (GET /users/{vendor.user_email}/)", "passed": False}
    try:
        t0 = time.time()
        resp = requests.get(
            vendor._user_url(),  # bare /users/{email}/
            headers=vendor._get_headers(),
            timeout=30,
        )
        test["elapsed_s"] = round(time.time() - t0, 3)
        test["status_code"] = resp.status_code
        if resp.status_code == 200:
            data = resp.json()
            test["response_keys"] = sorted(data.keys()) if isinstance(data, dict) else type(data).__name__
            test["passed"] = True
        else:
            test["body_preview"] = resp.text[:300]
    except Exception as e:
        test["error"] = str(e)
    results["tests"]["auth_probe"] = test

    if not test["passed"]:
        results["summary"] = "FAILED at auth probe — fix this before testing other endpoints."
        return results

    # --- Test 1: Codon Optimization Choices (lookup) ---
    test = {"name": "Codon Optimization Choices", "passed": False}
    try:
        t0 = time.time()
        choices = vendor.get_codon_optimization_choices()
        test["elapsed_s"] = round(time.time() - t0, 3)
        test["organism_count"] = len(choices.get("organism", []))
        test["igg_organism_count"] = len(choices.get("igg_organisms", []))
        test["avoid_introducing_count"] = len(choices.get("avoid_introducing", {}))
        test["sample_organisms"] = list(choices.get("organism", []))[:5]
        test["passed"] = test["organism_count"] > 0
    except Exception as e:
        test["error"] = str(e)
    results["tests"]["codon_opt_choices"] = test

    # --- Test 2: List Vectors ---
    test = {"name": "List Vectors", "passed": False}
    try:
        t0 = time.time()
        vectors = vendor.list_vectors()
        test["elapsed_s"] = round(time.time() - t0, 3)
        test["vector_count"] = len(vectors)
        test["sample_names"] = [v.get("name") for v in vectors[:5]]
        test["passed"] = True  # zero vectors is valid (account may have none onboarded)
    except Exception as e:
        test["error"] = str(e)
    results["tests"]["list_vectors"] = test

    # --- Test 3: Sequence Screening (POST /constructs/ + bulk-retrieve) ---
    # ~350 bp partial GFP CDS
    screen_seq = (
        "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGG"
        "ACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACC"
        "TACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCC"
        "CACCCTCGTGACCACCTTCGGCTACGGCCTGATGTGCTTCGCCCGCTACCCCGACCACA"
        "TGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACC"
        "ATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGA"
    )
    test = {"name": "Sequence Screening", "passed": False, "sequence_length": len(screen_seq)}
    try:
        t0 = time.time()
        r = vendor.screen(screen_seq, "GENE")
        test["elapsed_s"] = round(time.time() - t0, 3)
        test["feasible"] = r.feasible
        test["complexity_score"] = r.complexity_score
        test["turnaround_days"] = list(r.turnaround_days)
        test["warning_count"] = len(r.warnings)
        test["error_count"] = len(r.errors)
        test["warnings_preview"] = r.warnings[:3]
        test["errors_preview"] = r.errors[:3]
        # Mock fallback returns price = length * 0.07; real path returns 0.0.
        test["likely_mock"] = abs(r.estimated_price - len(screen_seq) * 0.07) < 1e-6
        test["passed"] = not test["likely_mock"]
    except Exception as e:
        test["error"] = str(e)
    results["tests"]["screening"] = test

    # --- Test 4: Codon Optimization (reverse-translate + codon-opt chain) ---
    test_protein = "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFS"
    test = {"name": "Codon Optimization (reverse-translate + optimize)", "passed": False,
            "protein_length": len(test_protein)}
    try:
        t0 = time.time()
        # Use a Twist-supported organism. "Escherichia coli" may not be in their
        # enum; "Homo sapiens" is shown explicitly in the docs.
        r = vendor.optimize_codons(test_protein, organism="Homo sapiens")
        test["elapsed_s"] = round(time.time() - t0, 3)
        test["optimized_length"] = len(r.optimized_sequence)
        test["optimized_preview"] = (r.optimized_sequence[:60] + "...") if r.optimized_sequence else "(empty)"
        test["gc_content"] = r.gc_content
        test["notes"] = r.notes
        test["passed"] = bool(r.optimized_sequence)
    except Exception as e:
        test["error"] = str(e)
    results["tests"]["codon_opt"] = test

    passed = sum(1 for t in results["tests"].values() if t.get("passed"))
    total = len(results["tests"])
    results["summary"] = f"{passed}/{total} tests passed"
    return results


def main():
    parser = argparse.ArgumentParser(description="Test Twist TAPI")
    parser.add_argument("--output", "-o", default="/tmp/twist_test_results.json")
    args = parser.parse_args()

    print("Running Twist API tests...")
    results = run_tests()

    with open(args.output, "w") as f:
        json.dump(results, f, indent=2)

    creds = results["credentials"]
    print(f"\n{'='*60}")
    print(f"Credentials: authenticated={creds['authenticated']} "
          f"(jwt={creds['jwt_token_chars']}c, eut={creds['end_user_token_chars']}c, "
          f"email={creds['user_email']})")
    print(f"Results: {results['summary']}")
    print(f"{'='*60}")
    for _, test in results.get("tests", {}).items():
        status = "PASS" if test.get("passed") else "FAIL"
        elapsed = f" ({test['elapsed_s']}s)" if "elapsed_s" in test else ""
        mock_flag = " [MOCK FALLBACK]" if test.get("likely_mock") else ""
        print(f"  [{status}] {test['name']}{elapsed}{mock_flag}")
        if not test.get("passed") and "error" in test:
            print(f"         {test['error'][:160]}")
    print(f"\nFull results: {args.output}")


if __name__ == "__main__":
    main()
