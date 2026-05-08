#!/usr/bin/env python3
"""
Annotate an eval JSONL with pricing — heuristic + live Twist quote, side by side.

Reads records produced by twist_codon_eval.py, recompiles each passing
spec to obtain the full assembled insert DNA, and writes a new JSONL
with two new fields per record:

  heuristic_price_usd   length(bp) * $0.07 (NON_CLONED_GENE rate),
                        scaled by Twist difficulty when known.
  twist_price_usd       real list_unit_price returned by a single
                        bundled Twist quote covering every spec.
  twist_difficulty      "STANDARD" / "MODERATE" / ... from screening.
  twist_turnaround_days [low, high] from quote.tat.business_days.

Pipeline:
  1. read input JSONL
  2. for each record where compile_error is None:
       - recompile (parse, resolve_parts, reverse_translate)
       - assembled insert = concat of non-Backbone parts in graph order
       - submit construct via _create_construct (NON_CLONED_GENE) and
         poll until scored, capturing difficulty + buildability
       - heuristic price = insert_length * BASE_RATE * difficulty_multiplier
  3. issue a single create_quote with all (BUILDABLE) construct UUIDs as
     containers, poll get_quote until status_info.status == "SUCCESS"
  4. parse line_items by index, attach twist_price_usd to each record
  5. write annotated JSONL

A global rate limiter (default 0.5s) gates every Twist HTTP call.

Requires the user to have a verified shipping address in their Twist
account (TwistVendor.list_addresses() must return at least one entry
with verification_status == "VERIFIED"). The script auto-picks the
first verified shipping address; pass --address-id to override.

Usage:
    uv run python scripts/twist_pricing_annotate.py \
        --input  evals/results/twist_codon_eval_v3.jsonl \
        --output evals/results/twist_codon_eval_v3_priced.jsonl
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

try:
    from dotenv import load_dotenv
    load_dotenv(ROOT / ".env")
except ImportError:
    pass

from Bio.Seq import Seq

from construct_compiler.frontend.parser import parse_spec
from construct_compiler.passes.part_resolution import resolve_parts
from construct_compiler.passes.reverse_translation import reverse_translate
from construct_compiler.core.parts import Backbone
from construct_compiler.vendors.twist import TwistVendor

# Reuse the rate-limited wrapper from the sibling script
sys.path.insert(0, str(ROOT / "scripts"))
from twist_codon_eval import RateLimitedTwist  # type: ignore

logger = logging.getLogger(__name__)

FALLBACK_USER_EMAIL = "REMOVED"
BASE_RATE_USD_PER_BP = 0.07          # NON_CLONED_GENE
DIFFICULTY_MULTIPLIER = {
    "STANDARD":     1.00,
    "MODERATE":     1.20,
    "DIFFICULT":    1.40,
    "COMPLEX":      1.60,
    "VERY_COMPLEX": 1.85,
}
QUOTE_POLL_INTERVAL = 5.0
QUOTE_POLL_TIMEOUT = 1800.0   # 30 min


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def assembled_insert(graph) -> str:
    """Concatenate every non-Backbone part's sequence in graph order
    (matches the genbank backend's insert assembly)."""
    parts: list[str] = []
    for part in graph.parts():
        if isinstance(part, Backbone):
            continue
        if part.sequence is None:
            continue
        parts.append(str(part.sequence))
    return "".join(parts)


def pick_shipping_address(vendor: TwistVendor,
                          allow_pending: bool = False) -> Optional[dict]:
    """Find the first usable shipping address on the user's Twist account.
    By default only VERIFIED/APPROVED addresses count; pass
    ``allow_pending=True`` to also accept PENDING_REVIEW."""
    accepted = {"VERIFIED", "APPROVED"}
    if allow_pending:
        accepted.add("PENDING_REVIEW")
    addrs = vendor.list_addresses()
    candidates = [a for a in addrs
                  if a.get("address_type") == "Shipping"
                  and a.get("verification_status") in accepted]
    if not candidates:
        return None
    default = [a for a in candidates if a.get("is_default")]
    return (default or candidates)[0]


def heuristic_price(length_bp: int, difficulty: str) -> float:
    mult = DIFFICULTY_MULTIPLIER.get(difficulty, 1.0)
    return round(length_bp * BASE_RATE_USD_PER_BP * mult, 2)


# ---------------------------------------------------------------------------
# Phases
# ---------------------------------------------------------------------------

NON_CLONED_GENE_MIN_BP = 300
NON_CLONED_GENE_MAX_BP = 5000
TWIST_NAME_MAX = 32


def submit_constructs(records: list[dict], vendor: Optional[RateLimitedTwist],
                      ) -> list[dict]:
    """For each compilable record, recompile and capture insert length +
    heuristic price. When ``vendor`` is provided, also submit a
    NON_CLONED_GENE construct (when length is within Twist's 300–5000 bp
    window) and poll for scoring to attach difficulty + buildability.
    Returns the list of records eligible for a Twist quote (BUILDABLE)."""
    submitted: list[dict] = []
    for i, rec in enumerate(records, 1):
        if rec.get("compile_error"):
            continue
        spec_path = rec.get("spec_path")
        if not spec_path:
            continue

        try:
            graph = parse_spec(spec_path)
            graph = resolve_parts(graph)
            graph = reverse_translate(graph)
        except Exception as exc:
            rec["pricing_error"] = f"recompile: {exc!r}"
            continue

        seq = assembled_insert(graph)
        if not seq:
            rec["pricing_error"] = "empty assembled insert"
            continue
        rec["insert_length_bp"] = len(seq)
        # Heuristic always available; difficulty defaults to STANDARD when
        # we don't have a Twist screening result.
        rec["heuristic_price_usd"] = heuristic_price(len(seq), "STANDARD")

        if vendor is None:
            print(f"  [{i}/{len(records)}] {Path(spec_path).name} "
                  f"len={len(seq)} heuristic=${rec['heuristic_price_usd']:.2f}",
                  flush=True)
            continue

        if len(seq) < NON_CLONED_GENE_MIN_BP or len(seq) > NON_CLONED_GENE_MAX_BP:
            rec["pricing_note"] = (
                f"length {len(seq)} outside NON_CLONED_GENE window "
                f"({NON_CLONED_GENE_MIN_BP}-{NON_CLONED_GENE_MAX_BP} bp); "
                f"heuristic only"
            )
            print(f"  [{i}/{len(records)}] {Path(spec_path).name} "
                  f"len={len(seq)} (skip-quote) "
                  f"heuristic=${rec['heuristic_price_usd']:.2f}", flush=True)
            continue

        # Submit + poll for scoring (gives difficulty rating).
        base_name = rec.get("construct_name") or Path(spec_path).stem
        try:
            cid = vendor._create_construct(
                seq,
                name=base_name[:TWIST_NAME_MAX],
                construct_type="NON_CLONED_GENE",
                external_id=f"price-{i:04d}",
            )
            scored = vendor._bulk_retrieve_construct(cid)
        except Exception as exc:
            rec["pricing_error"] = f"submit/score: {exc!r}"
            continue

        score_data = scored.get("score_data") or {}
        difficulty = score_data.get("difficulty", "STANDARD")
        score = scored.get("score", "")
        rec["twist_construct_id"] = cid
        rec["twist_difficulty"] = difficulty
        rec["twist_score"] = score
        # Refine heuristic with the real difficulty multiplier.
        rec["heuristic_price_usd"] = heuristic_price(len(seq), difficulty)

        if score == "BUILDABLE":
            submitted.append(rec)
        else:
            # Unbuildable constructs can't be quoted — leave heuristic only.
            rec["pricing_note"] = f"not buildable ({score}); excluded from quote"

        print(f"  [{i}/{len(records)}] {Path(spec_path).name} "
              f"len={len(seq)} difficulty={difficulty} score={score} "
              f"heuristic=${rec['heuristic_price_usd']:.2f}", flush=True)

    return submitted


DEFAULT_ORDER_SETTINGS = [{
    "name": "Delivery Format",
    "product_code": "SER_PKG_TUBE",
}]


def request_quote(vendor: RateLimitedTwist, submitted: list[dict],
                  recipient_address_id: str, first_name: str, last_name: str,
                  phone: str,
                  order_settings: Optional[list[dict]] = None) -> dict:
    """Issue ONE bundled quote covering every BUILDABLE construct."""
    containers = [{
        "constructs": [
            {"id": rec["twist_construct_id"], "index": i + 1}
            for i, rec in enumerate(submitted)
        ]
    }]
    external_id = f"price-batch-{int(time.time())}"
    print(f"\nFiling bundled quote with {len(submitted)} constructs "
          f"(external_id={external_id})...", flush=True)
    quote = vendor.create_quote(
        external_id=external_id,
        containers=containers,
        order_sub_product_type="NON_CLONAL_ADAPTERS_OFF",
        recipient_address_id=recipient_address_id,
        first_name=first_name,
        last_name=last_name,
        phone=phone,
        order_settings=order_settings or DEFAULT_ORDER_SETTINGS,
    )
    quote_id = quote.get("id")
    if not quote_id:
        raise RuntimeError(f"create_quote returned no id: {quote}")

    print(f"Quote {quote_id} filed; polling until SUCCESS "
          f"(timeout {QUOTE_POLL_TIMEOUT/60:.0f} min)...", flush=True)
    deadline = time.time() + QUOTE_POLL_TIMEOUT
    last_status = ""
    while time.time() < deadline:
        q = vendor.get_quote(quote_id) or {}
        status = (q.get("status_info") or {}).get("status", "")
        if status != last_status:
            print(f"  status: {status}", flush=True)
            last_status = status
        if status == "SUCCESS":
            return q
        if status == "ERROR":
            err = (q.get("status_info") or {}).get("error", "<no detail>")
            raise RuntimeError(f"quote {quote_id} ended in ERROR: {err}")
        time.sleep(QUOTE_POLL_INTERVAL)
    raise TimeoutError(f"quote {quote_id} did not reach SUCCESS within "
                       f"{QUOTE_POLL_TIMEOUT:.0f}s")


_SIZE_RANGE_RE = re.compile(
    r"(\d+(?:\.\d+)?)\s*(bp|kb)\s*-\s*(\d+(?:\.\d+)?)\s*(bp|kb)",
    re.IGNORECASE,
)


def _to_bp(n: float, unit: str) -> float:
    return n * 1000 if unit.lower() == "kb" else n


def parse_pricing_tiers(quote: dict) -> list[tuple[float, float, float, float, str]]:
    """Extract Twist gene-fragment pricing tiers from a SUCCESS quote.

    Returns a list of (lo_bp, hi_bp, setup_fee, per_bp_rate, tier_code).
    Twist's quote_lines split each tier into:
      - a 'GEN_FRG_*AF_*' setup line  (flat charge per construct in tier)
      - a 'GEN_FRG_AF_BP_*' line       (per-bp rate for that tier)
    Some tiers (notably 300-500bp) have only the setup line.
    """
    quote_lines = (quote.get("quote") or {}).get("quote_lines") or []
    setups: dict[tuple[float, float], tuple[float, str]] = {}
    perbps: dict[tuple[float, float], float] = {}

    for line in quote_lines:
        code = line.get("product_code") or ""
        desc = line.get("description") or ""
        if not code.startswith("GEN_FRG_"):
            continue
        m = _SIZE_RANGE_RE.search(desc)
        if not m:
            continue
        lo = _to_bp(float(m.group(1)), m.group(2))
        hi = _to_bp(float(m.group(3)), m.group(4))
        rate = float(line.get("list_unit_price") or 0)
        if "Base Pairs" in desc or "AF_BP_" in code:
            perbps[(lo, hi)] = rate
        else:
            setups[(lo, hi)] = (rate, code)

    tiers: list[tuple[float, float, float, float, str]] = []
    for (lo, hi), (setup, code) in setups.items():
        per_bp = perbps.get((lo, hi), 0.0)
        tiers.append((lo, hi, setup, per_bp, code))
    tiers.sort()
    return tiers


def price_from_tiers(length_bp: int,
                     tiers: list[tuple[float, float, float, float, str]]
                     ) -> Optional[tuple[float, str]]:
    for lo, hi, setup, per_bp, code in tiers:
        if lo <= length_bp <= hi:
            return setup + per_bp * length_bp, code
    return None


def attach_quote_prices(submitted: list[dict], quote: dict) -> None:
    """Compute per-construct prices from the quote's tier structure and
    set twist_price_usd / twist_tier_code / twist_turnaround_days on each
    submitted record. Matches by length (since Twist line-items don't
    carry a construct_id back; they carry gene_name + length)."""
    tiers = parse_pricing_tiers(quote)
    if not tiers:
        for rec in submitted:
            rec["pricing_error"] = "no pricing tiers parsed from quote"
        return

    for rec in submitted:
        length = rec.get("insert_length_bp") or 0
        result = price_from_tiers(length, tiers)
        if result is None:
            rec["pricing_error"] = f"no tier matches length {length}"
            continue
        price, code = result
        rec["twist_price_usd"] = round(price, 2)
        rec["twist_tier_code"] = code

    business_days = ((quote.get("tat") or {}).get("business_days"))
    if isinstance(business_days, (int, float)):
        for rec in submitted:
            rec["twist_turnaround_business_days"] = int(business_days)

    quote_total = (quote.get("quote") or {}).get("price")
    quote_subtotal = (quote.get("quote") or {}).get("subtotal")
    print(f"  Twist quote total: ${quote_total}  subtotal: ${quote_subtotal}")
    print(f"  Pricing tiers parsed:")
    for lo, hi, setup, per_bp, code in tiers:
        print(f"    {int(lo):>5}-{int(hi):>5} bp: ${setup:>5.2f} setup + ${per_bp:.2f}/bp  ({code})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True,
                   help="Eval JSONL produced by twist_codon_eval.py")
    p.add_argument("--output", type=Path, required=True,
                   help="Annotated JSONL to write")
    p.add_argument("--interval", type=float, default=0.5,
                   help="Min seconds between Twist HTTP calls")
    p.add_argument("--user-email", default=os.environ.get("TWIST_USER_EMAIL")
                                            or FALLBACK_USER_EMAIL)
    p.add_argument("--sandbox", action="store_true")
    p.add_argument("--address-id", default=None,
                   help="Override the picked shipping address ID")
    p.add_argument("--allow-pending-address", action="store_true",
                   help="Accept a PENDING_REVIEW shipping address "
                        "(Twist may still allow quoting before verification)")
    p.add_argument("--phone", default=os.environ.get("TWIST_USER_PHONE"),
                   help="Phone number for the quote. Defaults to "
                        "TWIST_USER_PHONE env var (read from .env). "
                        "Used when the user's Twist profile has none.")
    p.add_argument("--limit", type=int, default=None,
                   help="Process only the first N records (for testing)")
    p.add_argument("--dry-run-heuristic-only", action="store_true",
                   help="Compute heuristic price only; skip Twist screening + quote")
    p.add_argument("--quote-only", action="store_true",
                   help="Skip Phase 1; reuse twist_construct_id values "
                        "already in --input and just file the bundled quote")
    p.add_argument("--existing-quote-id", default=None,
                   help="Skip create_quote and parse this quote_id instead "
                        "(useful for retrying just the price-extraction step)")
    args = p.parse_args()

    logging.basicConfig(level=logging.WARNING, format="%(message)s")

    if not args.input.exists():
        print(f"input not found: {args.input}", file=sys.stderr)
        return 2

    records = [json.loads(line) for line in args.input.open()]
    if args.limit:
        records = records[:args.limit]
    print(f"Loaded {len(records)} records from {args.input}")

    vendor = RateLimitedTwist(
        user_email=args.user_email,
        sandbox=args.sandbox,
        min_interval=args.interval,
    )
    if not vendor.authenticated:
        print("Twist credentials missing.", file=sys.stderr)
        return 1

    first_name = last_name = phone = recipient_address_id = ""
    if not args.dry_run_heuristic_only:
        user = vendor.get_user()
        first_name = user.get("first_name") or ""
        last_name = user.get("last_name") or ""
        phone = user.get("phone_number") or ""
        if args.address_id:
            recipient_address_id = args.address_id
        else:
            addr = pick_shipping_address(vendor,
                                          allow_pending=args.allow_pending_address)
            if addr is None:
                print("No usable shipping address on the user's Twist "
                      "account. Add one (web UI or create_address); pass "
                      "--allow-pending-address to use a PENDING_REVIEW "
                      "address.", file=sys.stderr)
                return 1
            recipient_address_id = addr.get("id", "")
            print(f"Using shipping address id={recipient_address_id} "
                  f"({addr.get('city')}, {addr.get('state')}) "
                  f"[{addr.get('verification_status')}]")
        if args.phone:
            phone = args.phone
        if not phone:
            print("No phone_number available (profile has none and --phone "
                  "not given); Twist requires it for quotes.", file=sys.stderr)
            return 1

    # Phase 1: per-record heuristic pricing (and screening unless skipped)
    if args.quote_only:
        submitted = [r for r in records
                     if r.get("twist_construct_id")
                     and r.get("twist_score") == "BUILDABLE"]
        print(f"\nPhase 1 skipped (--quote-only). "
              f"Reusing {len(submitted)} buildable construct ids.")
    elif args.dry_run_heuristic_only:
        print(f"\nPhase 1 (heuristic only): computing for {len(records)} records...")
        submitted = submit_constructs(records, vendor=None)
    else:
        print(f"\nPhase 1: screening {len(records)} constructs...")
        submitted = submit_constructs(records, vendor)
        print(f"\n{len(submitted)}/{len(records)} buildable; "
              f"{len(records) - len(submitted)} skipped or unbuildable.")

    # Phase 2: bundled quote (skip in heuristic-only mode)
    if args.dry_run_heuristic_only or not submitted:
        print("\nSkipping bundled quote (heuristic-only or nothing buildable).")
    else:
        try:
            if args.existing_quote_id:
                print(f"\nFetching existing quote {args.existing_quote_id}...")
                quote = vendor.get_quote(args.existing_quote_id) or {}
                status = (quote.get("status_info") or {}).get("status", "")
                if status != "SUCCESS":
                    raise RuntimeError(
                        f"existing quote status is {status!r}, not SUCCESS"
                    )
            else:
                quote = request_quote(
                    vendor, submitted,
                    recipient_address_id=recipient_address_id,
                    first_name=first_name, last_name=last_name, phone=phone,
                )
            attach_quote_prices(submitted, quote)
            priced_n = sum(1 for r in submitted if "twist_price_usd" in r)
            print(f"\nQuote SUCCESS — {priced_n}/{len(submitted)} priced.")
        except Exception as exc:
            print(f"\nQuote failed: {exc!r}", file=sys.stderr)
            for rec in submitted:
                rec.setdefault("pricing_error", f"quote: {exc!r}")

    # Write annotated output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as f:
        for rec in records:
            f.write(json.dumps(rec) + "\n")
    print(f"\nWrote {len(records)} records to {args.output}")

    # Summary
    priced = [r for r in records if "twist_price_usd" in r]
    if priced:
        twist_total = sum(r["twist_price_usd"] for r in priced)
        heur_total = sum(r.get("heuristic_price_usd", 0.0) for r in priced)
        print(f"\nPriced: {len(priced)}")
        print(f"  Heuristic total: ${heur_total:,.2f}")
        print(f"  Twist total:     ${twist_total:,.2f}")
        if heur_total > 0:
            print(f"  Twist / Heuristic ratio: {twist_total/heur_total:.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main())
