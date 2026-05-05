"""
Twist Bioscience vendor plugin — full TAPI coverage.

Implements every endpoint in the Twist Application Programming Interface (TAPI)
v1 OpenAPI spec, organised by logical group:

    Account Info          get_user, list_addresses, create_address,
                          list_payment_methods, create_payment_method
    Vectors               list_vectors, get_vector
    Codon optimisation    get_codon_optimization_choices, reverse_translate,
                          optimize_codons  (chains reverse-translation → codon-fitting)
    Constructs            _create_construct, describe_constructs,
                          _bulk_retrieve_construct, screen
    Quotes                create_quote, get_quote
    Orders                create_order, list_orders, get_order_items
    Plate maps            list_plate_maps, get_shipment_plate_maps,
                          get_plate_map_by_barcode
    Certificate of Analysis  request_coa_download, get_coa_download_url

Order placement (quotes → orders) is fully supported but intentionally separated
from the higher-level synthesis workflow.  Use ordering methods only after
constructs have been screened (score == "BUILDABLE") and the quote has reached
status == "SUCCESS".

API credentials
---------------
Set the following environment variables, or pass them to the constructor:

    TWIST_JWT_TOKEN          Authorization header ("JWT " prefix added automatically)
    TWIST_END_USER_TOKEN     X-End-User-Token header
    TWIST_USER_EMAIL         Scopes all /v1/users/{email}/ paths

Twist also whitelists the requesting IP; contact b2b-support@twistbioscience.com
to register your static IPs before use.

Use ``sandbox=True`` (or set ``TWIST_SANDBOX=1``) to point at the staging
environment (https://twist-api.twistbioscience-staging.com).

Docs: https://developers.twistdna.com/docs/tapi/
"""

from __future__ import annotations

import os
import time
import logging
from typing import Optional

from .base import VendorPlugin, ScreeningResult, OptimizationResult

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Scoring issue codes (score_data.issues[].code)
# ---------------------------------------------------------------------------

ISSUE_MESSAGES: dict[int, str] = {
    4000: "General error",
    4001: "Problematic sequence",
    4002: "Repeats or extreme high/low GC in the highlighted region may have caused this problem",
    4003: "Provided sequence is invalid",
    4004: "Provided sequence contains invalid character(s)",
    4100: "Invalid sequence length",
    4101: "Sequence is too short",
    4102: "Sequence is too long",
    4103: "Sequences longer than 1,700 bp elevate risk marginally",
    4104: "Sequences longer than 3,100 bp elevate risk",
    4200: "Invalid GC content",
    4201: "Overall GC content must be under 65% (under 60% optimal)",
    4202: "Overall GC content must be over 25%",
    4203: "GC content delta between highest- and lowest-GC 50 bp windows exceeds 52%",
    4300: "Secondary structure",
    4301: "Hairpin detected",
    4303: "Long direct repeat (≥20 bp) detected",
    4304: "Direct repeat with high Tm (>60 °C) detected",
    4305: ">45% of sequence composed of small repeats (≥9 bp)",
    4306: "Long, low-homology repeat region detected",
    4401: "Hierarchical design failed: no acceptable split point",
    4402: "Cannot make sequence as-is — try codon optimisation or split",
    4403: "Hierarchical design failed: no acceptable split point (retry)",
    4404: "Fragment sizes <300 bp; retry with fewer fragments",
    4405: "Fragment sizes >1,800 bp; retry with more fragments",
    4501: "His-tag with 5+ identical codons increases complexity",
    4502: "CpG multimeric segments of 14+ bases increase complexity",
    4503: "Long homopolymer stretches increase complexity",
    4504: "Sequence contains Gateway cloning att site(s)",
    4505: "Sequence contains impermissible sub-sequence(s)",
    4506: "Clonal gene contains sub-sequence with high homology to ccdB",
    4507: "Methylation-sensitive enzyme site overlaps a known methylation site",
    4508: "Unable to design primers — increase GC in first/last 60 bases",
    5001: "Unable to split fragment for synthesis",
    5002: "Fragment design exception",
}


# ---------------------------------------------------------------------------
# TwistVendor
# ---------------------------------------------------------------------------

class TwistVendor(VendorPlugin):
    """Twist Bioscience TAPI integration — complete v1 endpoint coverage."""

    PRODUCTION_URL = "https://twist-api.twistdna.com/v1"
    SANDBOX_URL    = "https://twist-api.twistbioscience-staging.com/v1"

    # Polling / timeout defaults (also used as method-signature defaults below,
    # which is valid because class-body names are in scope for default exprs).
    DEFAULT_POLL_INTERVAL  = 3.0    # seconds between status polls
    DEFAULT_POLL_TIMEOUT   = 300.0  # 5 min — codon-opt jobs can be slow
    SCORE_POLL_TIMEOUT     = 180.0  # bulk-retrieve scoring is usually faster
    DEFAULT_REQUEST_TIMEOUT = 60.0

    DEFAULT_TURNAROUND_DAYS = (12, 18)

    # Map Twist difficulty strings → 0–1 complexity score for ScreeningResult.
    DIFFICULTY_SCORES: dict[str, float] = {
        "STANDARD":    0.2,
        "MODERATE":    0.4,
        "DIFFICULT":   0.6,
        "COMPLEX":     0.7,
        "VERY_COMPLEX": 0.85,
    }

    # Plugin-internal product_type → Twist construct type string.
    PRODUCT_TYPE_MAP: dict[str, str] = {
        "GENE":           "NON_CLONED_GENE",
        "NON_CLONED_GENE": "NON_CLONED_GENE",
        "CLONED_GENE":    "CLONED_GENE",
        "OLIGO_POOL":     "OLIGO_POOL",
    }

    # ---------------------------------------------------------------------------
    # Construction
    # ---------------------------------------------------------------------------

    def __init__(
        self,
        jwt_token: Optional[str] = None,
        end_user_token: Optional[str] = None,
        user_email: Optional[str] = None,
        sandbox: bool = False,
    ):
        self.jwt_token      = jwt_token      or os.environ.get("TWIST_JWT_TOKEN", "")
        self.end_user_token = end_user_token or os.environ.get("TWIST_END_USER_TOKEN", "")
        self.user_email     = user_email     or os.environ.get("TWIST_USER_EMAIL", "")
        use_sandbox         = sandbox or bool(os.environ.get("TWIST_SANDBOX", ""))
        self.BASE_URL       = self.SANDBOX_URL if use_sandbox else self.PRODUCTION_URL

    # ---------------------------------------------------------------------------
    # VendorPlugin interface
    # ---------------------------------------------------------------------------

    @property
    def name(self) -> str:
        return "Twist Bioscience"

    @property
    def authenticated(self) -> bool:
        return bool(self.jwt_token and self.end_user_token and self.user_email)

    # ---------------------------------------------------------------------------
    # HTTP plumbing
    # ---------------------------------------------------------------------------

    def _get_headers(self) -> dict:
        return {
            "Authorization":   f"JWT {self.jwt_token}",
            "X-End-User-Token": self.end_user_token,
            "Content-Type":    "application/json",
            "Accept":          "application/json",
        }

    def _user_url(self, *parts: str) -> str:
        """Build a /v1/users/{email}/… URL.

        >>> vendor._user_url("constructs", "bulk-retrieve")
        'https://.../v1/users/me@example.com/constructs/bulk-retrieve/'
        """
        if not self.user_email:
            raise RuntimeError(
                "TWIST_USER_EMAIL not configured — TAPI paths are scoped to "
                "/users/{email}/"
            )
        suffix = "/".join(parts)
        base   = f"{self.BASE_URL}/users/{self.user_email}"
        return f"{base}/{suffix}/" if suffix else f"{base}/"

    def _request(self, method: str, url: str, **kwargs) -> object:
        """Execute an HTTP request and return the parsed JSON body.

        Raises RuntimeError on 4xx/5xx.  Returns None for empty bodies.
        """
        import requests
        kwargs.setdefault("headers", self._get_headers())
        kwargs.setdefault("timeout", self.DEFAULT_REQUEST_TIMEOUT)
        resp = requests.request(method, url, **kwargs)
        if resp.status_code >= 400:
            raise RuntimeError(
                f"Twist {method} {url} → {resp.status_code}: {resp.text[:400]}"
            )
        return resp.json() if resp.content else None

    # ---------------------------------------------------------------------------
    # Account Info
    # ---------------------------------------------------------------------------

    def get_user(self) -> dict:
        """GET /users/{email}/ — fetch user profile and account ID."""
        return self._request("GET", self._user_url()) or {}

    def list_addresses(self) -> list[dict]:
        """GET /users/{email}/addresses/ — list all shipping and billing addresses.

        Useful for discovering ``recipient_address_id`` values to pass to
        ``create_quote()``.  Each entry includes ``id``, ``address_type``
        ("Shipping"/"Billing"), ``verification_status``, ``is_default``, and
        the full address fields.
        """
        result = self._request("GET", self._user_url("addresses"))
        return result if isinstance(result, list) else ([result] if result else [])

    def create_address(
        self,
        address_type: str,          # "Shipping" or "Billing"
        first_name: str = "",
        last_name: str = "",
        street_1: str = "",
        street_2: str = "",
        city: str = "",
        state: str = "",
        zip_code: str = "",
        country: str = "US",        # ISO-3166-1 alpha-2
        organization: str = "",
        phone_number: str = "",
        billing_email: str = "",
    ) -> dict:
        """POST /users/{email}/addresses/ — create a shipping or billing address.

        Notes
        -----
        * Addresses are validated via Google Maps; fake addresses will be rejected.
        * The created address starts with ``verification_status: PENDING_REVIEW``.
        * ``country`` must be a 2-letter ISO code (e.g. "US", "GB").
        """
        body: dict = {"address_type": address_type, "country": country}
        for key, val in [
            ("first_name",    first_name),
            ("last_name",     last_name),
            ("street_1",      street_1),
            ("street_2",      street_2),
            ("city",          city),
            ("state",         state),
            ("zip_code",      zip_code),
            ("organization",  organization),
            ("phone_number",  phone_number),
            ("billing_email", billing_email),
        ]:
            if val:
                body[key] = val
        return self._request("POST", self._user_url("addresses"), json=body) or {}

    def list_payment_methods(self) -> list[dict]:
        """GET /users/{email}/payments/ — list all PO payment methods."""
        result = self._request("GET", self._user_url("payments"))
        return result if isinstance(result, list) else ([result] if result else [])

    def create_payment_method(
        self,
        po_type: str,                           # "Regular PO" or "Blanket PO"
        purchase_order_reference: str,
        billing_address_id: str,
        currency_iso_code: str = "USD",
        available_po_balance: str = "0.00",
        starting_po_balance: str = "0.00",
        transient_po_document_url: str = "",
        sent_via: str = "EMAIL",                # "EMAIL" or "UPLOAD"
    ) -> dict:
        """POST /users/{email}/payments/ — register a new PO payment method.

        Credit cards are not supported via the API (use the eCommerce web app).
        ``sent_via`` should be "UPLOAD" when providing a ``transient_po_document_url``.
        """
        body: dict = {
            "type":                     po_type,
            "purchase_order_reference": purchase_order_reference,
            "billing_address":          billing_address_id,
            "currency_iso_code":        currency_iso_code,
            "available_po_balance":     available_po_balance,
            "starting_po_balance":      starting_po_balance,
            "sent_via":                 sent_via,
        }
        if transient_po_document_url:
            body["transient_po_document_url"] = transient_po_document_url
        return self._request("POST", self._user_url("payments"), json=body) or {}

    # ---------------------------------------------------------------------------
    # Vectors
    # ---------------------------------------------------------------------------

    def list_vectors(self) -> list[dict]:
        """GET /users/{email}/vectors/ — list catalog and custom vectors."""
        result = self._request("GET", self._user_url("vectors"))
        if isinstance(result, dict):
            return [result]
        return result or []

    def get_vector(self, vector_id: str) -> dict:
        """GET /users/{email}/vectors/{id}/ — fetch a single vector with insertion sites."""
        return self._request("GET", self._user_url("vectors", vector_id)) or {}

    # ---------------------------------------------------------------------------
    # Codon optimisation — choices
    # ---------------------------------------------------------------------------

    def get_codon_optimization_choices(self) -> dict:
        """GET /users/{email}/codon-optimizations/choices/ — valid organism names
        and restriction enzymes that can be passed to codon-optimisation jobs.
        """
        return self._request("GET", self._user_url("codon-optimizations", "choices")) or {}

    # ---------------------------------------------------------------------------
    # Async poll helper
    # ---------------------------------------------------------------------------

    def _poll_async_job(
        self,
        list_url: str,
        job_id: str,
        interval: float = DEFAULT_POLL_INTERVAL,
        timeout: float  = DEFAULT_POLL_TIMEOUT,
    ) -> dict:
        """Poll a Twist async-job list endpoint until the job with ``job_id``
        has ``completed: true``.

        ``id__in`` is passed as a required query parameter so the server filters
        the response to just our job rather than returning every job for the user.
        """
        deadline = time.time() + timeout
        while time.time() < deadline:
            jobs = self._request("GET", list_url, params={"id__in": job_id}) or []
            for job in jobs:
                if job.get("id") != job_id:
                    continue
                if job.get("error"):
                    raise RuntimeError(
                        f"Twist job {job_id} failed: {job['error']}"
                    )
                if job.get("completed"):
                    return job
                break   # job exists but not yet complete — keep polling
            time.sleep(interval)
        raise TimeoutError(
            f"Twist job {job_id} did not complete within {timeout:.0f} s"
        )

    # ---------------------------------------------------------------------------
    # Reverse translation  (protein → DNA)
    # ---------------------------------------------------------------------------

    def reverse_translate(
        self,
        protein_sequence: str,
        organism: str = "Homo sapiens",
        is_igg: bool = False,
        external_id: Optional[str] = None,
        callback_url: Optional[str] = None,
    ) -> str:
        """Submit a protein → DNA reverse-translation job and block until done.

        Returns the codon-optimised DNA string (upper-case ACTG).

        Parameters
        ----------
        protein_sequence:
            Single-letter amino-acid sequence.
        organism:
            Target expression host (use ``get_codon_optimization_choices()``
            for valid values).
        is_igg:
            Set True for IgG antibody sequences (enables Twist's antibody
            codon table).
        external_id:
            Optional caller-supplied ID echoed back in the response.
        callback_url:
            If provided the API POSTs results to this URL instead of requiring
            polling (the method still polls as a fallback safety net).
        """
        url = self._user_url("reverse-translations")
        body: list[dict] = [{
            "sequence": protein_sequence,
            "organism": organism,
            "is_igg":   is_igg,
        }]
        if external_id:
            body[0]["external_id"] = external_id

        params = {"callback_url": callback_url} if callback_url else {}
        resp = self._request("POST", url, json=body, params=params)
        if not resp or not isinstance(resp, list):
            raise RuntimeError(
                f"Unexpected reverse-translation POST response: {resp!r}"
            )
        job_id    = resp[0]["id"]
        completed = self._poll_async_job(url, job_id)
        dna = completed.get("translated_sequence")
        if not dna:
            raise RuntimeError(
                f"Reverse-translation job {job_id} completed without a sequence: "
                f"{completed}"
            )
        return dna

    # ---------------------------------------------------------------------------
    # Codon-fitting optimisation  (DNA → optimised DNA)
    # ---------------------------------------------------------------------------

    def _submit_codon_optimization(
        self,
        dna: str,
        organism: str,
        avoid_introducing: Optional[list[str]] = None,
        preserve: Optional[list] = None,
        optimization_start: int = 0,
        optimization_len: Optional[int] = None,
        adapters_on: bool = False,
        old_scoring: bool = False,
        external_id: Optional[str] = None,
    ) -> dict:
        """POST a codon-fitting optimisation job and poll until complete.

        Parameters
        ----------
        dna:
            Input DNA sequence (upper-case ACTG).  Automatically uppercased.
        organism:
            Target expression host (see ``get_codon_optimization_choices()``).
        avoid_introducing:
            List of restriction enzyme names whose sites must not be introduced
            (e.g. ``["EcoRI", "XbaI"]``).
        preserve:
            List of [start, end] pairs (0-based) marking regions whose codons
            must not be changed (e.g. ``[[0, 59], [300, 359]]``).
        optimization_start:
            0-based index of the first base to optimise.
        optimization_len:
            Number of bases to optimise; defaults to the full sequence length.
        adapters_on:
            Whether Twist should add synthesis adapters.
        old_scoring:
            Use the legacy scoring algorithm (usually False).
        external_id:
            Optional caller ID echoed in the response.

        Returns
        -------
        The completed job dict, which includes ``optimal_result`` with the
        optimised sequence and ``score_data``.
        """
        url = self._user_url("codon-optimizations")
        item: dict = {
            "sequence":           dna.upper(),
            "organism":           organism,
            "avoid_introducing":  avoid_introducing,        # None → server default
            "preserve":           preserve if preserve is not None else [],
            "optimization_start": optimization_start,
            "optimization_len":   optimization_len if optimization_len is not None else len(dna),
            "adapters_on":        adapters_on,
            "old_scoring":        old_scoring,
        }
        if external_id:
            item["external_id"] = external_id

        resp = self._request("POST", url, json=[item])
        if not resp or not isinstance(resp, list):
            raise RuntimeError(
                f"Unexpected codon-optimisation POST response: {resp!r}"
            )
        job_id = resp[0]["id"]
        return self._poll_async_job(url, job_id)

    def optimize_codons(
        self,
        protein_sequence: str,
        organism: str = "Homo sapiens",
        **kwargs,
    ) -> OptimizationResult:
        """Reverse-translate a protein sequence then codon-optimise the DNA.

        This chains two Twist async jobs:
          1. Reverse translation (protein → raw DNA, codon table for *organism*)
          2. Codon-fitting optimisation (DNA → manufacturability-optimised DNA)

        Falls back to a stub ``OptimizationResult`` (``optimized_sequence=""``)
        when credentials are missing or the API fails, so callers can detect
        the empty output and apply local optimisation as a fallback.

        Keyword arguments (all optional)
        ----------------------------------
        is_igg, external_id, avoid_introducing, preserve,
        optimization_start, optimization_len, adapters_on, old_scoring
            Passed through to the respective TAPI calls (see
            ``reverse_translate`` and ``_submit_codon_optimization``).
        """
        if not self.authenticated:
            return OptimizationResult(
                vendor="twist",
                original_sequence=protein_sequence,
                optimized_sequence="",
                notes=["Twist credentials not configured; use local codon optimisation."],
            )
        try:
            dna = self.reverse_translate(
                protein_sequence,
                organism=organism,
                is_igg=kwargs.get("is_igg", False),
                external_id=kwargs.get("external_id"),
            )
            optimised = self._submit_codon_optimization(
                dna,
                organism=organism,
                avoid_introducing=kwargs.get("avoid_introducing"),
                preserve=kwargs.get("preserve"),
                optimization_start=kwargs.get("optimization_start", 0),
                optimization_len=kwargs.get("optimization_len"),
                adapters_on=kwargs.get("adapters_on", False),
                old_scoring=kwargs.get("old_scoring", False),
                external_id=kwargs.get("external_id"),
            )
            optimal      = optimised.get("optimal_result") or {}
            optimal_seq  = optimal.get("sequence", "")
            score_data   = optimal.get("score_data") or {}
            metrics      = score_data.get("scoring_metrics") or {}
            gc           = float(metrics.get("overall_gc_percent", 0.0))
            part         = score_data.get("part_number", "unknown")
            difficulty   = score_data.get("difficulty", "unknown")
            return OptimizationResult(
                vendor="twist",
                original_sequence=protein_sequence,
                optimized_sequence=optimal_seq,
                changes_made=0,         # Twist does not expose a per-codon delta
                gc_content=gc,
                notes=self._build_optimization_notes(
                    organism, part, difficulty, metrics, score_data,
                ),
            )
        except Exception as exc:
            logger.warning("Twist optimize_codons failed: %s", exc)
            return OptimizationResult(
                vendor="twist",
                original_sequence=protein_sequence,
                optimized_sequence="",
                notes=[f"Twist API error: {exc}"],
            )

    @staticmethod
    def _build_optimization_notes(
        organism: str,
        part: str,
        difficulty: str,
        metrics: dict,
        score_data: dict,
    ) -> list[str]:
        """Render Twist score_data into a flat list of human-readable notes."""
        notes: list[str] = [
            f"Reverse-translated then codon-optimised for {organism}.",
            f"Twist part: {part} ({difficulty}).",
        ]
        if not metrics:
            return notes

        gc        = metrics.get("overall_gc_percent")
        gc_delta  = metrics.get("trimmed_gc_win_50_min_max_delta")
        max_homo  = metrics.get("max_homopolymer_length")
        max_rep   = metrics.get("max_long_repeat_length")
        rep_homo  = metrics.get("max_long_repeat_homology")
        banned    = metrics.get("has_banned_sequence")
        warn_seq  = metrics.get("has_warning_sequence")

        if gc is not None:
            notes.append(f"GC content: {float(gc) * 100:.1f}%.")
        if gc_delta is not None:
            notes.append(f"GC delta (50 bp windows): {float(gc_delta) * 100:.1f}%.")
        if max_homo is not None:
            notes.append(f"Max homopolymer run: {max_homo} bp.")
        if max_rep:
            notes.append(
                f"Max long repeat: {max_rep} bp (homology {float(rep_homo or 0) * 100:.0f}%)."
            )
        if banned:
            notes.append("Contains banned sub-sequence(s).")
        if warn_seq:
            notes.append("Contains warning sub-sequence(s).")

        for issue in score_data.get("issues") or []:
            code = issue.get("code") if isinstance(issue, dict) else None
            msg  = ISSUE_MESSAGES.get(code) if code else None
            if not msg and isinstance(issue, dict):
                msg = issue.get("message") or issue.get("text")
            notes.append(f"Issue [{code}] {msg}" if code else f"Issue: {issue}")

        return notes

    # ---------------------------------------------------------------------------
    # Constructs
    # ---------------------------------------------------------------------------

    def _create_construct(
        self,
        sequence: str,
        name: str,
        construct_type: str,
        vector_mes_uid: Optional[str] = None,
        insertion_point_mes_uid: Optional[str] = None,
        adapters_on: bool = False,
        external_id: Optional[str] = None,
    ) -> str:
        """POST /constructs/ — create a construct and return its UUID.

        For ``CLONED_GENE`` constructs both ``vector_mes_uid`` and
        ``insertion_point_mes_uid`` are required; obtain them from
        ``list_vectors()``.

        Scoring is asynchronous — poll with ``describe_constructs()`` or
        ``_bulk_retrieve_construct()``.
        """
        url  = self._user_url("constructs")
        body: dict = {
            "sequences":   [sequence.upper()],
            "name":        name,
            "type":        construct_type,
            "adapters_on": adapters_on,
        }
        if construct_type == "CLONED_GENE":
            if not (vector_mes_uid and insertion_point_mes_uid):
                raise ValueError(
                    "CLONED_GENE constructs require vector_mes_uid and "
                    "insertion_point_mes_uid (use list_vectors() to look these up)"
                )
            body["vector_mes_uid"]          = vector_mes_uid
            body["insertion_point_mes_uid"] = insertion_point_mes_uid
        if external_id:
            body["external_id"] = external_id

        resp = self._request("POST", url, json=body)
        if not isinstance(resp, dict) or "id" not in resp:
            raise RuntimeError(
                f"Unexpected POST /constructs/ response: {resp!r}"
            )
        return resp["id"]

    def describe_constructs(
        self,
        construct_ids: list[str],
        scored: Optional[bool] = None,
    ) -> list[dict]:
        """GET /constructs/describe/ — poll scoring results for one or more constructs.

        Parameters
        ----------
        construct_ids:
            List of construct UUIDs to retrieve.
        scored:
            Pass ``True`` to filter to only constructs whose scoring is complete.
            Omit to return all constructs regardless of scoring state.

        Returns
        -------
        List of construct dicts with ``score``, ``score_data``, ``scored``, etc.
        """
        params: dict = {"id__in": ",".join(construct_ids)}
        if scored is not None:
            params["scored"] = scored
        result = self._request(
            "GET", self._user_url("constructs", "describe"), params=params
        )
        return result if isinstance(result, list) else ([result] if result else [])

    def _bulk_retrieve_construct(
        self,
        construct_id: str,
        interval: float = DEFAULT_POLL_INTERVAL,
        timeout: float  = SCORE_POLL_TIMEOUT,
    ) -> dict:
        """POST /constructs/bulk-retrieve/ — poll until a construct has scored.

        Blocks until ``scored: true`` or *timeout* seconds elapse.
        """
        url      = self._user_url("constructs", "bulk-retrieve")
        deadline = time.time() + timeout
        while time.time() < deadline:
            results = self._request(
                "POST", url, json={"construct_ids": [construct_id]}
            ) or []
            if results and results[0].get("scored"):
                return results[0]
            time.sleep(interval)
        raise TimeoutError(
            f"Construct {construct_id} did not score within {timeout:.0f} s"
        )

    def screen(
        self,
        sequence: str,
        product_type: str = "GENE",
        **kwargs,
    ) -> ScreeningResult:
        """Submit a sequence as a Twist Construct and return its manufacturability score.

        Blocks until scoring is complete.  Falls back to a local heuristic
        ``ScreeningResult`` when credentials are absent or the API fails.

        Keyword arguments
        -----------------
        construct_type : str
            Explicit Twist type ("NON_CLONED_GENE", "CLONED_GENE", "OLIGO_POOL");
            takes precedence over *product_type*.
        name : str
            Human-readable construct name (auto-generated if omitted).
        adapters_on : bool
            Whether to request synthesis adapters.
        external_id : str
            Optional caller ID.
        vector_mes_uid, insertion_point_mes_uid : str
            Required for CLONED_GENE.
        """
        length = len(sequence)
        if not self.authenticated:
            return self._mock_screen(sequence, product_type)
        try:
            construct_type = (
                kwargs.get("construct_type")
                or self.PRODUCT_TYPE_MAP.get(product_type, "NON_CLONED_GENE")
            )
            name   = kwargs.get("name") or f"screen_{int(time.time())}"
            cid    = self._create_construct(
                sequence,
                name=name,
                construct_type=construct_type,
                vector_mes_uid=kwargs.get("vector_mes_uid"),
                insertion_point_mes_uid=kwargs.get("insertion_point_mes_uid"),
                adapters_on=kwargs.get("adapters_on", False),
                external_id=kwargs.get("external_id"),
            )
            scored     = self._bulk_retrieve_construct(cid)
            score_data = scored.get("score_data") or {}
            difficulty = score_data.get("difficulty", "")
            issues_raw = score_data.get("issues") or []

            warnings: list[str] = []
            errors:   list[str] = []
            for issue in issues_raw:
                code = issue.get("code") if isinstance(issue, dict) else None
                msg  = ISSUE_MESSAGES.get(code, "") if code else ""
                if not msg and isinstance(issue, dict):
                    msg = issue.get("message") or issue.get("text") or str(issue)
                rendered = f"[{code}] {msg}" if code else (msg or str(issue))
                # 5xxx codes = fatal synthesis errors
                (errors if isinstance(code, int) and code >= 5000 else warnings).append(rendered)

            score    = scored.get("score", "")
            feasible = (score == "BUILDABLE") and not errors
            return ScreeningResult(
                vendor="twist",
                feasible=feasible,
                sequence_length=length,
                estimated_price=0.0,    # real pricing requires a quote
                turnaround_days=self.DEFAULT_TURNAROUND_DAYS,
                complexity_score=self.DIFFICULTY_SCORES.get(difficulty, 0.5),
                warnings=warnings,
                errors=errors,
            )
        except Exception as exc:
            logger.warning(
                "Twist screen failed: %s — falling back to mock heuristics", exc
            )
            return self._mock_screen(sequence, product_type)

    # ---------------------------------------------------------------------------
    # Quotes
    # ---------------------------------------------------------------------------

    def create_quote(
        self,
        external_id: str,
        containers: list[dict],
        order_sub_product_type: str,
        recipient_address_id: str,
        first_name: str,
        last_name: str,
        phone: str,
        order_settings: Optional[list[dict]] = None,
        ecommerce_project_name: Optional[str] = None,
        cc_emails: Optional[str] = None,
    ) -> dict:
        """POST /users/{email}/quotes/ — create a quote for a set of constructs.

        Parameters
        ----------
        external_id:
            Must be unique per quote per owner.  Supplying the same
            ``external_id`` again deletes the previous quote and creates a new one.
        containers:
            List of container objects.  Each must contain ``constructs``, a list
            of ``{"id": construct_uuid, "index": int}`` dicts::

                containers=[{
                    "constructs": [
                        {"id": "uuid-1", "index": 1},
                        {"id": "uuid-2", "index": 2},
                    ]
                }]

        order_sub_product_type:
            One of:
              - ``"NON_CLONAL_ADAPTERS_OFF"``  — gene fragments without adapters
              - ``"NON_CLONAL_ADAPTERS_ON"``   — gene fragments with adapters
              - ``"OLIGO_POOLS_REGULAR"``       — oligo pools (max 6 constructs)
              - ``"CLONAL_GENES_SHORT"``        — clonal genes

        recipient_address_id:
            Shipping address ID from ``create_address()`` /
            ``get_user()["default_shipping_address"]``.
        order_settings:
            List of delivery-format / add-on service dicts, e.g.::

                [{"name": "Delivery Format", "product_code": "SER_PKG_TUBE",
                  "configuration": {"fill_method": "Vertical"}}]

        ecommerce_project_name:
            Optional project title shown in the quote PDF.
        cc_emails:
            Additional notification recipients, newline-separated
            (``"user1@co.com\\r\\nuser2@co.com"``).

        Returns
        -------
        ``{"id": quote_uuid}`` — use the ID with ``get_quote()`` to poll status.
        """
        body: dict = {
            "external_id": external_id,
            "shipment": {
                "first_name":           first_name,
                "last_name":            last_name,
                "phone":                phone,
                "recipient_address_id": recipient_address_id,
            },
            "containers":            containers,
            "order_sub_product_type": order_sub_product_type,
        }
        if order_settings:
            body["order_settings"] = order_settings
        if ecommerce_project_name:
            body["ecommerce_project_name"] = ecommerce_project_name
        if cc_emails:
            body["cc_emails"] = cc_emails

        return self._request("POST", self._user_url("quotes"), json=body) or {}

    def get_quote(self, quote_id: str) -> dict:
        """GET /users/{email}/quotes/{id}/ — retrieve quote details and status.

        The ``status_info.status`` field indicates progress:
          - ``"IN_PROGRESS"`` — quote is being generated
          - ``"SUCCESS"``     — ready to convert to an order
          - ``"ERROR"``       — generation failed; see ``status_info.error``

        A quote PDF download link is available at ``pdf_download_link`` once the
        status is ``"SUCCESS"``.  Turnaround time is in ``tat.business_days``.
        """
        return self._request("GET", self._user_url("quotes", quote_id)) or {}

    # ---------------------------------------------------------------------------
    # Orders
    # ---------------------------------------------------------------------------

    def create_order(
        self,
        quote_id: str,
        payment_method_id: Optional[str],
        po_reference: Optional[str] = None,
        payment_flow: Optional[str] = None,
    ) -> dict:
        """POST /users/{email}/orders/ — place an order from a completed quote.

        The quote must have ``status_info.status == "SUCCESS"`` before calling
        this method.

        Parameters
        ----------
        quote_id:
            UUID of the quote to order.
        payment_method_id:
            ID from ``list_payment_methods()``.  Pass ``None`` together with
            ``payment_flow="NO_PO"`` to submit without a PO (synthesis begins
            after a valid PO is later received).
        po_reference:
            Optional internal reference (max 29 chars) shown on invoices.
        payment_flow:
            Pass ``"NO_PO"`` to submit without an immediate purchase order.

        Returns
        -------
        ``{"id": sfdc_order_id, "currency_iso_code": ..., ...}``
        """
        body: dict = {
            "quote_id":          quote_id,
            "payment_method_id": payment_method_id,
        }
        if po_reference:
            body["po_reference"] = po_reference
        if payment_flow:
            body["payment_flow"] = payment_flow
        return self._request("POST", self._user_url("orders"), json=body) or {}

    def list_orders(
        self,
        order_status: Optional[str] = None,
        page: Optional[int] = None,
        page_size: Optional[int] = None,
        sort_by: Optional[str] = None,
        reverse: Optional[bool] = None,
    ) -> list[dict]:
        """GET /users/{email}/orders/ — list all orders with shipment summaries.

        Parameters
        ----------
        order_status:
            Filter by status: ``"pending"``, ``"open"``, ``"past"``,
            ``"cancelled"``, or ``"unknown"``.
        page, page_size:
            Pagination controls.
        sort_by:
            Field to sort by, e.g. ``"received_date"``.
        reverse:
            ``True`` for descending sort.

        Note
        ----
        The Twist API accepts these filters as *headers*, not query parameters.
        """
        extra_headers: dict = {}
        if order_status is not None:
            extra_headers["order_status"] = order_status
        if page is not None:
            extra_headers["page"] = str(page)
        if page_size is not None:
            extra_headers["page_size"] = str(page_size)
        if sort_by is not None:
            extra_headers["sort_by"] = sort_by
        if reverse is not None:
            extra_headers["reverse"] = str(reverse).lower()

        headers = {**self._get_headers(), **extra_headers}
        result  = self._request("GET", self._user_url("orders"), headers=headers)
        return result if isinstance(result, list) else ([result] if result else [])

    def get_order_items(self, sfdc_id: str) -> dict:
        """GET /users/{email}/orders/{sfdc_id}/items/ — items and status for one order.

        ``sfdc_id`` is the Salesforce order ID returned by ``list_orders()``
        (the ``sfdc_id`` field, not ``id``).

        Returns a rich dict with ``order_items``, ``shipments``, invoice details,
        item-level status counts, and shipping / billing addresses.
        """
        return self._request("GET", self._user_url("orders", sfdc_id, "items")) or {}

    # ---------------------------------------------------------------------------
    # Plate maps
    # ---------------------------------------------------------------------------

    def list_plate_maps(
        self,
        page: Optional[int] = None,
        page_size: Optional[int] = None,
    ) -> list[dict]:
        """GET /users/{email}/platemaps/ — plate-map index for all past shipments.

        Only relevant for orders delivered in plate format (96-well, 384-echo,
        96 PCR, or 96 deep-well).  Returns order/container/shipment IDs,
        delivery dates, product types, and item counts — use those IDs with
        ``get_shipment_plate_maps()`` or ``get_plate_map_by_barcode()``.

        Pagination is passed as headers.
        """
        extra_headers: dict = {}
        if page is not None:
            extra_headers["page"] = str(page)
        if page_size is not None:
            extra_headers["page_size"] = str(page_size)

        headers = {**self._get_headers(), **extra_headers}
        result  = self._request("GET", self._user_url("platemaps"), headers=headers)
        return result if isinstance(result, list) else ([result] if result else [])

    def get_shipment_plate_maps(
        self,
        order_id: str,
        shipment_id: str,
        file_name: Optional[str] = None,
        store_as_file: Optional[str] = None,
    ) -> str:
        """GET /users/{email}/orders/{id}/shipments/{shipment_id}/plate-maps/ —
        download link for a ZIP of CSV plate maps for an entire shipment.

        Returns the pre-signed S3 URL string (``platemaps_file_url``).
        """
        params: dict = {}
        if file_name:
            params["file_name"] = file_name
        if store_as_file:
            params["store_as_file"] = store_as_file

        url    = self._user_url("orders", order_id, "shipments", shipment_id, "plate-maps")
        result = self._request("GET", url, params=params or None) or {}
        return result.get("platemaps_file_url", "")

    def get_plate_map_by_barcode(
        self,
        order_id: str,
        barcode: str,
        file_name: Optional[str] = None,
        store_as_file: Optional[str] = None,
    ) -> str:
        """GET /users/{email}/orders/{id}/plate-maps/{barcode}/ —
        download link for a single-plate CSV identified by barcode.

        ``barcode`` (``container_id``) is available from ``list_plate_maps()``.
        Returns the pre-signed S3 URL string (``platemaps_file_url``).
        """
        params: dict = {}
        if file_name:
            params["file_name"] = file_name
        if store_as_file:
            params["store_as_file"] = store_as_file

        url    = self._user_url("orders", order_id, "plate-maps", barcode)
        result = self._request("GET", url, params=params or None) or {}
        return result.get("platemaps_file_url", "")

    # ---------------------------------------------------------------------------
    # Certificate of Analysis (CoA)
    # ---------------------------------------------------------------------------

    def request_coa_download(self, order_id: str, shipment_id: str) -> dict:
        """POST …/coa/request-download — trigger CoA document generation.

        This is an asynchronous request.  After calling this method, poll
        ``get_coa_download_url()`` until ``status`` indicates the document is ready.

        Returns the raw API response dict (``{"download_status": ...}``).
        """
        # Note: this endpoint has no trailing slash per the spec.
        url = (
            f"{self.BASE_URL}/users/{self.user_email}"
            f"/orders/{order_id}/shipments/{shipment_id}/coa/request-download"
        )
        return self._request("POST", url) or {}

    def get_coa_download_url(self, order_id: str, shipment_id: str) -> dict:
        """GET …/coa/download/ — retrieve the CoA download URL for a shipment.

        Call ``request_coa_download()`` first, then poll this endpoint until
        ``status`` signals completion.

        Returns
        -------
        dict with ``download_url``, ``status``, and ``file_size`` (bytes).
        """
        url = self._user_url(
            "orders", order_id, "shipments", shipment_id, "coa", "download"
        )
        return self._request("GET", url) or {}

    # ---------------------------------------------------------------------------
    # Mock fallback (used when unauthenticated or on API failure)
    # ---------------------------------------------------------------------------

    def _mock_screen(self, sequence: str, product_type: str) -> ScreeningResult:
        """Heuristic manufacturability estimate used when the API is unavailable."""
        length   = len(sequence)
        warnings: list[str] = []
        errors:   list[str] = []
        feasible = True

        if product_type == "GENE" and length > 5000:
            errors.append(
                f"Sequence length {length} bp exceeds gene synthesis limit (5000 bp)"
            )
            feasible = False

        seq_upper = sequence.upper()
        gc = (seq_upper.count("G") + seq_upper.count("C")) / max(length, 1)
        if gc < 0.25 or gc > 0.75:
            warnings.append(
                f"Extreme GC content ({gc * 100:.1f}%) may affect synthesis"
            )

        for base in "ATGC":
            if base * 11 in seq_upper:
                warnings.append(f"Long {base}-homopolymer (>10 bp) detected")

        price_per_bp = 0.07 if product_type == "GENE" else 0.09
        return ScreeningResult(
            vendor="twist",
            feasible=feasible,
            sequence_length=length,
            estimated_price=length * price_per_bp,
            turnaround_days=self.DEFAULT_TURNAROUND_DAYS,
            complexity_score=0.5 if warnings else 0.2,
            warnings=warnings,
            errors=errors,
        )
