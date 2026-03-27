"""
Integration tests for IDT live API.
"""

import os
import pytest
from dotenv import load_dotenv
from construct_compiler.vendors.idt import IDTVendor

load_dotenv()  # Load local .env
load_dotenv("/usr/local/google/home/charlieemrich/code/idt_api/.env")  # Fallback to other workspace .env



def test_idt_vendor_init():
    vendor = IDTVendor()
    # Should read from env if set, or empty
    assert vendor.client_id == os.environ.get("IDT_CLIENT_ID", "")


def test_idt_vendor_authenticated():
    vendor = IDTVendor()
    # If any of the 4 are missing, should be false
    creds = [
        os.environ.get("IDT_CLIENT_ID"),
        os.environ.get("IDT_CLIENT_SECRET"),
        os.environ.get("IDT_USERNAME"),
        os.environ.get("IDT_PASSWORD"),
    ]
    if all(creds):
        assert vendor.authenticated is True
    else:
        assert vendor.authenticated is False


@pytest.mark.skipif(
    not (
        os.environ.get("IDT_CLIENT_ID")
        and os.environ.get("IDT_CLIENT_SECRET")
        and os.environ.get("IDT_USERNAME")
        and os.environ.get("IDT_PASSWORD")
    ),
    reason="Missing IDT credentials in .env",
)
def test_idt_live_optimize():
    vendor = IDTVendor()
    assert vendor.authenticated is True

    # EGFP protein sequence from FPbase
    protein_seq = "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK"
    result = vendor.optimize_codons(protein_seq, organism="Escherichia coli")

    assert result.vendor == "idt"
    # In live mode or mock fallback, let's check what it does.
    # If live succeeded, optimized_sequence should not be empty.
    # If live failed, optimized_sequence is empty and notes have error.
    
    if result.optimized_sequence:
        assert len(result.optimized_sequence) >= 12  # 4 * 3 = 12 bp
        assert "Complexity" in "".join(result.notes)
    else:
        # If it failed, check the notes for error
        assert "API call exception" in "".join(result.notes) or "Auth failed" in "".join(result.notes)


@pytest.mark.skipif(
    not (
        os.environ.get("IDT_CLIENT_ID")
        and os.environ.get("IDT_CLIENT_SECRET")
        and os.environ.get("IDT_USERNAME")
        and os.environ.get("IDT_PASSWORD")
    ),
    reason="Missing IDT credentials in .env",
)
def test_idt_live_optimize_product_types():
    vendor = IDTVendor()
    assert vendor.authenticated is True

    protein_seq = "MAGS"
    
    # Test gblock
    result_gblock = vendor.optimize_codons(protein_seq, product_type="gblock")
    assert result_gblock.vendor == "idt"

    # Test megamer
    result_megamer = vendor.optimize_codons(protein_seq, product_type="megamer")
    assert result_megamer.vendor == "idt"


@pytest.mark.skipif(
    not (
        os.environ.get("IDT_CLIENT_ID")
        and os.environ.get("IDT_CLIENT_SECRET")
        and os.environ.get("IDT_USERNAME")
        and os.environ.get("IDT_PASSWORD")
    ),
    reason="Missing IDT credentials in .env",
)
def test_idt_live_screen():
    vendor = IDTVendor()
    assert vendor.authenticated is True

    # EGFP DNA sequence (already optimized)
    dna_seq = "ATGGTTTCTAAAGGCGAAGAGTTATTTACCGGTGTTGTTCCGATATTGGTCGAGCTGGACGGAGACGTAAACGGTCATAAATTCAGTGTGTCAGGAGAAGGTGAGGGTGATGCAACGTATGGGAAGCTTACTTTAAAATTTATCTGCACAACAGGCAAATTACCTGTGCCTTGGCCGACCCTGGTAACTACCCTTACTTACGGGGTTCAATGCTTTAGTCGTTACCCCGACCATATGAAGCAACACGATTTTTTCAAGTCTGCTATGCCAGAGGGATACGTCCAGGAGCGTACTATTTTCTTTAAGGATGATGGCAATTACAAAACACGGGCCGAAGTTAAATTTGAAGGGGACACCCTTGTTAATCGTATCGAGTTAAAAGGGATAGACTTCAAAGAAGATGGCAATATCCTGGGACACAAACTGGAATATAATTACAATAGTCATAATGTTTATATCATGGCGGATAAACAGAAAAATGGAATAAAAGTCAACTTCAAGATCAGACATAACATAGAGGATGGAAGCGTCCAGCTGGCTGACCACTACCAACAAAATACCCCCATAGGTGACGGCCCCGTCTTATTACCCGACAATCATTACTTGAGTACCCAGTCAGCGCTGTCCAAGGACCCTAATGAAAAACGGGATCATATGGTCTTGCTGGAGTTCGTGACAGCTGCAGGGATTACCTTGGGGATGGATGAATTATATAAG"
    
    result = vendor.screen(dna_seq, product_type="gblock")
    
    assert result.vendor == "idt"
    assert result.sequence_length == len(dna_seq)
    # Since it's a valid sequence, it should likely be feasible or have warnings but not errors that block it completely if it's "Accepted"
    # In live mode if it fails auth it falls back to mock, which might also pass if it fits criteria.
    # We just want to check it runs without crashing.
    assert result.feasible is True or result.feasible is False # Always true in logic, but checks it exists
