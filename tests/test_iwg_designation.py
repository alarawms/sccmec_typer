"""IWG-SCC designations, composite candidates, ccr copy number and assembly-edge flag."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.classifier import classify_sccmec, CONTIG_EDGE_BP


def hit(gene, start=0, end=1000, strand="+", contig="contig_1", contig_len=2_800_000):
    return {"gene": gene, "accession": "T", "scc_type": "T", "identity": 0.99, "coverage": 0.95,
            "contig": contig, "contig_len": contig_len, "start": start, "end": end,
            "strand": strand, "id_pct": 99.0, "cov_pct": 95.0, "len": end - start,
            "aln_len": end - start}


# mec class C2: mecA + mecR1 with IS431 copies in opposite orientation
C2 = [hit("mecA", 20000, 22000), hit("mecR1", 22000, 23000),
      hit("IS431_1", 17000, 17800, "+"), hit("IS431_2", 31000, 31800, "-")]


def test_single_type_designation():
    r = classify_sccmec(C2 + [hit("ccrC1", 10000, 11700)])
    assert r["sccmec_type"] == "Type V"
    assert r["iwg_type"] == "V(5C2)" and r["iwg_designation"] == "5C2"
    assert r["ccr_copies"] == {"Type 5": 1}


def test_overlapping_allele_hits_are_one_locus():
    r = classify_sccmec(C2 + [hit("ccrC1", 10000, 11700), hit("ccrC1", 10040, 11700)])
    assert r["iwg_type"] == "V(5C2)"


def test_second_ccrC_copy_is_reported():
    """Type V with two ccrC loci, e.g. ST672 / V(5C2&5)."""
    r = classify_sccmec(C2 + [hit("ccrC1", 10000, 11700), hit("ccrC1", 26600, 28300)])
    assert r["iwg_type"] == "V(5C2&5)"
    assert r["ccr_copies"] == {"Type 5": 2}


def test_composite_lists_all_fitting_types():
    """Class B + ccr1 + ccr2 is type I or IV; must not pick the first rule silently."""
    hits = [hit("mecA", 20000, 22000), hit("mecR1", 22000, 23000), hit("IS1272", 23000, 25000),
            hit("ccrA1", 50000, 51300), hit("ccrB1", 51300, 52800),
            hit("ccrA2", 60000, 61300), hit("ccrB2", 61300, 62800)]
    r = classify_sccmec(hits)
    assert r["type_candidates"] == ["Type I", "Type IV"]
    assert r["sccmec_type"] == "Composite (Type I / Type IV)"
    assert r["iwg_type"] == "I(1B&2) / IV(2B&1)"


def test_composite_with_extra_ccr1():
    """Class C1 + ccrC1 + ccrA1B1 (ST1-like): VII(5C1&1)."""
    hits = [hit("mecA", 20000, 22000), hit("mecR1", 22000, 23000),
            hit("IS431_1", 17000, 17800, "+"), hit("IS431_2", 31000, 31800, "+"),
            hit("ccrC1", 10000, 11700), hit("ccrA1", 50000, 51300), hit("ccrB1", 51300, 52800)]
    r = classify_sccmec(hits)
    assert r["mec_complex"] == "Class C1"
    assert r["iwg_type"] == "VII(5C1&1)"


def test_class_A_ccr3_ccr5_reports_both_defined_types():
    """Class A + ccr3 + ccrC is III(3A&5) or XIV(5A&3); IWG interpretation left open."""
    hits = [hit("mecA", 20000, 22000), hit("mecR1", 22000, 23000), hit("mecI", 23000, 23300),
            hit("ccrA3", 40000, 41300), hit("ccrB3", 41300, 42800), hit("ccrC1", 5000, 6700)]
    r = classify_sccmec(hits)
    assert r["type_candidates"] == ["Type III", "Type XIV"]
    assert r["iwg_type"] == "III(3A&5) / XIV(5A&3)"


def test_mec_and_ccr_without_defined_type():
    hits = [hit("mecA", 20000, 22000), hit("mecR1", 22000, 23000), hit("IS1272", 23000, 25000),
            hit("ccrC2", 10000, 11700)]          # class B + ccr9: no defined type
    r = classify_sccmec(hits)
    assert r["iwg_type"] == "nt(9B)" and r["type_candidates"] == []


def test_assembly_limited_when_mec_at_contig_end():
    """Short-read assemblies: mecA on a ~4.7 kb contig ending ~230 bp after the gene."""
    r = classify_sccmec([hit("mecA", 2487, 4494, contig="c19", contig_len=4724)])
    assert r["assembly_limited"] is True
    assert r["status"] == "Partial (Assembly-limited)"
    assert any("contig end" in w for w in r["warnings"])


def test_not_assembly_limited_when_mec_mid_contig():
    r = classify_sccmec([hit("mecA", 500000, 502000, contig_len=2_800_000)])
    assert r["assembly_limited"] is False
    assert r["status"] == "Partial (Unclassifiable)"


def test_complete_cassette_never_flagged_even_near_edge():
    """A fully typed cassette is not flagged even if mecA is near a contig end."""
    hits = [h for h in C2] + [hit("ccrC1", 10000, 11700)]
    hits[0] = hit("mecA", 100, 2100, contig_len=40000)
    r = classify_sccmec(hits)
    assert r["assembly_limited"] is False and r["iwg_type"] == "V(5C2)"


def test_edge_threshold_constant():
    assert 500 <= CONTIG_EDGE_BP <= 5000
