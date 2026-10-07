"""Estimator/classifier agreement, locus de-duplication and cassette (DR) segmentation."""
import os
import random
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.classifier import classify_sccmec, collapse_loci
from lib.estimator import estimate_closest_types
from lib.cassette import dr_candidates, segment_loci, interpret, summarize


def hit(gene, start, end, strand="+", contig="c1", contig_len=2_800_000, id_pct=99.0, cov_pct=95.0):
    return {"gene": gene, "accession": "T", "scc_type": "T", "identity": id_pct / 100, "coverage": cov_pct / 100,
            "contig": contig, "contig_len": contig_len, "start": start, "end": end, "strand": strand,
            "id_pct": id_pct, "cov_pct": cov_pct, "len": end - start, "aln_len": end - start}


def test_estimator_respects_classifier_mec_class():
    """Flanking IS431 opposite (C2) + ccrC + ccr1: VII (C1) must not be a full match."""
    hits = [hit("ccrC1", 5000, 6700), hit("IS431_1", 9000, 9790, "-"), hit("mecA", 10000, 12000),
            hit("IS431_2", 14400, 15190, "+"), hit("ccrA1", 50000, 51300), hit("ccrB1", 51300, 52800)]
    r = classify_sccmec(hits)
    assert r["mec_complex"] == "Class C2"
    est = estimate_closest_types(r, hits, [])
    by = {c["type"]: c for c in est["candidates"]}
    assert est["best_guess"] in ("Type V", "Type IX")          # both C2 types
    if "Type VII" in by:
        assert by["Type VII"]["matching_components"]["mec_complex"]["status"] == "conflict"
        assert by["Type VII"]["status"] != "full_match"
        assert by["Type VII"]["score"] < by["Type V"]["score"]


def test_estimator_ccr_bonus_works_for_composites():
    hits = [hit("ccrC1", 5000, 6700), hit("IS431_1", 9000, 9790, "-"), hit("mecA", 10000, 12000),
            hit("IS431_2", 14400, 15190, "+"), hit("ccrA1", 50000, 51300), hit("ccrB1", 51300, 52800)]
    r = classify_sccmec(hits)
    assert "/" in r["ccr_complex"]                              # "Type 1 / Type 5"
    est = estimate_closest_types(r, hits, [])
    assert est["candidates"][0]["score"] > 0.9


def test_collapse_loci_keeps_one_best_hit_per_locus():
    hits = [hit("ccrC1", 1000, 2700, id_pct=92), hit("ccrC1", 1040, 2700, id_pct=100),
            hit("ccrC1", 30000, 31700), hit("mecA", 5000, 7000)]
    out = collapse_loci(hits)
    cc = [h for h in out if h["gene"] == "ccrC1"]
    assert len(cc) == 2 and cc[0]["id_pct"] == 100 and cc[0]["start"] == 1000
    assert len(out) == 3


def test_collapse_loci_leaves_read_mode_hits():
    reads = [{"gene": "mecA", "id_pct": 99, "cov_pct": 100}, {"gene": "mecA", "id_pct": 98, "cov_pct": 100}]
    assert len(collapse_loci(reads)) == 2


def _genome_with_drs():
    random.seed(1)
    rnd = lambda n: "".join(random.choice("ACGT") for _ in range(n))
    core = "GAGGCTTATCATAAATAA"
    orfx_end = 1000
    seq = rnd(orfx_end - len(core)) + core               # orfX ends with the core (attL)
    seq += rnd(20000) + core                              # DR at +20 kb (0 mm)
    seq += rnd(15000) + core[:-2] + "GG"                  # DR at +35 kb (2 mm)
    seq += rnd(10000) + "TTTT" + core[4:]                 # +45 kb (4 mm, candidate only)
    seq += rnd(5000)
    return seq, orfx_end


def test_dr_candidates_and_segments():
    seq, end3 = _genome_with_drs()
    c = dr_candidates(seq, end3, "+")
    offs = {o: m for o, m in c}
    assert offs[0] == 0
    assert any(abs(o - 20018) < 30 and m == 0 for o, m in c)
    assert any(m == 2 for o, m in c) and any(m >= 3 for o, m in c)
    loci = [("ccrC1", 5000), ("mecA", 9000), ("ccrC1", 15000), ("ccrA1", 25000), ("ccrB1", 26000)]
    segs = segment_loci(c, loci)
    filled = [s for s in segs if s["genes"]]
    assert filled[0]["genes"] == ["ccrC1", "mecA", "ccrC1"]
    assert filled[1]["genes"] == ["ccrA1", "ccrB1"]
    assert "tandem" in interpret(segs)
    assert summarize(segs).startswith("seg1(0.0-20.0kb):ccrC1,mecA,ccrC1|seg2(20.0-35.0kb):ccrA1,ccrB1")


def test_dr_candidates_minus_strand():
    seq, end3 = _genome_with_drs()
    comp = str.maketrans("ACGT", "TGCA")
    rc = seq.translate(comp)[::-1]
    c_plus = dr_candidates(seq, end3, "+")
    c_minus = dr_candidates(rc, len(seq) - end3, "-")
    assert c_plus == c_minus


def test_single_segment_no_tandem_note():
    segs = segment_loci([(0, 0), (29719, 2)], [("ccrC1", 6700), ("mecA", 13000), ("ccrC1", 23100)])
    assert interpret(segs) == ""


def test_class_C_unresolved_names_both_types():
    """Read mode (no coordinates): IS431 orientation unknown -> V and VII both named."""
    reads = [{"gene": g, "id_pct": 99.0, "cov_pct": 100.0} for g in ("mecA", "IS431", "ccrC1")]
    r = classify_sccmec(reads)
    assert r["mec_complex"] == "Class C"
    assert r["sccmec_type"] == "Type V / Type VII" and r["iwg_type"] == "V(5C2) / VII(5C1)"


def test_read_mode_prefers_class_forming_a_defined_type():
    """Reads: IS1272 elsewhere in the genome must not leave an undefined 'nt(5B)'."""
    reads = [{"gene": g, "id_pct": 99.0, "cov_pct": 100.0} for g in ("mecA", "IS431", "IS1272", "ccrC1")]
    r = classify_sccmec(reads)
    assert r["mec_complex"] == "Class C"
    assert r["iwg_type"] == "V(5C2) / VII(5C1)"
    assert any("Read mode" in w for w in r["warnings"])


def test_read_mode_keeps_class_B_when_it_forms_a_type():
    reads = [{"gene": g, "id_pct": 99.0, "cov_pct": 100.0} for g in ("mecA", "IS431", "IS1272", "ccrA2", "ccrB2")]
    r = classify_sccmec(reads)
    assert r["mec_complex"] == "Class B" and r["iwg_type"] == "IV(2B)"
