import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.estimator import estimate_closest_types


def _make_hit(gene, id_pct=99.0, cov_pct=95.0, contig="contig_1", strand="+"):
    return {
        "gene": gene, "accession": "TEST", "scc_type": "TEST",
        "identity": id_pct / 100, "coverage": cov_pct / 100,
        "contig": contig, "start": 0, "end": 1000, "strand": strand,
        "id_pct": id_pct, "cov_pct": cov_pct, "len": 1000, "aln_len": 950,
    }


class TestSubThresholdEstimation:
    """Scenario: genes detected but below threshold — soft_hits provide evidence."""

    def test_sub_threshold_ccrB2_estimates_type_IV(self):
        """mecA + IS1272 above threshold, ccrA2 above, ccrB2 below → estimate Type IV."""
        hits = [_make_hit("mecA"), _make_hit("IS1272"), _make_hit("IS431"),
                _make_hit("ccrA2")]
        soft_hits = [_make_hit("ccrB2", id_pct=92.0, cov_pct=78.0)]
        result = {
            "status": "Partial (Unclassifiable)",
            "sccmec_type": "Unknown",
            "mec_complex": "Class B",
            "ccr_complex": "Negative",
            "genes_detected": ["mecA", "IS1272", "IS431", "ccrA2"],
        }
        estimation = estimate_closest_types(result, hits, soft_hits)
        assert estimation is not None
        assert estimation["best_guess"] == "Type IV"
        assert estimation["best_guess_score"] > 0.5
        assert "ccrB2" in estimation["sub_threshold_genes_used"]
        assert any("ccrB2" in c["ruling"] for c in estimation["candidates"])

    def test_sub_threshold_returns_top_3(self):
        hits = [_make_hit("mecA"), _make_hit("IS1272"), _make_hit("IS431")]
        soft_hits = [_make_hit("ccrA2", id_pct=88.0, cov_pct=75.0)]
        result = {
            "status": "Partial (Unclassifiable)",
            "sccmec_type": "Unknown",
            "mec_complex": "Class B",
            "ccr_complex": "Negative",
            "genes_detected": ["mecA", "IS1272", "IS431"],
        }
        estimation = estimate_closest_types(result, hits, soft_hits)
        assert estimation["n_candidates"] <= 3
        assert len(estimation["candidates"]) <= 3


class TestMissingGeneEstimation:
    """Scenario: gene completely absent."""

    def test_missing_ccr_estimates_closest(self):
        """Class B mec but no ccr at all → should still rank types that use Class B."""
        hits = [_make_hit("mecA"), _make_hit("IS1272"), _make_hit("IS431")]
        soft_hits = []
        result = {
            "status": "Partial (Unclassifiable)",
            "sccmec_type": "Unknown",
            "mec_complex": "Class B",
            "ccr_complex": "Negative",
            "genes_detected": ["mecA", "IS1272", "IS431"],
        }
        estimation = estimate_closest_types(result, hits, soft_hits)
        assert estimation is not None
        # All Class B types (I, IV, VI) should be candidates
        candidate_types = [c["type"] for c in estimation["candidates"]]
        assert any(t in candidate_types for t in ["Type I", "Type IV", "Type VI"])
        # Scores should be < 1.0 since ccr is missing
        assert all(c["score"] < 1.0 for c in estimation["candidates"])


class TestOrphanCcrEstimation:
    """Scenario: ccr found but no mec gene."""

    def test_orphan_ccr_estimates_type(self):
        hits = [_make_hit("ccrA2"), _make_hit("ccrB2")]
        soft_hits = []
        result = {
            "status": "Partial (Orphan ccr)",
            "sccmec_type": "Unknown",
            "mec_complex": "Negative",
            "ccr_complex": "Type 2",
            "genes_detected": ["ccrA2", "ccrB2"],
        }
        estimation = estimate_closest_types(result, hits, soft_hits)
        assert estimation is not None
        # Types with ccr Type 2: II, IV
        candidate_types = [c["type"] for c in estimation["candidates"]]
        assert any(t in candidate_types for t in ["Type II", "Type IV"])
        assert "ruling" in estimation["candidates"][0]


class TestWrongCombinationEstimation:
    """Scenario: mec + ccr found but combination doesn't match any approved type."""

    def test_novel_combination(self):
        """Class B + ccr Type 5 has no approved type → estimate closest."""
        hits = [_make_hit("mecA"), _make_hit("IS1272"), _make_hit("IS431"),
                _make_hit("ccrC1")]
        soft_hits = []
        result = {
            "status": "Positive",
            "sccmec_type": "Unknown",
            "mec_complex": "Class B",
            "ccr_complex": "Type 5",
            "genes_detected": ["mecA", "IS1272", "IS431", "ccrC1"],
        }
        estimation = estimate_closest_types(result, hits, soft_hits)
        assert estimation is not None
        assert estimation["n_candidates"] >= 1
        assert estimation["candidates"][0]["score"] > 0


class TestPositiveNoEstimation:
    """Positive results should not trigger estimation."""

    def test_positive_returns_none(self):
        hits = [_make_hit("mecA"), _make_hit("IS1272"), _make_hit("IS431"),
                _make_hit("ccrA1"), _make_hit("ccrB1")]
        soft_hits = []
        result = {
            "status": "Positive",
            "sccmec_type": "Type I",
            "mec_complex": "Class B",
            "ccr_complex": "Type 1",
            "genes_detected": ["mecA", "IS1272", "IS431", "ccrA1", "ccrB1"],
        }
        estimation = estimate_closest_types(result, hits, soft_hits)
        assert estimation is None


class TestNegativeNoEstimation:
    """Negative results should not trigger estimation."""

    def test_negative_returns_none(self):
        result = {
            "status": "Negative",
            "sccmec_type": "Negative",
            "mec_complex": "Negative",
            "ccr_complex": "Negative",
            "genes_detected": [],
        }
        estimation = estimate_closest_types(result, [], [])
        assert estimation is None
