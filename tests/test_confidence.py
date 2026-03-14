import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.confidence import (
    compute_gene_confidence,
    compute_component_confidence,
    compute_type_confidence,
    enrich_result_with_confidence,
)


class TestGeneConfidence:
    """Level 1: Per-gene confidence scoring."""

    def test_perfect_hit(self):
        score, tier = compute_gene_confidence(100.0, 100.0)
        assert score == 1.0
        assert tier == "High"

    def test_minimum_threshold(self):
        score, tier = compute_gene_confidence(90.0, 80.0)
        assert score == 0.0
        assert tier == "Low"

    def test_mid_range(self):
        score, tier = compute_gene_confidence(95.0, 90.0)
        # norm_id = (95-90)/10 = 0.5, norm_cov = (90-80)/20 = 0.5
        # score = 0.6*0.5 + 0.4*0.5 = 0.5
        assert abs(score - 0.5) < 0.001
        assert tier == "Medium"

    def test_high_threshold(self):
        score, tier = compute_gene_confidence(98.0, 96.0)
        # norm_id = 0.8, norm_cov = 0.8
        # score = 0.6*0.8 + 0.4*0.8 = 0.8
        assert abs(score - 0.8) < 0.001
        assert tier == "High"

    def test_clamp_above_100(self):
        """PAF can produce identity > 100% in edge cases."""
        score, tier = compute_gene_confidence(102.0, 105.0)
        assert score == 1.0

    def test_reads_mode_fallback(self):
        """When id_pct is 0 (reads mode), use coverage-only."""
        score, tier = compute_gene_confidence(0.0, 95.0)
        # coverage-only: (95-80)/20 = 0.75
        assert abs(score - 0.75) < 0.001
        assert tier == "High"

    def test_reads_mode_low(self):
        score, tier = compute_gene_confidence(0.0, 82.0)
        # (82-80)/20 = 0.1
        assert abs(score - 0.1) < 0.001
        assert tier == "Low"


class TestComponentConfidence:
    """Level 2: Component-level confidence."""

    def test_all_genes_present_high(self):
        gene_scores = {"mecA": 1.0, "mecR1": 0.98, "mecI": 0.99}
        expected_genes = ["mecA", "mecR1", "mecI"]
        score, tier, missing = compute_component_confidence(
            gene_scores, expected_genes
        )
        assert score > 0.9
        assert tier == "High"
        assert missing == []

    def test_missing_gene(self):
        gene_scores = {"mecA": 1.0, "mecR1": 0.98}
        expected_genes = ["mecA", "mecR1", "mecI"]
        score, tier, missing = compute_component_confidence(
            gene_scores, expected_genes
        )
        # (2/3) * mean(1.0, 0.98) = 0.667 * 0.99 = 0.660
        assert score < 0.7
        assert "mecI" in missing

    def test_orientation_penalty(self):
        gene_scores = {"mecA": 1.0}
        expected_genes = ["mecA"]
        score_no_penalty, _, _ = compute_component_confidence(
            gene_scores, expected_genes
        )
        score_penalty, _, _ = compute_component_confidence(
            gene_scores, expected_genes, orientation_undetermined=True,
            mec_class="Class C"
        )
        assert abs(score_penalty - score_no_penalty * 0.8) < 0.001

    def test_orientation_penalty_not_applied_to_class_E(self):
        gene_scores = {"mecC": 1.0}
        expected_genes = ["mecC"]
        score, _, _ = compute_component_confidence(
            gene_scores, expected_genes, orientation_undetermined=True,
            mec_class="Class E"
        )
        # No penalty for Class E
        assert score == 1.0

    def test_empty_expected(self):
        """ccr=ANY case: no genes expected."""
        score, tier, missing = compute_component_confidence({}, [])
        assert score == 1.0
        assert tier == "High"


class TestTypeConfidence:
    """Level 3: Type-level confidence."""

    def test_high_confidence(self):
        score, tier = compute_type_confidence(0.99, 0.93, is_split=False)
        assert score > 0.7
        assert tier == "High"

    def test_split_penalty(self):
        score_normal, _ = compute_type_confidence(0.99, 0.93, is_split=False)
        score_split, _ = compute_type_confidence(0.99, 0.93, is_split=True)
        assert abs(score_split - score_normal * 0.85) < 0.001

    def test_low_confidence(self):
        score, tier = compute_type_confidence(0.2, 0.3, is_split=True)
        assert score < 0.4
        assert tier == "Low"


class TestEnrichResult:
    """Integration: enrich_result_with_confidence adds all fields."""

    def test_positive_result_has_confidence(self):
        result = {
            "status": "Positive",
            "sccmec_type": "Type II",
            "mec_complex": "Class A",
            "ccr_complex": "Type 2",
            "genes_detected": ["mecA", "mecR1", "mecI", "ccrA2", "ccrB2"],
            "mecA_present": True,
            "warnings": [],
            "hits_summary": [
                {"gene": "mecA", "id_pct": 100, "cov_pct": 100, "contig": "c1",
                 "start": 0, "end": 1000, "strand": "+"},
                {"gene": "mecR1", "id_pct": 99, "cov_pct": 99, "contig": "c1",
                 "start": 1000, "end": 2000, "strand": "+"},
                {"gene": "mecI", "id_pct": 99, "cov_pct": 100, "contig": "c1",
                 "start": 2000, "end": 3000, "strand": "+"},
                {"gene": "ccrA2", "id_pct": 99, "cov_pct": 98, "contig": "c1",
                 "start": 5000, "end": 6000, "strand": "-"},
                {"gene": "ccrB2", "id_pct": 99, "cov_pct": 97, "contig": "c1",
                 "start": 6000, "end": 7000, "strand": "-"},
            ],
        }
        enriched = enrich_result_with_confidence(result)
        assert "confidence" in enriched
        assert "assembly" in enriched
        assert enriched["confidence"]["mode"] == "full"
        assert enriched["confidence"]["type_level"]["tier"] == "High"
        assert "mecA" in enriched["confidence"]["per_gene"]
        assert enriched["assembly"]["is_split"] is False

    def test_negative_result_has_confidence(self):
        result = {
            "status": "Negative",
            "reason": "No alignments found",
        }
        enriched = enrich_result_with_confidence(result)
        assert enriched["confidence"]["mode"] == "n/a"
        assert enriched["confidence"]["type_level"]["score"] == 0.0
        assert enriched["assembly"]["gene_locations"] == []
