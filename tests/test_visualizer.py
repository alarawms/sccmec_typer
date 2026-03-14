import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.visualizer import generate_cassette_svg, generate_report_html


class TestCassetteSVG:
    """Verify SVG generation for schematic diagram."""

    def test_positive_result_produces_svg(self, tmp_path):
        result = _make_positive_result()
        svg = generate_cassette_svg(result)
        assert svg.startswith("<?xml")
        assert "SCCmec Type II" in svg
        assert "mecA" in svg
        assert "ccrA2" in svg

    def test_negative_result_produces_svg(self, tmp_path):
        result = _make_negative_result()
        svg = generate_cassette_svg(result)
        assert "No SCCmec" in svg

    def test_split_assembly_shows_badges(self):
        result = _make_split_result()
        svg = generate_cassette_svg(result)
        assert "[c1]" in svg or "[c2]" in svg


class TestHTMLReport:
    """Verify HTML report generation."""

    def test_produces_html(self, tmp_path):
        result = _make_positive_result()
        html = generate_report_html(result, sample_name="test_sample")
        assert "<!DOCTYPE html>" in html
        assert "test_sample" in html
        assert "Gene Details" in html

    def test_negative_result_report(self):
        result = _make_negative_result()
        html = generate_report_html(result, sample_name="neg_sample")
        assert "No SCCmec Detected" in html


def _make_positive_result():
    """Helper: enriched positive Type II result."""
    return {
        "status": "Positive",
        "sccmec_type": "Type II",
        "designation": "2A",
        "mec_complex": "Class A",
        "ccr_complex": "Type 2",
        "ccr_genes": ["ccrA2", "ccrB2"],
        "genes_detected": ["IS431", "ccrA2", "ccrB2", "mecA", "mecI", "mecR1"],
        "mecA_present": True,
        "mecC_present": False,
        "mecB_present": False,
        "confidence": {
            "mode": "full",
            "type_level": {"score": 0.92, "tier": "High"},
            "mec_component": {
                "class": "Class A", "score": 0.99, "tier": "High",
                "genes_expected": ["mecA", "mecR1", "mecI"],
                "genes_detected": ["mecA", "mecR1", "mecI"],
                "genes_missing": [],
                "orientation": None, "orientation_resolved": None,
            },
            "ccr_component": {
                "type": "Type 2", "score": 0.93, "tier": "High",
                "genes_expected": ["ccrA2", "ccrB2"],
                "genes_detected": ["ccrA2", "ccrB2"],
                "genes_missing": [],
            },
            "per_gene": {
                "mecA": {"identity_pct": 100, "coverage_pct": 100, "score": 1.0, "tier": "High"},
                "mecR1": {"identity_pct": 99, "coverage_pct": 99, "score": 0.97, "tier": "High"},
                "mecI": {"identity_pct": 99, "coverage_pct": 100, "score": 0.99, "tier": "High"},
                "ccrA2": {"identity_pct": 99, "coverage_pct": 98, "score": 0.94, "tier": "High"},
                "ccrB2": {"identity_pct": 99, "coverage_pct": 97, "score": 0.92, "tier": "High"},
                "IS431": {"identity_pct": 98, "coverage_pct": 95, "score": 0.78, "tier": "High"},
            },
        },
        "assembly": {
            "contigs": ["contig_1"], "is_split": False, "split_penalty_applied": False,
            "gene_locations": [
                {"gene": "mecA", "contig": "contig_1", "start": 28100, "end": 30100, "strand": "+"},
                {"gene": "mecR1", "contig": "contig_1", "start": 30200, "end": 32000, "strand": "+"},
                {"gene": "mecI", "contig": "contig_1", "start": 32100, "end": 33000, "strand": "+"},
                {"gene": "IS431", "contig": "contig_1", "start": 27000, "end": 28000, "strand": "+"},
                {"gene": "ccrA2", "contig": "contig_1", "start": 12000, "end": 14000, "strand": "-"},
                {"gene": "ccrB2", "contig": "contig_1", "start": 10000, "end": 12000, "strand": "-"},
            ],
        },
        "warnings": [],
        "hits_summary": [],
    }


def _make_negative_result():
    return {
        "status": "Negative", "reason": "No alignments found",
        "sccmec_type": "Negative", "designation": None,
        "mec_complex": "Negative", "ccr_complex": "Negative",
        "ccr_genes": [], "genes_detected": [],
        "mecA_present": False, "mecC_present": False, "mecB_present": False,
        "confidence": {"mode": "n/a", "type_level": {"score": 0.0, "tier": "n/a"},
                       "mec_component": None, "ccr_component": None, "per_gene": {}},
        "assembly": {"contigs": [], "is_split": False, "split_penalty_applied": False,
                     "gene_locations": []},
        "warnings": [], "hits_summary": [],
    }


def _make_split_result():
    result = _make_positive_result()
    result["assembly"]["is_split"] = True
    result["assembly"]["contigs"] = ["contig_1", "contig_2"]
    result["assembly"]["gene_locations"][4]["contig"] = "contig_2"
    result["assembly"]["gene_locations"][5]["contig"] = "contig_2"
    result["warnings"] = ["Split Assembly: Components found on 2 different contigs"]
    return result
