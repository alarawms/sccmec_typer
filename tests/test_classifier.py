import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lib.classifier import classify_sccmec


class TestMecComplexClassification:
    """Verify correct mec complex class assignment."""

    def test_class_A(self, type_II_hits):
        result = classify_sccmec(type_II_hits)
        assert result["mec_complex"] == "Class A"

    def test_class_B(self, type_I_hits):
        result = classify_sccmec(type_I_hits)
        assert result["mec_complex"] == "Class B"

    def test_class_C1_same_orientation(self, type_VII_hits):
        result = classify_sccmec(type_VII_hits)
        assert result["mec_complex"] == "Class C1"

    def test_class_C2_opposite_orientation(self, type_V_hits):
        result = classify_sccmec(type_V_hits)
        assert result["mec_complex"] == "Class C2"

    def test_class_E(self, type_XI_hits):
        result = classify_sccmec(type_XI_hits)
        assert result["mec_complex"] == "Class E"

    def test_plasmid_mecB(self, plasmid_mecB_hits):
        result = classify_sccmec(plasmid_mecB_hits)
        assert result["mec_complex"] == "Plasmid-borne (mecB)"


class TestCcrComplexClassification:
    """Verify correct ccr complex type assignment via gene pairs."""

    def test_type_1_ccrA1B1(self, type_I_hits):
        result = classify_sccmec(type_I_hits)
        assert "Type 1" in result["ccr_complex"]

    def test_type_2_ccrA2B2(self, type_II_hits):
        result = classify_sccmec(type_II_hits)
        assert "Type 2" in result["ccr_complex"]

    def test_type_3_ccrA3B3(self, type_III_hits):
        result = classify_sccmec(type_III_hits)
        assert "Type 3" in result["ccr_complex"]

    def test_type_4_ccrA4B4(self, type_VI_hits):
        result = classify_sccmec(type_VI_hits)
        assert "Type 4" in result["ccr_complex"]

    def test_type_5_ccrC1(self, type_V_hits):
        result = classify_sccmec(type_V_hits)
        assert "Type 5" in result["ccr_complex"]

    def test_type_7_ccrA1B6(self, type_X_hits):
        """ccrA1 + ccrB6 = Type 7, NOT Type 1."""
        result = classify_sccmec(type_X_hits)
        assert "Type 7" in result["ccr_complex"]
        assert "Type 1" not in result["ccr_complex"]

    def test_type_8_ccrA1B3(self, type_XI_hits):
        """ccrA1 + ccrB3 = Type 8, NOT Type 1 or Type 3."""
        result = classify_sccmec(type_XI_hits)
        assert "Type 8" in result["ccr_complex"]
        assert "Type 1" not in result["ccr_complex"]
        assert "Type 3" not in result["ccr_complex"]

    def test_type_9_ccrC2(self, type_XII_hits):
        result = classify_sccmec(type_XII_hits)
        assert "Type 9" in result["ccr_complex"]


class TestSCCmecTypeAssignment:
    """Verify correct SCCmec type from mec+ccr combination."""

    def test_type_I(self, type_I_hits):
        result = classify_sccmec(type_I_hits)
        assert result["sccmec_type"] == "Type I"
        assert result["status"] == "Positive"

    def test_type_II(self, type_II_hits):
        result = classify_sccmec(type_II_hits)
        assert result["sccmec_type"] == "Type II"

    def test_type_III(self, type_III_hits):
        result = classify_sccmec(type_III_hits)
        assert result["sccmec_type"] == "Type III"

    def test_type_IV(self, type_IV_hits):
        result = classify_sccmec(type_IV_hits)
        assert result["sccmec_type"] == "Type IV"

    def test_type_V(self, type_V_hits):
        result = classify_sccmec(type_V_hits)
        assert result["sccmec_type"] == "Type V"

    def test_type_VI(self, type_VI_hits):
        result = classify_sccmec(type_VI_hits)
        assert result["sccmec_type"] == "Type VI"

    def test_type_VII(self, type_VII_hits):
        result = classify_sccmec(type_VII_hits)
        assert result["sccmec_type"] == "Type VII"

    def test_type_VIII(self, type_VIII_hits):
        result = classify_sccmec(type_VIII_hits)
        assert result["sccmec_type"] == "Type VIII"

    def test_type_IX(self, type_IX_hits):
        result = classify_sccmec(type_IX_hits)
        assert result["sccmec_type"] == "Type IX"

    def test_type_X(self, type_X_hits):
        result = classify_sccmec(type_X_hits)
        assert result["sccmec_type"] == "Type X"

    def test_type_XI(self, type_XI_hits):
        result = classify_sccmec(type_XI_hits)
        assert result["sccmec_type"] == "Type XI"

    def test_type_XII(self, type_XII_hits):
        result = classify_sccmec(type_XII_hits)
        assert result["sccmec_type"] == "Type XII"

    def test_type_XIII(self, type_XIII_hits):
        result = classify_sccmec(type_XIII_hits)
        assert result["sccmec_type"] == "Type XIII"

    def test_type_XIV(self, type_XIV_hits):
        result = classify_sccmec(type_XIV_hits)
        assert result["sccmec_type"] == "Type XIV"

    def test_type_XV(self, type_XV_hits):
        result = classify_sccmec(type_XV_hits)
        assert result["sccmec_type"] == "Type XV"

    def test_plasmid_mecB(self, plasmid_mecB_hits):
        result = classify_sccmec(plasmid_mecB_hits)
        assert result["sccmec_type"] == "n/a (Plasmid)"


class TestEdgeCases:
    """Verify handling of negative, partial, and composite results."""

    def test_negative(self, negative_hits):
        result = classify_sccmec(negative_hits)
        assert result["status"] == "Negative"

    def test_orphan_ccr(self, orphan_ccr_hits):
        result = classify_sccmec(orphan_ccr_hits)
        assert result["status"] == "Partial (Orphan ccr)"

    def test_composite_multiple_ccr(self, composite_hits):
        result = classify_sccmec(composite_hits)
        assert "Composite" in result["sccmec_type"]
        assert any("Multiple ccr" in w for w in result["warnings"])

    def test_mecA_present_flag(self, type_I_hits):
        result = classify_sccmec(type_I_hits)
        assert result["mecA_present"] is True

    def test_mecA_absent_for_type_XI(self, type_XI_hits):
        """Type XI uses mecC, not mecA."""
        result = classify_sccmec(type_XI_hits)
        assert result["mecA_present"] is False

    def test_split_assembly_warning(self):
        """Components on different contigs should warn."""
        from tests.conftest import make_hit
        hits = [
            make_hit("mecA", contig="contig_1"),
            make_hit("mecR1", contig="contig_1"),
            make_hit("IS1272", contig="contig_1"),
            make_hit("IS431", contig="contig_1"),
            make_hit("ccrA1", contig="contig_2"),
            make_hit("ccrB1", contig="contig_2"),
        ]
        result = classify_sccmec(hits)
        assert any("Split" in w for w in result["warnings"])
