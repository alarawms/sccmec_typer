"""Tests for --best-fit promotion logic (apply_best_fit)."""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from bin.sccmec_typer import apply_best_fit, _DEFAULT_BEST_FIT_THRESHOLD


def _make_unknown_result(best_guess="Type V", score=0.91, mec_complex="Class C"):
    """Return a minimal result dict that looks like a split-assembly Unknown."""
    return {
        "status": "Positive",
        "sccmec_type": "Unknown",
        "mec_complex": mec_complex,
        "ccr_complex": "Type 5",
        "genes_detected": ["mecA", "ccrC1", "IS431"],
        "mecA_present": True,
        "warnings": ["Split Assembly: Components found on 2 different contigs"],
        "estimation": {
            "best_guess": best_guess,
            "best_guess_designation": "5C2",
            "best_guess_score": score,
            "best_guess_ruling": "IS431 orientation undetermined",
            "n_candidates": 3,
            "method": "weighted_component_matching",
            "sub_threshold_genes_used": [],
            "note": "Estimation based on partial evidence",
            "candidates": [],
        },
    }


class TestApplyBestFit:

    # ------------------------------------------------------------------
    # Happy path
    # ------------------------------------------------------------------

    def test_promotes_unknown_above_threshold(self):
        result = _make_unknown_result(score=0.91)
        result, promoted = apply_best_fit(result)
        assert promoted is True
        assert result["sccmec_type"] == "Type V (est.)"

    def test_best_fit_applied_block_populated(self):
        result = _make_unknown_result(score=0.91)
        result, promoted = apply_best_fit(result)
        bfi = result["best_fit_applied"]
        assert bfi["original"] == "Unknown"
        assert bfi["promoted_to"] == "Type V (est.)"
        assert bfi["score"] == 0.91
        assert bfi["threshold"] == _DEFAULT_BEST_FIT_THRESHOLD

    def test_warning_appended(self):
        result = _make_unknown_result(score=0.91)
        original_warnings = len(result["warnings"])
        result, _ = apply_best_fit(result)
        assert len(result["warnings"]) == original_warnings + 1
        assert "Best-fit estimation promoted" in result["warnings"][-1]

    def test_custom_threshold_respected(self):
        """Score 0.80 promoted at threshold=0.75, not at threshold=0.85."""
        r1 = _make_unknown_result(score=0.80)
        r1, p1 = apply_best_fit(r1, threshold=0.75)
        assert p1 is True

        r2 = _make_unknown_result(score=0.80)
        r2, p2 = apply_best_fit(r2, threshold=0.85)
        assert p2 is False
        assert r2["sccmec_type"] == "Unknown"

    # ------------------------------------------------------------------
    # Boundary / threshold
    # ------------------------------------------------------------------

    def test_exactly_at_threshold_promotes(self):
        result = _make_unknown_result(score=_DEFAULT_BEST_FIT_THRESHOLD)
        result, promoted = apply_best_fit(result)
        assert promoted is True

    def test_just_below_threshold_does_not_promote(self):
        result = _make_unknown_result(score=_DEFAULT_BEST_FIT_THRESHOLD - 0.01)
        result, promoted = apply_best_fit(result)
        assert promoted is False
        assert result["sccmec_type"] == "Unknown"

    def test_very_low_score_does_not_promote(self):
        result = _make_unknown_result(score=0.37)
        result, promoted = apply_best_fit(result)
        assert promoted is False

    # ------------------------------------------------------------------
    # No-op conditions
    # ------------------------------------------------------------------

    def test_already_typed_not_touched(self):
        result = _make_unknown_result()
        result["sccmec_type"] = "Type V"  # already resolved
        result, promoted = apply_best_fit(result)
        assert promoted is False
        assert result["sccmec_type"] == "Type V"

    def test_negative_not_touched(self):
        result = {
            "status": "Negative",
            "sccmec_type": "Negative",
            "estimation": None,
        }
        result, promoted = apply_best_fit(result)
        assert promoted is False
        assert result["sccmec_type"] == "Negative"

    def test_no_estimation_data_not_touched(self):
        result = _make_unknown_result()
        result["estimation"] = None
        result, promoted = apply_best_fit(result)
        assert promoted is False

    def test_empty_best_guess_not_touched(self):
        result = _make_unknown_result()
        result["estimation"]["best_guess"] = ""
        result, promoted = apply_best_fit(result)
        assert promoted is False

    # ------------------------------------------------------------------
    # Idempotency
    # ------------------------------------------------------------------

    def test_idempotent_when_called_twice(self):
        result = _make_unknown_result(score=0.91)
        result, _ = apply_best_fit(result)
        # Second call: sccmec_type is now "Type V (est.)", not "Unknown"
        result, promoted2 = apply_best_fit(result)
        assert promoted2 is False
        assert result["sccmec_type"] == "Type V (est.)"
