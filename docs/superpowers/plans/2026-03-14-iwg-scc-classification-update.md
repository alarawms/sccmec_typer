# IWG-SCC Classification Update: SCCmec Types I-XV

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Update the SCCmec typer to support all 15 IWG-SCC approved types (I-XV) with scientifically accurate mec complex classes (A-E), pair-based ccr complex typing (types 1-9), and proper IS431 orientation detection for C1/C2 distinction.

**Architecture:** The rules.json schema is redesigned with explicit ccr gene pair rules replacing substring matching, proper IWG-SCC mec class nomenclature (Class E instead of "Class C2"), and C1/C2 mec class distinction via IS431 orientation analysis from PAF strand information. The classifier algorithm changes from substring-based ccr matching to pair-based matching with fallback logic.

**Tech Stack:** Python 3.9+, minimap2, pytest

---

## File Structure

| File | Action | Responsibility |
|------|--------|---------------|
| `db/rules.json` | **Rewrite** | Classification rules: mec complex classes A-E, ccr complex types 1-9, SCCmec types I-XV |
| `lib/classifier.py` | **Rewrite** | Classification algorithm: pair-based ccr matching, IS431 orientation analysis, mec class E support |
| `tests/test_classifier.py` | **Create** | Unit tests for classifier with all 15 type combinations |
| `tests/conftest.py` | **Create** | Shared test fixtures (mock hit lists for each reference strain) |

---

## Chunk 1: Rules Database Redesign

### Task 1: Rewrite `rules.json` with IWG-SCC compliant schema

**Files:**
- Modify: `db/rules.json`

The current rules.json has three problems:
1. ccr matching uses substring ("1" matches any gene with "1" in the name — ambiguous for ccrA1 which appears in types 1, 7, 8, 9, X, XI, XV)
2. mec Class E is misnamed "Class C2" and mec Classes C1/C2 are merged as "Class C"
3. Only 6 of 15 types are defined

The new schema uses explicit ccr gene pair matching and adds all IWG-SCC approved classes/types.

- [ ] **Step 1: Rewrite rules.json**

Replace the entire file with:

```json
{
    "mec_complex_rules": [
        {
            "name": "Class A",
            "required": ["mecA", "mecR1", "mecI"],
            "description": "IS431-mecA-mecR1(complete)-mecI",
            "note": "Full regulatory region intact"
        },
        {
            "name": "Class B",
            "required": ["mecA", "IS1272"],
            "description": "IS1272-mecA-mecR1(truncated)-IS431",
            "note": "IS1272 insertion truncates mecR1"
        },
        {
            "name": "Class C1",
            "required": ["mecA"],
            "any_of": ["IS431", "IS431_1", "IS431_2"],
            "excluded": ["IS1272", "mecI"],
            "orientation": "same",
            "description": "IS431-mecA-mecR1(truncated)-IS431 (same orientation)",
            "note": "Both IS431 copies in same orientation"
        },
        {
            "name": "Class C2",
            "required": ["mecA"],
            "any_of": ["IS431", "IS431_1", "IS431_2"],
            "excluded": ["IS1272", "mecI"],
            "orientation": "opposite",
            "description": "IS431-mecA-mecR1(truncated)-IS431 (opposite orientation)",
            "note": "IS431 copies in opposite orientations"
        },
        {
            "name": "Class C",
            "required": ["mecA"],
            "any_of": ["IS431", "IS431_1", "IS431_2"],
            "excluded": ["IS1272", "mecI"],
            "description": "IS431-mecA-mecR1(truncated)-IS431 (orientation undetermined)",
            "note": "Fallback when IS431 orientation cannot be resolved"
        },
        {
            "name": "Class D",
            "required": ["mecA"],
            "excluded": ["IS1272", "mecI", "IS431", "IS431_1", "IS431_2"],
            "requires_mecR1_truncated": true,
            "description": "mecA-mecR1(truncated), no flanking IS elements",
            "note": "mecR1 present but truncated, no IS431 or IS1272"
        },
        {
            "name": "Class E",
            "required": ["mecC"],
            "description": "blaZ-mecC-mecR1-mecI",
            "note": "mecC homolog with blaZ; formerly called 'mecC-associated'"
        },
        {
            "name": "Plasmid-borne (mecB)",
            "required": ["mecB"],
            "description": "Plasmid-encoded mecB",
            "note": "Not chromosomally integrated SCCmec"
        },
        {
            "name": "Unclassifiable (mecA)",
            "required": ["mecA"],
            "description": "mecA detected but mec complex class cannot be determined",
            "note": "Fallback when no IS or regulatory markers found"
        }
    ],
    "ccr_complex_rules": [
        {
            "name": "Type 1",
            "required_genes": ["ccrA1", "ccrB1"],
            "description": "ccrA1B1"
        },
        {
            "name": "Type 2",
            "required_genes": ["ccrA2", "ccrB2"],
            "description": "ccrA2B2"
        },
        {
            "name": "Type 3",
            "required_genes": ["ccrA3", "ccrB3"],
            "description": "ccrA3B3"
        },
        {
            "name": "Type 4",
            "required_genes": ["ccrA4", "ccrB4"],
            "description": "ccrA4B4"
        },
        {
            "name": "Type 5",
            "required_genes": ["ccrC1"],
            "description": "ccrC1 (single recombinase)"
        },
        {
            "name": "Type 6",
            "required_genes": ["ccrA4", "ccrB6"],
            "description": "ccrA4B6",
            "note": "Rare; reported in some CoNS"
        },
        {
            "name": "Type 7",
            "required_genes": ["ccrA1", "ccrB6"],
            "description": "ccrA1B6"
        },
        {
            "name": "Type 8",
            "required_genes": ["ccrA1", "ccrB3"],
            "description": "ccrA1B3"
        },
        {
            "name": "Type 9",
            "required_genes": ["ccrC2"],
            "description": "ccrC2 (single recombinase)"
        }
    ],
    "sccmec_type_rules": [
        {
            "name": "Type I",
            "designation": "1B",
            "mec": "Class B",
            "ccr": "Type 1",
            "reference_strain": "NCTC10442"
        },
        {
            "name": "Type II",
            "designation": "2A",
            "mec": "Class A",
            "ccr": "Type 2",
            "reference_strain": "N315"
        },
        {
            "name": "Type III",
            "designation": "3A",
            "mec": "Class A",
            "ccr": "Type 3",
            "reference_strain": "85/2082"
        },
        {
            "name": "Type IV",
            "designation": "2B",
            "mec": "Class B",
            "ccr": "Type 2",
            "reference_strain": "CA05"
        },
        {
            "name": "Type V",
            "designation": "5C2",
            "mec": "Class C2",
            "ccr": "Type 5",
            "reference_strain": "WIS"
        },
        {
            "name": "Type VI",
            "designation": "4B",
            "mec": "Class B",
            "ccr": "Type 4",
            "reference_strain": "HDE288"
        },
        {
            "name": "Type VII",
            "designation": "5C1",
            "mec": "Class C1",
            "ccr": "Type 5",
            "reference_strain": "JCSC6082"
        },
        {
            "name": "Type VIII",
            "designation": "4A",
            "mec": "Class A",
            "ccr": "Type 4",
            "reference_strain": "C10682"
        },
        {
            "name": "Type IX",
            "designation": "1C2",
            "mec": "Class C2",
            "ccr": "Type 1",
            "reference_strain": "JCSC6943"
        },
        {
            "name": "Type X",
            "designation": "7C1",
            "mec": "Class C1",
            "ccr": "Type 7",
            "reference_strain": "JCSC6945"
        },
        {
            "name": "Type XI",
            "designation": "8E",
            "mec": "Class E",
            "ccr": "Type 8",
            "reference_strain": "LGA251"
        },
        {
            "name": "Type XII",
            "designation": "9C2",
            "mec": "Class C2",
            "ccr": "Type 9",
            "reference_strain": "BA01611"
        },
        {
            "name": "Type XIII",
            "designation": "9A",
            "mec": "Class A",
            "ccr": "Type 9",
            "reference_strain": "55-99-44"
        },
        {
            "name": "Type XIV",
            "designation": "5A",
            "mec": "Class A",
            "ccr": "Type 5",
            "reference_strain": "SC792"
        },
        {
            "name": "Type XV",
            "designation": "7A",
            "mec": "Class A",
            "ccr": "Type 7",
            "reference_strain": "NV_1"
        },
        {
            "name": "n/a (Plasmid)",
            "designation": "n/a",
            "mec": "Plasmid-borne (mecB)",
            "ccr": "ANY",
            "reference_strain": "n/a"
        }
    ]
}
```

- [ ] **Step 2: Validate JSON syntax**

Run: `python3 -c "import json; json.load(open('db/rules.json')); print('Valid JSON')"`
Expected: `Valid JSON`

- [ ] **Step 3: Commit**

```bash
git add db/rules.json
git commit -m "feat: rewrite rules.json with IWG-SCC types I-XV, pair-based ccr, classes A-E"
```

---

## Chunk 2: Test Fixtures and Unit Tests

### Task 2: Create test infrastructure

**Files:**
- Create: `tests/__init__.py`
- Create: `tests/conftest.py`
- Create: `tests/test_classifier.py`

These tests use mock hit lists (no minimap2 needed) to verify that the classifier correctly assigns SCCmec types from gene presence/absence data.

- [ ] **Step 1: Create tests/__init__.py**

Empty file to make tests a package.

- [ ] **Step 2: Create tests/conftest.py with fixtures for all 15 types**

Each fixture returns a list of hit dicts mimicking what parser.py or coverage.py would produce for a reference strain. The key fields the classifier uses are: `gene`, `strand`, `accession`, `scc_type`.

```python
import pytest


def make_hit(gene, accession="TEST", scc_type="TEST", strand="+", contig="contig_1"):
    """Helper to create a minimal hit dict for classifier testing."""
    return {
        "gene": gene,
        "accession": accession,
        "scc_type": scc_type,
        "identity": 0.99,
        "coverage": 0.95,
        "contig": contig,
        "start": 0,
        "end": 1000,
        "strand": strand,
        "id_pct": 99.0,
        "cov_pct": 95.0,
        "len": 1000,
        "aln_len": 950,
    }


@pytest.fixture
def type_I_hits():
    """Type I (1B): Class B + ccr Type 1 (ccrA1B1). Ref: NCTC10442."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("IS1272"),
        make_hit("IS431"), make_hit("ccrA1"), make_hit("ccrB1"),
    ]


@pytest.fixture
def type_II_hits():
    """Type II (2A): Class A + ccr Type 2 (ccrA2B2). Ref: N315."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("IS431"), make_hit("ccrA2"), make_hit("ccrB2"),
    ]


@pytest.fixture
def type_III_hits():
    """Type III (3A): Class A + ccr Type 3 (ccrA3B3). Ref: 85/2082."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("IS431"), make_hit("ccrA3"), make_hit("ccrB3"),
    ]


@pytest.fixture
def type_IV_hits():
    """Type IV (2B): Class B + ccr Type 2 (ccrA2B2). Ref: CA05."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("IS1272"),
        make_hit("IS431"), make_hit("ccrA2"), make_hit("ccrB2"),
    ]


@pytest.fixture
def type_V_hits():
    """Type V (5C2): Class C2 + ccr Type 5 (ccrC1). Ref: WIS.
    IS431 copies in opposite orientations."""
    return [
        make_hit("mecA"), make_hit("mecR1"),
        make_hit("IS431_1", strand="+"), make_hit("IS431_2", strand="-"),
        make_hit("ccrC1"),
    ]


@pytest.fixture
def type_VI_hits():
    """Type VI (4B): Class B + ccr Type 4 (ccrA4B4). Ref: HDE288."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("IS1272"),
        make_hit("IS431"), make_hit("ccrA4"), make_hit("ccrB4"),
    ]


@pytest.fixture
def type_VII_hits():
    """Type VII (5C1): Class C1 + ccr Type 5 (ccrC1). Ref: JCSC6082.
    IS431 copies in same orientation."""
    return [
        make_hit("mecA"), make_hit("mecR1"),
        make_hit("IS431_1", strand="+"), make_hit("IS431_2", strand="+"),
        make_hit("ccrC1"),
    ]


@pytest.fixture
def type_VIII_hits():
    """Type VIII (4A): Class A + ccr Type 4 (ccrA4B4). Ref: C10682."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("IS431"), make_hit("ccrA4"), make_hit("ccrB4"),
    ]


@pytest.fixture
def type_IX_hits():
    """Type IX (1C2): Class C2 + ccr Type 1 (ccrA1B1). Ref: JCSC6943.
    IS431 copies in opposite orientations."""
    return [
        make_hit("mecA"), make_hit("mecR1"),
        make_hit("IS431_1", strand="+"), make_hit("IS431_2", strand="-"),
        make_hit("ccrA1"), make_hit("ccrB1"),
    ]


@pytest.fixture
def type_X_hits():
    """Type X (7C1): Class C1 + ccr Type 7 (ccrA1B6). Ref: JCSC6945.
    IS431 copies in same orientation."""
    return [
        make_hit("mecA"), make_hit("mecR1"),
        make_hit("IS431_1", strand="+"), make_hit("IS431_2", strand="+"),
        make_hit("ccrA1"), make_hit("ccrB6"),
    ]


@pytest.fixture
def type_XI_hits():
    """Type XI (8E): Class E + ccr Type 8 (ccrA1B3). Ref: LGA251."""
    return [
        make_hit("mecC"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("blaZ"),
        make_hit("ccrA1"), make_hit("ccrB3"),
    ]


@pytest.fixture
def type_XII_hits():
    """Type XII (9C2): Class C2 + ccr Type 9 (ccrC2). Ref: BA01611.
    IS431 copies in opposite orientations."""
    return [
        make_hit("mecA"), make_hit("mecR1"),
        make_hit("IS431_1", strand="+"), make_hit("IS431_2", strand="-"),
        make_hit("ccrC2"),
    ]


@pytest.fixture
def type_XIII_hits():
    """Type XIII (9A): Class A + ccr Type 9 (ccrC2). Ref: 55-99-44."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("IS431"), make_hit("ccrC2"),
    ]


@pytest.fixture
def type_XIV_hits():
    """Type XIV (5A): Class A + ccr Type 5 (ccrC1). Ref: SC792."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("IS431"), make_hit("ccrC1"),
    ]


@pytest.fixture
def type_XV_hits():
    """Type XV (7A): Class A + ccr Type 7 (ccrA1B6). Ref: NV_1."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("mecI"),
        make_hit("IS431"), make_hit("ccrA1"), make_hit("ccrB6"),
    ]


@pytest.fixture
def plasmid_mecB_hits():
    """Plasmid-borne mecB. Not a chromosomal SCCmec."""
    return [make_hit("mecB")]


@pytest.fixture
def negative_hits():
    """No mec or ccr genes — e.g., MSSA."""
    return []


@pytest.fixture
def orphan_ccr_hits():
    """ccr genes found but no mec gene."""
    return [make_hit("ccrA2"), make_hit("ccrB2")]


@pytest.fixture
def composite_hits():
    """Multiple ccr types detected (e.g., ccrA2B2 + ccrC1)."""
    return [
        make_hit("mecA"), make_hit("mecR1"), make_hit("IS1272"),
        make_hit("IS431"), make_hit("ccrA2"), make_hit("ccrB2"),
        make_hit("ccrC1"),
    ]
```

- [ ] **Step 3: Create tests/test_classifier.py**

```python
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
```

- [ ] **Step 4: Run tests to verify they FAIL (classifier not yet updated)**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/test_classifier.py -v 2>&1 | head -80`
Expected: Multiple FAILures (old classifier does not support new types/pair matching)

- [ ] **Step 5: Commit test infrastructure**

```bash
git add tests/__init__.py tests/conftest.py tests/test_classifier.py
git commit -m "test: add unit tests for all 15 IWG-SCC SCCmec types"
```

---

## Chunk 3: Classifier Rewrite

### Task 3: Rewrite classifier.py with pair-based ccr matching and IS431 orientation

**Files:**
- Modify: `lib/classifier.py`

Key algorithmic changes:
1. **ccr matching**: Check for specific gene pairs (ccrA1+ccrB6 = Type 7) instead of substring matching
2. **IS431 orientation**: Compare strand of IS431_1 vs IS431_2 to distinguish C1 (same) vs C2 (opposite)
3. **Excluded genes**: Class B requires IS1272 but must NOT have mecI (which would make it Class A)
4. **Class E**: mecC-based, replaces old "Class C2" naming
5. **All 15 types**: Lookup from rules.json

- [ ] **Step 1: Rewrite lib/classifier.py**

```python
import json
import os


def load_rules():
    """Load classification rules from JSON file."""
    rules_path = os.path.join(os.path.dirname(__file__), "../db/rules.json")
    try:
        with open(rules_path, "r") as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading rules: {e}")
        return None


def _determine_is431_orientation(hits_list):
    """Determine IS431 orientation from PAF strand information.

    Compares strands of IS431_1 and IS431_2 hits to distinguish
    mec complex Class C1 (same orientation) from C2 (opposite orientation).

    Returns:
        "same"     - both IS431 copies on same strand
        "opposite" - IS431 copies on different strands
        None       - cannot determine (only one copy, or no numbered variants)
    """
    is431_strands = {}
    for h in hits_list:
        gene = h["gene"]
        if gene in ("IS431_1", "IS431_2"):
            is431_strands[gene] = h.get("strand", "+")

    if "IS431_1" in is431_strands and "IS431_2" in is431_strands:
        if is431_strands["IS431_1"] == is431_strands["IS431_2"]:
            return "same"
        else:
            return "opposite"

    return None


def _classify_mec_complex(genes_found, hits_list, rules):
    """Determine mec complex class using priority-ordered rules.

    The rules are evaluated in order. The first matching rule wins.
    Rules may specify:
      - required: all these genes must be present
      - any_of: at least one of these genes must be present
      - excluded: none of these genes may be present
      - orientation: IS431 orientation must match ("same" or "opposite")
    """
    is431_orientation = _determine_is431_orientation(hits_list)

    for rule in rules["mec_complex_rules"]:
        # Check required genes
        if not all(g in genes_found for g in rule["required"]):
            continue

        # Check any_of genes
        if "any_of" in rule:
            if not any(g in genes_found for g in rule["any_of"]):
                continue

        # Check excluded genes
        if "excluded" in rule:
            if any(g in genes_found for g in rule["excluded"]):
                continue

        # Check IS431 orientation (for C1/C2 distinction)
        if "orientation" in rule:
            if is431_orientation is None:
                # Cannot determine orientation; skip this specific rule,
                # the generic "Class C" fallback will catch it
                continue
            if is431_orientation != rule["orientation"]:
                continue

        return rule["name"]

    return "Negative"


def _classify_ccr_complex(genes_found, rules):
    """Determine ccr complex type(s) using pair-based gene matching.

    Each ccr rule specifies required_genes (e.g., ["ccrA1", "ccrB6"] for Type 7).
    A ccr type matches only when ALL its required genes are present.

    Pair-based rules are checked first (ccrAx + ccrBx), then single-gene rules
    (ccrC1, ccrC2). When a pair matches, the individual genes are "consumed"
    so they don't also match a single-component rule.

    Returns a set of matched ccr type names.
    """
    ccr_genes = {g for g in genes_found if g.startswith("ccr")}
    if not ccr_genes:
        return set()

    ccr_types = set()
    consumed_genes = set()

    # Sort rules: multi-gene pairs first, then single-gene rules
    pair_rules = [r for r in rules["ccr_complex_rules"] if len(r["required_genes"]) > 1]
    single_rules = [r for r in rules["ccr_complex_rules"] if len(r["required_genes"]) == 1]

    # Match pair rules first
    for rule in pair_rules:
        required = set(rule["required_genes"])
        available = ccr_genes - consumed_genes
        if required.issubset(available):
            ccr_types.add(rule["name"])
            consumed_genes.update(required)

    # Match single-gene rules (ccrC1, ccrC2) only if gene not consumed
    for rule in single_rules:
        required = set(rule["required_genes"])
        available = ccr_genes - consumed_genes
        if required.issubset(available):
            ccr_types.add(rule["name"])
            consumed_genes.update(required)

    return ccr_types


def classify_sccmec(hits_list):
    """Classify SCCmec type from parsed alignment hits.

    Takes parsed hits (list of dicts with at minimum a "gene" key) and
    determines the SCCmec type using the IWG-SCC classification scheme
    defined in db/rules.json.

    Algorithm:
    1. Identify mec complex class (A, B, C1, C2, D, E) from gene presence
       and IS431 orientation
    2. Identify ccr complex type(s) (1-9) from ccr gene pairs
    3. Match mec+ccr combination to one of the 15 approved SCCmec types
    4. Handle edge cases: composite, orphan ccr, partial, negative
    """
    if not hits_list:
        return {"status": "Negative", "reason": "No alignments found"}

    rules = load_rules()
    if not rules:
        return {"status": "Error", "reason": "Could not load classification rules"}

    # Unique genes detected
    genes_found = set(h["gene"] for h in hits_list)

    # Split assembly check
    contigs_found = set(h.get("contig", "unknown") for h in hits_list)
    warnings = []
    if len(contigs_found) > 1:
        warnings.append(
            f"Split Assembly: Components found on {len(contigs_found)} "
            f"different contigs: {', '.join(sorted(contigs_found))}"
        )

    # mecA/mecC/mecB presence flags
    mecA_present = "mecA" in genes_found
    mecC_present = "mecC" in genes_found
    mecB_present = "mecB" in genes_found

    # Step 1: mec complex
    mec_complex = _classify_mec_complex(genes_found, hits_list, rules)

    # Step 2: ccr complex
    ccr_types = _classify_ccr_complex(genes_found, rules)
    ccr_complex = " / ".join(sorted(ccr_types)) if ccr_types else "Negative"

    # Step 3: SCCmec type assignment
    sccmec_type = "Unknown"
    status = "Positive"

    if mec_complex == "Negative" and not ccr_types:
        return {
            "status": "Negative",
            "reason": "No mec or ccr genes found",
            "genes_detected": sorted(genes_found),
            "mecA_present": mecA_present,
        }

    if mec_complex == "Negative":
        status = "Partial (Orphan ccr)"
        warnings.append("Found ccr genes but no mecA/mecB/mecC")
    elif not ccr_types:
        # Check if mec complex allows missing ccr (e.g., Plasmid-borne)
        matched_special = False
        for rule in rules["sccmec_type_rules"]:
            if rule["mec"] == mec_complex and rule.get("ccr") == "ANY":
                sccmec_type = rule["name"]
                status = "Positive"
                matched_special = True
                break
        if not matched_special:
            status = "Partial (Unclassifiable)"
            warnings.append("Found mec gene(s) but no ccr genes")

    # Match mec+ccr combination to SCCmec type
    if status == "Positive" and sccmec_type == "Unknown":
        # For mec classes that have C1/C2 fallback to generic C:
        # If mec_complex is "Class C" (orientation undetermined), try matching
        # both C1 and C2 type rules
        mec_candidates = [mec_complex]
        if mec_complex == "Class C":
            mec_candidates = ["Class C1", "Class C2", "Class C"]

        for mec_candidate in mec_candidates:
            for rule in rules["sccmec_type_rules"]:
                if rule["mec"] != mec_candidate:
                    continue
                target_ccr = rule["ccr"]
                if target_ccr == "ANY":
                    sccmec_type = rule["name"]
                    break
                if target_ccr in ccr_types:
                    sccmec_type = rule["name"]
                    # When falling back to generic Class C, report which
                    # specific class the matched type expects
                    if mec_complex == "Class C" and mec_candidate != "Class C":
                        warnings.append(
                            f"IS431 orientation undetermined; type {rule['name']} "
                            f"expects {mec_candidate}"
                        )
                    break
            if sccmec_type != "Unknown":
                break

    # Composite detection: multiple ccr types
    if len(ccr_types) > 1 and sccmec_type != "Unknown":
        sccmec_type = f"Composite ({sccmec_type})"
        warnings.append(f"Multiple ccr types detected: {', '.join(sorted(ccr_types))}")
    elif len(ccr_types) > 1:
        sccmec_type = "Composite"
        warnings.append(f"Multiple ccr types detected: {', '.join(sorted(ccr_types))}")

    return {
        "status": status,
        "sccmec_type": sccmec_type,
        "mec_complex": mec_complex,
        "ccr_complex": ccr_complex,
        "genes_detected": sorted(genes_found),
        "mecA_present": mecA_present,
        "warnings": warnings,
        "hits_summary": hits_list,
    }
```

- [ ] **Step 2: Run all tests**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 -m pytest tests/test_classifier.py -v`
Expected: All tests PASS

- [ ] **Step 3: If any tests fail, fix and re-run**

Debug by running individual failing test:
`python3 -m pytest tests/test_classifier.py::TestSCCmecTypeAssignment::test_type_X -v -s`

- [ ] **Step 4: Commit classifier**

```bash
git add lib/classifier.py
git commit -m "feat: rewrite classifier with pair-based ccr, IS431 orientation, classes A-E"
```

---

## Chunk 4: Integration Verification

### Task 4: Run against reference strains

**Files:**
- No new files; uses existing test data in `tests/data/`

Verify the updated classifier produces correct results against the real reference genome assemblies already in the test data directory.

- [ ] **Step 1: Run against COL (expected Type I)**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 bin/sccmec_typer.py --1 tests/data/COL.fna -d db/sccmec_targets.fasta -o /tmp/COL_iwg_test && cat /tmp/COL_iwg_test.json | python3 -m json.tool`
Expected: `sccmec_type: "Type I"`, `mec_complex: "Class B"`, `ccr_complex` contains `"Type 1"`

- [ ] **Step 2: Run against LGA251 (expected Type XI)**

Run: `cd /home/alarawms/dev/sccmec_typer && python3 bin/sccmec_typer.py --1 tests/data/LGA251.fna -d db/sccmec_targets.fasta -o /tmp/LGA251_iwg_test && cat /tmp/LGA251_iwg_test.json | python3 -m json.tool`
Expected: `sccmec_type: "Type XI"`, `mec_complex: "Class E"` (was previously "Class C2"), `ccr_complex` contains `"Type 8"`

- [ ] **Step 3: Run against all available reference strains and check results**

Run each available strain (COL, N315, Mu50, MW2, USA300, LGA251, JCSC5402) and verify types match expected.

- [ ] **Step 4: Commit any fixes needed**

```bash
git add -u
git commit -m "fix: integration test corrections for IWG-SCC classification"
```

---

## Summary of Changes

| Component | Before | After |
|-----------|--------|-------|
| SCCmec types | 6 (I-VI, XI) + Plasmid | **15 (I-XV)** + Plasmid |
| mec classes | A, B, C, "C2" (mecC) | **A, B, C1, C2, C, D, E** |
| ccr types | 1-5, 8 (substring match) | **1-9 (pair-based match)** |
| ccr matching | Substring in gene name | **Explicit gene pair rules** |
| IS431 orientation | Not analyzed | **Strand comparison for C1/C2** |
| Nomenclature | Non-standard | **IWG-SCC compliant** |
