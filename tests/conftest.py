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
