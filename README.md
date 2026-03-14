# SCCmec Typer
تحليل المحتوى الجيني للكاسيت الخاص بمقاومة الميثسيلين و البيتالاكتام عن طريق بروتين ب ب ب المعدل
  التطبيق  إلى فحص الجينوم او قراءات التسلسل الجينومي القصيرة او القراءات التسلسل الجينومي الطويلة لوجود اي جين من الجينات المكونة للكاسيت و من ثم تصنيفه الى الانواع المصفة عن طريق الجمعية الخاصة بتصنيف SCCmec

A robust, standalone bioinformatics tool for **Staphylococcal Cassette Chromosome mec (SCCmec)** typing in *Staphylococcus aureus* genomes. It uses `minimap2` for fast, accurate alignment and a component-based classification logic (identifying *mec* and *ccr* complexes) to determine the SCCmec type according to the **IWG-SCC** (International Working Group on the Classification of Staphylococcal Cassette Chromosome Elements) guidelines.

## Features
*   **Complete IWG-SCC Coverage**: Supports all **15 approved SCCmec types (I–XV)** plus plasmid-borne mecB.
*   **Scientifically Accurate**: Pair-based ccr complex matching (types 1–9), proper mec complex classes (A, B, C1, C2, D, E), and IS431 orientation analysis for C1/C2 distinction.
*   **Three-Level Confidence Scoring**: Per-gene, per-component, and per-type confidence scores with High/Medium/Low tiers.
*   **IWG-SCC Visualization**: Color-coded SVG cassette diagrams and genomic coordinate maps following IWG-SCC color conventions.
*   **Interactive HTML Reports**: Self-contained HTML reports with embedded SVGs, gene detail tables, and hover tooltips.
*   **Pipeline-Ready JSON**: Enhanced JSON output with confidence scores, assembly metadata, and IWG-SCC designations for downstream integration (e.g., Staphit pipeline).
*   **Fast**: Built on `minimap2` for rapid alignment.
*   **Versatile**: Supports **Assembly** (FASTA), **Long Reads** (Nanopore FASTQ), and **Short Reads** (Illumina Paired-End FASTQ).
*   **Configuration-Driven**: Classification logic defined in `db/rules.json` — add new types by editing a file, not code.
*   **Dockerized**: Fully containerized for easy deployment.

## Supported SCCmec Types

| Type | Designation | mec Class | ccr Type | ccr Genes | Reference Strain |
|------|-------------|-----------|----------|-----------|------------------|
| I | 1B | B | 1 | ccrA1B1 | NCTC10442 |
| II | 2A | A | 2 | ccrA2B2 | N315 |
| III | 3A | A | 3 | ccrA3B3 | 85/2082 |
| IV | 2B | B | 2 | ccrA2B2 | CA05 |
| V | 5C2 | C2 | 5 | ccrC1 | WIS |
| VI | 4B | B | 4 | ccrA4B4 | HDE288 |
| VII | 5C1 | C1 | 5 | ccrC1 | JCSC6082 |
| VIII | 4A | A | 4 | ccrA4B4 | C10682 |
| IX | 1C2 | C2 | 1 | ccrA1B1 | JCSC6943 |
| X | 7C1 | C1 | 7 | ccrA1B6 | JCSC6945 |
| XI | 8E | E | 8 | ccrA1B3 | LGA251 |
| XII | 9C2 | C2 | 9 | ccrC2 | BA01611 |
| XIII | 9A | A | 9 | ccrC2 | 55-99-44 |
| XIV | 5A | A | 5 | ccrC1 | SC792 |
| XV | 7A | A | 7 | ccrA1B6 | NV_1 |

## Installation

### Local Development
1.  Clone the repository:
    ```bash
    git clone https://github.com/alarawms/sccmec_typer.git
    cd sccmec_typer
    ```
2.  Set up the environment (creates `.venv`, installs `pytest`, downloads `minimap2`):
    ```bash
    ./setup_dev.sh
    source .venv/bin/activate
    ```
3.  Install Jinja2 (required for visualization):
    ```bash
    pip install jinja2
    ```

### Docker
Build the image locally:
```bash
docker build -t sccmec_typer .
```

## Usage

The tool accepts assemblies or reads (single or paired-end).

### 1. Assembly Mode
For finished genomes or contigs (FASTA).
```bash
python3 bin/sccmec_typer.py \
  --1 assembly.fasta \
  -d db/sccmec_targets.fasta \
  -o output_prefix
```

### 2. Read Mode (Long Reads / Single-End)
For Nanopore or PacBio reads (FASTQ). The tool automatically detects FASTQ input and switches to read mode.
```bash
python3 bin/sccmec_typer.py \
  --1 reads.fastq \
  -d db/sccmec_targets.fasta \
  -o output_prefix
```

### 3. Read Mode (Paired-End Short Reads)
For Illumina paired-end data.
```bash
python3 bin/sccmec_typer.py \
  --1 reads_R1.fastq \
  --2 reads_R2.fastq \
  -d db/sccmec_targets.fasta \
  -o output_prefix
```

### Pipeline Mode (JSON only, no visualization)
```bash
python3 bin/sccmec_typer.py \
  --1 assembly.fasta \
  -d db/sccmec_targets.fasta \
  -o output_prefix \
  --no-viz
```

### Running with Docker
```bash
docker run --rm -v $(pwd):/data sccmec_typer \
  --1 /data/input.fasta \
  -d /app/db/sccmec_targets.fasta \
  -o /data/output_prefix
```

## Outputs

The tool generates five output files per sample:

| File | Format | Description |
|------|--------|-------------|
| `{prefix}.json` | JSON | Full classification with confidence scores, assembly metadata, and IWG-SCC designation |
| `{prefix}.tsv` | TSV | One-line summary (Sample, Status, mecA, Type, mec Complex, ccr Complex, Genes, Warnings) |
| `{prefix}_elements.csv` | CSV | Per-gene alignment details (coordinates, identity, coverage) |
| `{prefix}_map.svg` | SVG | IWG-SCC color-coded schematic cassette diagram |
| `{prefix}_report.html` | HTML | Interactive report with diagrams, gene table, and confidence scores |

Use `--no-viz` to skip SVG/HTML generation (JSON/TSV/CSV still produced).

### JSON Output Structure

The JSON output includes all data needed for downstream pipeline integration:

```json
{
    "status": "Positive",
    "sccmec_type": "Type IV",
    "designation": "2B",
    "mec_complex": "Class B",
    "ccr_complex": "Type 2",
    "ccr_genes": ["ccrA2", "ccrB2"],
    "mecA_present": true,
    "confidence": {
        "mode": "full",
        "type_level": {"score": 0.92, "tier": "High"},
        "mec_component": {"class": "Class B", "score": 0.99, "tier": "High", "...": "..."},
        "ccr_component": {"type": "Type 2", "score": 0.93, "tier": "High", "...": "..."},
        "per_gene": {"mecA": {"score": 1.0, "tier": "High"}, "...": "..."}
    },
    "assembly": {
        "contigs": ["contig_1"],
        "is_split": false,
        "gene_locations": [{"gene": "mecA", "contig": "contig_1", "start": 28100, "end": 30100, "strand": "+"}]
    }
}
```

### Confidence Scoring

Three levels of confidence assessment:

| Level | Formula | What it measures |
|-------|---------|-----------------|
| **Per-gene** | `0.6 * norm_identity + 0.4 * norm_coverage` | Alignment quality of each detected gene |
| **Per-component** | `(detected/expected) * mean(gene_scores)` | Completeness and quality of mec/ccr complex |
| **Per-type** | `mec_conf * ccr_conf * assembly_penalty` | Overall trust in the SCCmec type assignment |

Tiers: **High** (>0.70), **Medium** (0.40–0.70), **Low** (<0.40)

See `docs/confidence_scoring_methods.md` for full methodology.

## Classification Algorithm

1. **Align** reference gene sequences against the input genome using minimap2
2. **Filter** hits by identity (≥90%) and coverage (≥80%)
3. **Identify mec complex** class (A, B, C1, C2, D, E) from gene presence, IS element markers, and IS431 orientation
4. **Identify ccr complex** type (1–9) from ccr gene pairs (e.g., ccrA1 + ccrB6 = Type 7)
5. **Match** mec + ccr combination to one of the 15 approved SCCmec types
6. **Score** confidence at gene, component, and type levels
7. **Generate** visualization and reports

## References

- IWG-SCC Guidelines: [Antimicrobial Agents and Chemotherapy, 2009](https://journals.asm.org/doi/full/10.1128/aac.00579-09)
- Current Status of SCCmec: [Antibiotics, 2022](https://www.mdpi.com/2079-6382/11/1/86)
- SCCmec Type XV: [J Antimicrob Chemother, 2022](https://academic.oup.com/jac/article/77/4/903/6510785)

## See Also

- [rpetit3/sccmec](https://github.com/rpetit3/sccmec) — Alternative SCCmec typing tool using BLAST+
