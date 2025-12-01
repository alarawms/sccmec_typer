# SCCmec Typer
تحليل المحتوى الجيني للكاسيت الخاص بمقاومة الميثسيلين و البيتالاكتام عن طريق بروتين ب ب ب المعدل 
  التطبيق  إلى فحص الجينوم او قراءات التسلسل الجينومي القصيرة او القراءات التسلسل الجينومي الطويلة لوجود اي جين من الجينات المكونة للكاسيت و من ثم تصنيفه الى الانواع المصفة عن طريق الجمعية الخاصة بتصنيف SCCmec
  
A robust, standalone bioinformatics tool for **Staphylococcal Cassette Chromosome mec (SCCmec)** typing in *Staphylococcus aureus* genomes. It uses `minimap2` for fast, accurate alignment and a component-based classification logic (identifying *mec* and *ccr* complexes) to determine the SCCmec type.

## Features
*   **Fast**: Built on `minimap2` for rapid alignment.
*   **Versatile**: Supports **Assembly** (FASTA), **Long Reads** (Nanopore FASTQ), and **Short Reads** (Illumina Paired-End FASTQ).
*   **Robust**: Component-based typing handles novel or composite elements better than whole-cassette matching.
*   **Accurate**: Verified against standard reference genomes (COL, N315, Mu50, MW2, USA300).
*   **Detailed Reporting**: Explicitly reports `mecA` presence and provides detailed CSV output for all detected elements.
*   **Configuration-Driven**: Classification logic is defined in `db/rules.json`, allowing for easy customization.
*   **Dockerized**: Fully containerized for easy deployment.

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

### Running with Docker
```bash
docker run --rm -v $(pwd):/data sccmec_typer \
  --1 /data/input.fasta \
  -d /app/db/sccmec_targets.fasta \
  -o /data/output_prefix
```

## Outputs
The tool generates three output files:
1.  **`{prefix}.tsv`**: A summary table containing the Sample Name, Status, **mecA Presence**, SCCmec Type, Mec Complex, Ccr Complex, and Detected Genes.
2.  **`{prefix}.json`**: A detailed JSON report with the full classification structure.
3.  **`{prefix}_elements.csv`**: A detailed list of all detected genetic elements (genes), including their coordinates, identity, and coverage.

also look at 
https://github.com/rpetit3/sccmec
