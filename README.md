# SCCmec Typer

A robust, standalone bioinformatics tool for **Staphylococcal Cassette Chromosome mec (SCCmec)** typing in *Staphylococcus aureus* genomes. It uses `minimap2` for fast, accurate alignment and a component-based classification logic (identifying *mec* and *ccr* complexes) to determine the SCCmec type.

## Features
*   **Fast**: Built on `minimap2` for rapid alignment of assemblies or long reads.
*   **Robust**: Component-based typing handles novel or composite elements better than whole-cassette matching.
*   **Accurate**: Verified against standard reference genomes (COL, N315, Mu50, MW2, USA300).
*   **Configuration-Driven**: Classification logic is defined in `db/rules.json`, allowing for easy customization and extension of SCCmec types.
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

### Running Locally
```bash
python3 bin/sccmec_typer.py \
  -i input.fasta \
  -d db/sccmec_targets.fasta \
  -o output_prefix
```

### Running with Docker
```bash
docker run --rm -v $(pwd):/data sccmec_typer \
  -i /data/input.fasta \
  -d /app/db/sccmec_targets.fasta \
  -o /data/output_prefix
```
