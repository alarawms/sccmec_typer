FROM continuumio/miniconda3

WORKDIR /app

# Copy environment file and create environment
COPY environment.yml .
RUN conda env create -f environment.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "sccmec_typer", "/bin/bash", "-c"]

# Copy source code
COPY . .

# Set entrypoint
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "sccmec_typer", "python", "bin/sccmec_typer.py"]
