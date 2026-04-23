# Detect potential pathogen genome pipeline

This repository provides a Python-based pathogen detection pipeline and a local
browser GUI for colleagues who are less familiar with command-line workflows.

## 1. Create the DPPP conda environment

Windows users should create this environment inside WSL. macOS/Linux users can
create it directly in their terminal.

```bash
conda env create -f environment.yml
conda activate DPPP
```

The environment name is fixed as `DPPP`. It includes Python packages and command
line tools used by the pipeline, including `seqkit`, `minimap2`, `samtools`,
`flye`, `kraken2`, `blast`, `recentrifuge`, `pandas`, `openpyxl`, `pysam`, and
`colorama`.

## 2. Start the browser GUI

Windows users should open WSL first, then move to the project folder. For
example, if the repository is stored on `D:`:

```bash
cd /mnt/d/YiLun/Detect_potential_pathogen_genome_for_publish
```

Then start the GUI from WSL, macOS, or Linux:

```bash
bash start_pipeline_web.sh
```

The shell script activates `DPPP` and creates it from `environment.yml` if the
environment is missing.

Then open:

```text
http://127.0.0.1:8080
```

Keep the terminal open while jobs are running.

## 3. Command-line usage

Advanced users can still run the original command-line pipeline:

```bash
conda activate DPPP
python tools/Detect_potential_pathogen_pipeline.py \
  -i reads.fastq \
  -sid Sample_001 \
  -o output \
  --human_index /path/to/hg38.mmi
```

## 4. Flye failure behavior

If Flye fails or does not produce a non-empty `assembly.fasta`, the pipeline no
longer stops. It continues with reads-based analyses and skips denovo-contig
Kraken/Recentrifuge and blastn analyses.

## 5. Output

GUI-submitted jobs write to:

```text
web_runs/results/<sample_id>/
```

The final report and collected files are in:

```text
web_runs/results/<sample_id>/final_results/
```

## 6. Windows notes

This project does not use a Windows `.bat` launcher. The recommended Windows
workflow is to run the server inside WSL:

```bash
cd /mnt/d/YiLun/Detect_potential_pathogen_genome_for_publish
bash start_pipeline_web.sh
```

After the server starts, open `http://127.0.0.1:8080` in the Windows browser.
