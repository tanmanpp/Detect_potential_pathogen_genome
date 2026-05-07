# Web GUI for the pathogen detection pipeline

This project now includes a local browser interface for users who are not
comfortable with command-line workflows. The web UI uploads a reads file to a
small local Python API server, then the server runs `tools/Detect_potential_pathogen_pipeline.py`.

## Start the web UI

Create both conda environments first:

```bash
conda env create -f environment.yml
conda env create -f environment.recentrifuge.yml
conda activate DPPP
```

Windows users should create and run these environments inside WSL because the pipeline
depends on Linux bioinformatics tools.

The main pipeline runs in `DPPP`. Recentrifuge runs in a separate environment
named `DPPP-rcf`, because it may require a newer Python version.

On Windows, open WSL first and move to the project folder. For example:

```bash
cd /mnt/d/YiLun/Detect_potential_pathogen_genome_for_publish
```

Then start the GUI from WSL, macOS, or Linux:

```bash
bash start_pipeline_web.sh
```

The shell script activates `DPPP` and creates both `DPPP` and `DPPP-rcf` from
their YAML files if needed.

Then open:

```text
http://127.0.0.1:8080
```

Keep the terminal open while the pipeline is running.

## Runtime behavior

The web app supports three runtime modes:

- `auto`: Windows uses WSL when `wsl.exe` is available; otherwise it uses native Python.
- `wsl`: always run the pipeline through WSL with `python3`.
- `native`: run the pipeline with the same Python that started the web server.

On Windows + WSL, uploaded files and project paths are automatically translated
from Windows paths such as `D:\...` to WSL paths such as `/mnt/d/...`.

## Optional environment variables

Set these before starting the server if needed:

```text
PATHOGEN_WEB_RUNTIME=wsl
PATHOGEN_WEB_WSL_DISTRO=Ubuntu
PATHOGEN_RCF_ENV=DPPP-rcf
```

`PATHOGEN_WEB_WSL_DISTRO` is optional. Leave it empty to use the default WSL
distribution.

`PATHOGEN_RCF_ENV` is optional. Use it only if your Recentrifuge environment is
named something other than `DPPP-rcf`.

## Flye behavior

If Flye fails or does not produce a usable `assembly.fasta`, the pipeline keeps
running. Reads-based Kraken/Recentrifuge, mapping, and VFDB blastn still run;
denovo-contig Kraken/Recentrifuge and denovo-contig blastn are skipped.

## Output

By default, web-submitted jobs write to:

```text
web_runs/results/<sample_id>/
```

The browser monitor shows the live pipeline log and links files from:

```text
web_runs/results/<sample_id>/final_results/
```

Users can also fill an output folder in the advanced settings.

## Windows troubleshooting

There is no Windows `.bat` launcher. Open WSL manually and run:

```bash
cd /mnt/d/YiLun/Detect_potential_pathogen_genome_for_publish
bash start_pipeline_web.sh
```

After the server starts, open `http://127.0.0.1:8080` in the Windows browser.
