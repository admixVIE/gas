[![license](https://img.shields.io/badge/license-GPL%20v3-black.svg?style=flat-square)](LICENSE)
[![build Status](https://img.shields.io/github/actions/workflow/status/admixVIE/gas/dry-run.yaml?branch=main&style=flat-square&label=dry-run)](https://github.com/admixVIE/gas/actions)

# Great Ape Selection

## Introduction

This repository contains a Snakemake workflow designed to reproduce the results for investigating genomic landscapes of natural selection in great apes. The workflow was tested on Rocky Linux 9.5 (Blue Onyx) at [the Life Science Compute Cluster](https://lisc.univie.ac.at/) and on Rocky Linux 9.4 (Blue Onyx) at [the Multi-Site Computer Austria](https://docs.vsc.ac.at/systems/musica.html).

## Usage

1. Install [miniforge](https://github.com/conda-forge/miniforge/releases). [Mambaforge (version 23.3.1)](https://github.com/conda-forge/miniforge/releases/download/23.3.1-1/Mambaforge-23.3.1-1-Linux-x86_64.sh) was used for analysis.

2. Clone this repository:

```
git clone https://github.com/admixVIE/gas
cd gas
```

3. Create the environment:

```
mamba env create -f workflow/envs/env.yaml
```

4. Activate the environment:

```
mamba activate gas
```

5. Run the analysis for the *Pan* genus locally:

```
snakemake -s workflow/within_species.smk -c 1 --configfile config/pan.yaml
```

6. Run the analysis for the *Pan* genus on HPC:

```
snakemake -s workflow/within_species.smk -c 1 --configfile config/pan.yaml --profile config/slurm
```

7. After analyzing each genus, run a cross-species summary:

```
snakemake -s workflow/cross_species.smk -c 1
```

Users should adjust the resource parameters in each Snakemake file to match their cluster settings and modify the `config.yaml` file in `config/slurm` according to their job scheduler.
Before analysis, users should manually download [ANNOVAR](https://annovar.openbioinformatics.org/en/latest/user-guide/download/) and place it in `resources/tools`.

Note that due to random factors in some tools such as `dadi` or `Gowinda`, and the regular update of annotation files such as `gene2go.gz`, the results may show slight differences between runs.
