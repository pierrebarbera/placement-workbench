# The Placement Workbench
A highly automated, configurable pair of pipelines, for creation and engineering of reference trees, and for their subsequent use with phylogenetic placement

## Structure of this repository
This repository contains two separate snakemake pipelines:
`place-pipe`:

`reftree-pipe`:

## Setup
### 1) install conda/mamba
The pipelines are built on top of `snakemake` and `conda`. Thankfully, we can get `snakemake` from `conda`, which means we only have one real dependency! In principle everything should work with just the normal `conda` installation, however I've had much better results using `mamba`.
You can find the instructions for installing `mamba` [here](https://github.com/mamba-org/mamba), though to summarize, the current reccomendation is:
1) install [miniconda](https://docs.conda.io/en/latest/miniconda.html) (by downloading the installer and running it)
2) use `conda` to install mamba into the base environment:
```
conda install mamba -n base -c conda-forge
```

### 2) create the environment
Next, we create our special "place-workbench" environment (in the  repository root directory):
```
mamba env create -v environment.yml
```
### 3) enable the environment
Then, we activate that environment.
```
mamba env plwb
```

This step has to be repeated every time you re-open your terminal or re-connect to a server. By default, mamba will add information about which is the active environment, to your command line:
```
(base) user@computer$ ... # by default, we are in the base environment
# now we enable our environment:
mamba env place-workbench
(place-workbench) user@computer$ ... # tada!
```

### 4) thats it! ready to go
No, seriously, thats it.

## Basic Usage

### First steps

## Advanced Usage and Configuration

