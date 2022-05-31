# The Placement Workbench
A highly automated, configurable pair of pipelines, for creation and engineering of reference trees, and for their subsequent use with phylogenetic placement

## Structure of this repository
This repository contains two separate snakemake pipelines:
`place-pipe`:

`reftree-pipe`:

## Setup
### 1) install conda/mamba


### 2) create the environment
Next, we create our special "place-workbench" environment (in the  repository root directory):
```
conda env create -v environment
```
### 3) enable the environment
```
conda env place-workbench
```

This step has to be repeated every time you re-open your terminal or re-connect to a server. By default, conda will add information about which is the active environment, to your command line:
```
(base) user@computer$ ... # by default, we are in the base environment
# now we enable our environment:
conda env place-workbench
(place-workbench) user@computer$ ... # tada!
```

### 4) thats it! ready to go
No, seriously, thats it.

## Basic Usage

### First steps

## Advanced Usage and Configuration

