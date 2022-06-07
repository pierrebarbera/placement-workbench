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
conda env plwb
```

This step has to be repeated every time you re-open your terminal or re-connect to a server. By default, mamba will add information about which is the active environment, to your command line:
```
(base) user@computer$ ... # by default, we are in the base environment
# now we enable our environment:
conda env plwb
(plwb) user@computer$ ... # tada!
```

### 4) thats it! ready to go
No, seriously, thats it.

## Basic Usage

### Examples
Both pipelines come with a convenient python script that should cover the majority of use-cases. 

By default, these scripts will output to a directory called `run-` followed by a timestamp. You can change the output directory by adding the `--out-dir ...` option. Use the `--help` command for a full list of options.

Some examples:
#### reftree-pipe
Infer phylogenetic trees from a set of unaligned sequences (`data/sequences.fasta`). It's DNA/Nucleotide data (`nt`), and we want to just use 4 threads:
```
./search.py --fasta-paths data/sequences.fasta  --align --datatype nt --threads 4
```

Same as the previous, except now we get some sequences by specifying their genbank accession labels in a file (`data/accessions.csv`), which are then automatically downloaded.
```
./search.py --csv-paths data/accessions.csv  --align --datatype nt --threads 4
```

You can also combine both commands, add multiple files each, and even add all fasta/csv files by the directory they're in:
```
./search.py --fasta-paths a.fasta b.fa data/ --csv-paths c.csv data/ ...
```

#### place-pipe
```
./place.py --fasta-paths data/query.fa --reference-tree data/reference/tree.newick --reference-msa data/reference/seqs.fa --model-file data/reference/model --threads 4
```

## Advanced Usage and Configuration

