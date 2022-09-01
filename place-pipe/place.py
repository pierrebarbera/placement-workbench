#!/usr/bin/env python3

import multiprocessing
import argparse, sys, os
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, os.path.join(script_dir, '..', 'common'))
import util
from datetime import datetime
import pandas as pd
import snakemake
import platform

# if we're on mac/arm, default back down to x86 compiled packages
# ideally this decision would be per-package availability
if platform.machine() == 'arm64':
  os.environ['CONDA_SUBDIR'] = "osx-64"

parser = argparse.ArgumentParser(description='Wrapper to run the pipeline with some default settings.')
###
#  Input Files
###
input_group = parser.add_argument_group('Input')
input_group.add_argument('--fasta-paths', dest='fasta_files', type=str, nargs='+',
                    help='input fasta files', required=True)
input_group.add_argument('--reference-tree', dest='ref_tree', type=str,
                    help='Reference tree, in newick format', required=True)
input_group.add_argument('--reference-msa', dest='ref_msa', type=str,
                    help='Reference MSA, in fasta format', required=True)
input_group.add_argument('--model-file', dest='model_file', type=str,
                    help='Reference tree, in newick format', required=True)
input_group.add_argument('--taxonomy-file', dest='taxon_file', type=str,
                    help='Tab-separated file mapping reference labels to their taxonomic paths', required=False)

###
#   Config Manip
###
pipeline_group = parser.add_argument_group('Pipeline Options')
pipeline_group.add_argument('--out-dir', dest='out_dir', type=str,
                    default=None,
                    help='optional output directory. By default, the script creates a timestamped output directory.')
pipeline_group.add_argument('-d', '--datatype', dest='datatype', type=str, nargs='?',
                    const='nt', default='nt', choices=['nt', 'aa'],
                    help="datatype, 'aa' for protein, 'nt' for DNA data")
pipeline_group.add_argument('-p', '--prefix', dest='prefix', type=str,
                    default=None,
                    help='prefix to fasta paths and output (useful to specify where data was mounted to in docker)')
pipeline_group.add_argument('--no-chunkify', dest='no_chunkify', action='store_true',
                    help='use the chunkify routine to split queries into blocks of more managable size')
pipeline_group.add_argument('--sequence-clustering', dest='sequence_clustering', type=str, nargs='?',
                    action='store', default=None, choices=['swarm'],
                    help="Which tools, if any, to use for query clustering.")

###
#   Cluster options
###
cluster_group = parser.add_argument_group('Cluster Execution')
cluster_group.add_argument('--cluster-exec', dest='on_cluster', action='store_true',
                    help="Starts the pipeline in computing-cluster (slurm/sge etc.) submission mode. "
                    "Highly recommended to do this from a screen/tmux session!")
cluster_group.add_argument('--cluster-env', dest='clust_env', type=str, nargs='?',
                    const='auto', default='auto', choices=['auto','slurm', 'sge'],
                    help="What job submission system we are on. 'auto' attempts to autodetect.")

###
#   Runtime Manip
###
misc_group = parser.add_argument_group('Misc. Options')
misc_group.add_argument('--threads', dest='threads', type=int,
                    default=multiprocessing.cpu_count(),
                    help='number of threads to use')
misc_group.add_argument('-v','--verbose', dest='verbose', action='store_true',
                    help='increase verbosity')
args = parser.parse_args()

if len(args.fasta_files) + len(args.unmerged_fastq_files) == 0:
  util.fail( "Must supply either query fasta or unmerged F/R fastq files" )
elif len(args.unmerged_fastq_files) % 2 != 0:
  util.fail( "When supplying unmerged F/R fastq files, must supply separate forward and reverse files." )

util.expect_file_exists( args.ref_tree )
util.expect_file_exists( args.ref_msa )
util.expect_file_exists( args.model_file )

if args.prefix:
  util.expect_dir_exists( args.prefix )

if args.taxon_file:
  util.expect_file_exists( args.taxon_file )

clustering_tool = args.sequence_clustering if args.sequence_clustering else "no_clustering"

use_chunkify = not args.no_chunkify


# make a unique output dir, labeled by date and time
out_dir = "run-{}".format(datetime.now().strftime("%Y-%m-%d-%H:%M:%S")) if( not args.out_dir ) else args.out_dir

# build up a samples.tsv file to be used in this execution
# 
# first read in a unique list of all the desired files
file_paths = []
if args.fasta_files:
  file_paths.extend( util.ingest_paths( args.fasta_files,
                                        extensions=['.fa', '.afa', '.fasta'],
                                        allow_gz=True ) )


# add file paths to the samples, giving it a sample name corresponding to the file name / directory
samples = pd.DataFrame({
  'sample':util.get_unique_names(file_paths),
  'input_file':file_paths
  }).set_index( 'sample' )

# finally, write the samples.tsv to the output folder
util.make_path( out_dir )
samples_file = os.path.join( out_dir, "samples.tsv" )
samples.to_csv( samples_file, sep='\t' )

# next, we set those config values that we wish to override from the defaults, by creating a dict
# of those values
config_overrrides = {
  'data':
  {
    'samples': samples_file,
    'reference-tree': args.ref_tree,
    'reference-alignment': args.ref_msa
  },
  'settings':
  {
    'datatype': args.datatype,
    'clustering-tool': clustering_tool,
    'use-chunkify': use_chunkify,
    'outdir': out_dir,
  },
  'params':
  {
    'threads': args.threads
  }
}

if args.taxon_file:
  config_overrrides['data']['taxonomy_file'] = args.taxon_file

calling_dir = os.path.dirname(os.path.abspath(__file__))

# check if mamba exists
conda_front = 'conda'
if util.is_tool('mamba'):
  conda_front = 'mamba'

# get cluster settings
cluster, cluster_config = (None, None) if not args.on_cluster else util.cluster_settings( args.clust_env, calling_dir )

snakemake.snakemake(
  snakefile=os.path.join( calling_dir, "Snakefile" ),
  use_conda=True,
  conda_frontend=conda_front,
  cores=args.threads,
  config=config_overrrides,
  cluster_config=cluster_config,
  cluster=cluster
  )
