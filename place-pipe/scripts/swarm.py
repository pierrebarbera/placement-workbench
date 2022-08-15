# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import snakeparser as sp
import os

shell.executable("bash")

# Get the output directory
sample_outdir = os.path.dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "swarm", snakemake )

# Output args
ps.add( snakemake.output.seeds,             "--seeds {}" )
ps.add( snakemake.output.statistics_file,   "--statistics-file {}" )

# General options
ps.add_threads()

# Clustering options
ps.add_opt( "differences",      "--differences {}",     sp.typ.UINT )
ps.add_opt( "no_otu_breaking",  "--no-otu-breaking",    sp.typ.FLAG )
# Fastidious options
ps.add_opt( "boundary",         "--boundary {}",        sp.typ.UINT )
ps.add_opt( "ceiling",          "--ceiling {}",         sp.typ.UINT )
ps.add_opt( "fastidious",       "--fastidious",         sp.typ.FLAG )
ps.add_opt( "bloom_bits",       "--bloom-bits {}",      sp.typ.UINT )

# (other) output options
ps.add_opt( "append_abundance",     "--append-abundance {}",    sp.typ.UINT )
ps.add_opt( "internal_structure",   "--internal-structure {}" )
ps.add_opt( "network_file",         "--network-file {}" )
ps.add_opt( "output_file",          "--output-file {}" )
ps.add_opt( "mothur",               "--mothur",                 sp.typ.FLAG )
ps.add_opt( "uclust_file",          "--uclust-file {}" )
ps.add_opt( "usearch_abundance",    "--usearch-abundance",      sp.typ.FLAG )


# Pairwise alignment advanced options
ps.add_opt( "match_reward",         "--match-reward {}",            sp.typ.UINT )
ps.add_opt( "mismatch_penalty",     "--mismatch-penalty {}",        sp.typ.UINT )
ps.add_opt( "gap_opening_penalty",  "--gap-opening-penalty {}",     sp.typ.UINT )
ps.add_opt( "gap_extension_penalty","--gap-extension-penalty {}",   sp.typ.UINT )

# Required fasta input file at the end
ps.add( snakemake.input[0], "{}",    sp.typ.FILE )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
