# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import util
import snakeparser as sp
import os

shell.executable("bash")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get the output directory
sample_outdir = util.dirname( snakemake.output )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "gappa prepare phat", snakemake )

# Required args
ps.add( snakemake.input.taxonomy_file,  "--taxonomy-file {}",   sp.typ.FILE )
ps.add( snakemake.input.sequence_file,  "--sequence-file {}",   sp.typ.FILE )
ps.add( snakemake.params.target_size,   "--target-size {}",     sp.typ.UINT )

# Optional args
ps.add_opt( "sub_taxonomy",          "--sub-taxonomy {}" )
ps.add_opt( "min_subclade_size",     "--min-subclade-size {}",  sp.typ.UINT )
ps.add_opt( "max_subclade_size",     "--max-subclade-size {}",  sp.typ.UINT )
ps.add_opt( "min_tax_level",         "--min-tax-level {}",      sp.typ.UINT )
ps.add_opt( "allow_approximation",   "--allow-approximation",   sp.typ.FLAG )
ps.add_opt( "no_taxa_selection",     "--no-taxa-selection",     sp.typ.FLAG )
ps.add_opt( "consensus_method",      "--consensus-method {}", 
    sp.typ.IN(['majorities','cavener','threshold']) )
ps.add_opt( "consensus_threshold",   "--consensus-threshold {}",
    sp.typ.FLOAT(0.0, 1.0) )

# Closing args
ps.add( sample_outdir, "--out-dir {}", sp.typ.DIR )
if snakemake.threads: ps.add( snakemake.threads, "--threads {}", sp.typ.UINT )
ps.add_opt( "extra" )
ps.add( log )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.shell_string )

# also, move the result file to what is expected by the output
shell( "mv {} {}".format( 
    os.path.join( sample_outdir, "consensus_sequences.fa" ), 
    snakemake.output )
)
