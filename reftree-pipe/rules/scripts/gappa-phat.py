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

ps = sp.Parser("gappa prepare phat")

# Parse arguments
# Required args
ps.add( snakemake.input.taxonomy_file, "--taxonomy-file {}", sp.typ.FILE )
ps.add( snakemake.input.sequence_file, "--sequence-file {}", sp.typ.FILE )
ps.add( snakemake.params.target_size, "--target-size {}" )

# Optional args
if snakemake.input.sub_taxonomy:
    ps.add( snakemake.params.sub_taxonomy, "--sub-taxonomy {}" )
if snakemake.params.min_subclade_size:
    ps.add( snakemake.params.min_subclade_size, "--min-subclade-size {}" )
if snakemake.params.max_subclade_size:
    ps.add( snakemake.params.max_subclade_size, "--max-subclade-size {}" )
if snakemake.params.min_tax_level:
    ps.add( snakemake.params.min_tax_level, "--min-tax-level {}" )
if snakemake.params.allow_approximation:
    ps.add( snakemake.params.allow_approximation, "--allow-approximation", sp.typ.FLAG )
if snakemake.params.no_taxa_selection:
    ps.add( snakemake.params.no_taxa_selection, "--no-taxa-selection", sp.typ.FLAG )
if snakemake.params.consensus_method:
    ps.add( snakemake.params.consensus_method, "--consensus-method" )
if snakemake.params.consensus_threshold:
    ps.add( snakemake.params.consensus_threshold, "--consensus-threshold" )

# Closing args
ps.add( sample_outdir, "--out-dir {}" , sp.typ.DIR )
ps.add( snakemake.threads, "--threads {}" )
if snakemake.params.extra:
    ps.add( snakemake.params.extra )
ps.add( log )

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.shell_string +
        # also, move the result file to what is expected by the output
        " && mv {} {}".format( os.path.join( sample_outdir, "consensus_sequences.fa" ), snakemake.output ) 
)
