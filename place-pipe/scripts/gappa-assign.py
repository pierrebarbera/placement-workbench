# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import util

shell.executable("bash")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Get the output directory
sample_outdir = util.dirname( snakemake.output )

# =================================================================================================
#     Run
# =================================================================================================

shell(
    "gappa examine assign"
    " --jplace-path {snakemake.input.jplace}"
    " --taxon-file {snakemake.input.taxon_file}"
    " --out-dir {sample_outdir}"
    " --threads {snakemake.threads}"
    " {snakemake.params.extra}"
    " {log}"
)
