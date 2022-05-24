# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell

shell.executable("bash")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# =================================================================================================
#     Rule Variables
# =================================================================================================

# We except the following variables to be set:
#     snakemake.input.msa
#     snakemake.input.hmmprofile
#     snakemake.input.sample
#     snakemake.output.sequences

# =================================================================================================
#     Run
# =================================================================================================

shell(
    "hmmalign "
    "--{snakemake.params.states} "
    "--outformat afa "
    "-o {snakemake.output.sequences} "
    "--mapali {snakemake.input.msa} "
    "{snakemake.params.extra} {snakemake.input.hmmprofile} {snakemake.input.sample} "
    "{log}"
)
