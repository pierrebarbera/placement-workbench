# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
from tempfile import TemporaryDirectory

shell.executable("bash")

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# =================================================================================================
#     Rule Variables
# =================================================================================================

# We except the following variables to be set:
#     snakemake.input.sequences - containing all aligned sequences (MSA + queries)
#     snakemake.input.msa
#     snakemake.input.tree
#     snakemake.output.jplace
# and one of
#     snakemake.input.model
#     snakemake.params.model
# for the model - either a RAxML model file (in the former case), or a model string as expected by
# RAxML and epa-ng (in the latter case).

# Excatly one of the ways of model input has to be given.
if bool(snakemake.input.model) == bool(snakemake.params.model):
    raise Exception(
        "Exactly one of input:model and params:model has to be given in the epa-ng snakemake rule."
    )

# =================================================================================================
#     Run
# =================================================================================================

# epa-ng always uses the same output file names, so for parallel instances of this rule,
# we have to use temp directories where we can do our work.
with TemporaryDirectory() as tempdir:

    # We need to split the combined query+ref msa into individual (temp) files.
    shell(
        "epa-ng "
        '--split "{snakemake.input.msa}" {snakemake.input.sequences} '
        "--out-dir {tempdir:q} "
        "{log}"
    )
    # The output of this step are two files:
    # query.fasta
    # reference.fasta
    # which we use in the following

    # Now run the actual placement. We here do not use the original reference alignment,
    # but the one resulting from the above split step instead, as this is guaranteed to have
    # the same width as the query alignment. Unfortunately, some aligners (such as version 3 of
    # hmmalign) mess around with gaps, so that the produced alignment does not fit with the
    # original any more. Hence, we cannot use the original here, and use the split one instead.
    shell(
        "epa-ng "
        # "--redo "
        "--query {tempdir:q}/query.fasta "
        "--ref-msa {tempdir:q}/reference.fasta "
        '--tree "{snakemake.input.tree}" '
        "--outdir {tempdir:q} "
        '--model "{snakemake.input.model}{snakemake.params.model}" '
        "--threads {snakemake.threads} "
        "{snakemake.params.extra} "
        "{log}"
    )

    # Finally, move the result files that we care about to our actual output,
    # and remove the rest of the temp directory, out of courtesy.
    shell( "mv {tempdir:q}/epa_result.jplace {snakemake.output.jplace}" )
    # "mv {tempdir:q}/epa_info.log logs/.../epa_info.log"
    shell( "rm -rf {tempdir:q}" )
