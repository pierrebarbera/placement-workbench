# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
from tempfile import TemporaryDirectory
import os

shell.executable("bash")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# =================================================================================================
#     Run
# =================================================================================================

# Get all files types of trees that we want to produce by their extensions in the output files.
# Also, get the output directory of the first out file, assuming that all others need to go into
# the same output dir. This works for now, and is good enough.
exts = []
out = None
for file in snakemake.output:
    ext = os.path.splitext( file )[1][1:]
    exts.append(ext)
    if out:
        if out != file[: -(len(ext) + 1)]:
            raise Exception(
                "Output file paths for different file types in gappa heat-tree rule have to "
                "be identical except for their file extension."
            )
    else:
        out = file[: -(len(ext) + 1)]

# Create command line arguments and check that the extensions are valid.
write_trees = []
for ext in exts:
    if ext not in [ "newick", "nexus", "phyloxml", "svg" ]:
        raise Exception("Unknown gappa heat-tree file extension: " + ext)
    write_trees.append( "--write-" + ext + "-tree" )

# gappa always uses the same output file names, so for independent instances of this rule,
# we have to use temp directories where we can do our work.
with TemporaryDirectory() as tempdir:

    # We need to split the combined query+ref msa into individual (temp) files.
    shell(
        "gappa examine heat-tree"
        " --jplace-path {snakemake.input} "
        " {write_trees} "
        " --out-dir {tempdir:q} "
        " {snakemake.params.extra} "
        " {log}"
    )

    # Move the result files that we care about to our actual output,
    # and remove the rest of the temp directory, out of courtesy.
    for ext in exts:
        shell( "mv {tempdir:q}/tree.{ext} {out}.{ext}" )
    shell( "rm -rf {tempdir:q}" )
