# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import util
import snakeparser as sp


shell.executable("bash")

# Get the output directory
outdir = util.dirname( snakemake.output[0] )

# =================================================================================================
#     File prep / aggregation
# =================================================================================================
# If there is only one file on the input, we assume its a file containing multiple newick trees
if len(snakemake.input) == 1:
    ml_trees = snakemake.input[0]
# If there are multiple, we assume that multiple newick files with single trees are specified
# In this case we combine them into one file first
else:
    util.make_path( outdir )
    ml_trees = os.path.join( outdir, "ml_trees.newick" )
    with open( ml_trees, 'w' ) as outfile:
        for in_path in snakemake.input:
            with open( in_path, 'r' ) as infile:
                outfile.write( infile.read() )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser("raxml-ng", snakemake)

# select the run mode
ps.add( "--rfdist" )

# Required args
ps.add( ml_trees, "--tree {}", sp.typ.FILE )

# Optional args


# Closing args
ps.add( os.path.join( outdir, "rf_calc" ), "--prefix {}" )
ps.add_threads()


# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )

result_file = os.path.join( outdir, "rf_calc.raxml.rfDistances" )
util.expect_file_exists( result_file )

# rename the result file
shell( f"mv {result_file} {snakemake.output[0]}" )

