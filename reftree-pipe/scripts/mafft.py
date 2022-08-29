# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
import sys, os
common_dir = os.path.abspath(os.path.join( os.path.dirname(__file__), "..", "..", "common" ))
sys.path.insert(0, common_dir)
import snakeparser as sp

shell.executable("bash")

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "mafft", snakemake )

# as per https://mafft.cbrc.jp/alignment/software/manual/manual.html
# OPTIONS

# Algorithm
ps.add_opt( "auto",         "--auto",           sp.typ.FLAG )
ps.add_opt( "6merpair",     "--6merpair",       sp.typ.FLAG )
ps.add_opt( "globalpair",   "--globalpair",     sp.typ.FLAG )
ps.add_opt( "localpair",    "--localpair",      sp.typ.FLAG )
ps.add_opt( "genafpair",    "--genafpair",      sp.typ.FLAG )
ps.add_opt( "fastapair",    "--fastapair",      sp.typ.FLAG )
ps.add_opt( "weighti",      "--weighti {}",     sp.typ.UINT )
ps.add_opt( "retree",       "--retree {}",      sp.typ.UINT )
ps.add_opt( "maxiterate",   "--maxiterate {}",  sp.typ.UINT )
ps.add_opt( "fft",          "--fft",            sp.typ.FLAG )
ps.add_opt( "nofft",        "--nofft",          sp.typ.FLAG )
ps.add_opt( "noscore",      "--noscore",        sp.typ.FLAG )
ps.add_opt( "memsave",      "--memsave",        sp.typ.FLAG )
ps.add_opt( "parttree",     "--parttree",       sp.typ.FLAG )
ps.add_opt( "dpparttree",   "--dpparttree",     sp.typ.FLAG )
ps.add_opt( "fastaparttree","--fastaparttree",  sp.typ.FLAG )
ps.add_opt( "partsize",     "--partsize {}",    sp.typ.UINT )
ps.add_opt( "groupsize",    "--groupsize {}",   sp.typ.UINT )

# Parameter
ps.add_opt( "op",           "--op {}",          sp.typ.FLOAT )
ps.add_opt( "ep",           "--ep {}",          sp.typ.FLOAT )
ps.add_opt( "lop",          "--lop {}",         sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "lep",          "--lep {}",         sp.typ.FLOAT )
ps.add_opt( "lexp",         "--lexp {}",        sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "LOP",          "--LOP {}",         sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "LEXP",         "--LEXP {}",        sp.typ.FLOAT(float("-inf"),0.0) )
ps.add_opt( "bl",           "--bl {}",          sp.typ.UINT )
ps.add_opt( "jtt",          "--jtt {}",         sp.typ.UINT )
ps.add_opt( "tm",           "--tm {}",          sp.typ.UINT )
ps.add_opt( "aamatrix",     "--aamatrix {}",    sp.typ.FILE )
ps.add_opt( "fmodel",       "--fmodel",         sp.typ.FLAG )

# Output
ps.add_opt( "clustalout",   "--clustalout",     sp.typ.FLAG )
ps.add_opt( "inputorder",   "--inputorder",     sp.typ.FLAG )
ps.add_opt( "reorder",      "--reorder",        sp.typ.FLAG )
ps.add_opt( "treeout",      "--treeout",        sp.typ.FLAG )
ps.add_opt( "quiet",        "--quiet",          sp.typ.FLAG )

# Input
ps.add_opt( "nuc",          "--nuc",            sp.typ.FLAG )
ps.add_opt( "amino",        "--amino",          sp.typ.FLAG )
ps.add_opt( "seed",         "--seed",           sp.typ.FLAG )

# Runtime / Other (undocumented?)
ps.add_opt( "dash",         "--dash",           sp.typ.FLAG )
ps.add_threads("--thread {}")

# Input/Output (they're positional)
ps.add( snakemake.input[0], "{}",   sp.typ.FILE )
ps.add( ">", "{}" )
ps.add( snakemake.output[0], "{}" )

# =================================================================================================
#     Run
# =================================================================================================



shell( ps.get_shell_string( log_stdout=False ) )
