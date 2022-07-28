# =================================================================================================
#     Dependencies and Setup
# =================================================================================================

from snakemake.shell import shell
sys.path.insert(0, '../common')
import snakeparser as sp
import os

shell.executable("bash")

# Get the output directory
sample_outdir = os.path.dirname( snakemake.output[0] )

# =================================================================================================
#     Parse arguments
# =================================================================================================
ps = sp.Parser( "gappa examine assign", snakemake )

# Required args
ps.add( snakemake.input.jplace,     "--jplace-path {}", sp.typ.FILE )
ps.add( snakemake.input.taxon_file, "--taxon-file {}",  sp.typ.FILE )

# Optional args
ps.add_opt( "root_outgroup",        "--root-outgroup {}",           sp.typ.FILE )
ps.add_opt( "taxonomy",             "--taxonomy {}",                sp.typ.FILE )
ps.add_opt( "ranks_string",         "--ranks-string {}" )
ps.add_opt( "sub_taxopath",         "--sub-taxopath {}" )
ps.add_opt( "max_level",            "--max-level {}",               sp.typ.UINT )
ps.add_opt( "distribution_ratio",   "--distribution-ratio {}",
    sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "consensus_thresh",     "--consensus-thresh {}",
    sp.typ.FLOAT(0.0,1.0) )
ps.add_opt( "resolve_missing_paths","--resolve-missing-paths",   sp.typ.FLAG )
ps.add_opt( "distant_label",        "--distant-label",           sp.typ.FLAG )

# Output args
ps.add_opt( "file_prefix",              "--file-prefix {}" )
ps.add_opt( "file_suffix",              "--file-suffix {}" )
ps.add_opt( "cami",                     "--cami",                    sp.typ.FLAG )
ps.add_opt( "krona",                    "--krona",                   sp.typ.FLAG )
ps.add_opt( "sativa",                   "--sativa",                  sp.typ.FLAG )
ps.add_opt( "sample_id",                "--sample-id {}" )
ps.add_opt( "per_query_results",        "--per-query-results",       sp.typ.FLAG )
ps.add_opt( "best_hit",                 "--best-hit",                sp.typ.FLAG )
ps.add_opt( "allow_file_overwriting",   "--allow-file-overwriting",  sp.typ.FLAG )
ps.add_opt( "verbose",                  "--verbose",                 sp.typ.FLAG )

# Closing args
ps.add( sample_outdir, "--out-dir {}" , sp.typ.DIR )
ps.add_threads()

# =================================================================================================
#     Run
# =================================================================================================

shell( ps.get_shell_string() )
