from snakemake.shell import shell
sys.path.insert(0, '../common')
import util
import os

log = snakemake.log_fmt_shell(stdout=True, stderr=True)

shell.executable("bash")

"""
Symlinks the input to the output, with correct relative paths
"""

rel_input = os.path.relpath( str(snakemake.input[0]), os.path.dirname( str(snakemake.output[0]) ) )

shell( f"ln -s {rel_input} {snakemake.output[0]} {log}" )
