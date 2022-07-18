import util
from enum import Enum
from dataclasses import dataclass

class typ(Enum):
    NONE  = 0
    FILE  = 1
    EXEC  = 2
    DIR   = 3
    FLAG  = 4

@dataclass
class Parser:
    """Class for progressively building a validated shell command."""
    shell_string: str
    snakemake: object

    def add( self, arg: str, format_string: str = "{}", validation_type: typ = typ.NONE ):
        
        # Validate input
        if validation_type is typ.FILE:
            util.expect_file_exists( arg )
        elif validation_type is typ.EXEC:
            util.expect_executable_exists( arg )
        elif validation_type is typ.DIR:
            util.expect_dir_exists( arg )

        # add to the shell string
        if not validation_type is typ.FLAG:
            format_string = " " + format_string.format( arg )
        elif not arg:
            format_string = ""

        self.shell_string = self.shell_string + format_string

        return format_string
    
    def add_opt( self, arg: str, format_string: str = "{}", validation_type: typ = typ.NONE ):
        if arg in self.snakemake.params:
            return self.add( self.snakemake.params[arg], format_string, validation_type )

