from argparse import OPTIONAL
from email.errors import ObsoleteHeaderDefect
import util
from enum import Enum
from dataclasses import dataclass

class req(Enum):
    REQUIRED = 1
    OPTIONAL = 2

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

    def add(self, arg: str, format_string: str = "{}", validation_type: typ = typ.NONE ):
        
        # Validate input
        if validation_type is typ.FILE:
            util.expect_file_exists( arg )
        elif validation_type is typ.EXECUTABLE:
            util.expect_executable_exists( arg )
        elif validation_type is typ.DIRECTORY:
            util.expect_dir_exists( arg )

        # add to the shell string
        if not validation_type is typ.FLAG:
            format_string = " " + format_string.format( arg )
        elif not arg:
            format_string = ""

        self.shell_string = self.shell_string + format_string

        return format_string
