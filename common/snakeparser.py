import util
from dataclasses import dataclass
from collections.abc import Callable
from sys import maxsize

class typ:
    @staticmethod
    def NONE( arg ):
        return True
    @staticmethod
    def FILE( arg ):
        util.expect_file_exists( arg )
    @staticmethod
    def EXEC( arg ):
        util.expect_executable_exists( arg )
    @staticmethod
    def DIR( arg ):
        util.expect_dir_exists( arg )
    @staticmethod
    def FLAG( arg ):
        if not arg in ['False', 'True']:
            util.fail("expected flag (True or False), but got '{}'".format(arg))
    @staticmethod
    def IN( set: list[str] ):
        def func( arg: str ):
            if not arg in set:
                util.fail( "{} not in ".format( arg, set ) )
        return func
    @staticmethod
    def FLOAT( lower: float = 0.0, upper: float = float("inf") ):
        """
        Returns validation function that checks if arg is within [lower,upper]
        """
        def func( arg ):
            f = float(arg)
            if not ( (f >= lower) and (f <= upper) ):
                util.fail( "{} not within [{},{}] ".format( arg, lower, upper ) )
        return func
    @staticmethod
    def UINT( lower: int = 0, upper: int = maxsize*2+1 ):
        """
        Returns validation function that checks if arg is within [lower,upper]
        """
        def func( arg ):
            f = int(arg)
            if not ( (f >= lower) and (f <= upper) ):
                util.fail( "{} not within [{},{}] ".format( arg, lower, upper ) )
        return func

@dataclass
class Parser:
    """Class for progressively building a validated shell command.
    Basically a fancy StringBuilder class with some input validation.
    """
    shell_string: str
    snakemake: object

    def add( self, arg: str, format_string: str = "{}", valid_func: Callable[[str], bool] = typ.NONE ):
        
        # Validate input
        valid_func( arg )

        # add to the shell string
        if not valid_func is typ.FLAG:
            format_string = " " + format_string.format( arg )
        elif arg == 'False':
            format_string = ""

        self.shell_string = self.shell_string + format_string

        return format_string
    
    def add_opt( self, arg: str, format_string: str = "{}", valid_func: Callable[[str], bool] = typ.NONE ):
        if arg in self.snakemake.params:
            return self.add( self.snakemake.params[arg], format_string, valid_func )

