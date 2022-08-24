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
    def FLAG( arg: str ):
        if not isinstance( arg, bool ):
            util.fail("expected flag (True or False), but got '{}'".format(arg))
    @staticmethod
    def IN( set: list ):
        def func( arg ):
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
    do_log: bool = True

    def add( self, arg, format_string: str = "{}", valid_func: Callable[[str], bool] = typ.NONE ):
        
        # Validate input
        valid_func( arg )

        # add to the shell string
        if valid_func is typ.FLAG:
            # differentiate between flags (they don't have arg values attached)
            if arg:
                format_string = " " + format_string
            else:
                format_string = ""
        else:
            # ...and normal CLI arguments
            format_string = " " + format_string.format( arg )

        # add to the complete shell string
        self.shell_string = self.shell_string + format_string

        # return the format string incase we want to check it
        return format_string
    
    def add_opt( self, arg: str, format_string: str = "{}", valid_func: Callable[[str], bool] = typ.NONE ):
        if arg in self.snakemake.params.keys():
            return self.add( self.snakemake.params[arg], format_string, valid_func )

    def add_threads( self, format_string: str = "--threads {}", valid_func: Callable[[str], bool] = typ.UINT ):
        return self.add( self.snakemake.threads, format_string, valid_func )

    def get_shell_string( self ):
        if self.do_log:
            return self.shell_string + " {}".format( self.snakemake.log_fmt_shell(stdout=True, stderr=True) )
        else:
            return self.shell_string 
