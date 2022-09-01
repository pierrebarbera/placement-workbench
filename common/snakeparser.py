from os import access
import util
from collections.abc import Callable
from sys import maxsize
from operator import getitem
from functools import reduce

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

class Parser:
    """Class for progressively building a validated shell command.
    Basically a fancy StringBuilder class with some input validation.
    """
    _shell_string: str
    _snakemake: object
    # list of accessors (lists of progressively deeper entries into config). In descending order
    # of preference (i.e., first try snakemake["config"]["params"][tool][key], then more if specified)
    # note that the params set in the calling rule always takes precedence
    _accessors: list[list] = []

    # custom constructor, to clarify inputs and to give direct control over accessor behaviour
    def __init__(   self,
                    tool: str,
                    snakemake: object,
                    accessors: list[list] = None ):
        self._shell_string = tool
        self._snakemake = snakemake
        if not accessors and tool in snakemake.config['params'].keys():
            self._accessors = [['params',tool]] + self._accessors
        if accessors:
            if type(accessors) is not list[list]: accessors = [ accessors ]
            self._accessors = accessors + self._accessors

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
        self._shell_string = self._shell_string + format_string

        # return the format string incase we want to check it
        return format_string
    
    def add_opt( self, key: str, format_string: str = "{}", valid_func: Callable[[str], bool] = typ.NONE ):
        # give absolute priority to the rule params, as they may have critical settings
        if key in self._snakemake.params.keys():
            return self.add( self._snakemake.params[key], format_string, valid_func )

        # if the key isn't set in rule/params, look through the config accessors
        # and return the first hit
        for accessor in self._accessors:
            entry = reduce( getitem, accessor, self._snakemake.config )
            if key in entry.keys():
                return self.add( entry[key], format_string, valid_func )

    def add_threads( self, format_string: str = "--threads {}", valid_func: Callable[[str], bool] = typ.UINT ):
        return self.add( self._snakemake.threads, format_string, valid_func )

    def get_shell_string( self, do_log: bool = True, log_stdout: bool = True ):
        if do_log:
            return self._shell_string + " {}".format( self._snakemake.log_fmt_shell(stdout=log_stdout, stderr=True) )
        else:
            return self._shell_string 
