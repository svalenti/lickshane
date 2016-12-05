import util 
import lickshane2Ddef
import lickshane1Ddef
import cosmics

__version__ = "unknown"
try:
    from _version import __version__
except ImportError:
    # We're running in a tree that doesn't have a _version.py, so we don't know what our version is.
    pass
