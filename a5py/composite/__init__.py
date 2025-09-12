"""Collection of classes that the acts as the user interface to ASCOT.

Since these classes act as the user interface, the focus here is on the
documentation while the actual implementation is usually done elsewhere. Another
motivation is to avoid repetition, because one class usually has methods that
are very similar to the methods of other classes (e.g. both AscotRun and
BeamRun have getstate(), but the arguments and where the data is read differ).

So try to keep the actual logic here as simple as possible and focus on the
interface. Logic that is tightly knitted with libascot.so should be in engine,
and higher level logic should be in the routines.
"""
from .ascot import Ascot
from .run import Run
from .bbnbi import BeamRun
from .afsi import AfsiRun

__all__ = [
    "Ascot",
    "Run",
    "BeamRun",
    "AfsiRun",
    ]