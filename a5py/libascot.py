"""Imports libascot.so.

The libascot.so should be found via relative path "../build/libascot.so"
(if the code was installed from source) or in ``LD_LIBRARY_PATH`` (if installed
as a package).
"""
import ctypes
from pathlib import Path

LIBASCOT = None
"""The ctypes.CDLL object for libascot.so."""

def _find_libascot():
    global LIBASCOT
    err = 0
    libpath = (str(Path(__file__).absolute().parent.parent)
               + "/build/libascot.so")
    try:
        LIBASCOT = ctypes.CDLL(libpath)
    except OSError as error:
        # Either libascot.so was not in "../build/" or it could not be loaded
        # due to missing dependencies.
        err = error
        # Missing dependencies error would contain the name of the missing
        # library, not libascot.so.
        unresolved_dependencies = not "libascot.so" in str(err)
        # Undefined symbol error should only concern developers, but it can be
        # a major headache.
        undefined_symbol = "undefined symbol" in str(err)

    if err:
        if unresolved_dependencies:
            raise ImportError(str(err))
        elif undefined_symbol:
            raise ImportError("Error in the source: " + str(err))
        else:
            try:
                # This looks for libascot.so in LD_LIBRARY_PATH.
                LIBASCOT = ctypes.CDLL("libascot.so")
            except OSError as error:
                final_err = error
            if final_err:
                raise ImportError(
                    "Failed to load libascot.so: make sure it is compiled "
                    "(if installed from source) or that it is included in "
                    "LD_LIBRARY_PATH."
                )

_find_libascot()
