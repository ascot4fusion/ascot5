"""Modifies ascot2py.py so that it first try to find libascot.so using relative
path (if installing from the source) and only then from LD_LIBRARY_PATH.
"""
import fileinput
import sys

for line in fileinput.input("src/ascot2py.py", inplace=True):
    if line.strip() == "_libraries['libascot.so'] = ctypes.CDLL('libascot.so')":
        sys.stdout.write(
            "# Try to locate libascot.so from ../../build/ or LD_LIBRARY_PATH\n"
            "from pathlib import Path\n"
            "err = 0\n"
            "libpath = str(Path(__file__).absolute().parent.parent.parent) \\\n"
            "    + \"/build/libascot.so\"\n"
            "try:\n"
            "    _libraries['libascot.so'] = ctypes.CDLL(libpath)\n"
            "except OSError as error:\n"
            "    err = error\n"
            "if err:\n"
            "    if not 'libascot.so' in str(err): raise ImportError(str(err))\n"
            "if 'libascot.so' not in _libraries:\n"
            "    try:\n"
            "        _libraries['libascot.so'] = ctypes.CDLL('libascot.so')\n"
            "    except OSError as error:\n"
            "        err = error\n"
            "    if err: raise ImportError(str(err))\n"
            "\n"
            )
    else:
        sys.stdout.write(line)
