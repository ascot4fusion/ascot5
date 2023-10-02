"""Modifies ascot2py.py so that is has the correct CDLL path.
"""
import fileinput
import sys

for line in fileinput.input("src/ascot2py.py", inplace=True):
    if line.strip() == "_libraries['libascot.so'] = ctypes.CDLL('libascot.so')":
        sys.stdout.write(
            "# Ugly solution to find libascot.so\n"
            "from pathlib import Path\n"
            "libpath = str(Path(__file__).absolute().parent.parent.parent) \\\n"
            "+ \"/build/libascot.so\"\n"
            "_libraries['libascot.so'] = ctypes.CDLL(libpath)\n"
        )
    else:
        sys.stdout.write(line)
