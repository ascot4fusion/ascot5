"""
Python interface to ascot library.

This module defines Ascotpy which is a class that inherits all other
classes found in ascotpy package. Thus having an instance of Ascotpy is
enough to gain access to all available tools.

File: ascotpy.py
"""
from a5py.ascotpy.libbfield import LibBfield

class Ascotpy(LibBfield):
    """
    One class to rule them all.
    """
    pass
