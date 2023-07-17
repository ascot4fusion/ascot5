"""Contains definitions of this package's exceptions.
"""

class AscotIOException(Exception):
    """Raised when there is an internal error in `Ascot5IO`.

    This exception should be raised in cases where the issue is most certainly
    a bug in code or inconsistent HDF5 file.
    """
    pass

class AscotNoDataException(Exception):
    """Raised when required input or output is not present.

    This exception should be raised when user tries to query input or results
    that is not present. Other case is when using methods that require the
    presence of certain input or output data (e.g. when plotting orbits but
    orbit data is not present).
    """
    pass

class AscotInitException(Exception):
    """Raised when there is an issue with data initialized by `Ascot`.

    This exception should be raised whenever there is an issue with `libascot`
    or any related data.
    """
    pass
