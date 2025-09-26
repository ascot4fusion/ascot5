"""Contains definitions of this package's exceptions.

.. autoclass:: a5py.exceptions.AscotDataException

.. autoclass:: a5py.exceptions.AscotMeltdownError

"""

class AscotMeltdownError(RuntimeError):
    """Indicates an internal bug in Ascot itself."""
    def __init__(self, message: str):
        super().__init__(
            f"{message}\n\n"
            "This appears to be an internal bug in Ascot. "
            "Please report it at https://github.com/ascot4fusion/ascot5/issues"
        )

class AscotDataException(Exception):
    """Invalid operation related to data.

    This exception should be raised when user tries to do something with data
    that is not valid or if the required data is not available.
    """
    pass

class AscotIOException(Exception):
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

class AscotUnitWarning(UserWarning):
    """Warning raised when quantities are provided without units.
    """
    pass
