"""Templates for... templates!"""
from a5py import Ascot

class InputTemplate():
    """Template for input templates.

    The data is read and processed when this object is initialized. The methods
    provide access to the data, it's type, and generates the actual input.

    The motivation for this class is that this way the user can review the data,
    or even modify it, before the actual input is created (which is immutable).
    On the other hand, some templates have the option to return different kinds
    of inputs, and the user can't know the input type beforehand.

    For these reasons we separate the data processing and input creation, and
    put the relevant factory method(s) to a single anonymous method, so that the
    user don't have to bother about the input type.
    """


    def __init__(self, ascot: Ascot, input_type: str, data: dict):
        self._ascot: Ascot = ascot
        self._type: str = input_type
        self.data: dict = data

    @property
    def input_type(self) -> str:
        return self._type

    def create_input(self) -> object:
        """Uses the Master to create the correct input object."""
        method_name = f"create_{self.input_type}"
        factory_method = getattr(self._ascot.data, method_name, None)
        if factory_method is None:
            raise ValueError(f"Unsupported input type: {self.input_type}")
        return factory_method(**self.data)
