import ctypes

import unyt
import numpy as np

from a5py.exceptions import AscotIOException

from a5py.data.access import Leaf, InputLeaf, Tree, OutputLeaf
from a5py.data.access.dataholder import DataStruct, Status
from a5py.data.access import variants

input_prototype = "inputprototype"

output_prototype = "runprototype"


def free(ptr):
    libc = ctypes.CDLL(None)
    free = libc.free
    free.argtypes = [ctypes.c_void_p]
    free(ptr)


def malloc(n, ctype):
    libc = ctypes.CDLL(None)
    malloc = libc.malloc
    malloc.restype = ctypes.c_void_p
    malloc.argtypes = [ctypes.c_size_t]
    arr = malloc(n * ctypes.sizeof(ctype))
    return ctypes.cast(arr, ctypes.POINTER(ctype))


class StructInput(DataStruct):
    """Data stored in C."""
    _fields_ = [
        ("n", ctypes.c_size_t),
        ("bval", ctypes.POINTER(ctypes.c_double)),
        ]


class StructDiag(DataStruct):
    """Data stored in C."""
    _fields_ = [
        ("n", ctypes.c_size_t),
        ("dval", ctypes.POINTER(ctypes.c_double)),
        ]


@Leaf.register(input_prototype)
class InputPrototype(InputLeaf):
    """A prototype class for testing purposes."""

    @property
    def npoint(self):
        """Number of points."""
        if self.status is Status.SAVED:
            return self._file.read("bval").size
        elif self.status & Status.STAGED:
            return self._cdata.n

    @property
    def bval(self):
        """Number of points."""
        if self.status is Status.SAVED:
            return self._file.read("bval")
        elif self.status & Status.STAGED:
            return self._cdata.readonly_carray("bval", (self.npoint,), "T")

    def save(self):
        super().save()
        self._file.write("bval", self.bval)

    def export(self):
        return {"bval": self.bval}

    def stage(self):
        super().stage()
        self._stage(bval=self.bval)

    def unstage(self):
        super().unstage()
        free(self._cdata.bval)
        self._cdata = None

    def _stage(self, bval):
        n = bval.size
        cdata = StructInput()
        cdata.n = n
        cdata.bval = malloc(n, ctypes.c_double)
        for i in range(n):
            cdata.bval[i] = bval[i]

        self._cdata = cdata


class PrototypeMixin():

    #pylint: disable=protected-access
    def create_inputprototype(self, bval, note=None, activate=False, dryrun=False, save=True) -> InputPrototype:
        """Creates a point cloud instance.
        """
        parameters = variants.parse_parameters(bval)
        n = 3 if parameters["bval"] is None else parameters["bval"].size
        variants.validate_required_parameters(
            parameters, names=["bval",], units=["T",], shape=[(n,),], dtype="f8",
            default=[np.array([0., 1., 2.]),],
        )
        #variants.validate_optional_parameters(
        #    parameters,
        #    ["psival"], ["Wb/m"], [()], "f8", [parameters["rhoval"]],
        #)
        meta = variants.new_metadata(input_prototype, note=note)
        leaf = InputPrototype(*meta)
        leaf._stage(bval=parameters["bval"])
        self._treemanager.enter_input(
            leaf, activate=activate, dryrun=dryrun, save=save, category="group",
            )
        return leaf


@Leaf.register(output_prototype)
class OutputPrototype(OutputLeaf):
    """A prototype class for testing purposes."""

    def __init__(self, qid, date, note, variant, inputs):
        super().__init__(qid=qid, date=date, note=note, variant=variant, inputs=inputs)

    def _setup(self, params):
        if params["diag1"]:
            self._diagnostics["diag1"] = (
                DiagnosticPrototype.from_params(params["dval"])
                )

    def _load(self, file):
        super()._load(file)
        for name in self._file.get_children():
            if name == "diag1":
                self._diagnostics[name] = DiagnosticPrototype()

    def get_diag1(self):
        """Get values from 'diag1'.

        This is just a test method to mock real methods such as 'getstate'.
        """
        if self._file is None:
            return self._diagnostics["diag1"].get()
        return self._diagnostics["diag1"].get(self._file.get_minimanager("diag1"))

    def save(self):
        super().save()
        for name, diag in self._diagnostics.items():
            if name == "diag1":
                file = self._file.get_minimanager("diag1")
                diag.save(file)


class DiagnosticPrototype():
    """A prototype class for testing purposes."""

    def __init__(self):
        super().__init__()
        self.name = "DataArray"
        self._cdata: StructInput | None = None

    def save(self, file):
        """"""
        file.write("dval", self.get())
        free(self._cdata.dval)
        del self._cdata

    def combine(self, status):
        """"""

    def get(self, file=None):
        """"""
        if file is None:
            n = self._cdata.n
            return self._cdata.readonly_carray("dval", (n,), "T")
        return file.read("dval")

    @classmethod
    def from_params(cls, dval):
        """"""
        obj = cls()

        # Normally the diagnostic is initialized to clean state but here we set
        # values to mimic the state at the end of the simulation.
        n = dval.size

        cdata = StructDiag()
        cdata.n = n
        cdata.dval = malloc(n, ctypes.c_double)
        for i in range(n):
            cdata.dval[i] = dval[i]

        obj._cdata = cdata
        return obj



class TreePrototype(Tree, PrototypeMixin):
    """A prototype class for testing purposes."""
