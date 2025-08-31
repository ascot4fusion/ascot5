"""Interface for accessing data in ASCOT5 HDF5 files."""
from typing import Dict, List, Any, Optional

from .access import Tree, InputVariant, RunVariant
from .access.treeparts import Leaf
from .access.metadata import MetaData
from . import (bfield, efield, plasma, neutral, wall, marker, mhd, boozer, nbi,
               asigma, options)

#from a5py.routines.runmixin import RunMixin
#from a5py.routines.afsi5 import AfsiMixin
#from a5py.routines.bbnbi5 import BBNBIMixin


class AscotData(
    Tree, bfield.CreateBfieldMixin, efield.CreateEfieldMixin,
    plasma.CreatePlasmaMixin, neutral.CreateNeutralMixin,
    wall.CreateWallMixin, marker.CreateMarkerMixin,
    mhd.CreateMhdMixin, boozer.CreateBoozerMixin,
    nbi.CreateNbiMixin, asigma.CreateAsigmaLocMixin,
    options.CreateOptionsMixin,
    ):
    """A class through which all data is accessed and managed.

    This class forms a tree-like structure for accessing different data
    containers. Creating new data sets is managed via this class to make sure
    that the tree-like structure stays functional.
    """

    @classmethod
    def _leaf_factory(
        cls,
        meta: MetaData,
        inputs: Optional[Dict[str, InputVariant]] = None,
        diagnostics: Optional[List] = None,
        ) -> Leaf:
        """Create `Leaf` instances of different variety.

        Override this method to create specific `Leaf` instances that are stored
        in this tree. This function is automatically called when using the
        `_add_input_dataset` and `_add_simulation_output` methods. This base
        implementation simply creates a generic `Leaf` instance.
        """
        if(meta.variant in ["run"]
            and (inputs is not None and diagnostics is not None) ):
            return RunVariant(
                **meta._asdict(), inputs=inputs, diagnostics=diagnostics,
                #**kwargs,
            )
        for input in [bfield, efield, plasma, neutral, wall, mhd, boozer,
                      asigma, nbi, marker, options]:
            if hasattr(input, meta.variant):
                Variant = getattr(input, meta.variant)
                break
        else:
            raise ValueError(
                f"Unknown variant {meta.variant} for {meta.qid}"
            )

        return Variant(meta.qid, meta.date, meta.note)

    def _create_inputgroup(self, path, h5):
        """Create an input group to be added to the treeview.

        Parameters
        ----------
        path : str
            Path to the input group within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.

        Returns
        -------
        group : `InputGroup`
            The input group that was created
        """
        return InputGroup(self, path, h5)

    def _create_resultgroup(self, path, h5, inputqids):
        """Create a result group to be added to the treeview.

        Parameters
        ----------
        path : str
            Path to the result node within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.
        inputqids : dict [str, str]
            Dictionary containing the name of the input parent group
            (e.g. "bfield") and the QID of the input used for this result.

        Returns
        -------
        group : `ResultGroup`
            The result group that was created.
        """
        if path.split("/")[-1].split("_")[0] == "run":
            return RunGroup(self, path, h5, inputqids)
        elif path.split("/")[-1].split("_")[0] == "afsi":
            return AfsiGroup(self, path, h5, inputqids)
        elif path.split("/")[-1].split("_")[0] == "bbnbi":
            return BBNBIGroup(self, path, h5, inputqids)
        else:
            raise ValueError("Unknown")


class InputGroup(InputVariant):
    """Node containing input data groups.
    """


#class RunGroup(RunMixin):
    """Node containing results and methods to process them.
    """


#class AfsiGroup(AfsiMixin):
    """Node containing AFSI results and methods to process them.
    """


#class BBNBIGroup(BBNBIMixin):
    """Node containing BBNBI results and methods to process them.
    """
