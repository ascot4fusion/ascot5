"""Interface for accessing data in ASCOT5 HDF5 files.
"""
from .bfield  import B_TC, B_GS, B_2DS, B_3DS, B_3DST, B_STS
from .efield  import E_TC, E_1DS, E_3D, E_3DS, E_3DST
from .marker  import Marker, Prt, GC, FL
from .plasma  import plasma_1D, plasma_1DS, plasma_1Dt
from .wall    import wall_2D, wall_3D
from .neutral import N0_1D, N0_3D
from .boozer  import Boozer
from .mhd     import MHD_STAT, MHD_NONSTAT
from .asigma  import Asigma_loc
from .options import Opt
from .nbi     import NBI

from .state      import State
from .orbits     import Orbits
from .transcoef  import Transcoef
from .dist import Dist_5D, Dist_6D, Dist_rho5D, Dist_rho6D, Dist_COM, Dist
from .reaction import Reaction

from .coreio.fileapi import INPUTGROUPS
from .coreio.treeview import RootNode, InputNode, ResultNode
from .coreio.treedata import DataGroup
from a5py.routines.runmixin import RunMixin
from a5py.routines.afsi5 import AfsiMixin
from a5py.routines.bbnbi5 import BBNBIMixin
from a5py.templates import Template

HDF5TOOBJ = {
    "B_TC" : B_TC, "B_GS" : B_GS, "B_2DS" : B_2DS, "B_3DS" : B_3DS,
    "B_3DST" : B_3DST, "B_STS" : B_STS,
    "E_TC" : E_TC, "E_1DS" : E_1DS,
    #"E_3D" : E_3D, "E_3DS" : E_3DS, "E_3DST" : E_3DST,
    "prt" : Prt, "gc" : GC, "fl" : FL,
    "wall_2D" : wall_2D, "wall_3D" : wall_3D,
    "plasma_1D" : plasma_1D, "plasma_1DS" : plasma_1DS,
    "plasma_1Dt" : plasma_1Dt,
    "N0_1D" : N0_1D, "N0_3D" : N0_3D,
    "Boozer" : Boozer, "MHD_STAT" : MHD_STAT, "MHD_NONSTAT" : MHD_NONSTAT,
    "asigma_loc" : Asigma_loc,
    "opt" : Opt,
    "nbi" : NBI,
    "inistate" : State, "endstate" : State, "state" : State,
    "orbit" : Orbits,
    "reaction" : Reaction,
    "dist5d" : Dist_5D, "prod1dist5d" : Dist_5D, "prod2dist5d" : Dist_5D,
    "dist6d" : Dist_6D,
    "distrho5d" : Dist_rho5D,
    "distrho6d" : Dist_rho6D,
    "distcom" : Dist_COM,
    "transcoef": Transcoef
}
"""Dictionary connecting group names in HDF5 file to corresponding data objects.
"""

class Ascot5IO(RootNode):
    """Entry node for accessing data in the HDF5 file.

    Initializing this node builds rest of the tree. This object and its child
    nodes act as container objects for the data in the HDF5 file. At top level
    (this node) run nodes containing simulation results can be accessed as well
    as the parent nodes (bfield, efield, etc.) that in turn contain the actual
    input groups.

    The data can be accessed as

    .. code-block:: python

       data.bfield.B_2DS_1234567890

    or, equivalently,

    .. code-block:: python

       data["bfield"]["B_2DS_1234567890"]

    In each input group, one input is always set as "active" (meaning it would
    be used for the next simulation) and it can be accessed as

    .. code-block:: python

       data.bfield.active

    QID can be used as a reference as well

    .. code-block:: python

       data.bfield.q1234567890

    Run groups are accessed in a similar fashion, e.g.

    .. code-block:: python

       data.run_1234567890 .

    However, the easiest way to access the simulation output is via the methods
    in the RunNode. The active run (by default the most recent simulation) can
    be accessed with

    .. code-block:: python

       data.active

    and its inputs as

    .. code-block:: python

       data.active.bfield

    Most of these examples also work with dictionary-like reference but here we
    use only the attribute-like referencing for brevity.

    You can even use user defined tag taken from description to refer to it

    .. code-block:: python

       data.THATPRLRUN

    However, there are few rules to this:

    - Tag is the first word in description converted to all caps for brevity.
    - Maximum length is ten characters.
    - Only letters and numbers are allowed (all special characters are removed).
    - First character must be a letter.
    - If two fields would have identical tags, they are changed to format
      <tag>_<i>, where i is running index starting from zero (corresponding
      to most recent data).

    Finally, you can print the contents of a node with

    .. code-block:: python

       data.ls()
    """

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

    def _create_datagroup(self, grouptype, path):
        """Create data group based on the given type.

        Parameters
        ----------
        grouptype : str
            Type of the group as it appears in `HDF5TOOBJ`.
        path : str
            Path to the data in the HDF5 file that corresponds to the group.

        Returns
        -------
        `DataContainer`
            The data group that was created.
        """
        return HDF5TOOBJ[grouptype](self, path)

    def create_input(self, inp, desc=None, activate=None, dryrun=False,
                     **kwargs):
        """Create input and write the data to the HDF5 file.

        This method can be used in two ways to write input data.

        1. From :mod:`ascot5io` choose the input type and find what is requested
           by the corresponding ``write_hdf5`` function to write inputs
           explicitly.

        2. From :mod:`template` choose a template or import, and find what
           parameters are requested to write premade inputs or imported data.

        Parameters
        ----------
        inp : str
            Type of the input e.g. "B_2DS" or name of the template/import
            e.g. "analytical ITER-like".
        desc : str, optional
            Input description.
        activate : bool, optional
            Set created input as active.
        dryrun : bool, optional
            If True, the data is not written to HDF5 but instead it is returned
            as a dictionary.

            Only works for templates.
        **kwargs : dict
            The parameters that are needed to create the requested input type
            or template.

            If input type is given but kwargs is empty, dummy input of given
            type is created.

        Returns
        -------
        data : str or dict
            Name, i.e. "<type>_<qid>", of the new input that was written or the
            data that was created but not written if ``dryrun`` is True.
        """
        if inp in HDF5TOOBJ.keys():
            if len(kwargs) == 0:
                data = HDF5TOOBJ[inp].create_dummy()
                if dryrun:
                    return data
                else:
                    name = HDF5TOOBJ[inp].write_hdf5(
                        self._ascot.file_getpath(), **data)
            else:
                name = HDF5TOOBJ[inp].write_hdf5(
                    self._ascot.file_getpath(), **kwargs)
        else:
            gtype, data = Template(self._ascot).usetemplate(inp, **kwargs)
            if dryrun:
                return data
            else:
                name = HDF5TOOBJ[gtype].write_hdf5(
                    self._ascot.file_getpath(), **data)

        self._build(self._ascot.file_getpath())
        if activate:
            self._activate_group(name)
        if desc is not None:
            self._get_group(name).set_desc(desc)
        return name

    def templates(self, template):
        """Show information about a template.

        Parameters
        ----------
        template : str
            Name of the template.

            If not given, all templates are listed.
        """
        Template(self._ascot).showtemplate(template)

class RunGroup(ResultNode, RunMixin):
    """Node containing results and methods to process them.
    """
    pass

class AfsiGroup(ResultNode, AfsiMixin):
    """Node containing AFSI results and methods to process them.
    """
    pass

class BBNBIGroup(ResultNode, BBNBIMixin):
    """Node containing BBNBI results and methods to process them.
    """
    pass
