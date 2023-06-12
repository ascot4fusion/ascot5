"""Interface for accessing data in Ascot5 HDF5 files.
"""
from a5py.ascot5io.B_TC       import B_TC
from a5py.ascot5io.B_GS       import B_GS
from a5py.ascot5io.B_2DS      import B_2DS
from a5py.ascot5io.B_3DS      import B_3DS
from a5py.ascot5io.B_3DST     import B_3DST
from a5py.ascot5io.B_STS      import B_STS
from a5py.ascot5io.E_TC       import E_TC
from a5py.ascot5io.E_1DS      import E_1DS
from a5py.ascot5io.E_3D       import E_3D
from a5py.ascot5io.E_3DS      import E_3DS
from a5py.ascot5io.E_3DST     import E_3DST
from a5py.ascot5io.mrk_prt    import mrk_prt
from a5py.ascot5io.mrk_prt_shined    import mrk_prt_shined
from a5py.ascot5io.mrk_gc     import mrk_gc
from a5py.ascot5io.mrk_fl     import mrk_fl
from a5py.ascot5io.wall_2D    import wall_2D
from a5py.ascot5io.wall_3D    import wall_3D
from a5py.ascot5io.plasma_1D  import plasma_1D
from a5py.ascot5io.plasma_1DS import plasma_1DS
from a5py.ascot5io.N0_3D      import N0_3D
from a5py.ascot5io.boozer     import Boozer
from a5py.ascot5io.mhd        import MHD
from a5py.ascot5io.options    import Opt
from a5py.ascot5io.nbi        import nbi

from a5py.ascot5io.state      import State
from a5py.ascot5io.orbits     import Orbits
from a5py.ascot5io.transcoef  import Transcoef
from a5py.ascot5io.dist_5D    import Dist_5D
from a5py.ascot5io.dist_6D    import Dist_6D
from a5py.ascot5io.dist_rho5D import Dist_rho5D
from a5py.ascot5io.dist_rho6D import Dist_rho6D

from ._iohelpers.fileapi import INPUTGROUPS
from ._iohelpers.treeview import RootNode, InputNode, ResultNode
from ._iohelpers.treedata import DataGroup
from a5py.routines.runmixin import RunMixin
import a5py.premade as premade

HDF5TOOBJ = {
    "B_TC" : B_TC, "B_GS" : B_GS, "B_2DS" : B_2DS, "B_3DS" : B_3DS,
    "B_3DST" : B_3DST, "B_STS" : B_STS,
    "E_TC" : E_TC, "E_1DS" : E_1DS, "E_3D" : E_3D, "E_3DS" : E_3DS,
    "E_3DST" : E_3DST,
    "prt" : mrk_prt, "prt_shined" : mrk_prt_shined, "gc" : mrk_gc,
    "fl" : mrk_fl,
    "wall_2D" : wall_2D, "wall_3D" : wall_3D,
    "plasma_1D" : plasma_1D, "plasma_1DS" : plasma_1DS,
    "N0_3D" : N0_3D,
    "Boozer" : Boozer, "MHD_STAT" : MHD, "MHD_NONSTAT" : MHD,
    "opt" : Opt,
    "nbi" : nbi,
    "inistate" : State,
    "endstate" : State,
    "orbit" : Orbits,
    "dist5d" : Dist_5D,
    "dist6d" : Dist_6D,
    "distrho5d" : Dist_rho5D,
    "distrho6d" : Dist_rho6D,
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
        path : `str`
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
        path : `str`
            Path to the result node within the HDF5 file.
        h5 : `h5py.File`
            The HDF5 file from which the tree is constructed.
        inputqids : `dict` [`str`, `str`]
            Dictionary containing the name of the input parent group
            (e.g. "bfield") and the QID of the input used for this result.

        Returns
        -------
        group : `ResultGroup`
            The result group that was created.
        """
        return RunGroup(self, path, h5, inputqids)

    def _create_datagroup(self, grouptype, path):
        """Create data group based on the given type.

        Parameters
        ----------
        grouptype : `str`
            Type of the group as it appears in `HDF5TOOBJ`.
        path : `str`
            Path to the data in the HDF5 file that corresponds to the group.

        Returns
        -------
        `DataContainer`
            The data group that was created.
        """
        return HDF5TOOBJ[grouptype](self, path)

    def create_input(self, inputtype, inputdata):
        """Create input and write the data to the HDF5 file.

        Parameters
        ----------
        inputtype : `str`
            Type of the input e.g. "B_2DS" or "options".
        inputdata : `dict`
            Dictionary containing all the data that is needed to create the
            requested input type.

            If None, dummy input of given type is created.

        Returns
        -------
        qid : `str`
            QID of the created input.
        """
        if inputdata is None:
            qid = HDF5TOOBJ[inputtype](self, None).write_dummy(
                self._ascot.file_getpath())
        else:
            qid = HDF5TOOBJ[inputtype](self, None).write(
                self._ascot.file_getpath(), data=inputdata)
        self._build(self._ascot.file_getpath())
        return qid

    def create_premade(self, predef, write=True, activate=True, **kwargs):
        """Create inputs based on predefined simulations or existing interfaces
        that import data to Ascot5.

        Predefined simulations are simulations whose inputs have alredy been
        prepared for a given purpose e.g. to create Poincaré plots.

        This method also acts as an interface to import data from external
        sources using existing interfaces e.g. to convert EQDSK file to magnetic
        field input.

        See `a5py.premade` for details on what this method can produce.

        Parameters
        ----------
        predef : `str`
            Type of the premade data.

            Available premades are listed in `a5py.premade`.
        write : `str`
            If `True`, the created inputs are written to the HDF5 file.
        activate : `bool`, optional
            If `True`, the created inputs are set as active.

            This option is ignored if `write` is `False`.
        **kwargs
            Any required or optional parameters that are passed to
            the constructor of the requested premade type.

            See `a5py.premade` for more details.

        Returns
        -------
        qids : `list` [`str`] or `dict` [`str`, `dict`]
            QIDs of the inputs that were created and written or, if `write` is
            `False`, dictionary containing input type and data that can be
            passed to the corresponding `write_hdf5` function.
        """
        simdata = premade.construct(predef, **kwargs)
        if not write:
            return simdata

        qids = []
        for name, data in simdata.items():
            qids.append(self.create_input(name, data))

        if activate:
            for q in qids:
                self._activate_group(q)

        return qids

class InputGroup(InputNode):
    """Node containing input data groups.
    """
    pass

class RunGroup(ResultNode, RunMixin):
    """Node containing results and methods to process them.
    """
    pass
