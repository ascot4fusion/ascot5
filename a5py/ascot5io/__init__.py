"""Interface for accessing data in ASCOT5 HDF5 files.
"""
from typing import Any

from a5py.ascot5io.coreio.treestructure import Leaf
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

from .coreio.treestructure import Root
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

class Ascot5IO(Root):
    """Tree structure for accessing ASCOT5 data.

    This class is an interface to ASCOT5 data wherever it is stored (in HDF5
    file, IMAS IDS, memory). The data is organized into a tree where the top
    level has one node for each simulation category (bfield, efield, etc.) and
    groups for each simulation run. Each input category holds all inputs of that
    kind in separate groups. Multiple simulations and multiple inputs of same
    category are supported.

    Groups can be accessed as

    .. code-block:: python

       data.bfield.B_2DS_1234567890

    or, equivalently,

    .. code-block:: python

       data["bfield"]["B_2DS_1234567890"]

    Groups can also be accessed through their QID (quasi-unique identifier) or
    user specified tag (first word in the description, capitalized, without
    special characters, and starting with a letter).

    .. code-block:: python

       # These are equivalent references
       data.bfield.B_2DS_1234567890
       data.bfield.q1234567890
       data.bfield.TAG

    In each input category, one input is always set as "active" indicating that
    it will be used in the next simulation or in post-processing. The active
    group is accessed with

    .. code-block:: python

       data.bfield.active

    Simulation outputs are located on the top level and they are accessed in
    identical manner, including that one run is always set as active (by default
    the most recent simulation):

    .. code-block:: python

       data.run_1234567890
       data.active

    To quickly show the contents of the tree, input category, or simulation
    output, use the `show_contents` methods:

    .. code-block:: python

       data.show_contents()
       data.bfield.show_contents()
       data.run_1234567890.show_contents()

    """

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

    def create_BTC(
            self,
            bxyz,
            jacobian,
            rhoval,
            description=None,
            activate=None,
            dryrun=False,
            store_hdf5=True,
            ):
        """Create a magnetic field input which is defined on a Cartesian basis.

        This input represents a magnetic field input, where the field vector is
        defined at the origo on a Cartesian basis and the Jacobian is constant.

        The purpose of this field is to validate the orbit-integrators.

        Parameters
        ----------
        bxyz : array_like (3,1)
            Magnetic field in cartesian coordinates at origo.
        jacobian : array_like (3,3)
            Magnetic field Jacobian, jacobian[i,j] = dB_i/dx_j.
        rhoval: float
            Constant rho value.
        psival: float, optional
            Constant psi value.

            If None, same as rhoval.
        axisr: float, optional
            Magnetic axis R coordinate.
        axisz: real, optional
            Magnetic axis z coordinate.
        description : str, optional
            User-defined description of the data.

            Use this to document the contents of this input. The first word of
            the description is capitalized and used as a "tag" to reference the
            input as a5.data.bfield.MYTAG.
        activate : bool, optional
            Set this input as active on creation.
        dryrun : bool, optional
            Do not add this input to the data structure or store it on disk.

            Use this flag to modify the input manually after it has been
            created, but before it is stored.
        store_hdf5 : bool, optional
            Write this input to the HDF5 file if one has been specified when
            `Ascot` was initialized.

        Returns
        -------
        """
        obj = B_TC(bxyz, jacobian, rhoval)
        self._add_input(
            obj, parent="bfield", dryrun=dryrun, description=description,
            store_hdf5=store_hdf5
            )
        return obj

class RunGroup(RunMixin):
    """Node containing results and methods to process them.
    """
    pass

class AfsiGroup(AfsiMixin):
    """Node containing AFSI results and methods to process them.
    """
    pass

class BBNBIGroup(BBNBIMixin):
    """Node containing BBNBI results and methods to process them.
    """
    pass
