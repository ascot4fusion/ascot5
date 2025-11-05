import numpy as np
import unyt

from scipy.interpolate import interpn
from scipy.integrate   import quad

from a5py import Ascot
from a5py.physlib import aeq
from a5py.engine.interpolate import evaluate
from a5py.data.bfield import BfieldAnalytical

from .template import InputTemplate


class TransmuteAnalyticalBfieldToSplines(InputTemplate):


    def __init__(self, ascot: Ascot, field: BfieldAnalytical, nr:int, rlim, nz:int, zlim):
        """Transmute analytical magnetic field to spline-interpolated one.

        The result is either BfieldSpline2D or BfieldSpline3D depending on
        whether the analytical field has ripple or not.

        Parameters
        ----------
        ascot : :class:`Ascot`
            The Ascot object where the input will be created.
        field : BfieldAnalytical
            The analytical field to transmute.
        """
        rgrid = np.linspace(rlim[0], rlim[1], nr)
        zgrid = np.linspace(zlim[0], zlim[1], nz)

        psi, bphi = evaluate(
            rgrid, 0*unyt.deg, zgrid,
            0.*unyt.s, "psi", "bphi", grid=True, bfield=field
        )
        input_type = "bfieldspline3d" if field.nripple else "bfieldspline2d"
        data = {
            "rgrid": rgrid,
            "zgrid": zgrid,
            "psi": np.squeeze(psi),
            "bphi": np.squeeze(bphi),
            "psilimits": field.psilimits,
            "axisrz": field.axisrz,
        }
        super().__init__(ascot, input_type, data)


class ToroidallyAveragedBfield(InputTemplate):

    def convert_B_3DS(**kwargs):
        """Convert :class:`B_3DS` input to `B_2DS` input.

        This function takes a toroidal average of Bphi and sets BR and Bz to
        zero.

        Parameters
        ----------
        **kwargs
            Arguments passed to :meth:`B_3DS.write_hdf5` excluding ``fn`` and
            ``desc``.

        Returns
        -------
        out : dict
            :class:`B_3DS` converted as an input for :meth:`write_hdf5`.
        """
        bphi = np.mean(kwargs["bphi"], axis=1)
        br   = np.mean(kwargs["bphi"], axis=1) * 0
        bz   = np.mean(kwargs["bphi"], axis=1) * 0
        return {
            "rmin":kwargs["b_rmin"], "rmax":kwargs["b_rmax"],
            "nr":kwargs["b_nr"], "zmin":kwargs["b_zmin"],
            "zmax":kwargs["b_zmax"], "nz":kwargs["b_nz"],
            "axisr":kwargs["axisr"], "axisz":kwargs["axisz"],
            "psi0":kwargs["psi0"], "psi1":kwargs["psi1"],
            "psi":kwargs["psi"], "br":br, "bz":bz,
            "bphi":bphi}
