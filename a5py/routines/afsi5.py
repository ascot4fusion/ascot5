"""AFSI5: Versatile fusion source integrator AFSI for fast ion and neutron
studies in fusion devices
"""
import ctypes
import numpy as np
import unyt
import numpy.ctypeslib as npctypes

from a5py.ascotpy.libascot import _LIBASCOT, STRUCT_HIST, STRUCT_AFSIDATA, \
    STRUCT_AFSITHERMAL, PTR_REAL, AFSI_REACTIONS
from a5py.exceptions import AscotNoDataException
from a5py.routines.distmixin import DistMixin

class Afsi():
    """ASCOT Fusion Source Integrator AFSI.

    ASCOT Fusion Source Integrator (AFSI) is a tool for calculating fusion rates
    and fusion products for arbitrary particle populations. It can be used
    either as a preprocessing tool e.g. to generate source distribution
    of alphas or as a postprocessing tool e.g. to evaluate neutron source
    from beam-thermal fusion.

    AFSI can be used in three ways:

    - thermal: Calculates fusion rates between two Maxwellian populations.
    - beam-thermal: Calculates fusion rates between a Maxwellian and arbitrary
      population e.g. beam ions as obtained from ASCOT5 simulation.
    - beam-beam: Calculates fusion rates between two arbitrary populations.

    Reference:
    https://iopscience.iop.org/article/10.1088/1741-4326/aa92e9

    Attributes
    ----------
    _ascot : :class:`.Ascot`
        Ascot object used to run AFSI.
    """

    def __init__(self, ascot):
        self._ascot = ascot

    def thermal(self, reaction, minr, maxr, nr, minz, maxz, nz,
                minphi=0, maxphi=360, nphi=1, nmc=1000,
                minppara=-1.3e-19, maxppara=1.3e-19, nppara=80,
                minpperp=0, maxpperp=1.3e-19, npperp=40):
        """Calculate thermonuclear fusion between two thermal (Maxwellian)
        species.

        Parameters
        ----------
        reaction : int or str
            Fusion reaction index or name.
        minr : float
            Minimum R value in the output distribution.
        maxr : float
            Maximum R value in the output distribution.
        nr : int
            Number of R bins in the output distribution.
        minz : float
            Minimum z value in the output distribution.
        maxz : float
            Maximum z value in the output distribution.
        nz : int
            Number of z bins in the output distribution.
        minphi : float, optional
            Minimum phi value in the output distribution.
        maxphi : float, optional
            Maximum phi value in the output distribution.
        nphi : int, optional
            Number of phi bins in the output distribution.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        minppara : float, optional
            Minimum ppara value in the output distribution.
        maxppara : float, optional
            Maximum ppara value in the output distribution.
        nppara : int, optional
            Number of pperp bins in the output distribution.
        minpperp : float, optional
            Minimum pperp value in the output distribution.
        maxpperp : float, optional
            Maximum ppara value in the output distribution.
        npperp : int, optional
            Number of pperp bins in the output distribution.

        Returns
        -------
        prod1 : array_like
            Fusion product 1 distribution.
        prod2 : array_like
            Fusion product 2 distribution.
        """
        m1, q1, m2, q2, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        anum1 = np.round(m1.to("amu").v)
        anum2 = np.round(m2.to("amu").v)
        znum1 = np.round(q1.to("e").v)
        znum2 = np.round(q2.to("e").v)
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        time = 0*unyt.s
        r_edges   = np.linspace(minr, maxr, nr+1)*unyt.m
        phi_edges = np.linspace(minphi, maxphi, nphi+1)*unyt.deg
        z_edges   = np.linspace(minz, maxz, nz+1)*unyt.m
        r   = 0.5 * (r_edges[1:]   + r_edges[:-1])
        phi = 0.5 * (phi_edges[1:] + phi_edges[:-1])
        z   = 0.5 * (z_edges[1:]   + z_edges[:-1])

        self._ascot.input_init(bfield=True, plasma=True)
        nspec, _, _, anums, znums = self._ascot.input_getplasmaspecies()
        ispecies1 = np.nan; ispecies2 = np.nan
        for i in range(nspec):
            if( anum1 == anums[i] and znum1 == znums[i] ):
                ispecies1 = i
            if( anum2 == anums[i] and znum2 == znums[i] ):
                ispecies2 = i
        if np.isnan(ispecies1) or np.isnan(ispecies2):
            self._ascot.input_free(bfield=True, plasma=True)
            raise ValueError("Reactant species not present in plasma input.")
        mult = 0.5 if ispecies1 == ispecies2 else 1.0

        temp  = np.zeros((nr, nphi, nz))
        dens1 = np.zeros((nr, nphi, nz))
        dens2 = np.zeros((nr, nphi, nz))
        for j in range(nphi):
            for i in range(nr):
                for k in range(nz):
                    temp[i,j,k]  = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ti1").to("J")
                    dens1[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ni"+str(ispecies1))
                    dens2[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time,"ni"+str(ispecies2))

        thermal1 = self._init_thermal_data(
            minr, maxr, nr, minphi, maxphi, nphi, minz, maxz, nz, temp, dens1)
        thermal2 = self._init_thermal_data(
            minr, maxr, nr, minphi, maxphi, nphi, minz, maxz, nz, temp, dens2)

        react1 = self._init_afsi_data(dist_thermal=thermal1)
        react2 = self._init_afsi_data(dist_thermal=thermal2)

        class DummyDist(object):
            pass

        react = DummyDist()
        react.n_r      = nr
        react.min_r    = r_edges[0]
        react.max_r    = r_edges[-1]
        react.n_phi    = nphi
        react.min_phi  = phi_edges.to('rad')[0]
        react.max_phi  = phi_edges.to('rad')[-1]
        react.n_z      = nz
        react.min_z    = z_edges[0]
        react.max_z    = z_edges[-1]
        react.n_time   = 1
        react.min_time = 0.0
        react.max_time = 1.0

        prod1 = self._init_product_dist(
            react, q1, minppara, maxppara, nppara, minpperp, maxpperp, npperp)
        prod2 = self._init_product_dist(
            react, q2, minppara, maxppara, nppara, minpperp, maxpperp, npperp)

        _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                           react1, react2, mult, prod1, prod2,
                           )
        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2

    def thermal_exi(
            self, reaction, minr, maxr, nr, minz, maxz, nz,
            minphi=0, maxphi=360, nphi=1, nmc=1000,
            minekin=1.0e6, maxekin=6.0e6, nekin=40,
            minpitch=-1.0, maxpitch=1.0, npitch=20):
        """Calculate thermonuclear fusion between two thermal (Maxwellian)
        species in (E,xi) basis.

        Parameters
        ----------
        reaction : int or str
            Fusion reaction index or name.
        minr : float
            Minimum R value in the output distribution.
        maxr : float
            Maximum R value in the output distribution.
        nr : int
            Number of R bins in the output distribution.
        minz : float
            Minimum z value in the output distribution.
        maxz : float
            Maximum z value in the output distribution.
        nz : int
            Number of z bins in the output distribution.
        minphi : float, optional
            Minimum phi value in the output distribution.
        maxphi : float, optional
            Maximum phi value in the output distribution.
        nphi : int, optional
            Number of phi bins in the output distribution.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        minekin : float, optional
            Minimum ekin value in the output distribution.
        maxekin : float, optional
            Maximum ekin value in the output distribution.
        nekin : int, optional
            Number of ekin bins in the output distribution.
        minpitch : float, optional
            Minimum pitch value in the output distribution.
        maxpitch : float, optional
            Maximum pitch value in the output distribution.
        npitch : int, optional
            Number of pitch bins in the output distribution.

        Returns
        -------
        prod1 : array_like
            Fusion product 1 distribution.
        prod2 : array_like
            Fusion product 2 distribution.
        """
        m1, q1, m2, q2, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        anum1 = np.round(m1.to("amu").v)
        anum2 = np.round(m2.to("amu").v)
        znum1 = np.round(q1.to("e").v)
        znum2 = np.round(q2.to("e").v)
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        time = 0*unyt.s
        r_edges   = np.linspace(minr, maxr, nr+1)*unyt.m
        phi_edges = np.linspace(minphi, maxphi, nphi+1)*unyt.deg
        z_edges   = np.linspace(minz, maxz, nz+1)*unyt.m
        r   = 0.5 * (r_edges[1:]   + r_edges[:-1])
        phi = 0.5 * (phi_edges[1:] + phi_edges[:-1])
        z   = 0.5 * (z_edges[1:]   + z_edges[:-1])

        self._ascot.input_init(bfield=True, plasma=True)
        nspec, _, _, anums, znums = self._ascot.input_getplasmaspecies()
        ispecies1 = np.nan; ispecies2 = np.nan
        for i in range(nspec):
            if( anum1 == anums[i] and znum1 == znums[i] ):
                ispecies1 = i
            if( anum2 == anums[i] and znum2 == znums[i] ):
                ispecies2 = i
        if np.isnan(ispecies1) or np.isnan(ispecies2):
            self._ascot.input_free(bfield=True, plasma=True)
            raise ValueError("Reactant species not present in plasma input.")
        mult = 0.5 if ispecies1 == ispecies2 else 1.0

        temp  = np.zeros((nr, nphi, nz))
        dens1 = np.zeros((nr, nphi, nz))
        dens2 = np.zeros((nr, nphi, nz))
        for j in range(nphi):
            for i in range(nr):
                for k in range(nz):
                    temp[i,j,k]  = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ti1").to("J")
                    dens1[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ni"+str(ispecies1))
                    dens2[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time,"ni"+str(ispecies2))

        thermal1 = self._init_thermal_data(
            minr, maxr, nr, minphi, maxphi, nphi, minz, maxz, nz, temp, dens1)
        thermal2 = self._init_thermal_data(
            minr, maxr, nr, minphi, maxphi, nphi, minz, maxz, nz, temp, dens2)

        react1 = self._init_afsi_data(dist_thermal=thermal1)
        react2 = self._init_afsi_data(dist_thermal=thermal2)

        class DummyDist(object):
            pass

        react = DummyDist()
        react.n_r      = nr
        react.min_r    = r_edges[0]
        react.max_r    = r_edges[-1]
        react.n_phi    = nphi
        react.min_phi  = phi_edges.to('rad')[0]
        react.max_phi  = phi_edges.to('rad')[-1]
        react.n_z      = nz
        react.min_z    = z_edges[0]
        react.max_z    = z_edges[-1]
        react.n_time   = 1
        react.min_time = 0.0
        react.max_time = 1.0

        minekin = (minekin * unyt.eV).to("J").v
        maxekin = (maxekin * unyt.eV).to("J").v
        prod1 = self._init_product_dist(
            react, q1, minekin, maxekin, nekin, minpitch, maxpitch, npitch,
            exi=True)
        prod2 = self._init_product_dist(
            react, q2, minekin, maxekin, nekin, minpitch, maxpitch, npitch,
            exi=True)

        _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                           react1, react2, mult, prod1, prod2,
                           )
        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2

    def beamthermal(self, reaction, beam, swap=False, nmc=1000,
                    minppara=-1.3e-19, maxppara=1.3e-19, nppara=80,
                    minpperp=0, maxpperp=1.3e-19, npperp=40):
        """Calculate beam-thermal fusion.

        Parameters
        ----------
        reaction : int
            Fusion reaction index
        beam : dict
            Beam distribution that acts as the first reactant.
        swap : bool, optional
            If True, beam distribution acts as the second reactant and
            the first reactant is a background species.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        minppara : float, optional
            Minimum ppara value in the output distribution.
        maxppara : float, optional
            Maximum ppara value in the output distribution.
        nppara : int, optional
            Number of pperp bins in the output distribution.
        minpperp : float, optional
            Minimum pperp value in the output distribution.
        maxpperp : float, optional
            Maximum ppara value in the output distribution.
        npperp : int, optional
            Number of pperp bins in the output distribution.

        Returns
        -------
        prod1 : array_like
            Fusion product 1 distribution.
        prod2 : array_like
            Fusion product 2 distribution.
        """
        m1, q1, m2, q2, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        anum1 = np.round(m1.to("amu").v)
        anum2 = np.round(m2.to("amu").v)
        znum1 = np.round(q1.to("e").v)
        znum2 = np.round(q2.to("e").v)
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        time   = 0*unyt.s
        nr     = beam.abscissa("r").size
        nphi   = beam.abscissa("phi").size
        nz     = beam.abscissa("z").size
        r      = beam.abscissa("r")
        phi    = beam.abscissa("phi")
        z      = beam.abscissa("z")
        minr   = beam.abscissa_edges("r")[0]
        maxr   = beam.abscissa_edges("r")[-1]
        minphi = beam.abscissa_edges("phi")[0]
        maxphi = beam.abscissa_edges("phi")[-1]
        minz   = beam.abscissa_edges("z")[0]
        maxz   = beam.abscissa_edges("r")[-1]

        self._ascot.input_init(bfield=True, plasma=True)
        nspec, _, _, anums, znums = self._ascot.input_getplasmaspecies()
        ispecies = np.nan
        for i in range(nspec):
            if( swap and anum1 == anums[i] and znum1 == znums[i] ):
                ispecies = i
            if( not swap and anum2 == anums[i] and znum2 == znums[i] ):
                ispecies = i
        if np.isnan(ispecies):
            self._ascot.input_free(bfield=True, plasma=True)
            raise ValueError("Reactant species not present in plasma input.")

        temp = np.zeros((nr, nphi, nz))
        dens = np.zeros((nr, nphi, nz))
        for j in range(nphi):
            for i in range(nr):
                for k in range(nz):
                    temp[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ti1").to("J")
                    dens[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ni"+str(ispecies))

        thermal = self._init_thermal_data(
            minr, maxr, nr, minphi, maxphi, nphi, minz, maxz, nz, temp, dens)

        dist = self._init_dist_5d(beam)
        react1 = self._init_afsi_data(dist_thermal=thermal)
        react2 = self._init_afsi_data(dist_5D=dist)

        prod1 = self._init_product_dist(
            dist, q1, minppara, maxppara, nppara, minpperp, maxpperp, npperp)
        prod2 = self._init_product_dist(
            dist, q2, minppara, maxppara, nppara, minpperp, maxpperp, npperp)

        mult = 1.0
        if swap:
            _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                               react1, react2, mult, prod1, prod2,
                               )
        else:
            _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                               react2, react1, mult, prod1, prod2,
                               )
        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2

    def beambeam(self, reaction, beam1, beam2=None, nmc=1000,
                 minppara=-1.3e-19, maxppara=1.3e-19, nppara=80,
                 minpperp=0, maxpperp=1.3e-19, npperp=40):
        """Calculate beam-beam fusion.

        Parameters
        ----------
        reaction : int
            Fusion reaction index.
        beam1 : dict
            Beam1 distribution.
        beam2 : dict, optional
            Beam2 distribution or None to calculate fusion generation with
            beam1 itself.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        minppara : float, optional
            Minimum ppara value in the output distribution.
        maxppara : float, optional
            Maximum ppara value in the output distribution.
        nppara : int, optional
            Number of pperp bins in the output distribution.
        minpperp : float, optional
            Minimum pperp value in the output distribution.
        maxpperp : float, optional
            Maximum ppara value in the output distribution.
        npperp : int, optional
            Number of pperp bins in the output distribution.

        Returns
        -------
        prod1 : array_like
            Fusion product 1 distribution.
        prod2 : array_like
            Fusion product 2 distribution.
        """
        _, _, _, _, _, qprod1, _, qprod2, _ = self.reactions(reaction)
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        reaction = reactions[reaction]
        q1 = np.round(qprod1.to("e").v)
        q2 = np.round(qprod2.to("e").v)

        self._ascot.input_init(bfield=True, plasma=True)
        dist1 = self._init_dist_5d(beam1)

        react1 = self._init_afsi_data(dist_5D=dist1)
        if beam2 is not None:
            dist2  = self._init_dist_5d(beam2)
            react2 = self._init_afsi_data(dist_5D=dist2)
            mult = 1.0
        else:
            react2 = react1
            mult = 0.5

        prod1 = self._init_product_dist(
            dist1, q1, minppara, maxppara, nppara, minpperp, maxpperp, npperp)
        prod2 = self._init_product_dist(
            dist1, q2, minppara, maxppara, nppara, minpperp, maxpperp, npperp)

        _LIBASCOT.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                           react1, react2, mult, prod1, prod2,
                           )

        self._ascot.input_free(bfield=True, plasma=True)

        # Reload Ascot
        self._ascot.file_load(self._ascot.file_getpath())
        return prod1, prod2

    def _init_afsi_data(self, dist_5D=None, dist_thermal=None):
        afsidata = STRUCT_AFSIDATA()
        if dist_5D is not None:
            afsidata.type = 1
            afsidata.beam = ctypes.pointer(dist_5D)
        elif dist_thermal is not None:
            afsidata.type = 2
            afsidata.thermal = ctypes.pointer(dist_thermal)
        else:
            afsidata.type = 0
        return afsidata

    def _init_thermal_data(self, minr, maxr, nr, minphi, maxphi, nphi, minz,
                           maxz, nz, temperature, density):
        thermaldata         = STRUCT_AFSITHERMAL()
        thermaldata.n_r     = nr
        thermaldata.min_r   = minr
        thermaldata.max_r   = maxr
        thermaldata.n_phi   = nphi
        thermaldata.min_phi = minphi * np.pi/180
        thermaldata.max_phi = maxphi * np.pi/180
        thermaldata.n_z     = nz
        thermaldata.min_z   = minz
        thermaldata.max_z   = maxz
        thermaldata.temperature = npctypes.as_ctypes(
            np.ascontiguousarray(temperature.ravel(), dtype="f8"))
        thermaldata.density = npctypes.as_ctypes(
            np.ascontiguousarray(density.ravel(), dtype="f8"))
        return thermaldata

    def _init_dist_5d(self, dist):
        data = STRUCT_HIST()
        coordinates = np.array([0, 1, 2, 5, 6, 14, 15], dtype="uint32")
        nbin = np.array([
            dist.abscissa("r").size, dist.abscissa("phi").size,
            dist.abscissa("z").size, dist.abscissa("ppar").size,
            dist.abscissa("pperp").size, dist.abscissa("time").size,
            dist.abscissa("charge").size], dtype="u8")
        binmin = np.array([
            dist.abscissa_edges("r")[0], dist.abscissa_edges("phi")[0],
            dist.abscissa_edges("z")[0], dist.abscissa_edges("ppar")[0],
            dist.abscissa_edges("pperp")[0], dist.abscissa_edges("time")[0],
            dist.abscissa_edges("charge")[0]])
        binmax = np.array([
            dist.abscissa_edges("r")[-1], dist.abscissa_edges("phi")[-1],
            dist.abscissa_edges("z")[-1], dist.abscissa_edges("ppar")[-1],
            dist.abscissa_edges("pperp")[-1], dist.abscissa_edges("time")[-1],
            dist.abscissa_edges("charge")[-1]])
        _LIBASCOT.hist_init(
            ctypes.byref(data),
            coordinates.size,
            coordinates.ctypes.data_as(ctypes.POINTER(ctypes.c_uint)),
            binmin.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            binmax.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            nbin.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t))
            )
        return data

    def _init_product_dist(self, react, charge, minp1, maxp1, np1,
                           minp2, maxp2, np2, exi=False):
        prod = STRUCT_HIST()
        if(exi):
            coordinates = np.array([0, 1, 2, 10, 11, 14, 15], dtype="uint32")
        else:
            coordinates = np.array([0, 1, 2, 5, 6, 14, 15], dtype="uint32")
        nbin = np.array([
            react.n_r, react.n_phi, react.n_z, np1, np2, react.n_time, 1
            ], dtype="u8")
        binmin = np.array([
            react.min_r, react.min_phi, react.min_z, minp1, minp2,
            react.min_time, charge[0] - 1])
        binmax = np.array([
            react.max_r, react.max_phi, react.max_z, maxp1, maxp2,
            react.max_time, charge[0] + 1])
        _LIBASCOT.hist_init(
            ctypes.byref(prod),
            coordinates.size,
            coordinates.ctypes.data_as(ctypes.POINTER(ctypes.c_uint)),
            binmin.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            binmax.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            nbin.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t))
            )
        return prod

    def reactions(self, reaction=None):
        """Return reaction data for a given reaction.

        Parameters
        ----------
        reaction : str, optional
            Reaction or None to return list of available reactions.

        Returns
        -------
        reactions : [str]
            List of available reactions if ``reaction=None``.
        m1 : float
            Mass of reactant 1.
        q1 : float
            Charge of reactant 1.
        m2 : float
            Mass of reactant 2.
        q2 : float
            Charge of reactant 2.
        mprod1 : float
            Mass of product 1.
        qprod1 : float
            Charge of product 1.
        mprod2 : float
            Mass of product 2.
        qprod2 : float
            Charge of product 2.
        q : float
            Energy released in the reaction.
        """
        reactions = {v: k for k, v in AFSI_REACTIONS.items()}
        if reaction is None:
            return reactions.keys()
        if not reaction in reactions:
            raise ValueError("Unknown reaction")

        m1     = np.zeros((1,), dtype="f8")
        q1     = np.zeros((1,), dtype="f8")
        m2     = np.zeros((1,), dtype="f8")
        q2     = np.zeros((1,), dtype="f8")
        mprod1 = np.zeros((1,), dtype="f8")
        qprod1 = np.zeros((1,), dtype="f8")
        mprod2 = np.zeros((1,), dtype="f8")
        qprod2 = np.zeros((1,), dtype="f8")
        q      = np.zeros((1,), dtype="f8")
        fun = _LIBASCOT.boschhale_reaction
        fun.restype  = None
        fun.argtypes = [ctypes.c_int, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL,
                        PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL, PTR_REAL]
        fun(reactions[reaction], m1, q1, m2, q2, mprod1, qprod1, mprod2,
            qprod2, q)

        return m1*unyt.kg, q1*unyt.C, m2*unyt.kg, q2*unyt.C, mprod1*unyt.kg, \
            qprod1*unyt.C, mprod2*unyt.kg, qprod2*unyt.C, q*unyt.eV

class AfsiMixin(DistMixin):
    """Mixin class with post-processing results related to AFSI.
    """

    def _require(self, *args):
        """Check if required data is present and raise exception if not.

        This is a helper function to quickly check that the data is available.

        Parameters
        ----------
        *args : `str`
            Name(s) of the required data.

        Raises
        ------
        AscotNoDataException
            Raised if the required data is not present.
        """
        for arg in args:
            if not hasattr(self, arg):
                raise AscotNoDataException(
                    "Data for \"" +  arg + "\" is required but not present.")

    def getdist(self, dist, exi=False, ekin_edges=None, pitch_edges=None,
                plotexi=False):
        """Return 5D distribution function of one of the fusion products.

        Parameters
        ----------
        dist : {"prod1", "prod2"}
            Which product to return.
        exi : bool, optional
            Convert the momentum space to energy-pitch.

            The distribution is normalized to conserve the particle number.
        ekin_edges : int or array_like, optional
            Number of bins or bin edges in the energy abscissa if ``exi=True``.
        pitch_edges : int or array_like, optional
            Number of bins or bin edges in the pitch abscissa if ``exi=True``.
        plotexi : bool, optional
            Visualize the transformation from ppar-perp to energy-pitch if
            if ``exi=True``.

            Use this option to adjust energy and pitch abscissae to your liking.

        Returns
        -------
        data : :class:`DistData`
            The distribution data object.
        """
        if dist == "prod1":
            self._require("_prod1dist5d", "_reaction")
            distout = self._prod1dist5d.get()
            mass = self._reaction.get()[2]
        elif dist == "prod2":
            self._require("_prod2dist5d", "_reaction")
            distout = self._prod2dist5d.get()
            mass = self._reaction.get()[3]
        else:
            raise ValueError("dist must be either 'prod1' or 'prod2'")

        return self._getdist(distout, mass, exi=exi, ekin_edges=ekin_edges,
                             pitch_edges=pitch_edges, plotexi=plotexi)
