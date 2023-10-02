"""AFSI5: Versatile fusion source integrator AFSI for fast ion and neutron
studies in fusion devices
"""
import ctypes
import copy
import numpy as np
import numpy.ctypeslib as npctypes

from a5py.ascotpy import ascot2py

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

    Possible fusion reactions are listed below in a format
    reactant 1 + reactant 2 -> product1 + product2:

    1. D + T   -> He4 + n
    2. D + He3 -> He4 + p
    3. D + D   -> T   + p
    4. D + D   -> He3 + n

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
                minphi=0, maxphi=2*np.pi, nphi=1, nmc=1000, mult=1.0,
                ispecies1=1, ispecies2=1,
                minppara=-1.3e-19, maxppara=1.3e-19, nppara=80,
                minpperp=0, maxpperp=1.3e-19, npperp=40):
        """Calculate thermonuclear fusion between two thermal (Maxwellian)
        species.

        Parameters
        ----------
        reaction : int
            Fusion reaction index.
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
        mult : float, optional
            Multiplier for the fusion rate.

            Use 0.5 if the population is interacting with itself to avoid double
            counting.
        ispecies1 : int, optional
            Ion species index (as they are listed in the plasma input) for
            the first reactant.
        ispecies2 : int, optional
            Ion species index (as they are listed in the plasma input) for
            the second reactant.
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
        time = 0
        r_edges   = np.linspace(minr, maxr, nr+1)
        phi_edges = np.linspace(minphi, maxphi, nphi+1)
        z_edges   = np.linspace(minz, maxz, nz+1)
        r   = (r_edges[1:]   + r_edges[:-1])   / 2
        phi = (phi_edges[1:] + phi_edges[:-1]) / 2
        z   = (z_edges[1:]   + z_edges[:-1])   / 2
        temp  = np.zeros((nr, nphi, nz))
        dens1 = np.zeros((nr, nphi, nz))
        dens2 = np.zeros((nr, nphi, nz))

        for j in range(nphi):
            for i in range(nr):
                for k in range(nz):
                    temp[i,j,k]  = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ti1")
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

        d = {
            "nr" : nr, "nphi" : nphi, "nz" : nz, "nppar" : 1, "npperp" : 1,
            "ntime" : 1, "ncharge" : 1,
            "r_edges"      : r_edges,
            "phi_edges"    : phi_edges,
            "z_edges"      : z_edges,
            "ppar_edges"   : np.linspace(0, 1, 2),
            "pperp_edges"  : np.linspace(0, 1, 2),
            "charge_edges" : np.linspace(0, 2, 2),
            "time_edges"   : np.linspace(0, 1, 2),
            "histogram"    : np.zeros((nr, nphi, nz, 1, 1, 1, 1))
        }
        temp = self._init_dist_5d(d)

        q1, q2 = self._product_charge(reaction)
        prod1, prod1_data = self._init_product_dist(
            temp, q1, minppara, maxppara, nppara, minpperp, maxpperp, npperp)
        prod2, prod2_data = self._init_product_dist(
            temp, q2, minppara, maxppara, nppara, minpperp, maxpperp, npperp)

        ascot2py.afsi_run(reaction, nmc, react1, react2, mult, prod1, prod2,
                          prod1_data, prod2_data)

        return prod1, prod2

    def beamthermal(self, reaction, beam, it=1, nmc=1000, mult=1.0, ispecies=1,
                    swap=False, minppara=-1.3e-19, maxppara=1.3e-19, nppara=80,
                    minpperp=0, maxpperp=1.3e-19, npperp=40):
        """Calculate beam-thermal fusion.

        Parameters
        ----------
        reaction : int
            Fusion reaction index
        beam1 : dict
            Beam distribution that acts as reactant 1.
        it : int, optional
            Time index from 5D distribution
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        mult : float, optional
            Multiplier for the fusion rate.

            Use 0.5 if the population is interacting with itself to avoid double
            counting.
        ispecies : int, optional
            Ion species index (as they are listed in the plasma input) for
            thermal reactant.
        swap : bool, optional
            Swap reactants so that the beam becomes reactant 2 and thermal
            species reactant 1.
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
        time   = 0
        nr     = beam["nr"]
        nphi   = beam["nphi"]
        nz     = beam["nz"]
        r      = beam["r"]
        phi    = beam["phi"]
        z      = beam["z"]
        minr   = beam["r_edges"][0]
        maxr   = beam["r_edges"][-1]
        minphi = beam["phi_edges"][0]
        maxphi = beam["phi_edges"][-1]
        minz   = beam["z_edges"][0]
        maxz   = beam["z_edges"][-1]

        temp = np.zeros((nr, nphi, nz))
        dens = np.zeros((nr, nphi, nz))
        for j in range(nphi):
            for i in range(nr):
                for k in range(nz):
                    temp[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ti1")
                    dens[i,j,k] = self._ascot.input_eval(
                        r[i], phi[j], z[k], time, "ni"+str(ispecies))

        thermal = self._init_thermal_data(
            minr, maxr, nr, minphi, maxphi, nphi, minz, maxz, nz, temp, dens)

        dist = self._init_dist_5d(beam)
        react1 = self._init_afsi_data(dist_thermal=thermal)
        react2 = self._init_afsi_data(dist_5D=dist)

        q1, q2 = self._product_charge(reaction)
        prod1, prod1_data = self._init_product_dist(
            dist, q1, minppara, maxppara, nppara, minpperp, maxpperp, npperp)
        prod2, prod2_data = self._init_product_dist(
            dist, q2, minppara, maxppara, nppara, minpperp, maxpperp, npperp)

        if swap:
            ascot2py.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                              react1, react2, mult, prod1, prod2,
                              prod1_data, prod2_data)
        else:
            ascot2py.afsi_run(ctypes.byref(self._ascot._sim), reaction, nmc,
                              react2, react1, mult, prod1, prod2,
                              prod1_data, prod2_data)
        return prod1, prod2

    def beambeam(self, reaction, beam1, beam2=None, it=0, nmc=1000, mult=1.0,
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
        it : int, optional
            Time index from 5D distribution.
        nmc : int, optional
            Number of MC samples used in each (R, phi, z) bin.
        mult : float, optional
            Multiplier for the fusion rate.

            Use 0.5 if the population is interacting with itself to avoid double
            counting.
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
        beam1["ntime"] = 1
        beam1["time_edges"] = beam1["time_edges"][it:it+1]
        beam1["histogram"]  = beam1["histogram"][:,:,:,:,:,[it],:]
        dist1 = self._init_dist_5d(beam1)

        react1 = self._init_afsi_data(dist_5D=dist1)
        if beam2 is not None:
            beam2["ntime"] = 1
            beam2["time_edges"] = beam2["time_edges"][it:it+1]
            beam2["histogram"]  = beam2["histogram"][:,:,:,:,:,[it],:]
            dist2  = self._init_dist_5d(beam2)
            react2 = self._init_afsi_data(dist_5D=dist2)
        else:
            react2 = react1

        q1, q2 = self._product_charge(reaction)
        prod1, prod1_data = self._init_product_dist(
            dist1, q1, minppara, maxppara, nppara, minpperp, maxpperp, npperp)
        prod2, prod2_data = self._init_product_dist(
            dist1, q2, minppara, maxppara, nppara, minpperp, maxpperp, npperp)

        ascot2py.afsi_run(reaction, nmc, react1, react2, mult, prod1, prod2,
                          prod1_data, prod2_data)

        return prod1, prod2

    def _init_afsi_data(self, dist_5D=None, dist_thermal=None):
        afsidata = ascot2py.afsi_data()
        if dist_5D is not None:
            afsidata.type = 1
            afsidata.dist_5D = ctypes.pointer(dist_5D)
        elif dist_thermal is not None:
            afsidata.type = 2
            afsidata.dist_thermal = ctypes.pointer(dist_thermal)
        else:
            afsidata.type = 0
        return afsidata

    def _init_thermal_data(self, minr, maxr, nr, minphi, maxphi, nphi, minz,
                           maxz, nz, temperature, density):
        thermaldata         = ascot2py.afsi_thermal_data()
        thermaldata.n_r     = nr
        thermaldata.min_r   = minr
        thermaldata.max_r   = maxr
        thermaldata.n_phi   = nphi
        thermaldata.min_phi = minphi
        thermaldata.max_phi = maxphi
        thermaldata.n_z     = nz
        thermaldata.min_z   = minz
        thermaldata.max_z   = maxz
        thermaldata.temperature = npctypes.as_ctypes(
            np.ascontiguousarray(temperature.ravel(), dtype="f8"))
        thermaldata.density = npctypes.as_ctypes(
            np.ascontiguousarray(density.ravel(), dtype="f8"))
        return thermaldata

    def _init_dist_5d(self, dist):
        data           = ascot2py.struct_c__SA_dist_5D_data()
        data.n_r       = dist["nr"]
        data.min_r     = dist["r_edges"][0]
        data.max_r     = dist["r_edges"][-1]
        data.n_phi     = dist["nphi"]
        data.min_phi   = dist["phi_edges"][0]
        data.max_phi   = dist["phi_edges"][-1]
        data.n_z       = dist["nz"]
        data.min_z     = dist["z_edges"][0]
        data.max_z     = dist["z_edges"][-1]
        data.n_ppara   = dist["nppar"]
        data.min_ppara = dist["ppar_edges"][0]
        data.max_ppara = dist["ppar_edges"][-1]
        data.n_pperp   = dist["npperp"]
        data.min_pperp = dist["pperp_edges"][0]
        data.max_pperp = dist["pperp_edges"][-1]
        data.n_time    = dist["ntime"]
        data.min_time  = dist["time_edges"][0]
        data.max_time  = dist["time_edges"][-1]
        data.n_q       = dist["ncharge"]
        data.min_q     = dist["charge_edges"][0]
        data.max_q     = dist["charge_edges"][-1]
        data.histogram = npctypes.as_ctypes(np.ascontiguousarray(
            dist["histogram"].ravel(), dtype="f8"))
        return data

    def _init_product_dist(self, react, charge, minppara, maxppara, nppara,
                           minpperp, maxpperp, npperp):
        prod            = ascot2py.struct_c__SA_dist_5D_offload_data()
        prod.n_r        = react.n_r
        prod.min_r      = react.min_r
        prod.max_r      = react.max_r
        prod.n_phi      = react.n_phi
        prod.min_phi    = react.min_phi
        prod.max_phi    = react.max_phi
        prod.n_z        = react.n_z
        prod.min_z      = react.min_z
        prod.max_z      = react.max_z
        prod.n_ppara    = nppara
        prod.min_ppara  = minppara
        prod.max_ppara  = maxppara
        prod.n_pperp    = npperp
        prod.min_pperp  = minpperp
        prod.max_pperp  = maxpperp
        prod.n_time     = react.n_time
        prod.min_time   = react.min_time
        prod.max_time   = react.max_time
        prod.n_charge   = 1
        prod.min_charge = charge - 1
        prod.max_charge = charge + 1

        distsize = prod.n_r * prod.n_phi * prod.n_z * prod.n_ppara \
            * prod.n_pperp
        #prod.histogram  = npctypes.as_ctypes(
        #    np.ascontiguousarray(np.zeros(distsize), dtype="f8") )
        #shape = [prod.n_r, prod.n_phi, prod.n_z, prod.n_ppara, prod.n_pperp]
        #prod.dist = npctypes.as_array(prod.histogram, shape=shape)
        #offloadarray = libascot_allocate_reals()
        offload_array = npctypes.as_ctypes(
            np.ascontiguousarray(np.zeros(distsize), dtype="f8") )
        return prod, offload_array

    def _product_charge(self, reaction):
        if reaction == 1:
            return 2, 0
        elif reaction == 2:
            return 2, 1
        elif reaction == 3:
            return 1, 1
        elif reaction == 4:
            return 2, 0
