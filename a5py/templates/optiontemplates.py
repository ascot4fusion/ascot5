"""Options templates for various simulations
"""
from a5py.ascot5io.options import Opt

class OptionTemplates():
    """All options templates.
    """

    def options_tutorial(self):
        """Create fast slowing-down simulation options with orbits and
        distributions enabled.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        out = Opt.get_default()
        out.update(
            {"SIM_MODE":2, "ENABLE_ADAPTIVE":1,
             "ENDCOND_SIMTIMELIM":1,      "ENDCOND_CPUTIMELIM":1,
             "ENDCOND_ENERGYLIM":1,       "ENDCOND_WALLHIT":1,
             "ENDCOND_LIM_SIMTIME":0.5,   "ENDCOND_MAX_MILEAGE":0.5,
             "ENDCOND_MAX_CPUTIME":1.0e1, "ENDCOND_MIN_ENERGY":2.0e3,
             "ENDCOND_MIN_THERMAL":2.0,
             "ENABLE_ORBIT_FOLLOWING":1,  "ENABLE_COULOMB_COLLISIONS":1,
             "ENABLE_DIST_5D":1,          "ENABLE_DIST_6D":1,
             "ENABLE_DIST_RHO5D":1,       "ENABLE_DIST_RHO6D":1,
             "ENABLE_DIST_COM":1,
             "DIST_MIN_R":3.5, "DIST_MAX_R":8.5, "DIST_NBIN_R":12,
             "DIST_MIN_PHI":0, "DIST_MAX_PHI":360, "DIST_NBIN_PHI":20,
             "DIST_MIN_Z":-2.45, "DIST_MAX_Z":2.45, "DIST_NBIN_Z":24,
             "DIST_MIN_RHO":0, "DIST_MAX_RHO":1, "DIST_NBIN_RHO":11,
             "DIST_MIN_THETA":0, "DIST_MAX_THETA":360, "DIST_NBIN_THETA":13,
             "DIST_MIN_PPA":-10e-20,"DIST_MAX_PPA":10e-20,"DIST_NBIN_PPA":36,
             "DIST_MIN_PPE":0, "DIST_MAX_PPE":10e-20, "DIST_NBIN_PPE":18,
             "DIST_MIN_PR":-10e-20, "DIST_MAX_PR":10e-20,
             "DIST_NBIN_PR":14,
             "DIST_MIN_PPHI":-10e-20, "DIST_MAX_PPHI":10e-20,
             "DIST_NBIN_PPHI":15,
             "DIST_MIN_PZ":-10e-20, "DIST_MAX_PZ":10e-20, "DIST_NBIN_PZ":16,
             "DIST_MIN_TIME":0, "DIST_MAX_TIME":3e-2, "DIST_NBIN_TIME":2,
             "DIST_MIN_CHARGE":-1.5, "DIST_MAX_CHARGE":2.5,
             "DIST_NBIN_CHARGE":1,
             "ENABLE_ORBITWRITE":1, "ORBITWRITE_MODE":1,
             "ORBITWRITE_NPOINT":100, "ORBITWRITE_INTERVAL":0,}
        )
        return ("opt", out)
    
    def options_gctracer(self, tmax: float, npoints=200):
        """Generate options to trace markers for a fixed number of orbits.

        Collisionless orbits are traced only for a single poloidal transit,
        and their orbits are recorded.

        Parameters
        ----------
        tmax : float
            Maximum simulation time in seconds.
        npoints : int, optional
            Number of points to record in the orbit. Default is 200.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        out = Opt.get_default()
        out.update({
            "ENABLE_ORBIT_FOLLOWING":1,
            "ENABLE_ORBITWRITE": 1,
            "ORBITWRITE_MODE":1,
            "ORBITWRITE_INTERVAL": tmax / npoints,
            "ORBITWRITE_TOROIDALANGLES":0.0,
            "ORBITWRITE_POLOIDALANGLES":0.0,
            "ORBITWRITE_NPOINT": npoints,
            "ENDCOND_MAX_POLOIDALORBS":1,
            "ENDCOND_MAX_TOROIDALORBS":1000,
            "ENDCOND_MAXORBS":1,
            "SIM_MODE":2,
            "ENABLE_ADAPTIVE":1,
            "ENDCOND_MAX_MILEAGE": tmax,
        })
        return ("opt", out)
