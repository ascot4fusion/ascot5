import numpy as np
import unyt

from a5py.ascot5io import State, Orbits
from .runmixin import RunMixin

class VirtualRun(RunMixin):

    def __init__(self, ascot, nmrk, inistate, endstate, npoint, mode, orbmode,
                 orbits):
        self.inistate = VirtualState(ascot, nmrk, inistate)
        self.endstate = VirtualState(ascot, nmrk, endstate)
        self.orbit = VirtualOrbits(ascot, nmrk, npoint, mode, orbmode, orbits)

class VirtualState():

    def __init__(self, ascot, nmrk, state):
        self.ascot = ascot
        self.state = state
        self.nmrk = nmrk.value

    def get(self, *qnt, mode="gc"):

        def _val(q):
            """Read quantity from HDF5.
            """
            arr = lambda q : np.array([getattr(self.state[i], q)
                                       for i in range(self.nmrk)])

            m = unyt.m; rad = unyt.rad; kg = unyt.kg; s = unyt.s; J = unyt.J
            T = unyt.T; amu = unyt.amu; C = unyt.C; nodim = unyt.dimensionless
            if   q == "ids":       return arr("id") * nodim
            elif q == "rprt":      return arr(q) * m
            elif q == "phiprt":    return arr(q) * rad
            elif q == "zprt":      return arr(q) * m
            elif q == "prprt":     return arr("p_r") * kg*m/s
            elif q == "pphiprt":   return arr("p_phi") * kg*m/s
            elif q == "pzprt":     return arr("p_z") * kg*m/s
            elif q == "r":         return arr(q) * m
            elif q == "phi":       return arr(q) * rad
            elif q == "z":         return arr(q) * m
            elif q == "ppar":      return arr(q) * kg*m/s
            elif q == "mu":        return arr(q) * J/T
            elif q == "zeta":      return arr(q) * rad
            elif q == "weight":    return arr(q) * 1
            elif q == "time":      return arr(q) * s
            elif q == "mileage":   return arr(q) * s
            elif q == "cputime":   return arr(q) * s
            elif q == "rho":       return arr(q) * nodim
            elif q == "theta":     return arr(q) * rad
            elif q == "mass":      return arr(q) * amu
            elif q == "br":        return arr("B_r") * T
            elif q == "bphi":      return arr("B_phi") * T
            elif q == "bz":        return arr("B_z") * T
            elif q == "ids":       return arr(q) * 1
            elif q == "endcond":   return arr(q) * 1
            elif q == "walltile":  return arr(q) * 1
            elif q == "charge":    return arr(q) * C
            elif q == "anum":      return arr(q) * 1
            elif q == "znum":      return arr(q) * 1
            elif q == "errormsg":  return arr("err") * 1
            elif q == "errorline": return arr("err") * 1
            elif q == "errormod":  return arr("err") * 1
            return arr

        def _eval(r, phi, z, t, *q):
            """Evaluate input quantities at marker position.
            """
            return self.ascot.input_eval(r, phi, z, t, *q)

        return State._getactual(mode, _val, _eval, *qnt)

class VirtualOrbits():

    def __init__(self, ascot, nmrk, npoint, mode, orbmode, orbits):
        self.ascot = ascot
        self.nmrk = nmrk.value
        self.npnt = npoint
        self.mode = mode
        self.orbmode = orbmode
        self.arr  = orbits

    def get(self, inistate, endstate, *qnt):
        """
        """
        # Prepare helper variables and functions
        arrlen = self.nmrk * self.npnt
        def _val(q, mask=None):
            """Read quantity from the orbit array.
            """
            idx = np.array(self.arr[arrlen*0:arrlen*1]) > 0
            m = unyt.m; rad = unyt.rad; kg = unyt.kg; s = unyt.s; J = unyt.J
            T = unyt.T; amu = unyt.amu; C = unyt.C; nodim = unyt.dimensionless
            if q == "ids":
                val = np.array(self.arr[arrlen*0:arrlen*1]) * nodim
            elif q == "mileage":
                val = np.array(self.arr[arrlen*1:arrlen*2]) * s
            elif q == "r":
                val = np.array(self.arr[arrlen*2:arrlen*3]) * m
            elif q == "phi":
                val = np.array(self.arr[arrlen*3:arrlen*4]) * rad
            elif q == "z":
                val = np.array(self.arr[arrlen*4:arrlen*5]) * m

            if self.mode == Orbits.GYROORBIT:
                if q == "pr":
                    val = np.array(self.arr[arrlen*5:arrlen*6]) * kg*m/s
                elif q == "pphi":
                    val = np.array(self.arr[arrlen*6:arrlen*7]) * kg*m/s
                elif q == "pz":
                    val = np.array(self.arr[arrlen*7:arrlen*8]) * kg*m/s
                elif q == "weight":
                    val = np.array(self.arr[arrlen*8:arrlen*9]) * nodim
                elif q == "charge":
                    val = np.array(self.arr[arrlen*9:arrlen*10]) * C
                elif q == "rho":
                    val = np.array(self.arr[arrlen*10:arrlen*11]) * nodim
                elif q == "theta":
                    val = np.array(self.arr[arrlen*11:arrlen*12]) * rad
                elif q == "br":
                    val = np.array(self.arr[arrlen*12:arrlen*13]) * T
                elif q == "bphi":
                    val = np.array(self.arr[arrlen*13:arrlen*14]) * T
                elif q == "bz":
                    val = np.array(self.arr[arrlen*14:arrlen*15]) * T
                elif q == "simmode":
                    val = np.array(self.arr[arrlen*15:arrlen*16]) * nodim
                if q == "pncrid" and self.orbmode:
                    val = np.array(self.arr[arrlen*16:arrlen*17]) * nodim
                if q == "pncrdir" and self.orbmode:
                    val = np.array(self.arr[arrlen*17:arrlen*18]) * nodim
            elif self.mode == Orbits.GUIDINGCENTER:
                if q == "ppar":
                    val = np.array(self.arr[arrlen*5:arrlen*6]) * kg*m/s
                elif q == "mu":
                    val = np.array(self.arr[arrlen*6:arrlen*7]) * J/T
                elif q == "zeta":
                    val = np.array(self.arr[arrlen*7:arrlen*8]) * rad
                elif q == "weight":
                    val = np.array(self.arr[arrlen*8:arrlen*9]) * nodim
                elif q == "charge":
                    val = np.array(self.arr[arrlen*9:arrlen*10]) * C
                elif q == "rho":
                    val = np.array(self.arr[arrlen*10:arrlen*11]) * nodim
                elif q == "theta":
                    val = np.array(self.arr[arrlen*11:arrlen*12]) * rad
                elif q == "br":
                    val = np.array(self.arr[arrlen*12:arrlen*13]) * T
                elif q == "bphi":
                    val = np.array(self.arr[arrlen*13:arrlen*14]) * T
                elif q == "bz":
                    val = np.array(self.arr[arrlen*14:arrlen*15]) * T
                elif q == "simmode":
                    val = np.array(self.arr[arrlen*15:arrlen*16]) * nodim
                if q == "pncrid" and self.orbmode:
                    val = np.array(self.arr[arrlen*16:arrlen*17]) * nodim
                if q == "pncrdir" and self.orbmode:
                    val = np.array(self.arr[arrlen*17:arrlen*18]) * nodim
            elif self.mode == Orbits.FIELDLINE:
                if q == "rho":
                    val = np.array(self.arr[arrlen*5:arrlen*6]) * nodim
                elif q == "theta":
                    val = np.array(self.arr[arrlen*6:arrlen*7]) * rad
                elif q == "br":
                    val = np.array(self.arr[arrlen*7:arrlen*8]) * T
                elif q == "bphi":
                    val = np.array(self.arr[arrlen*8:arrlen*9]) * T
                elif q == "bz":
                    val = np.array(self.arr[arrlen*9:arrlen*10]) * T
                elif q == "simmode":
                    val = np.array(self.arr[arrlen*10:arrlen*11]) * nodim
                if q == "pncrid" and self.orbmode:
                    val = np.array(self.arr[arrlen*11:arrlen*12]) * nodim
                if q == "pncrdir" and self.orbmode:
                    val = np.array(self.arr[arrlen*12:arrlen*13]) * nodim
            elif self.mode == Orbits.HYBRID:
                if q == "pr":
                    val = np.array(self.arr[arrlen*5:arrlen*6]) * kg*m/s
                elif q == "pphi":
                    val = np.array(self.arr[arrlen*6:arrlen*7]) * kg*m/s
                elif q == "pz":
                    val = np.array(self.arr[arrlen*7:arrlen*8]) * kg*m/s
                elif q == "ppar":
                    val = np.array(self.arr[arrlen*8:arrlen*9]) * kg*m/s
                elif q == "mu":
                    val = np.array(self.arr[arrlen*9:arrlen*10]) * J/T
                elif q == "zeta":
                    val = np.array(self.arr[arrlen*10:arrlen*11]) * rad
                elif q == "weight":
                    val = np.array(self.arr[arrlen*11:arrlen*12]) * nodim
                elif q == "charge":
                    val = np.array(self.arr[arrlen*12:arrlen*13]) * C
                elif q == "rho":
                    val = np.array(self.arr[arrlen*13:arrlen*14]) * nodim
                elif q == "theta":
                    val = np.array(self.arr[arrlen*14:arrlen*15]) * rad
                elif q == "br":
                    val = np.array(self.arr[arrlen*15:arrlen*16]) * T
                elif q == "bphi":
                    val = np.array(self.arr[arrlen*16:arrlen*17]) * T
                elif q == "bz":
                    val = np.array(self.arr[arrlen*17:arrlen*18]) * T
                elif q == "simmode":
                    val = np.array(self.arr[arrlen*18:arrlen*19]) * nodim
                if q == "pncrid" and self.orbmode:
                    val = np.array(self.arr[arrlen*19:arrlen*20]) * nodim
                if q == "pncrdir" and self.orbmode:
                    val = np.array(self.arr[arrlen*20:arrlen*21]) * nodim
            return val[idx][mask].ravel()

        mode = _val("simmode")
        ids  = inistate.get("ids")[0]
        val  = inistate.get("mass")[0]
        _, idx = np.unique(_val("ids").v, return_inverse=True)
        mass = val[idx]
        val     = inistate.get("time")[0]
        time = val[idx]
        val     = inistate.get("mileage")[0]
        connlen = val[idx]
        connlen -= _val("mileage")
        if not Orbits.FIELDLINE in mode:
            time = time + _val("mileage")

        def _eval(q, mask=None):
            """Evaluate input quantities at marker position.
            """
            return self.ascot.input_eval(
                _val("r", mask=mask), _val("phi",  mask=mask),
                _val("z", mask=mask), time[mask], *[q])

        return Orbits.getactual(mass, time, connlen, mode, _val, _eval, *qnt)
