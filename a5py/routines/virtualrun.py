"""Defines virtual run which reads data from the C arrays instead of HDF5.
"""
import numpy as np
import unyt
import ctypes

from a5py.ascot5io import State, Orbits, Dist
from a5py.ascot5io.dist import DistData
from .runmixin import RunMixin
from .bbnbi5 import BBNBIMixin

from a5py.ascotpy.libascot import _LIBASCOT
if _LIBASCOT:
    import a5py.ascotpy.ascot2py as ascot2py

class VirtualRun(RunMixin):
    """Virtual :class:`RunGroup` whose data exists solely in the memory.
    """

    def __init__(self, ascot, nmrk, inistate, endstate, options, markers,
                 diag_offload_array, diagorb=None, dist5d=None, dist5drho=None):
        """Initialize fields that allow this instance to replicate
        :class:`RunGroup` behavior.

        Sets state and orbit attributes which are used by the methods in
        :class:`RunMixin`.

        Parameters
        ----------
        ascot : :class:`Ascot`
            Ascot instance for input interpolation.
        nmrk : int
            Number of markers in the simulation.
        inistate :
            Pointer to marker inistate array.
        endstate :
            Pointer to marker endstate array.
        orbits : array_like, optional
            The diagnostics array containing the orbit data, if present.
        dist5d : array_like, optional
            The 5d dist data struct if present.
        dist5drho : array_like, optional
            The diagnostics array containing the 5D rhodist data, if present.
        """
        self.options = options
        self.markers = markers
        # There's a need for better solution here in case inputs are provided
        for inp in ["bfield", "efield", "plasma", "neutral", "wall", "boozer",
                   "mhd", "asigma"]:
            try:
                grp = ascot.data[inp].active
                setattr(self, inp, grp)
            except:
                pass

        self._inistate = VirtualState(ascot, ascot._nmrk, inistate)
        self._endstate = VirtualState(ascot, ascot._nmrk, endstate)

        def pointer_increment(idx):
            """Returns pointer to the diagnostics array on index idx.
            """
            array_at_idx = ctypes.cast(diag_offload_array,
                                       ctypes.POINTER(ctypes.c_double))
            ptr  = ctypes.cast(ctypes.pointer(array_at_idx),
                               ctypes.POINTER(ctypes.c_void_p))
            ptr.contents.value += idx*ctypes.sizeof(ctypes.c_double)
            return array_at_idx

        if diagorb is not None:
            data = pointer_increment(
                ascot._sim.diag_offload_data.offload_diagorb_index)

            self._orbit = VirtualOrbits(
                ascot, ascot._nmrk, diagorb, data)
        if dist5d is not None:
            data = pointer_increment(
                ascot._sim.diag_offload_data.offload_dist5D_index)

            self._dist5d = VirtualDist("5d", dist5d, diag_offload_array)

class VirtualBBNBIRun(BBNBIMixin):
    """Virtual :class:`BBNBIGroup` whose data exists solely in the memory.
    """

    def __init__(self, ascot, nmrk, state, options, diag_offload_array,
                 dist5d=None, dist5drho=None):
        """Initialize fields that allow this instance to replicate
        :class:`BBNBIGroup` behavior.

        Sets state and orbit attributes which are used by the methods in
        :class:`BBNBIMixin`.

        Parameters
        ----------
        ascot : :class:`Ascot`
            Ascot instance for input interpolation.
        nmrk : int
            Number of markers in the simulation.
        state :
            Pointer to marker state array.
        dist5d : array_like, optional
            The 5d dist data struct if present.
        dist5drho : array_like, optional
            The diagnostics array containing the 5D rhodist data, if present.
        """
        self.options = options
        # There's a need for better solution here in case inputs are provided
        #for inp in ["bfield", "efield", "plasma", "neutral", "wall", "boozer",
        #            "mhd", "asigma"]:
        #    grp = ascot.data[inp].active
        #    setattr(self, inp, grp)
        self._state = VirtualState(ascot, ascot._nmrk, state)

        def pointer_increment(idx):
            """Returns pointer to the diagnostics array on index idx.
            """
            array_at_idx = ctypes.cast(diag_offload_array,
                                       ctypes.POINTER(ctypes.c_double))
            ptr  = ctypes.cast(ctypes.pointer(array_at_idx),
                               ctypes.POINTER(ctypes.c_void_p))
            ptr.contents.value += idx*ctypes.sizeof(ctypes.c_double)
            return array_at_idx

        if dist5d is not None:
            data = pointer_increment(
                ascot._sim.diag_offload_data.offload_dist5D_index)

            self._dist5d = VirtualDist("5d", dist5d, diag_offload_array)

class VirtualState():
    """Like :class:`State` but the data is in C array.
    """

    def __init__(self, ascot, nmrk, state):
        """Initialize fields that allow this instance to replicate
        :class:`RunGroup` behavior.

        Parameters
        ----------
        ascot : :class:`Ascot`
            Ascot instance for input interpolation.
        nmrk : int
            Number of markers in the simulation.
        state :
            Pointer to marker state array.
        """
        self.ascot = ascot
        self.state = state
        self.nmrk  = nmrk.value

    def get(self, *qnt, mode="gc"):

        def _val(q):
            """Read quantity from HDF5.
            """
            arr = lambda q : np.array([getattr(self.state[i], q)
                                       for i in range(self.nmrk)])

            m = unyt.m; rad = unyt.rad; kg = unyt.kg; s = unyt.s; J = unyt.J
            T = unyt.T; C = unyt.C; nodim = unyt.dimensionless
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
            elif q == "weight":    return arr(q) * nodim
            elif q == "time":      return arr(q) * s
            elif q == "mileage":   return arr(q) * s
            elif q == "cputime":   return arr(q) * s
            elif q == "rho":       return arr(q) * nodim
            elif q == "theta":     return arr(q) * rad
            elif q == "mass":      return arr(q) * kg
            elif q == "br":        return arr("B_r") * T
            elif q == "bphi":      return arr("B_phi") * T
            elif q == "bz":        return arr("B_z") * T
            elif q == "ids":       return arr(q) * nodim
            elif q == "endcond":   return arr(q) * nodim
            elif q == "walltile":  return arr(q) * nodim
            elif q == "charge":    return arr(q) * C
            elif q == "anum":      return arr(q) * nodim
            elif q == "znum":      return arr(q) * nodim
            elif q == "errormsg":  return arr("err") * nodim
            elif q == "errorline": return arr("err") * nodim
            elif q == "errormod":  return arr("err") * nodim
            return arr

        def _eval(r, phi, z, t, *q):
            """Evaluate input quantities at marker position.
            """
            return self.ascot.input_eval(r, phi, z, t, *q)

        return State._getactual(mode, _val, _eval, *qnt)

class VirtualOrbits():
    """Like :class:`Orbits` but the data is in C array.
    """

    def __init__(self, ascot, nmrk, diagorb, diag_offload_array):
        """Initialize fields that allow this instance to replicate
        :class:`Orbits` behavior.

        Parameters
        ----------
        ascot : :class:`Ascot`
            Ascot instance for input interpolation.
        nmrk : int
            Number of markers in the simulation.
        diagorb :
            Diag_orb offload struct.
        diag_offload_array : array_like
            The offload array where the orbit data begins at the first index.
        """
        self.ascot   = ascot
        self.nmrk    = nmrk.value
        self.npnt    = diagorb.Npnt
        self.mode    = diagorb.record_mode
        self.orbmode = diagorb.mode
        self.arr     = diag_offload_array

    def get(self, inistate, endstate, *qnt):
        """Return marker quantity.

        This function accesses the orbit data within the C array and uses that
        to evaluate the queried quantity. The evaluated quantity at a given
        position corresponds to that mode which was active at the simulation:
        GO simulations return particle phase-space, GC guiding-center
        phase-space and hybrid depends on whether marker was GO or GC at that
        moment when data point was written.

        Parameters
        ----------
        inistate : State
            Inistate is needed to evaluate some orbit quantities.
        endstate : State
            Endstate is needed to evaluate some orbit quantities.
        *qnt : str
            Names of the quantities.

        Returns
        -------
        *value : array_like
            The quantities as an array ordered by marker ID (major) and mileage
            (minor).
        """
        # Prepare helper variables and functions
        arrlen = self.nmrk * self.npnt
        def _val(q, mask=None):
            """Read quantity from the orbit array.
            """
            # C array contains dummy data where ids == 0
            idx = np.array(self.arr[arrlen*0:arrlen*1]) > 0
            m = unyt.m; rad = unyt.rad; kg = unyt.kg; s = unyt.s; J = unyt.J
            T = unyt.T; C = unyt.C; nodim = unyt.dimensionless
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

            if self.mode == ascot2py.simulate_mode_fo:
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
                if q == "pncrid" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*16:arrlen*17]) * nodim
                if q == "pncrdir" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*17:arrlen*18]) * nodim
            elif self.mode == ascot2py.simulate_mode_gc:
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
                if q == "pncrid" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*16:arrlen*17]) * nodim
                if q == "pncrdir" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*17:arrlen*18]) * nodim
            elif self.mode == ascot2py.simulate_mode_ml:
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
                if q == "pncrid" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*11:arrlen*12]) * nodim
                if q == "pncrdir" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*12:arrlen*13]) * nodim
            elif self.mode == ascot2py.simulate_mode_hybrid:
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
                if q == "pncrid" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*19:arrlen*20]) * nodim
                if q == "pncrdir" and self.orbmode == 0:
                    val = np.array(self.arr[arrlen*20:arrlen*21]) * nodim
            return val[idx][mask].ravel()

        # Sort using the fact that inistate.get return values ordered by ID
        # and also np.unique returns indices that produce a sorted array.
        mode    = _val("simmode")
        _, idx  = np.unique(_val("ids").v, return_inverse=True)
        mass    = inistate.get("mass")[0][idx]
        time    = inistate.get("time")[0][idx]
        connlen = inistate.get("mileage")[0][idx] - _val("mileage")

        # Only field lines are constant in time
        if not Orbits.FIELDLINE in mode: time = time + _val("mileage")

        def _eval(q, mask=None):
            """Evaluate input quantities at marker position.
            """
            return self.ascot.input_eval(
                _val("r", mask=mask), _val("phi",  mask=mask),
                _val("z", mask=mask), time[mask], *[q])

        return Orbits._getactual(mass, time, connlen, mode, _val, _eval, *qnt)

class VirtualDist(Dist):
    """Distribution shared by C and Python.
    """

    def __init__(self, disttype, dist, data):
        self._histogram = data
        if disttype == "5d":
            names = ["r", "phi", "z", "ppar", "pperp", "time", "charge"]

        self._abscissa_edges = {}
        for n in names:
            match n:
                case "r":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_r, dist.max_r, dist.n_r+1) * unyt.m
                case "phi":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_phi, dist.max_phi, dist.n_phi+1) \
                            * unyt.rad
                    self._abscissa_edges[n].convert_to_units('deg')
                case "z":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_z, dist.max_z, dist.n_z+1) * unyt.m
                case "ppar":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_ppara, dist.max_ppara,
                                    dist.n_ppara+1) * unyt.kg*unyt.m**2/unyt.s
                case "pperp":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_pperp, dist.max_pperp,
                                    dist.n_pperp+1) * unyt.kg*unyt.m**2/unyt.s
                case "time":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_time, dist.max_time,
                                    dist.n_time+1) * unyt.s
                case "charge":
                    self._abscissa_edges[n] = \
                        np.linspace(dist.min_q, dist.max_q, dist.n_q+1) * unyt.e

    def get(self):
        """Return the distribution data.

        Returns
        -------
        dist : class:`DistData`
            Distribution data.
        """
        arrlen = 1
        dim = np.array([], dtype='i8')
        for d in self._abscissa_edges.values():
            dim = np.append(dim, d.size-1)
            arrlen *= d.size-1

        histogram = np.ctypeslib.as_array(
            (ctypes.c_double * arrlen).from_address(ctypes.addressof(self._histogram.contents)))
        histogram = np.reshape(histogram, dim) * unyt.particles

        return DistData(histogram, **self._abscissa_edges)