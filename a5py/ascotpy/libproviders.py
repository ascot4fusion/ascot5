"""Methods to inject input dependencies to the C-structures directly
from python.

This is in contrast to reading inputs from HDF5 files.
"""
import copy
import inspect
import ctypes
import unyt
import numpy as np
import numpy.ctypeslib as npctypes

import a5py.physlib.analyticequilibrium as psifun

from a5py.ascot5io.coreio import fileapi
from .libascot import _LIBASCOT, PTR_ARR, PTR_INT
if _LIBASCOT:
    from a5py.ascotpy import ascot2py

class LibProviders():
    """Mixin class to provide dependency injectors/providers.
    """

    def _find_input_based_on_kwargs(self, inputs, **kwargs):
        """Find the input type that corresponds to the given parameters.

        This function looks through all ``inputs`` ``write_hdf5`` methods and
        tries to find the one that takes the provided ``**kwargs`` as input
        parameters.

        Parameters
        ----------
        inputs : [str]
            List of possible input types e.g. all bfield inputs.
        **kwargs
            Parameters provided to the input's ``write_hdf5`` function excluding
            ``fn`` and ``desc``.

        Returns
        -------
        inp : str
            Name of the input type that matched.
        data : dict
            Same as ``**kwargs`` but with all the optional arguments included
            with default values if those were not present.

        Raises
        ------
        ValueError
            If no match was found.
        """
        from a5py.ascot5io import HDF5TOOBJ
        leastmissing = [None] * 100
        leastunknown = [None] * 100
        bestmatch = None

        # Loop to determine what magnetic field input is provided
        for inp in inputs:

            # Get names of all arguments and which of those are optional
            argspec = inspect.getfullargspec(HDF5TOOBJ[inp].write_hdf5)
            args = list(argspec.args)

            args.remove("fn")   # Remove filename argument
            args.remove("desc") # Remove description which is optional argument

            # In args, the last len(defaults) ar the optional arguments. Note -1
            # here since we have removed one optional argument (desc)
            required = args[:len(args)-(len(argspec.defaults) - 1)]

            # Check that kwargs contain all required arguments. It may also
            # contain optional arguments but nothing else
            missing = []; unknown = []
            for k in required:
                if k not in kwargs:
                    missing.append(k)
            for k in kwargs:
                if k not in args:
                    unknown.append(k)
            if len(missing) > 0:
                if len(missing) < len(leastmissing):
                    leastmissing = missing
                    bestmatch = inp
                continue
            leastmissing = []
            bestmatch = inp
            if len(unknown) > 0:
                if len(unknown) < len(leastunknown):
                    leastunknown = unknown
                    bestmatch = inp
                continue

            # Everything checks out; this is our input type. Add possible
            # missing optional arguments
            i = 0
            data = copy.copy(kwargs)
            for k in args:
                if k in required: continue
                if k not in kwargs.keys():
                    data[k] = argspec.defaults[i]
                i += 1

            # Make sure all data is given as numpy arrays
            for k in data.keys():
                if not isinstance(data[k], np.ndarray):
                    if isinstance(data[k], list):
                        data[k] = np.array(data[k])
                    else:
                        data[k] = np.array([data[k]])

            return inp, data

        # Matching input was not found. Produce error message that shows
        # what input was the best match
        if len(leastmissing) > 0:
            missing = ""
            for k in leastmissing:
                missing += "'%s', " % k
            raise ValueError(
                ("Input was not recognized. Best match was %s " +
                 "but you did not provide parameters %s")
                % (bestmatch, missing[:-2]))
        else:
            unknown = ""
            for k in leastunknown:
                unknown += "'%s', " % k
            raise ValueError(
                ("Input was not recognized. Best match was %s " +
                 "but you provided unknown parameters %s") %
                (bestmatch, unknown[:-2]))

    def _init_offload_array(self, *args, intdata=False):
        """Convert and combine numpy arrays into ctypes offload data array.

        This is a helper function to initialize input data from dictionaries.

        Parameters
        ----------
        args
            Data arrays in the order they are appended to the offload array.
        intdata : bool, optional
            Whether the double or int offload array is assumed.

        Returns
        -------
        offload_array_length : int
            Size of the offload array.
        offload_array : array_like
            Allocated and initialized offload array.
        """
        size = 0
        for a in args:
            size += a.size
        offload_array_length = size

        offload_array = ascot2py.libascot_allocate_reals(size)
        if intdata:
            array = ctypes.cast(offload_array,
                                ctypes.POINTER(ctypes.c_int*size) )[0]
            offload_array = ctypes.cast(offload_array,
                        ctypes.POINTER(ctypes.c_int) )
        else:
            array = ctypes.cast(offload_array,
                                ctypes.POINTER(ctypes.c_double*size) )[0]

        pos = 0
        for a in args:
            n = a.size
            array[(pos):(pos+n)] = a.flatten(order='F')
            pos += n

        return offload_array_length, offload_array

    def _provide_bfield(self, **kwargs):
        """Use the provided input parameters to initialize a magnetic field
        input.

        Parameters
        ----------
        **kwargs
            Dictionary with the magnetic field data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["B_TC", "B_GS", "B_2DS", "B_3DS", "B_STS"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_bfield = bytes(qid, "utf-8")

    def _provide_B_TC(self, **kwargs):
        """Initialize :class:`B_TC` straight from dictionary bypassing HDF5.
        """
        _LIBASCOT.B_TC_init(
            ctypes.byref(self._sim.B_data.BTC), kwargs["axisr"][0],
            kwargs["axisz"][0], kwargs["psival"][0], kwargs["rhoval"][0],
            (ctypes.c_double * 3)(*kwargs["bxyz"].flatten()),
            (ctypes.c_double * 9)(*kwargs["jacobian"].T.flatten()),
            )
        self._sim.B_data.type = ascot2py.B_field_type_TC

    def _provide_B_GS(self, **kwargs):
        """Initialize :class:`B_GS` straight from dictionary bypassing HDF5.
        """
        # Search for magnetic axis psi
        if kwargs["psi0"][0] is None:
            r0 = kwargs["r0"][0]
            z0 = kwargs["z0"][0]
            c = kwargs["coefficients"]
            x = psifun.find_axis(r0, z0, c[0], c[1], c[2], c[3], c[4], c[5],
                                 c[6], c[7], c[8], c[9], c[10], c[11], c[12])
            psi0 = psifun.psi0(x[0], x[1], c[0], c[1], c[2], c[3], c[4],
                               c[5], c[6], c[7], c[8], c[9], c[10], c[11],
                               c[12]) * kwargs["psimult"][0]
            kwargs["psi1"]  = np.array([0])
            kwargs["raxis"] = np.array([x[0]*r0])
            kwargs["zaxis"] = np.array([x[1]*r0])
            if psi0 < kwargs["psi1"][0]:
                kwargs["psi0"] = np.array([psi0 - 1e-8])
            else:
                kwargs["psi0"] = np.array([psi0 + 1e-8])

        coeff = (ctypes.c_double * 14)(*kwargs["coefficients"].flatten())
        _LIBASCOT.B_GS_init(
            ctypes.byref(self._sim.B_data.BGS), kwargs["r0"][0], kwargs["z0"][0],
            kwargs["raxis"][0], kwargs["zaxis"][0], kwargs["bphi0"][0],
            kwargs["psi0"][0], kwargs["psi1"][0], kwargs["psimult"][0],
            coeff, int(kwargs["nripple"][0]), kwargs["a0"][0],
            kwargs["alpha0"][0], kwargs["delta0"][0]
            )
        self._sim.B_data.type = ascot2py.B_field_type_GS

    def _provide_B_2DS(self, **kwargs):
        """Initialize :class:`B_2DS` from dictionary.
        """
        _LIBASCOT.B_2DS_init(
            ctypes.byref(self._sim.B_data.B2DS),
            int(kwargs["nr"][0]), kwargs["rmin"][0], kwargs["rmax"][0],
            int(kwargs["nz"][0]), kwargs["zmin"][0], kwargs["zmax"][0],
            kwargs["axisr"][0], kwargs["axisz"][0], kwargs["psi0"][0],
            kwargs["psi1"][0],
            kwargs["psi"].ctypes.data_as(PTR_ARR),
            kwargs["br"].ctypes.data_as(PTR_ARR),
            kwargs["bphi"].ctypes.data_as(PTR_ARR),
            kwargs["bz"].ctypes.data_as(PTR_ARR)
            )
        self._sim.B_data.type = ascot2py.B_field_type_2DS

    def _provide_B_3DS(self, **kwargs):
        """Initialize :class:`B_3DS` from dictionary.
        """
        if(kwargs["psi_rmin"] is None or kwargs["psi_rmax"] is None or
           kwargs["psi_nr"] is None or kwargs["psi_zmin"] is None or
           kwargs["psi_zmax"] is None or kwargs["psi_nz"] is None):
            kwargs["psi_rmin"] = kwargs["b_rmin"]
            kwargs["psi_rmax"] = kwargs["b_rmax"]
            kwargs["psi_nr"]   = kwargs["b_nr"]
            kwargs["psi_zmin"] = kwargs["b_zmin"]
            kwargs["psi_zmax"] = kwargs["b_zmax"]
            kwargs["psi_nz"]   = kwargs["b_nz"]

        _LIBASCOT.B_3DS_init(
            ctypes.byref(self._sim.B_data.B3DS),
            int(kwargs["psi_nr"][0]), kwargs["psi_rmin"][0],
            kwargs["psi_rmax"][0],
            int(kwargs["psi_nz"][0]), kwargs["psi_zmin"][0],
            kwargs["psi_zmax"][0],
            int(kwargs["b_nr"][0]), kwargs["b_rmin"][0], kwargs["b_rmax"][0],
            int(kwargs["b_nphi"][0]), kwargs["b_phimin"][0] * np.pi / 180,
            kwargs["b_phimax"][0] * np.pi / 180,
            int(kwargs["b_nz"][0]), kwargs["b_zmin"][0], kwargs["b_zmax"][0],
            kwargs["axisr"][0], kwargs["axisz"][0], kwargs["psi0"][0],
            kwargs["psi1"][0],
            kwargs["psi"].ctypes.data_as(PTR_ARR),
            kwargs["br"].ctypes.data_as(PTR_ARR),
            kwargs["bphi"].ctypes.data_as(PTR_ARR),
            kwargs["bz"].ctypes.data_as(PTR_ARR)
        )
        self._sim.B_data.type = ascot2py.B_field_type_3DS

    def _provide_B_STS(self, **kwargs):
        """Initialize :class:`B_STS` from dictionary.
        """
        _LIBASCOT.B_STS_init(
            ctypes.byref(self._sim.B_data.BSTS),
            int(kwargs["psi_nr"][0]), kwargs["psi_rmin"][0],
            kwargs["psi_rmax"][0],
            int(kwargs["psi_nphi"][0]), kwargs["psi_phimin"][0] * np.pi / 180,
            kwargs["psi_phimax"][0] * np.pi / 180,
            int(kwargs["psi_nz"][0]), kwargs["psi_zmin"][0],
            kwargs["psi_zmax"][0],
            int(kwargs["b_nr"][0]), kwargs["b_rmin"][0], kwargs["b_rmax"][0],
            int(kwargs["b_nphi"][0]), kwargs["b_phimin"][0] * np.pi / 180,
            kwargs["b_phimax"][0] * np.pi / 180,
            int(kwargs["b_nz"][0]), kwargs["b_zmin"][0], kwargs["b_zmax"][0],
            int(kwargs['axis_nphi'][0]), kwargs["axis_phimin"][0] * np.pi / 180,
            kwargs["axis_phimax"][0] * np.pi / 180,
            kwargs['axisr'].ctypes.data_as(PTR_ARR),
            kwargs['axisz'].ctypes.data_as(PTR_ARR),
            kwargs["psi0"][0], kwargs["psi1"][0],
            kwargs["psi"].ctypes.data_as(PTR_ARR),
            kwargs["br"].ctypes.data_as(PTR_ARR),
            kwargs["bphi"].ctypes.data_as(PTR_ARR),
            kwargs["bz"].ctypes.data_as(PTR_ARR)
            )
        self._sim.B_data.type = ascot2py.B_field_type_STS

    def _provide_efield(self, **kwargs):
        """Use the provided input parameters to initialize an electric field
        input.

        Parameters
        ----------
        **kwargs
            Dictionary with the electric field data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["E_TC", "E_1DS"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)

        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_efield = bytes(qid, "utf-8")

    def _provide_E_TC(self, **kwargs):
        """Initialize :class:`E_TC` from dictionary.
        """
        exyz = (ctypes.c_double * 3)(*kwargs["exyz"].flatten())
        _LIBASCOT.E_TC_init(ctypes.byref(self._sim.E_data.ETC), exyz)
        self._sim.E_data.type = ascot2py.E_field_type_TC

    def _provide_E_1DS(self, **kwargs):
        """Initialize :class:`E_1DS` from dictionary.
        """
        _LIBASCOT.E_1DS_init(
            ctypes.byref(self._sim.E_data.E1DS), int(kwargs["nrho"]),
            kwargs["rhomin"][0,0], kwargs["rhomax"][0,0], kwargs["reff"][0],
            kwargs["dvdrho"].flatten().ctypes.data_as(PTR_ARR))
        self._sim.E_data.type = ascot2py.E_field_type_1DS

    def _provide_plasma(self, **kwargs):
        """Use the provided input parameters to initialize a plasma input.

        Parameters
        ----------
        **kwargs
            Dictionary with the plasma data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["plasma_1D", "plasma_1DS", "plasma_1Dt"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_plasma = bytes(qid, "utf-8")

    def _provide_plasma_1D(self, **kwargs):
        """Initialize :class:`plasma_1D` from dictionary.
        """
        Te = kwargs["etemperature"] * unyt.elementary_charge
        Ti = kwargs["itemperature"] * unyt.elementary_charge
        mass = np.append(
            (1.0 * unyt.electron_mass),
             (kwargs["mass"].flatten() * unyt.atomic_mass_unit).to("kg")
            )
        charge = np.append(
            (-1.0 * unyt.elementary_charge),
             (kwargs["charge"].flatten() * unyt.elementary_charge).to("C")
            )
        anum = np.ascontiguousarray(kwargs["anum"], dtype=np.int32)
        znum = np.ascontiguousarray(kwargs["znum"], dtype=np.int32)
        _LIBASCOT.plasma_1D_init(
            ctypes.byref(self._sim.plasma_data.plasma_1D), int(kwargs["nrho"]),
            int(kwargs["nion"]), kwargs["rho"].ctypes.data_as(PTR_ARR),
            anum.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            znum.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            mass.ctypes.data_as(PTR_ARR), charge.ctypes.data_as(PTR_ARR),
            Te.ctypes.data_as(PTR_ARR), Ti.ctypes.data_as(PTR_ARR),
            kwargs["edensity"].ctypes.data_as(PTR_ARR),
            kwargs["idensity"].T.ctypes.data_as(PTR_ARR),
            kwargs["vtor"].ctypes.data_as(PTR_ARR),
            )
        self._sim.plasma_data.type = ascot2py.plasma_type_1D

    def _provide_plasma_1DS(self, **kwargs):
        """Initialize :class:`plasma_1DS` from dictionary.
        """
        Te = kwargs["etemperature"] * unyt.elementary_charge
        Ti = kwargs["itemperature"] * unyt.elementary_charge
        mass = np.append(
            (1.0 * unyt.electron_mass),
             (kwargs["mass"].flatten() * unyt.atomic_mass_unit).to("kg")
            )
        charge = np.append(
            (-1.0 * unyt.elementary_charge),
             (kwargs["charge"].flatten() * unyt.elementary_charge).to("C")
            )
        _LIBASCOT.plasma_1DS_init(
            ctypes.byref(self._sim.plasma_data.plasma_1DS), int(kwargs["nrho"]),
            kwargs["rhomin"][0,0], kwargs["rhomax"][0,0], int(kwargs["nion"]),
            kwargs["anum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["znum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            mass.ctypes.data_as(PTR_ARR), charge.ctypes.data_as(PTR_ARR),
            Te.ctypes.data_as(PTR_ARR), Ti.ctypes.data_as(PTR_ARR),
            kwargs["edensity"].ctypes.data_as(PTR_ARR),
            kwargs["idensity"].T.ctypes.data_as(PTR_ARR),
            kwargs["vtor"].ctypes.data_as(PTR_ARR),
        )
        self._sim.plasma_data.type = ascot2py.plasma_type_1DS

    def _provide_plasma_1Dt(self, **kwargs):
        """Initialize :class:`plasma_1Dt` from dictionary.
        """
        Te  = kwargs["etemperature"].ravel() * unyt.elementary_charge
        Ti  = kwargs["itemperature"].ravel() * unyt.elementary_charge
        mass = np.append(
            (1.0 * unyt.electron_mass),
             (kwargs["mass"].flatten() * unyt.atomic_mass_unit).to("kg")
            )
        charge = np.append(
            (-1.0 * unyt.elementary_charge),
             (kwargs["charge"].flatten() * unyt.elementary_charge).to("C")
            )
        _LIBASCOT.plasma_1Dt_init(
            ctypes.byref(self._sim.plasma_data.plasma_1Dt),
            int(kwargs["nrho"]), int(kwargs["ntime"]), int(kwargs["nion"]),
            kwargs["rho"].ctypes.data_as(PTR_ARR),
            kwargs["time"].ctypes.data_as(PTR_ARR),
            kwargs["anum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["znum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            mass.ctypes.data_as(PTR_ARR), charge.ctypes.data_as(PTR_ARR),
            Te.ctypes.data_as(PTR_ARR), Ti.ctypes.data_as(PTR_ARR),
            kwargs["edensity"].ctypes.data_as(PTR_ARR),
            kwargs["idensity"].T.ctypes.data_as(PTR_ARR),
            kwargs["vtor"].ctypes.data_as(PTR_ARR),
        )
        self._sim.plasma_data.type = ascot2py.plasma_type_1Dt

    def _provide_wall(self, **kwargs):
        """Use the provided input parameters to initialize a wall input.

        Parameters
        ----------
        **kwargs
            Dictionary with the wall data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["wall_2D", "wall_3D"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_wall = bytes(qid, "utf-8")

    def _provide_wall_2D(self, **kwargs):
        """Initialize :class:`wall_2D` from dictionary.
        """
        _LIBASCOT.wall_2d_init(
            ctypes.byref(self._sim.wall_data.w2d), int(kwargs["nelements"][0]),
            kwargs["r"].ctypes.data_as(PTR_ARR),
            kwargs["z"].ctypes.data_as(PTR_ARR),
            kwargs["flag"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
        )
        self._sim.wall_data.type = ascot2py.wall_type_2D

    def _provide_wall_3D(self, **kwargs):
        """Initialize :class:`wall_3D` from dictionary.
        """
        try:
            n = int(kwargs["nelements"].ravel()[0])
        except:
            n = int(kwargs["nelements"])
        _LIBASCOT.wall_3d_init(
            ctypes.byref(self._sim.wall_data.w3d), n,
            kwargs["x1x2x3"].ctypes.data_as(PTR_ARR),
            kwargs["y1y2y3"].ctypes.data_as(PTR_ARR),
            kwargs["z1z2z3"].ctypes.data_as(PTR_ARR),
            kwargs["flag"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            )
        self._sim.wall_data.type = ascot2py.wall_type_3D

    def _provide_neutral(self, **kwargs):
        """Use the provided input parameters to initialize a neutral input.

        Parameters
        ----------
        **kwargs
            Dictionary with the neutral data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["N0_1D", "N0_3D"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_neutral = bytes(qid, "utf-8")

    def _provide_N0_1D(self, **kwargs):
        """Initialize :class:`N0_1D` from dictionary.
        """
        _LIBASCOT.N0_1D_init(
            ctypes.byref(self._sim.neutral_data.N01D), int(kwargs["nrho"][0]),
            kwargs["rhomin"][0], kwargs["rhomax"][0],
            int(kwargs["nspecies"][0]),
            kwargs["anum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["znum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["maxwellian"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["density"].ctypes.data_as(PTR_ARR),
            kwargs["temperature"].ctypes.data_as(PTR_ARR)
            )
        self._sim.neutral_data.type = ascot2py.neutral_type_1D

    def _provide_N0_3D(self, **kwargs):
        """Initialize :class:`N0_3D` from dictionary.
        """
        _LIBASCOT.N0_3D_init(
            ctypes.byref(self._sim.neutral_data.N03D),
            int(kwargs["nr"][0]), kwargs["rmin"][0], kwargs["rmax"][0],
            int(kwargs["nphi"][0]), kwargs["phimin"][0] * np.pi / 180,
            kwargs["phimax"][0] * np.pi / 180,
            int(kwargs["nz"][0]), kwargs["zmin"][0], kwargs["zmax"][0],
            int(kwargs["nspecies"][0]),
            kwargs["anum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["znum"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["maxwellian"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["density"].ctypes.data_as(PTR_ARR),
            kwargs["temperature"].ctypes.data_as(PTR_ARR)
            )
        self._sim.neutral_data.type = ascot2py.neutral_type_3D

    def _provide_boozer(self, **kwargs):
        """Use the provided input parameters to initialize a boozer input.

        Parameters
        ----------
        **kwargs
            Dictionary with the boozer data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["Boozer"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_boozer = bytes(qid, "utf-8")

    def _provide_Boozer(self, **kwargs):
        """Initialize :class:`Boozer` from dictionary.
        """
        nthetag = int(kwargs["nthetag"])
        padding = 4
        data = np.copy(kwargs["theta_psithetageom"])
        theta_psithetageom = np.concatenate(
            (data, data[-1,:] + data[1:padding+1,:]) )
        theta_psithetageom = np.concatenate(
            (data[int(nthetag-padding-1):-1,:] - data[-1,:],
             theta_psithetageom) )
        nthetag += padding*2
        _LIBASCOT.boozer_init(
            ctypes.byref(self._sim.boozer_data), int(kwargs["npsi"]),
            kwargs["psimin"][0], kwargs["psimax"][0],
            int(kwargs["ntheta"]), nthetag,
            kwargs["nu_psitheta"].ctypes.data_as(PTR_ARR),
            theta_psithetageom.ctypes.data_as(PTR_ARR), int(kwargs["nrzs"]),
            kwargs["rs"].ctypes.data_as(PTR_ARR),
            kwargs["zs"].ctypes.data_as(PTR_ARR)
            )

    def _provide_mhd(self, **kwargs):
        """Use the provided input parameters to initialize a plasma input.

        Parameters
        ----------
        **kwargs
            Dictionary with the plasma data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["MHD_STAT", "MHD_NONSTAT"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_mhd = bytes(qid, "utf-8")

    def _provide_MHD_STAT(self, **kwargs):
        """Initialize :class:`MHD_STAT` from dictionary.
        """
        nmodes = np.ascontiguousarray(kwargs["nmodes"], dtype=np.int32)
        mmodes = np.ascontiguousarray(kwargs["mmodes"], dtype=np.int32)
        _LIBASCOT.mhd_stat_init(
            ctypes.byref(self._sim.mhd_data.stat), int(kwargs["nmode"][0]),
            int(kwargs["nrho"][0]), kwargs["rhomin"][0], kwargs["rhomax"][0],
            nmodes.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            mmodes.ctypes.data_as(ctypes.POINTER(ctypes.c_int32)),
            kwargs["amplitude"].ctypes.data_as(PTR_ARR),
            kwargs["omega"].ctypes.data_as(PTR_ARR),
            kwargs["phase"].ctypes.data_as(PTR_ARR),
            kwargs["alpha"].ctypes.data_as(PTR_ARR),
            kwargs["phi"].ctypes.data_as(PTR_ARR),
            )
        self._sim.mhd_data.type = ascot2py.mhd_type_stat

    def _provide_MHD_NONSTAT(self, **kwargs):
        """Initialize :class:`MHD_NONSTAT` from dictionary.
        """
        _LIBASCOT.mhd_nonstat_init(
            ctypes.byref(self._sim.mhd_data.nonstat), int(kwargs["nmode"][0]),
            int(kwargs["nrho"][0]), int(kwargs["ntime"][0]),
            kwargs["rhomin"][0], kwargs["rhomax"][0],
            kwargs["tmin"][0], kwargs["tmax"][0],
            kwargs["nmodes"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["mmodes"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["amplitude"].ctypes.data_as(PTR_ARR),
            kwargs["omega"].ctypes.data_as(PTR_ARR),
            kwargs["phase"].ctypes.data_as(PTR_ARR),
            kwargs["alpha"].ctypes.data_as(PTR_ARR),
            kwargs["phi"].ctypes.data_as(PTR_ARR),
            )
        self._sim.mhd_data.type = ascot2py.mhd_type_nonstat

    def _provide_asigma(self, **kwargs):
        """Use the provided input parameters to initialize an atomic data input.

        Parameters
        ----------
        **kwargs
            Dictionary with the atomic data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["asigma_loc"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_asigma = bytes(qid, "utf-8")

    def _provide_asigma_loc(self, **kwargs):
        """Initialize :class:`Asigma_loc` from dictionary.
        """
        _LIBASCOT.asigma_loc_init(
            ctypes.byref(self._sim.asigma_data.asigma_loc),
            int(kwargs["nreac"]),
            kwargs["z1"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["a1"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["z2"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["a2"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["reactype"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["nenergy"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["energymin"].ctypes.data_as(PTR_ARR),
            kwargs["energymax"].ctypes.data_as(PTR_ARR),
            kwargs["ndensity"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["densitymin"].ctypes.data_as(PTR_ARR),
            kwargs["densitymax"].ctypes.data_as(PTR_ARR),
            kwargs["ntemperature"].ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            kwargs["temperaturemin"].ctypes.data_as(PTR_ARR),
            kwargs["temperaturemax"].ctypes.data_as(PTR_ARR),
            kwargs["sigma"].ctypes.data_as(PTR_ARR),
            )
        self._sim.asigma_data.type = ascot2py.asigma_type_loc

    def _provide_nbi(self, **kwargs):
        """Use the provided input parameters to initialize a NBI input.

        Parameters
        ----------
        **kwargs
            Dictionary with the NBI data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["nbi"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)
        qid, _, _ = fileapi._generate_meta()
        self._sim.qid_nbi = bytes(qid, "utf-8")

    def _provide_NBI(self, **kwargs):
        """Initialize :class:`NBI` from dictionary.
        """
        _LIBASCOT_nbi_init()
