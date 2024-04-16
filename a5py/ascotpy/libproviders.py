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
from a5py.ascot5io.coreio import fileapi
from .libascot import _LIBASCOT
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

        ascot2py.B_field_init_offload(
            ctypes.byref(self._sim.B_offload_data),
            ctypes.byref(self._bfield_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_bfield = bytes(qid, "utf-8")

    def _provide_B_TC(self, **kwargs):
        """Initialize :class:`B_TC` straight from dictionary bypassing HDF5.
        """
        BTC = self._sim.B_offload_data.BTC
        BTC.axisr  = kwargs["axisr"][0]
        BTC.axisz  = kwargs["axisz"][0]
        BTC.psival = kwargs["psival"][0]
        BTC.rhoval = kwargs["rhoval"][0]
        BTC.B      = npctypes.as_ctypes(
            np.ascontiguousarray(kwargs["bxyz"].flatten(), dtype="f8") )
        BTC.dB     = npctypes.as_ctypes(
            np.ascontiguousarray(kwargs["jacobian"].flatten(), dtype="f8") )
        BTC.offload_array_length = 0

        self._sim.B_offload_data.type = ascot2py.B_field_type_TC

    def _provide_B_GS(self, **kwargs):
        """Initialize :class:`B_GS` straight from dictionary bypassing HDF5.
        """
        BGS = self._sim.B_offload_data.BGS
        BGS.R0        = kwargs["r0"][0]
        BGS.z0        = kwargs["z0"][0]
        BGS.raxis     = kwargs["raxis"][0]
        BGS.zaxis     = kwargs["zaxis"][0]
        BGS.B_phi0    = kwargs["bphi0"][0]
        BGS.psi0      = kwargs["psi0"][0]
        BGS.psi1      = kwargs["psi1"][0]
        BGS.psi_mult  = kwargs["psimult"][0]
        BGS.Nripple   = int(kwargs["nripple"][0])
        BGS.a0        = kwargs["a0"][0]
        BGS.alpha0    = kwargs["alpha0"][0]
        BGS.delta0    = kwargs["delta0"][0]
        BGS.psi_coeff = npctypes.as_ctypes(
            np.ascontiguousarray(kwargs["coefficients"].flatten(), dtype="f8") )
        BGS.offload_array_length = 0

        self._sim.B_offload_data.type = ascot2py.B_field_type_GS

    def _provide_B_2DS(self, **kwargs):
        """Initialize :class:`B_2DS` from dictionary.
        """
        B2DS = self._sim.B_offload_data.B2DS
        B2DS.n_r    = int(kwargs["nr"][0])
        B2DS.n_z    = int(kwargs["nz"][0])
        B2DS.r_min  = kwargs["rmin"][0]
        B2DS.r_max  = kwargs["rmax"][0]
        B2DS.z_min  = kwargs["zmin"][0]
        B2DS.z_max  = kwargs["zmax"][0]
        B2DS.psi0   = kwargs["psi0"][0]
        B2DS.psi1   = kwargs["psi1"][0]
        B2DS.axis_r = kwargs["axisr"][0]
        B2DS.axis_z = kwargs["axisz"][0]

        B2DS.offload_array_length, self._bfield_offload_array = \
            self._init_offload_array(
                kwargs["psi"], kwargs["br"], kwargs["bphi"], kwargs["bz"]
            )

        self._sim.B_offload_data.type = ascot2py.B_field_type_2DS

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

        B3DS = self._sim.B_offload_data.B3DS
        B3DS.psigrid_n_r   = int(kwargs["psi_nr"][0])
        B3DS.psigrid_n_z   = int(kwargs["psi_nz"][0])
        B3DS.psigrid_r_min = kwargs["psi_rmin"][0]
        B3DS.psigrid_r_max = kwargs["psi_rmax"][0]
        B3DS.psigrid_z_min = kwargs["psi_zmin"][0]
        B3DS.psigrid_z_max = kwargs["psi_zmax"][0]
        B3DS.Bgrid_n_r     = int(kwargs["b_nr"][0])
        B3DS.Bgrid_n_z     = int(kwargs["b_nz"][0])
        B3DS.Bgrid_r_min   = kwargs["b_rmin"][0]
        B3DS.Bgrid_r_max   = kwargs["b_rmax"][0]
        B3DS.Bgrid_z_min   = kwargs["b_zmin"][0]
        B3DS.Bgrid_z_max   = kwargs["b_zmax"][0]
        B3DS.Bgrid_n_phi   = int(kwargs["b_nphi"][0])
        B3DS.Bgrid_phi_min = kwargs["b_phimin"][0] * np.pi / 180
        B3DS.Bgrid_phi_max = kwargs["b_phimax"][0] * np.pi / 180
        B3DS.psi0          = kwargs["psi0"][0]
        B3DS.psi1          = kwargs["psi1"][0]
        B3DS.axis_r        = kwargs["axisr"][0]
        B3DS.axis_z        = kwargs["axisz"][0]

        B3DS.offload_array_length, self._bfield_offload_array = \
            self._init_offload_array(
                kwargs["br"], kwargs["bphi"], kwargs["bz"],kwargs["psi"]
            )

        self._sim.B_offload_data.type = ascot2py.B_field_type_3DS

    def _provide_B_STS(self, **kwargs):
        """Initialize :class:`B_STS` from dictionary.
        """
        BSTS = self._sim.B_offload_data.BSTS

        BSTS.psigrid_n_r     = kwargs['psi_nr'][0]
        BSTS.psigrid_n_z     = kwargs['psi_nz'][0]
        BSTS.psigrid_n_phi   = kwargs['psi_nphi'][0]
        BSTS.psigrid_r_min   = kwargs['psi_rmin'][0]
        BSTS.psigrid_r_max   = kwargs['psi_rmax'][0]
        BSTS.psigrid_z_min   = kwargs['psi_zmin'][0]
        BSTS.psigrid_z_max   = kwargs['psi_zmax'][0]
        BSTS.psigrid_phi_min = kwargs['psi_phimin'][0] * np.pi * 2.0 / 360.0
        BSTS.psigrid_phi_max = kwargs['psi_phimax'][0] * np.pi * 2.0 / 360.0
        BSTS.Bgrid_n_r       = kwargs['b_nr'][0]
        BSTS.Bgrid_n_z       = kwargs['b_nz'][0]
        BSTS.Bgrid_n_phi     = kwargs['b_nphi'][0]
        BSTS.Bgrid_r_min     = kwargs['b_rmin'][0]
        BSTS.Bgrid_r_max     = kwargs['b_rmax'][0]
        BSTS.Bgrid_z_min     = kwargs['b_zmin'][0]
        BSTS.Bgrid_z_max     = kwargs['b_zmax'][0]
        BSTS.Bgrid_phi_min   = kwargs['b_phimin'][0] * np.pi * 2.0 / 360.0
        BSTS.Bgrid_phi_max   = kwargs['b_phimax'][0] * np.pi * 2.0 / 360.0
        BSTS.psi0            = kwargs['psi0'][0]
        BSTS.psi1            = kwargs['psi1'][0]
        BSTS.n_axis          = kwargs['axis_nphi'][0]
        BSTS.axis_min        = kwargs['axis_phimin'][0] * np.pi * 2.0 / 360.0
        BSTS.axis_max        = kwargs['axis_phimax'][0] * np.pi * 2.0 / 360.0

        BSTS.offload_array_length, self._bfield_offload_array = \
            self._init_offload_array(
                kwargs['br'],  kwargs['bphi'],  kwargs['bz'],  kwargs['psi'],
                kwargs['axisr'], kwargs['axisz']
            )

        self._sim.B_offload_data.type = ascot2py.B_field_type_STS

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

        ascot2py.E_field_init_offload(
            ctypes.byref(self._sim.E_offload_data),
            ctypes.byref(self._efield_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_efield = bytes(qid, "utf-8")

    def _provide_E_TC(self, **kwargs):
        """Initialize :class:`E_TC` from dictionary.
        """
        self._sim.E_offload_data.ETC.Exyz[0] = kwargs["exyz"][0]
        self._sim.E_offload_data.ETC.Exyz[1] = kwargs["exyz"][1]
        self._sim.E_offload_data.ETC.Exyz[2] = kwargs["exyz"][2]
        self._sim.E_offload_data.type = ascot2py.E_field_type_TC

    def _provide_E_1DS(self, **kwargs):
        """Initialize :class:`E_1DS` from dictionary.
        """
        E1DS = self._sim.E_offload_data.E1DS
        E1DS.n_rho = int(kwargs["nrho"])
        E1DS.rho_min = kwargs["rhomin"][0,0]
        E1DS.rho_max = kwargs["rhomax"][0,0]
        reff = kwargs["reff"][0]

        E1DS.offload_array_length, self._efield_offload_array = \
            self._init_offload_array(kwargs["dvdrho"])
        self._sim.E_offload_data.type = ascot2py.E_field_type_1DS

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

        ascot2py.plasma_init_offload(
            ctypes.byref(self._sim.plasma_offload_data),
            ctypes.byref(self._plasma_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_plasma = bytes(qid, "utf-8")

    def _provide_plasma_1D(self, **kwargs):
        """Initialize :class:`plasma_1D` from dictionary.
        """
        nion = int(kwargs["nion"])
        P1D = self._sim.plasma_offload_data.plasma_1D
        P1D.n_species = nion + 1
        P1D.n_rho = int(kwargs["nrho"])

        P1D.mass[0]   =  1.0 * unyt.electron_mass
        P1D.charge[0] = -1.0 * unyt.elementary_charge
        for i in range(nion):
            P1D.mass[i+1]   = kwargs["mass"][i,0] * unyt.atomic_mass_unit
            P1D.charge[i+1] = 1.0 * kwargs["charge"][i,0] * unyt.elementary_charge
            P1D.anum[i]     = int(kwargs["anum"][i,0])
            P1D.znum[i]     = int(kwargs["znum"][i,0])

        Te = kwargs["etemperature"] * unyt.elementary_charge
        Ti = kwargs["itemperature"] * unyt.elementary_charge

        P1D.offload_array_length, self._plasma_offload_array = \
            self._init_offload_array(kwargs["rho"], Te, Ti, kwargs["edensity"],
                                     kwargs["idensity"].T)
        self._sim.plasma_offload_data.type = ascot2py.plasma_type_1D

    def _provide_plasma_1DS(self, **kwargs):
        """Initialize :class:`plasma_1DS` from dictionary.
        """
        nion = int(kwargs["nion"])
        P1DS = self._sim.plasma_offload_data.plasma_1DS
        P1DS.n_species = nion + 1
        P1DS.n_rho = int(kwargs["nrho"])
        P1DS.rho_min = kwargs["rhomin"][0,0]
        P1DS.rho_max = kwargs["rhomax"][0,0]

        P1DS.mass[0]   =  1.0 * unyt.electron_mass
        P1DS.charge[0] = -1.0 * unyt.elementary_charge
        for i in range(nion):
            P1DS.mass[i+1]   = kwargs["mass"][i,0] * unyt.atomic_mass_unit
            P1DS.charge[i+1] = 1.0 * kwargs["charge"][i,0] * unyt.elementary_charge
            P1DS.anum[i]     = int(kwargs["anum"][i,0])
            P1DS.znum[i]     = int(kwargs["znum"][i,0])

        Te = kwargs["etemperature"] * unyt.elementary_charge
        Ti = kwargs["itemperature"] * unyt.elementary_charge

        P1DS.offload_array_length, self._plasma_offload_array = \
            self._init_offload_array(Te, Ti, kwargs["edensity"],
                                     kwargs["idensity"].flatten())

        self._sim.plasma_offload_data.type = ascot2py.plasma_type_1DS

    def _provide_plasma_1Dt(self, **kwargs):
        """Initialize :class:`plasma_1Dt` from dictionary.
        """
        nion  = int(kwargs["nion"])
        nrho  = int(kwargs["nrho"])
        ntime = int(kwargs["ntime"])
        nspec = nion + 1
        P1Dt = self._sim.plasma_offload_data.plasma_1Dt
        P1Dt.n_species = nspec
        P1Dt.n_rho  = nrho
        P1Dt.n_time = ntime

        P1Dt.mass[0]   =  1.0 * unyt.electron_mass
        P1Dt.charge[0] = -1.0 * unyt.elementary_charge
        for i in range(nion):
            P1Dt.mass[i+1]   = kwargs["mass"][i,0] * unyt.atomic_mass_unit
            P1Dt.charge[i+1] = 1.0 * kwargs["charge"][i,0] * unyt.elementary_charge
            P1Dt.anum[i]     = int(kwargs["anum"][i,0])
            P1Dt.znum[i]     = int(kwargs["znum"][i,0])

        Te  = kwargs["etemperature"].ravel() * unyt.elementary_charge
        Ti  = kwargs["itemperature"].ravel() * unyt.elementary_charge
        Tei = np.zeros((nrho, ntime, 2)).ravel()
        nei = np.zeros((nrho, ntime, nspec)).ravel()

        # Rearrange the input arrays
        for irho in range(nrho):
            for itime in range(ntime):
                Tei[itime*2*nrho + irho]        = Te[itime*nrho + irho]
                Tei[itime*2*nrho + nrho + irho] = Ti[itime*nrho + irho]

                for ispec in range(nspec):
                    if ispec == 0:
                        nei[itime*nspec*nrho + ispec*nrho + irho] = \
                            kwargs["edensity"][itime,irho]
                    else:
                        nei[itime*nspec*nrho + ispec*nrho + irho] = \
                            kwargs["idensity"][itime,ispec-1,irho]

        P1Dt.offload_array_length, self._plasma_offload_array = \
            self._init_offload_array(kwargs["rho"], kwargs["time"], Tei, nei)
        self._sim.plasma_offload_data.type = ascot2py.plasma_type_1Dt

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

        ascot2py.wall_init_offload(
            ctypes.byref(self._sim.wall_offload_data),
            ctypes.byref(self._wall_offload_array),
            ctypes.byref(self._wall_int_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_wall = bytes(qid, "utf-8")

    def _provide_wall_2D(self, **kwargs):
        """Initialize :class:`wall_2D` from dictionary.
        """
        W2D = self._sim.wall_offload_data.w2d
        W2D.n = kwargs["nelements"]

        W2D.offload_array_length, self._wall_offload_array = \
            self._init_offload_array(kwargs["r"], kwargs["z"])

        self._sim.wall_offload_data.type = ascot2py.wall_type_2D

    def _provide_wall_3D(self, **kwargs):
        """Initialize :class:`wall_3D` from dictionary.
        """
        W3D = self._sim.wall_offload_data.w3d
        try:
            W3D.n = int(kwargs["nelements"].ravel()[0])
        except:
            W3D.n = int(kwargs["nelements"])

        # The data in the offload array is to be in the format
        #  [x1 y1 z1 x2 y2 z2 x3 y3 z3; ... ]
        xyz = np.c_[kwargs["x1x2x3"].flatten(), kwargs["y1y2y3"].flatten(),
                    kwargs["z1z2z3"].flatten()].flatten()

        W3D.offload_array_length, self._wall_offload_array = \
            self._init_offload_array(xyz)

        self._sim.wall_offload_data.type = ascot2py.wall_type_3D

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

        ascot2py.neutral_init_offload(
            ctypes.byref(self._sim.neutral_offload_data),
            ctypes.byref(self._neutral_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_neutral = bytes(qid, "utf-8")

    def _provide_N0_1D(self, **kwargs):
        """Initialize :class:`N0_1D` from dictionary.
        """
        nspecies = int(kwargs["nspecies"][0])
        N1D = self._sim.neutral_offload_data.N01D
        N1D.n_rho     = int(kwargs["nrho"][0])
        N1D.rho_min   = kwargs["rhomin"][0]
        N1D.rho_max   = kwargs["rhomax"][0]
        N1D.n_species = nspecies

        for i in range(nspecies):
            N1D.anum[i]       = int(kwargs["anum"][i])
            N1D.znum[i]       = int(kwargs["znum"][i])
            N1D.maxwellian[i] = int(kwargs["maxwellian"][i])

        N1D.offload_array_length, self._neutral_offload_array = \
            self._init_offload_array(kwargs["density"], kwargs["temperature"])
        self._sim.neutral_offload_data.type = ascot2py.neutral_type_1D

    def _provide_N0_3D(self, **kwargs):
        """Initialize :class:`N0_3D` from dictionary.
        """
        nspecies = int(kwargs["nspecies"][0])
        N3D = self._sim.neutral_offload_data.N03D
        N3D.n_r     = int(kwargs["nr"][0])
        N3D.r_min   = kwargs["rmin"][0]
        N3D.r_max   = kwargs["rmax"][0]
        N3D.n_phi   = int(kwargs["nphi"][0])
        N3D.phi_min = kwargs["phimin"][0] * np.pi / 180
        N3D.phi_max = kwargs["phimax"][0] * np.pi / 180
        N3D.n_z     = int(kwargs["nz"][0])
        N3D.z_min   = kwargs["zmin"][0]
        N3D.z_max   = kwargs["zmax"][0]
        N3D.n_species = nspecies

        for i in range(nspecies):
            N3D.anum[i]       = int(kwargs["anum"][i])
            N3D.znum[i]       = int(kwargs["znum"][i])
            N3D.maxwellian[i] = int(kwargs["maxwellian"][i])

        N3D.offload_array_length, self._neutral_offload_array = \
            self._init_offload_array(kwargs["density"], kwargs["temperature"])
        self._sim.neutral_offload_data.type = ascot2py.neutral_type_3D

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

        ascot2py.boozer_init_offload(
            ctypes.byref(self._sim.boozer_offload_data),
            ctypes.byref(self._boozer_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_boozer = bytes(qid, "utf-8")

    def _provide_Boozer(self, **kwargs):
        """Initialize :class:`Boozer` from dictionary.
        """
        BZ2 = self._sim.boozer_offload_data
        BZ2.npsi    = int(kwargs["npsi"])
        BZ2.psi_min = int(kwargs["psimin"][0])
        BZ2.psi_max = int(kwargs["psimax"][0])
        BZ2.ntheta  = int(kwargs["ntheta"])
        BZ2.nthetag = int(kwargs["nthetag"])
        BZ2.nrzs    = int(kwargs["nrzs"])

        padding = 4
        data = np.copy(kwargs["theta_psithetageom"])
        theta_psithetageom = np.concatenate(
            (data, data[-1,:] + data[1:padding+1,:]) )
        theta_psithetageom = np.concatenate(
            (data[int(BZ2.nthetag-padding-1):-1,:] - data[-1,:],
             theta_psithetageom) )
        BZ2.nthetag += padding*2

        BZ2.offload_array_length, self._boozer_offload_array = \
            self._init_offload_array(kwargs["nu_psitheta"], theta_psithetageom,
                                     kwargs["rs"], kwargs["zs"])

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

        ascot2py.mhd_init_offload(
            ctypes.byref(self._sim.mhd_offload_data),
            ctypes.byref(self._mhd_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_mhd = bytes(qid, "utf-8")

    def _provide_MHD_STAT(self, **kwargs):
        """Initialize :class:`MHD_STAT` from dictionary.
        """
        MHD = self._sim.mhd_offload_data.stat
        MHD.n_modes = int(kwargs["nmode"][0])
        MHD.nrho    = int(kwargs["nrho"][0])
        MHD.rho_min = kwargs["rhomin"][0]
        MHD.rho_max = kwargs["rhomax"][0]

        for i in range(MHD.n_modes):
            MHD.nmode[i]        = int(kwargs["nmodes"][i])
            MHD.mmode[i]        = int(kwargs["mmodes"][i])
            MHD.amplitude_nm[i] = kwargs["amplitude"][i]
            MHD.omega_nm[i]     = kwargs["omega"][i]
            MHD.phase_nm[i]     = kwargs["phase"][i]

        MHD.offload_array_length, self._mhd_offload_array = \
            self._init_offload_array(kwargs["alpha"], kwargs["phi"])
        self._sim.mhd_offload_data.type = ascot2py.mhd_type_stat

    def _provide_MHD_NONSTAT(self, **kwargs):
        """Initialize :class:`MHD_NONSTAT` from dictionary.
        """
        MHD = self._sim.mhd_offload_data.nonstat
        MHD.n_modes = int(kwargs["nmode"][0])
        MHD.nrho    = int(kwargs["nrho"][0])
        MHD.rho_min = kwargs["rhomin"][0]
        MHD.rho_max = kwargs["rhomax"][0]
        MHD.ntime   = int(kwargs["ntime"][0])
        MHD.t_min   = kwargs["tmin"][0]
        MHD.t_max   = kwargs["tmax"][0]

        for i in range(MHD.n_modes):
            MHD.nmode[i]        = int(kwargs["nmodes"][i])
            MHD.mmode[i]        = int(kwargs["mmodes"][i])
            MHD.amplitude_nm[i] = kwargs["amplitude"][i]
            MHD.omega_nm[i]     = kwargs["omega"][i]
            MHD.phase_nm[i]     = kwargs["phase"][i]

        MHD.offload_array_length, self._mhd_offload_array = \
            self._init_offload_array(kwargs["alpha"], kwargs["phi"])
        self._sim.mhd_offload_data.type = ascot2py.mhd_type_nonstat

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

        ascot2py.asigma_init_offload(
            ctypes.byref(self._sim.asigma_offload_data),
            ctypes.byref(self._asigma_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_asigma = bytes(qid, "utf-8")

    def _provide_asigma_loc(self, **kwargs):
        """Initialize :class:`Asigma_loc` from dictionary.
        """
        nreac = int(kwargs["nreac"])
        LOC = self._sim.asigma_offload_data.asigma_loc
        LOC.N_reac = nreac
        for i in range(nreac):
            LOC.N_E[i] = int(kwargs["nenergy"][i])
            LOC.N_n[i] = int(kwargs["ndensity"][i])
            LOC.N_T[i] = int(kwargs["ntemperature"][i])
            LOC.z_1[i] = int(kwargs["z1"][i])
            LOC.a_1[i] = int(kwargs["a1"][i])
            LOC.z_2[i] = int(kwargs["z2"][i])
            LOC.a_2[i] = int(kwargs["a2"][i])

        LOC.offload_array_length, self._asigma_offload_array = \
            self._init_offload_array(
                kwargs["energymin"], kwargs["energymax"],
                kwargs["densitymin"], kwargs["densitymax"],
                kwargs["temperaturemin"], kwargs["temperaturemax"],
                kwargs["sigma"])
        self._sim.asigma_offload_data.type = ascot2py.asigma_type_loc

    def _provide_nbi(self, **kwargs):
        """Use the provided input parameters to initialize a NBI input.

        Parameters
        ----------
        **kwargs
            Dictionary with the NBI data.
        """
        inp, data = self._find_input_based_on_kwargs(
            ["NBI"], **kwargs)
        getattr(self, "_provide_" + inp)(**data)

        ascot2py.nbi_init_offload(
            ctypes.byref(self._sim.nbi_offload_data),
            ctypes.byref(self._nbi_offload_array)
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_nbi = bytes(qid, "utf-8")

    def _provide_NBI(self, **kwargs):
        """Initialize :class:`NBI` from dictionary.
        """
        #NBI. = kwargs[""]
        ascot2py.nbi_init_offload(
            ctypes.byref(self._sim.plasma_offload_data),
            self._plasma_offload_array
            )
