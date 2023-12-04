"""Methods to inject input dependencies to the C-structures directly
from python.

This is in contrast to reading inputs from HDF5 files.
"""
import copy
import inspect
import ctypes
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

    def _provide_wall_2d(self,r,z):

        if( len(r) != len(z) ):
            raise ValueError("R and z must be of equal length.")

        nelements = len(r)

        # Create temporary variable for the wall
        R_c = (ctypes.c_double * nelements)()
        z_c = (ctypes.c_double * nelements)()

        for i in range(nelements):
            R_c[i] = r[i]
            z_c[i] = z[i]

        ascot2py.hdf5_wall_2d_to_offload(
            ctypes.byref(self._sim.wall_offload_data.w2d),
            ctypes.byref(self._wall_offload_array),
            nelements,
            R_c, z_c
            )

        self._sim.wall_offload_data.type = ascot2py.wall_type_2D

        ascot2py.wall_init_offload(
            ctypes.byref(self._sim.wall_offload_data),
            self._wall_offload_array,
            self._wall_int_offload_array
            )

    def _provide_wall_3d(self,x1x2x3,y1y2y3,z1z2z3):

        nelements = int(x1x2x3.shape[0])

        # Create temporary variable for the wall
        x1x2x3_c = (ctypes.c_double * (3*nelements) )()
        y1y2y3_c = (ctypes.c_double * (3*nelements) )()
        z1z2z3_c = (ctypes.c_double * (3*nelements) )()

        # Let's hope the ordering is correct...
        x1x2x3_c[:] = x1x2x3.flatten()[:]
        y1y2y3_c[:] = y1y2y3.flatten()[:]
        z1z2z3_c[:] = z1z2z3.flatten()[:]


        ascot2py.hdf5_wall_3d_to_offload(
            ctypes.byref(self._sim.wall_offload_data.w3d),
            ctypes.byref(self._wall_offload_array),
            nelements,
            x1x2x3_c, y1y2y3_c, z1z2z3_c,
            )

        self._sim.wall_offload_data.type = ascot2py.wall_type_3D

        ascot2py.wall_init_offload(
            ctypes.byref(self._sim.wall_offload_data),
            self._wall_offload_array,
            self._wall_int_offload_array
            )

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
            self._bfield_offload_array
        )

        qid, _, _, _ = fileapi._generate_meta()
        self._sim.qid_bfield = bytes(qid, "utf-8")

    def _provide_B_GS(self, **kwargs):
        """Initialize :class:`B_GS` straight from dictionary bypassing HDF5.
        """
        BGS = self._sim.B_offload_data.BGS
        BGS.R0        = kwargs["r0"]
        BGS.z0        = kwargs["z0"]
        BGS.raxis     = kwargs["raxis"]
        BGS.zaxis     = kwargs["zaxis"]
        BGS.B_phi0    = kwargs["bphi0"]
        BGS.psi0      = kwargs["psi0"]
        BGS.psi1      = kwargs["psi1"]
        BGS.psi_mult  = kwargs["psimult"]
        BGS.Nripple   = int(kwargs["nripple"])
        BGS.a0        = kwargs["a0"]
        BGS.alpha0    = kwargs["alpha0"]
        BGS.delta0    = kwargs["delta0"]
        BGS.psi_coeff = npctypes.as_ctypes(
            np.ascontiguousarray(kwargs["coefficients"].flatten(), dtype="f8") )
        BGS.offload_array_length = 0

        self._sim.B_offload_data.type = ascot2py.B_field_type_GS

    def _provide_BSTS(self,
                     b_rmin, b_rmax, b_nr, b_zmin, b_zmax, b_nz,
                     b_phimin, b_phimax, b_nphi, psi0, psi1,
                     br, bphi, bz, psi,
                     axis_phimin, axis_phimax, axis_nphi, axisr, axisz,
                     psi_rmin=None,   psi_rmax=None, psi_nr=None,
                     psi_zmin=None,   psi_zmax=None, psi_nz=None,
                     psi_phimin=None, psi_phimax=None, psi_nphi=None ):
        # bsts is the dictionary that comes from reading the hdf5
        bsts=locals()

        # 1. First fill in the meta-data struct
        #--------------------------------------
        # Mimic the C-function    hdf5_bfield_read_STS()

        #phimin/max deg2rad

        #B_STS_offload_data
        #sts = struct_c__SA_B_STS_offload_data()
        BSTS = self._sim.B_offload_data.BSTS

        BSTS.psigrid_n_r     = bsts['psi_nr'][0]
        BSTS.psigrid_n_z     = bsts['psi_nz'][0]
        BSTS.psigrid_n_phi   = bsts['psi_nphi'][0]
        BSTS.psigrid_r_min   = bsts['psi_rmin'][0]
        BSTS.psigrid_r_max   = bsts['psi_rmax'][0]
        BSTS.psigrid_z_min   = bsts['psi_zmin'][0]
        BSTS.psigrid_z_max   = bsts['psi_zmax'][0]
        BSTS.psigrid_phi_min = bsts['psi_phimin'][0] * np.pi * 2.0 / 360.0
        BSTS.psigrid_phi_max = bsts['psi_phimax'][0] * np.pi * 2.0 / 360.0

        BSTS.Bgrid_n_r       = bsts['b_nr'][0]
        BSTS.Bgrid_n_z       = bsts['b_nz'][0]
        BSTS.Bgrid_n_phi     = bsts['b_nphi'][0]
        BSTS.Bgrid_r_min     = bsts['b_rmin'][0]
        BSTS.Bgrid_r_max     = bsts['b_rmax'][0]
        BSTS.Bgrid_z_min     = bsts['b_zmin'][0]
        BSTS.Bgrid_z_max     = bsts['b_zmax'][0]
        BSTS.Bgrid_phi_min   = bsts['b_phimin'][0] * np.pi * 2.0 / 360.0
        BSTS.Bgrid_phi_max   = bsts['b_phimax'][0] * np.pi * 2.0 / 360.0

        BSTS.psi0            = bsts['psi0'][0]
        BSTS.psi1            = bsts['psi1'][0]


        BSTS.n_axis          = bsts['axis_nphi'][0]
        BSTS.axis_min        = bsts['axis_phimin'][0] * np.pi * 2.0 / 360.0
        BSTS.axis_max        = bsts['axis_phimax'][0] * np.pi * 2.0 / 360.0
        # BSTS.axis_grid       = bsts['']   # Not really used

        # 2. Get the right sized offload array
        #--------------------------------------
        # Does this deed to be allocated by C code or can we do it with a python array?

        '''
        /* Allocate offload_array storing psi and the three components of B */
        int psi_size = offload_data->psigrid_n_r*offload_data->psigrid_n_z
           * offload_data->psigrid_n_phi;
        int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
           * offload_data->Bgrid_n_phi;
        int axis_size = offload_data->n_axis;

        *offload_array = (real*) malloc((psi_size + 3 * B_size + 2 * axis_size)
                                    * sizeof(real));
        offload_data->offload_array_length = psi_size + 3 * B_size + 2 * axis_size;
        '''

        offload_size = 0

        npsi = bsts['psi_nr'][0] * bsts['psi_nz'][0] * bsts[ 'psi_nphi'][0]
        nB   = bsts[  'b_nr'][0] * bsts[  'b_nz'][0] * bsts[   'b_nphi'][0]
        naxis=                                         bsts['axis_nphi'][0]

        # psi_size
        offload_size +=  1 * npsi
        # B_size
        offload_size +=  3 * nB
        # axis_size
        offload_size +=  2 * naxis

        BSTS.offload_array_length = offload_size

        # Python side allocation
        #-----------------------
        # B_offload_array = (ctypes.c_double * (offload_size) )()
        #self._B_offload_array.contents = B_offload_array

        # C side allocation
        #------------------
        self._bfield_offload_array = \
            ascot2py.libascot_allocate_reals(offload_size)

        # Cast the pointer into an array
        B_offload_array = ctypes.cast(
            self._bfield_offload_array,
            ctypes.POINTER(ctypes.c_double*offload_size) )[0]

        # 3. copy the large arrays to the offload array
        #----------------------------------------------

        offset = 0
        order = 'F'

        B_offload_array[(offset):(offset+nB  )] = \
            bsts['br'].flatten(order=order)
        offset += nB

        B_offload_array[(offset):(offset+nB  )] = \
            bsts['bphi'].flatten(order=order)
        offset += nB

        B_offload_array[(offset):(offset+nB  )] = \
            bsts['bz'].flatten(order=order)
        offset += nB

        B_offload_array[(offset):(offset+npsi)] = \
            bsts['psi'].flatten(order=order)
        offset += npsi

        B_offload_array[(offset):(offset+naxis)] = bsts['axisr']
        offset += naxis

        B_offload_array[(offset):(offset+naxis)] = bsts['axisz']
        offset += naxis

        # 4. Set the correct data type
        #------------------------------
        self._sim.B_offload_data.type = ascot2py.B_field_type_STS

        # 5. Do the init offload
        #-----------------------

        #B_field_init_offload.argtypes = \
        #    [ctypes.POINTER(struct_c__SA_B_field_offload_data),
        #     ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]

        ascot2py.B_field_init_offload(
            ctypes.byref(self._sim.B_offload_data),
            ctypes.byref(self._bfield_offload_array) )
