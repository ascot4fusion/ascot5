import numpy as np
import warnings
import re
import unyt
import copy

from scipy.interpolate import RegularGridInterpolator

from a5py.physlib import cocos as cocosmod
from a5py.physlib import species as physlibspecies

try:
    import adas
except ImportError:
   adas = None

try:
    from freeqdsk import geqdsk
except ImportError:
    geqdsk = None

class ImportData():

    def import_adas(self,
                    z1cx=1, a1cx=2, m1cx=2.0135, z2cx=1, a2cx=2, m2cx=2.0135,
                    ekinmincx=5e1, ekinmaxcx=1e5, nekincx=1000,
                    tempmincx=1e0, tempmaxcx=2e4, ntempcx=100,
                    mltpresekinbms=1, mltpresdensbms=1, mltprestempbms=1,
                    show_progress=False, **kwargs):
        """Import data from ADAS files.

        This function may take several minutes.

        Parameters
        ----------
        z1cx : int, optional
            Atomic number (znum) of fast ion (receiver) in CX reaction.
        a1cx : int, optional
            Atomic mass number (anum) of fast ion in CX reaction.
        m1cx : float, optional
            Mass (amu) of fast ion in CX reaction.
        z2cx : int, optional
            Atomic number of background neutral (donor) in CX reaction.
        a2cx : int, optional
            Atomic mass number of background neutral in CX reaction.
        m2cx : float, optional
            Mass (amu) of background neutral in CX reaction.
        ekinmincx : float, optional
            Minimum of energy (eV) abscissa for CX rate coefficients.
        ekinmaxcx : float, optional
            Maximum of energy (eV) abscissa for CX rate coefficients.
        nekincx : float, optional
            Number of points in energy abscissa for CX rate coefficients.
        tempmincx : float, optional
            Minimum of temperature (eV) abscissa for CX rate coefficients.
        tempmaxcx : float, optional
            Maximum of temperature (eV) abscissa for CX rate coefficients.
        ntempcx : float, optional
            Number of points in temperature abscissa for CX rate coefficients.
        mltpresekin : int, optional
            Resolution multiplier for energy abscissa for BMS coefficients.
        mltpresdens : int, optional
            Resolution multiplier for density abscissa for BMS coefficients.
        mltprestemp : int, optional
            Resolution multiplier for temperature abscissa for
            BMS coefficients.
        show_progress : bool, optional
            Flag to determine if progress in converting cross-sections into
            rate coefficients should be printed to terminal.
        **kwargs
            ADAS data files in format: ``reaction``="/path/to/reaction/data".

            The key ``reaction`` is used to interpret the specific reaction
            and reactant species (charge state included) the data corresponds
            to, and it must follow the format ``"input_"<reaction>_<fast
            particle species><bulk particle species>``. Examples of valid
            key-value pairs are
            ``input_CX_H1H0="/home/adas/adas/adf24/scx#h0/scx#h0_ornl#h1.dat"``
            and
            ``input_BMS_H0H1="/home/adas/adas/adf21/bms10#h/bms10#h_h1.dat"``.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        if adas is None:
            raise Exception("Could not import adas package.")

        # Define relevant reaction types
        reac_type_sigma_cx   = 3
        reac_type_sigmav_cx  = 6
        reac_type_sigmav_bms = 7

        # Count number of reactions based on kwargs and initialize data list
        nreac = len(kwargs)
        sigmalist = [None]*nreac

        # Initialize arrays for reaction identifiers and abscissae
        z1       = np.zeros(nreac, dtype=int)
        a1       = np.zeros(nreac, dtype=int)
        z2       = np.zeros(nreac, dtype=int)
        a2       = np.zeros(nreac, dtype=int)
        reactype = np.zeros(nreac, dtype=int)
        nekin    = np.zeros(nreac, dtype=int)
        ekinmin  = np.zeros(nreac, dtype=float)
        ekinmax  = np.zeros(nreac, dtype=float)
        ndens    = np.zeros(nreac, dtype=int)
        densmin  = np.zeros(nreac, dtype=float)
        densmax  = np.zeros(nreac, dtype=float)
        ntemp    = np.zeros(nreac, dtype=int)
        tempmin  = np.zeros(nreac, dtype=float)
        tempmax  = np.zeros(nreac, dtype=float)
        n        = np.zeros(nreac, dtype=int)

        # Loop through kwargs
        ireac = 0
        for key,val in kwargs.items():

            if "CX" in key:
                # Reaction type for CX cross-section data is 3
                reactype[ireac] = reac_type_sigma_cx
                # Determine znum and anum values based on key
                match = re.match(r"([a-z]+)([0-9]+)([a-z]+)([0-9]+)",
                                 key.split('_')[-1],re.I)
                items = match.groups()
                if(items[0] == 'H'):
                    z1[ireac] = 1
                    a1[ireac] = 1
                else:
                    raise(Exception("Unsupported receiver species in CX."))
                if(items[2] == 'H'):
                    z2[ireac] = 1
                    a2[ireac] = 1
                else:
                    raise(Exception("Unsupported donor species in CX."))
                # Read data using ADAS function
                dat = adas.xxdata_24(val)
                ecol = dat['eea'][:,0]
                sigma = dat['scx'][:,0]

                ## Convert CX reaction data from cross-sections (sigma) to
                ## rate coefficients (sigmav).
                # Construct abscissae for rate coefficients.
                # First, check that wanted fast-ion energy range is possible
                # with given collision energy range of cross-section data.
                if((ekinmincx/m1cx >= ecol[0]) &
                   (ekinmaxcx/m1cx <= ecol[-1])):
                    ekin = np.linspace(ekinmincx,ekinmaxcx,nekincx)
                    ekinmin[ireac] = ekin[0]
                    ekinmax[ireac] = ekin[-1]
                    nekin[ireac] = nekincx
                else:
                    raise(Exception(
                        "Energy span goes outside input data domain."))
                # No density abscissa for this data so density dimension is 1
                ndens[ireac] = 1
                temp = np.linspace(tempmincx,tempmaxcx,ntempcx)
                tempmin[ireac] = temp[0]
                tempmax[ireac] = temp[-1]
                ntemp[ireac] = ntempcx
                n[ireac] = ntemp[ireac]*ndens[ireac]*nekin[ireac]
                # Initialize an array for the rate coefficient data, which is
                # a function of fast-ion energy and neutral temperature.
                sigmav = np.zeros(ntemp[ireac]*nekin[ireac])
                # Temporarily convert to 2D array to take advantage of
                # vectorization in ADAS function ceevth, which, at a given
                # energy, allows evaluation for several temperatures at once.
                sigmav2d = sigmav.reshape(ntemp[ireac],nekin[ireac])
                # Convert cross-sections into rate coefficients using
                # ADAS function ceevth.
                if(show_progress):
                    print('Converting cross-sections into rate coefficients:')
                for iekin in range(nekin[ireac]):
                    sigmav2d[:,iekin] = adas.ceevth(amdon=m2cx,amrec=m1cx,
                                                    catyp='td',
                                                    dren=ekin[iekin],iextyp=1,
                                                    log=False,
                                                    enin=ecol,sgin=sigma,
                                                    enout=temp)
                    if(show_progress):
                        quot,rem = np.divmod(iekin,int(nekin[ireac]/10))
                        if(rem == 0):
                            print('%d%%' % int(quot*10))
                if(show_progress):
                    print('Completed.')
                # Convert back to long 1D array where energy abscissa
                # runs fastest.
                sigmav = sigmav2d.flatten()
                # Convert units from cm3/s to m3/s
                sigmav = 1e-6*sigmav
                # Store data in data list
                sigmalist[ireac] = sigmav
                # Cross-section data (sigma(collision energy)) was converted
                # into rate coefficient data (sigmav(fast-ion energy,
                # background neutral temperature)). Hence, we update the
                # reaction type and znum and anum values accordingly.
                # Reaction type for CX rate coefficient data is 6
                reactype[ireac] = 6
                z1[ireac] = z1cx
                a1[ireac] = a1cx
                z2[ireac] = z2cx
                a2[ireac] = a2cx
                # Increment reaction index
                ireac += 1

            elif "BMS" in key:
                # Reaction type for BMS coefficient data is 7
                reactype[ireac] = 7
                # Determine znum and anum values based on key
                match = re.match(r"([a-z]+)([0-9]+)([a-z]+)([0-9]+)",
                                 key.split('_')[-1],re.I)
                items = match.groups()
                if(items[0] == 'H'):
                    z1[ireac] = 1
                    a1[ireac] = 1
                else:
                    raise(Exception("Unsupported beam species in BMS."))
                if(items[2] == 'H'):
                    z2[ireac] = 1
                    a2[ireac] = 1
                else:
                    raise(Exception("Unsupported target species in BMS."))
                # Read abscissa data using ADAS function
                dat = adas.xxdata_21(val)
                ekinmin[ireac] = dat['be'][0]
                ekinmax[ireac] = dat['be'][-1]
                nekin[ireac] = len(dat['be'])
                densmin[ireac] = dat['tdens'][0]
                densmax[ireac] = dat['tdens'][-1]
                ndens[ireac] = len(dat['tdens'])
                tempmin[ireac] = dat['ttemp'][0]
                tempmax[ireac] = dat['ttemp'][-1]
                ntemp[ireac] = len(dat['ttemp'])
                # Apply resolution multipliers for abscissa linearization
                nekin[ireac] = (nekin[ireac]-1)*mltpresekinbms + 1;
                ndens[ireac] = (ndens[ireac]-1)*mltpresdensbms + 1;
                ntemp[ireac] = (ntemp[ireac]-1)*mltprestempbms + 1;
                # Construct new, uniformly spaced (linear) abscissae
                ekin = np.linspace(ekinmin[ireac],ekinmax[ireac],nekin[ireac])
                dens = np.linspace(densmin[ireac],densmax[ireac],ndens[ireac])
                temp = np.linspace(tempmin[ireac],tempmax[ireac],ntemp[ireac])
                n[ireac] = ntemp[ireac]*ndens[ireac]*nekin[ireac]
                # Organize abscissae into a series of 3D points
                # for ADAS function read_adf21.
                pts = np.zeros((ntemp[ireac]*ndens[ireac]*nekin[ireac],3))
                for itemp in range(ntemp[ireac]):
                    for idens in range(ndens[ireac]):
                        for iekin in range(nekin[ireac]):
                            pts[itemp*ndens[ireac]*nekin[ireac] +
                                idens*nekin[ireac] + iekin,:] \
                                = np.array([temp[itemp],
                                            dens[idens],ekin[iekin]])
                # BMS effective rate coefficients
                sigmav = adas.read_adf21(files=val,fraction=[1.0],
                                         energy=pts[:,2],te=pts[:,0],
                                         dens=pts[:,1])
                # Convert units 1/cm3 to 1/m3
                densmin[ireac] = 1e6*densmin[ireac]
                densmax[ireac] = 1e6*densmax[ireac]
                # Convert units cm3/s to m3/s
                sigmav = 1e-6*sigmav
                # Store data in data list
                sigmalist[ireac] = sigmav
                # Increment reaction index
                ireac += 1
            else:
                raise(Exception("Unsupported kwarg."))

        ## Write data to HDF5
        # Place reaction data for all reactions in one long array
        ntot = 0
        for ireac in range(nreac):
            ntot += n[ireac]
        sigma = np.zeros((1,ntot))
        isigma0 = 0
        for ireac in range(nreac):
            sigma[0,isigma0:isigma0+n[ireac]] = sigmalist[ireac]
            isigma0 += n[ireac]

        out = {
            "nreac" : nreac,
            "z1" : z1, "a1" : a1, "z2" : z2, "a2" : a2, "reactype" : reactype,
            "nenergy" : nekin, "energymin" : ekinmin, "energymax" : ekinmax,
            "ndensity" : ndens, "densitymin" : densmin, "densitymax" : densmax,
            "ntemperature" : ntemp,
            "temperaturemin" : tempmin, "temperaturemax" : tempmax,
            "sigma" : sigma
        }
        return ("asigma_loc", out)

    def import_geqdsk(self, fn="input.eqdsk", cocos=None, phiclockwise=None,
                      weberperrad=None, verbose=True, interpolate_psi0=False):
        """Import axisymmetric magnetic field from EQDSK.

        Parameters
        ----------
        fn : str
            Filename of the G-EQDSK to read.
        cocos : int, optional
            Expected COCOS or None to deduce from the data.
        interpolate_psi0 : bool, optional
            Instead of using the psi on-axis value in EQDSK, interpolate it
            using libascot.

            Enabling this setting is recommended since having an incorrect value
            for psi0 could make rho imaginary near the axis, which leads to all
            kinds of trouble.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        if not geqdsk: raise ImportError("Package freeqdsk not found")
        with open(fn, "r") as f:
            eqd = geqdsk.read(f)

        if cocos:
            verbose and print("Ignoring COCOS check and using the give value")
        else:
            cocos = cocosmod.assign(
                eqd["qpsi"][0], eqd["cpasma"], eqd["bcentr"], eqd["simagx"],
                eqd["sibdry"], phiclockwise, weberperrad)

        verbose and print("Eqdsk cocos: "+str(cocos))
        verbose and print("ASCOT cocos: "+str(cocosmod.COCOS_ASCOT))

        if cocos != cocosmod.COCOS_ASCOT:
            warnings.warn(
                "EQDSK COCOS is %d while ASCOT5 expects 3. Transforming COCOS"
                % cocos)
            eqd = cocosmod.fromCocosNtoCocosM(eqd, cocosmod.COCOS_ASCOT,
                                              cocos_n=cocos)

        b2d = {
            "nr" : eqd["nx"],
            "rmin" : eqd["rleft"], "rmax" : eqd["rleft"]+eqd["rdim"],
            "nz" : eqd["ny"],
            "zmin" : eqd["zmid"] - 0.5*eqd["zdim"],
            "zmax" : eqd["zmid"] + 0.5*eqd["zdim"],
            "axisr" : eqd["rmagx"], "axisz" : eqd["zmagx"],
            "psi" : eqd["psi"], "psi0" : eqd["simagx"], "psi1" : eqd["sibdry"],
            "br" : eqd["psi"]*0, "bz" : eqd["psi"]*0
        }

        # Make sure the data grid does not extend to R=0
        if b2d["rmin"] == 0:
            b2d["nr"]  -= 1
            b2d["rmin"] = b2d["rmax"] / (b2d["nr"]+1)
            b2d["psi"]  = b2d["psi"][1:,:]
            b2d["br"]   = b2d["br"][1:,:]
            b2d["bz"]   = b2d["bz"][1:,:]

        # Toroidal component is more complicated for it can be evaluated from
        # Btor = F/R but we need to map F(psi) to F(R,z) first. However, F(psi)
        # is given only inside the plasma.
        psigrid = np.linspace(eqd["simagx"], eqd["sibdry"], eqd["nx"])
        if eqd["simagx"] < eqd["sibdry"]:
            fpolrz  = np.interp(eqd["psi"], psigrid, eqd["fpol"],
                                right=eqd["fpol"][-1])
        else:
            fpolrz  = np.interp(eqd["psi"], psigrid[::-1], eqd["fpol"][::-1],
                                right=eqd["fpol"][-1])
        rvec   = np.linspace(b2d["rmin"], b2d["rmax"], b2d["nr"])
        zvec   = np.linspace(b2d["zmin"], b2d["zmax"], b2d["nz"])
        R, Z   = np.meshgrid(rvec, zvec, indexing="ij")

        if eqd["nx"] != b2d["nr"]: fpolrz = fpolrz[1:,:] # If we had rmin=0
        b2d["bphi"] = fpolrz/R

        # Interpolate psi if needed
        if interpolate_psi0:
            self._ascot.input_init(bfield=b2d)
            b2d["axisr"], b2d["axisz"], b2d["psi0"] = \
                self._ascot.input_findpsi0(b2d["psi0"])
            self._ascot.input_free()

        return ("B_2DS", b2d)

    def import_wall_vtk(self, fn=None):
        """Import 3D wall from VTK file.

        Parameters
        ----------
        fn : str
            Path to the VTK file.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        import pyvista as pv
        mesh = pv.read(fn).extract_surface()
        cell_indices = mesh.faces

        # Ensure all elements are triangles
        if not ( not cell_indices.size % 4 and \
                np.all(cell_indices.reshape((-1, 4))[:,0] == 3) ):
            raise ValueError("Mesh contains non-triangular elements.")

        # cell_indices has format [Nvertices, v1_idx, v2_idx, v3_idx]
        cell_node_ids = cell_indices.reshape((-1, 4))[:,1:]
        cell_nodes = mesh.points[cell_node_ids]

        # cell_nodes has now format N x 3 x 3 but it is in millimeters
        x1x2x3 = cell_nodes[:,:,0] / 1000
        y1y2y3 = cell_nodes[:,:,1] / 1000
        z1z2z3 = cell_nodes[:,:,2] / 1000
        wall = {"nelements":x1x2x3.shape[0], "x1x2x3":x1x2x3, "y1y2y3":y1y2y3,
                "z1z2z3":z1z2z3}
        return ("wall_3D", wall)

    def import_plasma_profiles(self, fn=None, ne=None, ni=None, Te=None,
        Ti=None, anum=None, znum=None, charge=None, mass=None, species=None,
        pls=None, nrho=100, ionfrac=None, inputgrid="rho", extrapolate=None,
        extrapolate_len=0.1, **kwargs):
        """Import 1D plasma profiles.

        This function interpolates and extrapolates plasma profiles on a given
        grid. The input data can be in one of the three formats:

        1. Plasma profiles are in separate files:

           Specify the filenames containing the corresponding data in ``ne``,
           ``ni``, ``Te``, and ``Te``. The data is assumed to have a format
           where the first column is the radial coordinate and the second column
           is the value.

        2. Plasma profiles are in a single file:

           The data is assumed to be in a format where the first column is the
           radial coordinate and the following columns contain the values. Use
           ``ne``, ``ni``, ``Te``, and ``Te`` to specify the number of the
           column where the corresponding data is located (indexing starts from
           zero which corresponds to the first column with the radial
           coordinate).

        3. Plasma profiles are already read in a dictionary and they are only
           extrapolated here.

           Use this if the data you have doesn't conform with 1 or 2 and you
           have to read it separately. The dictionary ``pls`` must be in format
           conforming to :class:`~a5py.ascot5io.plasma.plasma_1D`.

        Parameters
        ----------
        fn : str, optional
            Filename if all plasma profiles are in a single file.
        ne : str or int, optional
            Filename containing electron density if profiles are in separate
            files, or column index if the profiles are in a single file.
        ni : str or int or list, optional
            Filename containing ion density if profiles are in separate
            files, or column index if the profiles are in a single file.

            If there are several ion species with separate profiles, this
            argument can be a list.
        Te : str or int, optional
            Filename containing electron temperature if profiles are in separate
            files, or column index if the profiles are in a single file.
        Ti : str or int, optional
            Filename containing ion temperature if profiles are in separate
            files, or column index if the profiles are in a single file.
        species : [str], (nion,), optional
            Names of the ion species if anum, znum, charge, and mass are not
            explicitly given.
        anum : [int]
            Ion species' atomic mass number if ``species`` is not given.
        znum : [int]
            Ion species' charge number if ``species`` is not given.
        mass : [int]
            Ion species' mass number if ``species`` is not given.
        charge : [int]
            Ion species' charge if ``species`` is not given.
        pls : dict, optional
            Dictionary containing the data in
            :class:`~a5py.ascot5io.plasma.plasma_1D` format.
        nrho : int, optional
            Number of radial grid points where the data is interpolated or
            extrapolated.

            The range is [0,1] by default and [0,``extrapolate``] if the data
            is to be extrapolated beyond rho=1.

            Interpolation is done linearly but if the data doesn't span the
            whole interval [0,1], nearest neighbour is used to extrapolate. Note
            that this is separate from the extrapolation done at rho > 1.
        ionfrac : list[float], (nion-i,), optional
            In case ion densities are not given independently, this list of
            fractions is used to divide ``ni`` in to individual densities.

            The order of the species should be same as in ``species``. Note that
            the fraction for the last species is not given. Instead, it is
            chosen so that the plasma is quasi-neutral.
        inputgrid : {"rho", "rhosquared"}, optional
            Specifies the format of the radial grid where the input values are
            given.

            - "rho": Square root of the normalized poloidal flux.
            - "rhosquared": Normalized poloidal flux.
        extrapolate : float, optional
            Extrapolate profiles up to this rho value.

            The extrapolated profiles are in form :math:`exp[-(\rho-1)/\lambda]`
            where :math:`\lambda` is ``extrapolate_len``.
        extrapolate_len : float, optional
            The ``e``-fold length (in rho) for the extrapolated profiles.
        **kwargs
            Arguments passed to :obj:`numpy.loadtxt` when reading data from
            file(s).

            For example, use ``skiprows`` if the input file(s) contain headers.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        # Check that inputs are valid
        if pls is None and species is None:
            if anum is None:   raise ValueError("species or anum is required")
            if znum is None:   raise ValueError("species or znum is required")
            if mass is None:   raise ValueError("species or mass is required")
            if charge is None: raise ValueError("species or charge is required")
            if np.ndim(anum)   == 0: anum   = [anum]
            if np.ndim(znum)   == 0: znum   = [znum]
            if np.ndim(mass)   == 0: mass   = [mass]
            if np.ndim(charge) == 0: charge = [charge]
        if pls is None and species:
            if not isinstance(species, list): species = [species]
            anum   = [None] * len(species)
            znum   = [None] * len(species)
            charge = [None] * len(species)
            mass   = [None] * len(species)
            for i in range(len(species)):
                d = physlibspecies.species(species[i])
                anum[i]   = d["anum"]
                znum[i]   = d["znum"]
                mass[i]   = d["mass"]
                charge[i] = d["charge"]

        # Set rho grid (include extrapolation)
        rho = np.linspace(0.0, 1.0, nrho)
        if extrapolate:
            rho = np.linspace(0.0, extrapolate, nrho)

        def interp(rad, val):
            """Function to interpolate and extrapolate values to rho grid"""
            # Mapping to rho from different radial coordinates
            if inputgrid == "rho":
                radgrid = rad
            elif inputgrid == "rhosquared":
                radgrid = np.sqrt(rad)

            # Extrapolate & interpolate (nearest-neighbour) in interval [0,1]
            val = np.interp(rho, radgrid, val, left=val[0], right=val[-1])
            if not extrapolate: return val

            # Exponential extrapolation at rho > 1
            extrp = rho >= 1.0
            val[extrp] = val[extrp][0] * np.exp( -( rho[extrp] - 1.0 ) /
                                                    extrapolate_len)
            return val

        # Read the data if dictionary was not provided
        if not pls:
            if not isinstance(ni, list): ni = [ni]
            if fn:
                # Single file containing all the data
                ints = [isinstance(x, int) for x in ni]
                if not (isinstance(ne, int) and isinstance(Te, int) and
                        isinstance(Ti, int) and all(ints)):
                    raise ValueError(
                        "ni, ne, Te, and Ti should contain the column numbers "
                        + "where the data is located.")
                data = np.loadtxt(fn, **kwargs)
                rad = data[:,0]
                ne  = interp(rad, data[:,ne])
                Te  = interp(rad, data[:,Te])
                Ti  = interp(rad, data[:,Ti])
                idx = ni
                ni  = np.zeros((nrho, len(idx)), dtype="f8")
                for i, j in enumerate(idx):
                    ni[:,i] = interp(rad, data[:,j])
            else:
                # Data is in separate files
                strs = [isinstance(x, str) for x in ni]
                if not (isinstance(ne, str) and isinstance(Te, str) and
                        isinstance(Ti, str) and all(strs)):
                    raise ValueError(
                        "ni, ne, Te, and Ti should contain the filenames "
                        "where the data is located.")
                data = np.loadtxt(ne, **kwargs)
                ne = interp(data[:,0], data[:,1])
                data = np.loadtxt(Te, **kwargs)
                Te = interp(data[:,0], data[:,1])
                data = np.loadtxt(Ti, **kwargs)
                Ti = interp(data[:,0], data[:,1])
                fns = ni
                ni  = np.zeros((nrho, len(fns)), dtype="f8")
                for i, fnin in enumerate(fns):
                    data = np.loadtxt(fnin, **kwargs)
                    ni[:,i] = interp(data[:,0], data[:,1])

            # Make the plasma input
            nion = ni.shape[1]
            if ionfrac:
                if nion > 1:
                    raise ValueError("Cannot use 'ionfrac' when input data "
                                     "already has multiple ion species")
                nion = len(ionfrac)+1
                if nion != len(anum) or nion != len(anum) or \
                   nion != len(mass) or nion != len(charge):
                    raise ValueError("Sizes of ionfrac and anum/znum/charge/"
                                     "mass are not consistent")
                ni_new = np.zeros((rho.size,nion))
                qdens  = np.zeros((rho.size,1)) * unyt.e
                for i in range(nion-1):
                    ni_new[:,i] = ni[:,0] * ionfrac[i]
                    qdens[:,0] += ni_new[:,i] * charge[i]

                ni_new[:,-1] = (ne[:] * unyt.e - qdens[:,0]) / charge[-1]
                if (ne[:] * unyt.e - qdens[:,0] < 0).any():
                    raise ValueError("Given ionfrac does not yield quasineutral"
                                     "plasma")
                ni = ni_new

            pls = {"nrho":rho.size, "rho":rho, "mass":mass, "charge":charge,
                   "anum":anum, "znum":znum, "nion":nion,
                   "etemperature":Te, "itemperature":Ti,
                   "edensity":ne, "idensity":ni}
        else:
            # Data is read already and only needs to be extrapolated
            pls["ne"] = interp(pls["rho"], pls["ne"])
            interp(pls["rho"], pls["Te"])
            interp(pls["rho"], pls["Ti"])
            interp(pls["rho"], pls["ni"])
            pls["rho"]  = rho
            pls["nrho"] = rho.size

        return ("plasma_1D", pls)

    def import_marsf(self, fn=None, n=None, b2d=None, b3d=None, phigrid=None):
        """Import toroidal harmonics calculated with MARS-F.

        This function converts the toroidal harmonics given in (R,z) grid to
        magnetic field components in (R,phi,z).

        The output consists of the input field (either 2D or 3D) with the 3D
        perurbation added. If the input is 2D, the toroidal grid must be
        specified in ``phigrid``.

        Parameters
        ----------
        fn : str
            Name of the MARS-F file.

            The file should contain the data in ASCII format organized in eight
            columns: R, z, Re(Br), Im(Br), Re(Bz), Im(Bz), Re(Bphi), Im(Bphi)
        n : int
            The toroidal harmonic the data corresponds to.
        b2d : dict, optional
            Dictionary containing 2D magnetic field data.

            If given, ``phigrid`` must also be specified.
        b3d : dict, optional
            Dictionary containing the 3D magnetic field data if ``b2d`` is not
            used.
        phigrid : array_like, optional
            Toroidal grid where the output field will be given if ``b2d`` is
            used [rad].

            Note that this grid must not contain the last point if it is the
            duplicate of the first point due to periodicity. In other words,
            if the data has 2pi periodicity then this argument should be
            ``np.linspace(0, 2pi, nphi+1)[:-1]``.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        if (b2d is None and b3d is None) or \
           (b2d is not None and b3d is not None):
            raise ValueError("Provide either 'b2d' or 'b3d'")
        if b2d is not None and phigrid is None:
            raise ValueError("'b2d' requires that 'phigrid' is also given")

        if b2d is not None:
            rgrid = np.linspace(b2d["rmin"], b2d["rmax"], b2d["nr"][0]).ravel()
            zgrid = np.linspace(b2d["zmin"], b2d["zmax"], b2d["nz"][0]).ravel()
        if b3d is not None:
            rgrid = np.linspace(b3d["b_rmin"], b3d["b_rmax"],
                                b3d["b_nr"][0]).ravel()
            zgrid = np.linspace(b3d["b_zmin"], b3d["b_zmax"],
                                b3d["b_nz"][0]).ravel()
            phigrid = np.linspace(b3d["b_phimin"], b3d["b_phimax"],
                                  b3d["b_nphi"][0]+1)[:-1].ravel()

        d = np.loadtxt(fn)
        marsf_rgrid = np.unique(d[:,0])
        marsf_zgrid = np.unique(d[:,1])
        nr = marsf_rgrid.size
        nz = marsf_zgrid.size

        def ifft(real, imag):
            """Inverse Fourier transform"""
            # Reshape the input 1D arrays to 2D (R,z) and then interpolate them
            # on given Rz grid.
            r,z = np.meshgrid(rgrid, zgrid, indexing='ij')
            real = RegularGridInterpolator(
                (marsf_rgrid, marsf_zgrid), np.reshape(real, (nr,nz)),
                method="cubic", bounds_error=False, fill_value=0)((r,z))
            imag = RegularGridInterpolator(
                (marsf_rgrid, marsf_zgrid), np.reshape(imag, (nr,nz)),
                method="cubic", bounds_error=False, fill_value=0)((r,z))

            # Calculate the inverse Fourier transform as
            # B(phi) = Re(B)*cos(phi) + Im(B)*sin(phi)
            return np.multiply.outer(real, np.cos(phigrid/n)) \
                 + np.multiply.outer(imag, np.sin(phigrid/n))

        br   = np.transpose(ifft(d[:,2], d[:,3]), (0,2,1))
        bz   = np.transpose(ifft(d[:,4], d[:,5]), (0,2,1))
        bphi = np.transpose(ifft(d[:,6], d[:,7]), (0,2,1))

        # Construct the magnetic field input
        if b3d is None:
            bphi += np.transpose(np.tile(b2d["bphi"], (phigrid.size,1,1)),
                                 (1,0,2))
            b3d = {
                "b_rmin":b2d["rmin"], "b_rmax":b2d["rmax"], "b_nr":b2d["nr"],
                "b_zmin":b2d["zmin"], "b_zmax":b2d["zmax"], "b_nz":b2d["nz"],
                "b_phimin":phigrid[0],
                "b_phimax":phigrid[-1]+phigrid[1]-phigrid[1],
                "b_nphi":phigrid.size,
                "axisr":b2d["axisr"], "axisz":b2d["axisz"], "psi":b2d["psi"],
                "psi0":b2d["psi0"], "psi1":b2d["psi1"],
                "br":br, "bphi":bphi, "bz":bz
            }
        else:
            b3d = copy.deepcopy(b3d)
            b3d["br"]   += br
            b3d["bphi"] += bphi
            b3d["bz"]   += bz

        return ("B_3DS", b3d)