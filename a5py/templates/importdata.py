import numpy as np
import warnings
import re
import unyt
import copy

from scipy.interpolate import RegularGridInterpolator,griddata,NearestNDInterpolator

from a5py.physlib import cocos as cocosmod
from a5py.physlib import species as physlibspecies
from a5py.ascot5io.wall import wall_3D

try:
    import adas
except ImportError:
   adas = None

try:
    from freeqdsk import geqdsk
except ImportError:
    geqdsk = None

#for desc_field
try:
    import desc.io as dscio
    import desc.grid as dscg
except:
    dscio = None
    dsccg = None

#for vmec_field and extender_field
try:
    import netCDF4 as nc
except:
    nc = None


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
                      weberperrad=True, verbose=True, interpolate_psi0=False):
        """Import axisymmetric magnetic field from EQDSK.

        Parameters
        ----------
        fn : str
            Filename of the G-EQDSK to read.
        cocos : int, optional
            Expected COCOS or None to deduce from the data.
        phiclockwise : Boolean
            If true, the phi-coordinate direction of the G-EQDSK file is assumed
            clockwise from above
        weberperrad : Boolean
            If true, the flux function is assumed to have been divided
            by 2*pi (COCOS ID 1-8)(default)
        verbose : Boolean
            If true, the function will talk a lot!
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

        verbose and print("G-EQDSK COCOS ID: "+str(cocos))
        verbose and print("ASCOT COCOS ID: "+str(cocosmod.COCOS_ASCOT))

        if cocos != cocosmod.COCOS_ASCOT:
            warnings.warn(
                "G-EQDSK HAS COCOS ID "+str(cocos) + " while ASCOT5 expects "
                + str(cocosmod.COCOS_ASCOT) + ". Transforming COCOS... "
            )
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
            print(b2d["psi0"])
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

    def wall_geqdsk(self, fn="input.eqdsk", phiclockwise=None, weberperrad=True,
                    return3Dwall=True, verbose=False):
        """Import wall data from EQDSK.

        Parameters
        ----------
        fn : str
            Filename of the G-EQDSK to read.
        phiclockwise : bool
            If True, the phi-coordinate direction of the G-EQDSK file is
            assumed clockwise from above.
        weberperrad : bool, optional
            If True, the flux function is assumed to have been divided by 2*pi
            (COCOS ID 1-8).
        return3Dwall : bool, optional
            Returns 3D wall object instead of 2D.
        verbose : bool, optional
            If True, make the function verbose.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        verbose and print("Loading wall data from G-EQDSK file... ")
        with open(fn, "r") as f:
            eqd = geqdsk.read(f)
        nwall  = eqd["nlim"] # Number of (R,z) points in the G-EQDSK wall data
        wall_r = eqd["rlim"]
        wall_z = eqd["zlim"]

        verbose and print(
            "Determining COCOS ID from G-EQDSK file data, phiclockwise and "
            "weberperrad keyword arguments... ")
        cocos_n = cocosmod.assign(eqd["qpsi"][0], eqd["cpasma"], eqd["bcentr"],
                                  eqd["simagx"], eqd["sibdry"], phiclockwise,
                                  weberperrad)

        if cocos_n!=cocosmod.COCOS_ASCOT:
            warnings.warn("G-EQDSK file has COCOS ID "+str(cocos_n)+
                          " while ASCOT5 expects "+str(cocosmod.COCOS_ASCOT)+
                          ". Transforming wall data... ")
            transform_dict = cocosmod.transform_cocos(
                cocosmod.cocos(cocos_n), cocosmod.cocos(cocosmod.COCOS_ASCOT))
            wall_r = wall_r*transform_dict["R"]
            wall_z = wall_z*transform_dict["Z"]

        # Remove last element since wall data from G-EQDSK are closed loops
        wall_r = wall_r[:-1]
        wall_z = wall_z[:-1]
        nwall = int(nwall-1)
        w2d = {"nelements":nwall, "r":wall_r, "z":wall_z}
        if return3Dwall:
            w3d = wall_3D.convert_wall_2D(180, **w2d)
            return ("wall_3D", w3d)
        return ("wall_2D", w2d)

    def import_plasma_profiles(self, fn=None, ne=None, ni=None, Te=None,
        Ti=None, anum=None, znum=None, charge=None, mass=None, species=None,
        pls=None, nrho=100, ionfrac=None, inputgrid="rho", extrapolate=None,
        extrapolate_len=0.1, Tmin=10, nmin=1.0e0, **kwargs):
        """Import 1D plasma profiles.

        This function interpolates and extrapolates plasma profiles on a given
        grid. Density is assumed to be in units of m^-3 and temperature in eV.

        The input data can be in one of the three formats:

        1. Plasma profiles are in separate files:

           Specify the filenames containing the corresponding data in arguments
           ``ne``, ``ni``, ``Te``, and ``Te``. The data is assumed to have
           a format where the first column is the radial coordinate and
           the second column is the value.

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

            The extrapolated profiles are in form
            :math:`exp[-(\\rho-1)/\\lambda]` where
            :math:`\\lambda` is ``extrapolate_len``.

            If the extrapolated values were to fall below ```Tlim`` or ``nlim``,
            these constant values are to be used instead.
        extrapolate_len : float, optional
            The ``e``-fold length (in rho) for the extrapolated profiles.
        Tmin : float, optional
            Minimum temperature in the extrapolated profile.
        nmin : float, optional
            Minimum density in the extrapolated profile.
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
                mass[i]   = d["mass"].to("amu")
                charge[i] = d["charge"]

        # Set rho grid (include extrapolation)
        rho = np.linspace(0.0, 1.0, nrho)
        if extrapolate:
            rho = np.linspace(0.0, extrapolate, nrho)

        def interp(rad, val, lim):
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

            # Make sure extrapolated values are not below the limits
            val[extrp] = np.maximum(val[extrp], lim)
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
                ne0 = data[-1,ne]
                ne  = interp(rad, data[:,ne], nmin)
                Te  = interp(rad, data[:,Te], Tmin)
                Ti  = interp(rad, data[:,Ti], Tmin)
                idx = ni
                ni  = np.zeros((nrho, len(idx)), dtype="f8")
                for i, j in enumerate(idx):
                    # The limit value is nmin for electrons, and the modified
                    # species-wise for ions to maintain quasi-neutrality
                    ni[:,i] = interp(rad, data[:,j], nmin * data[-1,j] / ne0)
            else:
                # Data is in separate files
                strs = [isinstance(x, str) for x in ni]
                if not (isinstance(ne, str) and isinstance(Te, str) and
                        isinstance(Ti, str) and all(strs)):
                    raise ValueError(
                        "ni, ne, Te, and Ti should contain the filenames "
                        "where the data is located.")
                data = np.loadtxt(ne, **kwargs)
                ne0 = data[-1,1]
                ne = interp(data[:,0], data[:,1], nmin)
                data = np.loadtxt(Te, **kwargs)
                Te = interp(data[:,0], data[:,1], Tmin)
                data = np.loadtxt(Ti, **kwargs)
                Ti = interp(data[:,0], data[:,1], Tmin)
                fns = ni
                ni  = np.zeros((nrho, len(fns)), dtype="f8")
                for i, fnin in enumerate(fns):
                    data = np.loadtxt(fnin, **kwargs)
                    # The limit value is nmin for electrons, and the modified
                    # species-wise for ions to maintain quasi-neutrality
                    ni[:,i] = interp(data[:,0], data[:,1],
                                     nmin * data[-1,1] / ne0)

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
                   "vtor":np.zeros((rho.size,1)),
                   "etemperature":Te, "itemperature":Ti,
                   "edensity":ne, "idensity":ni}
        else:
            # Data is read already and only needs to be extrapolated
            pls["ne"] = interp(pls["rho"], pls["ne"], nmin)
            pls["Te"] = interp(pls["rho"], pls["Te"], Tmin)
            pls["Ti"] = interp(pls["rho"], pls["Ti"], Tmin)
            pls["ni"] = interp(pls["rho"], pls["ni"], nmin)
            pls["vtor"] = interp(pls["rho"], pls["vtor"], pls["vtor"][-1])
            pls["rho"]  = rho
            pls["nrho"] = rho.size

        return ("plasma_1D", pls)

    def import_marsf(self, fn=None, n=None, scale=1.0, b2d=None, b3d=None,
                     phigrid=None):
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
        scale : float
            The value used to scale the MARS-F perturbation.

            Usually MARS-F data contains the perturbation caused by external
            coils. If the coil current used in the MARS-F simulation was 1 kAt,
            then the correct value for the coil current (e.g. 100 kAt) can be
            used here to obtain properly scaled fields.
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
        if n is None:
            raise ValueError("Toroidal number 'n' is required")

        if b2d is not None:
            rgrid = np.linspace(b2d["rmin"], b2d["rmax"], b2d["nr"][0]).ravel()
            zgrid = np.linspace(b2d["zmin"], b2d["zmax"], b2d["nz"][0]).ravel()
        if b3d is not None:
            rgrid = np.linspace(b3d["b_rmin"], b3d["b_rmax"],
                                b3d["b_nr"][0]).ravel()
            zgrid = np.linspace(b3d["b_zmin"], b3d["b_zmax"],
                                b3d["b_nz"][0]).ravel()
            phigrid = np.linspace(b3d["b_phimin"], b3d["b_phimax"],
                                  b3d["b_nphi"][0]+1)[:-1].ravel() * np.pi/180

        d = np.loadtxt(fn)
        marsf_rgrid = np.unique(d[:,0])
        marsf_zgrid = np.unique(d[:,1])
        nr = marsf_rgrid.size
        nz = marsf_zgrid.size

        # Check if the data is F ordered (R is incremented first) or C ordered
        # (z is incremented first). The first column is R so we can sniff that
        order = 'F' if d[0,0] < d[1,0] else 'C'

        def ifft(real, imag):
            """Inverse Fourier transform"""
            # Reshape the input 1D arrays to 2D (R,z) and then interpolate them
            # on given Rz grid.
            r,z = np.meshgrid(rgrid, zgrid, indexing='ij')
            real = RegularGridInterpolator(
                (marsf_rgrid, marsf_zgrid),
                np.reshape(real, (nr,nz), order=order),
                method="cubic", bounds_error=False, fill_value=0)((r,z))
            imag = RegularGridInterpolator(
                (marsf_rgrid, marsf_zgrid),
                np.reshape(imag, (nr,nz), order=order),
                method="cubic", bounds_error=False, fill_value=0)((r,z))

            # Calculate the inverse Fourier transform as
            # B(phi) = Re(B)*cos(n*phi) - Im(B)*sin(n*phi)
            return np.multiply.outer(real, np.cos(phigrid*n)) \
                 - np.multiply.outer(imag, np.sin(phigrid*n))

        br   = np.transpose(ifft(d[:,2], d[:,3]), (0,2,1)) * scale
        bz   = np.transpose(ifft(d[:,4], d[:,5]), (0,2,1)) * scale
        bphi = np.transpose(ifft(d[:,6], d[:,7]), (0,2,1)) * scale

        # Construct the magnetic field input
        if b3d is None:
            bphi += np.transpose(np.tile(b2d["bphi"], (phigrid.size,1,1)),
                                 (1,0,2))
            b3d = {
                "b_rmin":b2d["rmin"], "b_rmax":b2d["rmax"], "b_nr":b2d["nr"],
                "b_zmin":b2d["zmin"], "b_zmax":b2d["zmax"], "b_nz":b2d["nz"],
                "b_phimin":phigrid[0]*180/np.pi,
                "b_phimax":(phigrid[-1] + phigrid[1] - phigrid[0])*180/np.pi,
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

    def import_marker_locust(self, fn=None, power=None):
        """Convert LOCUST marker input to ASCOT5 markers.

        Parameters
        ----------
        fn : str
            Name of the input file.
        power : float
            The total power this marker population represents.

            The marker weight in LOCUST input is power/volume, so the easiest
            way to convert them to ASCOT5 weight (particles/s) is to just
            renormalize the weights.
        """
        if fn is None:
            raise ValueError("Input filename 'fn' is required")
        if power is None:
            raise ValueError("Total power 'power' is required")
        species = np.loadtxt(fn, max_rows=2)
        anum   = 4#int(species[0])
        znum   = 2#int(species[1])
        mass   = anum * unyt.amu
        charge = znum * unyt.e

        data = np.loadtxt(fn, skiprows=2)
        r    = data[:,0] * unyt.m
        phi  = data[:,1] * unyt.rad
        z    = data[:,2] * unyt.m
        vr   = data[:,3] * unyt.m/unyt.s
        vphi = data[:,4] * unyt.m/unyt.s
        vz   = data[:,5] * unyt.m/unyt.s
        w    = data[:,6]

        ekin = (0.5 * mass * (vr**2 + vphi**2 + vz**2)).to("J")
        rate = (power / ekin).v * unyt.particles / unyt.s
        weight = rate * w / np.sum(w)

        nmrk = r.size
        prt = {"n":nmrk, "ids":np.arange(nmrk)+1,
               "anum":anum*np.ones((nmrk,)), "znum":znum*np.ones((nmrk,)),
               "mass":mass*np.ones((nmrk,)), "charge":charge*np.ones((nmrk,)),
               "r":r, "phi":phi, "z":z, "vr":vr, "vphi":vphi, "vz":vz,
               "weight":weight, "time":np.zeros((nmrk,))}

        return ("prt", prt)
    @staticmethod
    def costransform(theta, phi, rmnc, xm, xn):
        """Cosine transform for VMEC coefficients

        Returns
        -------
        f : array
            Returns cosine expansion of rmnc coefficients as function of poloidal
            angle theta, toroidal angle phi, poloidal mode number xm, and
            toroidal mode number xn
        """
        ns  = rmnc.shape[0]
        xmt = xm*theta[:,None]
        xnz = xn*phi[:,None]
        f   = np.zeros((ns,theta.size,phi.size))
        for i in range(ns):
            f[i,:,:] = (np.matmul(rmnc[i,:]*np.cos(xmt), np.cos(xnz).T)
                        + np.matmul(rmnc[i,:]*np.sin(xmt), np.sin(xnz).T))
        return f
    @staticmethod
    def sintransform(theta, phi, rmnc, xm, xn):
        """Sine transform for VMEC coefficients

        Returns
        -------
        f : array
            Returns sine expansion of rmnc coefficients as function of poloidal
            angle theta, toroidal angle phi, poloidal mode number xm, and
            toroidal mode number xn
        """
        
        ns  = rmnc.shape[0]
        xmt = xm*theta[:,None]
        xnz = xn*phi[:,None]
        f   = np.zeros((ns,theta.size,phi.size))
        for i in range(ns):
            f[i,:,:] = (np.matmul(rmnc[i,:]*np.sin(xmt), np.cos(xnz).T)
                        - np.matmul(rmnc[i,:]*np.cos(xmt), np.sin(xnz).T))
        return f
    @staticmethod
    def vmec_field(ncfile,gtype="B_STS",phimin=0,phimax=361,nphi=361,ntheta=120,
                   nr=100,nz=100,psipad=0.0, ):
        """Load magnetic field data from a VMEC equilibrium.

        Notes
        -----
        VMEC coordinates:
           s = psi/psi1 (normalized toroidal flux)
           u = theta (poloidal angle)
           v = phi (toroidal angle)

        VMEC outputs quantities on two different radial meshes:
           `s_full = np.linspace(0, 1, ns) # full mesh`
           `s_half = np.insert(s_full[0:-1] + 0.5 / (ns - 1),0,np.nan) # half mesh`

        The interpolation/extrapolation scheme to convert from the half mesh 
        to the full meshis adapted from Hirshman et al. 1990: 
           https://doi.org/10.1016/0021-9991(90)90259-4.

        The toroidal magnetic flux is saved in the VMEC output as the variable `phi`
    
        In B_STS B_STS_eval_rho defines psi as:
        `rho[0] = sqrt( (psi - Bdata->psi0) / delta );`
        So psi is the toroidal magnetic flux.

        Values outside the LCFS are interpolated to the nearest values. Simulations
        should not extend past rho>1.

        Parameters
        ----------
        ncfile : str
           File path to VMEC NetCDF output.
        phimin : float, optional
           Minimum toroidal angle phi (deg). Default = 0.
        phimax : float, optional
           Maximum toroidal angle phi (deg). Default = 360.
        nphi : int, optional
           Number of toroidal angle phi grid points. Default = 360.
        ntheta : int, optional
           Number of poloidal angle theta grid points. Default = 120.
        nr : int, optional
           Number of radial coordinate R grid points. Default = 100.
        nz : int, optional
           Number of vertical coordinate Z grid points. Default = 100.
        psipad : float, optional
           Padding to slightly alter flux on axis. 

        Returns
        -------
        out : dict
        Dictionary with the following items:
           - `'axis_nphi'`, `'b_nphi'`, `'psi_nphi'`: nphi phi bins
           - `'b_nr'`, `'psi_nr'`: nr radial bins
           - `'b_nz'`, `'psi_nz'`: nz z-bins
           - `'axis_phimin'`, `'b_phimin'`, `'psi_phimin'`: phimin (deg)
           - `'axis_phimax'`, `'b_phimax'`, `'psi_phimax'`: (phimax-phimin)*(nphi-1)/nphi (deg)
           - `'b_rmin'`, `'psi_rmin'`: minimum radial coordinate R of output grids (m)
           - `'b_rmax'`, `'psi_rmax'`: maximum radial coordinate R of output grids (m)
           - `'b_zmin'`, `'psi_zmin'`: minimum vertical coordinate Z of output grids (m)
           - `'b_zmax'`, `'psi_zmax'`: maximum vertical coordinate Z of output grids (m)
           - `'axis_r'`: R(phi) on the magnetic axis (m)
           - `'axis_z'`: Z(phi) on the magnetic axis (m)
           - `'rlcfs'`: R(phi,len(u)) for the LCFS (m)
           - `'zlcfs'`: Z(phi,len(u)) for the LCFS (m)
           - `'psi0'`: toroidal magnetic flux on the magnetic axis (Wb)
           - `'psi1'`: toroidal magnetic flux through the last closed flux surface (Wb)
           - `'psi'`: toroidal magnetic flux psi(R,phi,Z) (Wb)
           - `'br'`: radial magnetic field B_R(R,phi,Z) (T)
           - `'bphi'`: toroidal magnetic field B_phi(R,phi,Z) (T)
           - `'bz'`: vertical magnetic field B_Z(R,phi,Z) (T)
        """
        
        # read VMEC NetCDF file
        if not nc: raise ImportError("Package netCDF4 not found")
        data = nc.Dataset(ncfile)
        xm = np.array(data.variables["xm"])  # poloidal mode numbers
        xn = np.array(data.variables["xn"])  # toroidal mode numbers
        xn_nyq = np.array(data.variables["xn_nyq"])  # tor mode numbers (Nyquist)
        xm_nyq = np.array(data.variables["xm_nyq"])  # pol mode numbers (Nyquist)
        psi_1d = np.array(data.variables["phi"])  # toroidal flux, full mesh
        rmnc = np.array(data.variables["rmnc"])  # cos(mn) comp of cyl R, full mesh
        zmns = np.array(data.variables["zmns"])  # sin(mn) comp of cyl Z, full mesh
        bsupumnc = np.array(data.variables["bsupumnc"])  #cos(mn) comp B^u,half mesh
        bsupvmnc = np.array(data.variables["bsupvmnc"])  #cos(mn) comp B^v,half mesh

        # poloidal angle array
        theta = np.linspace(0, 2.0 * np.pi, ntheta)  # rad

        # toroidal angle array
        # note phi should start at 0 and end on 360, inclusive
        phi = np.deg2rad(np.linspace(phimin, phimax, nphi, endpoint=False))  # rad
   
        # derivatives
        rumns = rmnc * (-1 * xm)  # drmn*cos(m*u-n*v)/du = -m*rmn*sin(m*u-n*v)
        zumnc = zmns * (xm)  # dzmn*sin(m*u-n*v)/du = m*zmn*cos(m*u-n*v)
        rvmns = rmnc * (xn)  # drmn*cos(m*u-n*v)/dv = n*rmn*sin(m*u-n*v)
        zvmnc = zmns * (-1 * xn)  # dzmn*sin(m*u-n*v)/dv = -n*zmn*cos(m*u-n*v)

        # convert from half mesh to full mesh
        ns = len(psi_1d)  # number of flux surfaces
        midx_e = np.nonzero(np.mod(xm_nyq, 2) == 0)[0]  # indices of even m modes
        midx_o = np.nonzero(np.mod(xm_nyq, 2) != 0)[0]  # indices of odd m modes
        w1 = np.atleast_2d(  # interpolation weights
            0.5 * np.sqrt((np.arange(ns-2) + 2 - 1)/(np.arange(ns-2) + 2 - 1.5))).T
        w2 = np.atleast_2d(  # interpolation weights
            0.5 * np.sqrt((np.arange(ns-2) + 2 - 1)/(np.arange(ns-2) + 2 - 0.5))).T
        w3 = 2 * np.sqrt((ns - 2) / (ns - 1.5))  # extrapolation weight

        # interpolation on intermediate surfaces
        bsupumnc[1:-1, midx_e] = 0.5*bsupumnc[1:-1,midx_e] + 0.5*bsupumnc[2:,midx_e]
        bsupvmnc[1:-1, midx_e] = 0.5*bsupvmnc[1:-1,midx_e] + 0.5*bsupvmnc[2:,midx_e]
        bsupumnc[1:-1, midx_o] = w1*bsupumnc[1:-1, midx_o] + w2*bsupumnc[2:, midx_o]
        bsupvmnc[1:-1, midx_o] = w1*bsupvmnc[1:-1, midx_o] + w2*bsupvmnc[2:, midx_o]

        # extrapolation to magnetic axis
        bsupumnc[0, midx_e] = 2 * bsupumnc[1, midx_e] - bsupumnc[2, midx_e]
        bsupvmnc[0, midx_e] = 2 * bsupvmnc[1, midx_e] - bsupvmnc[2, midx_e]
        bsupumnc[0, midx_o] = 0.0
        bsupvmnc[0, midx_o] = 0.0

        # extrapolation to boundary surface
        bsupumnc[-1, midx_e] = 2 * bsupumnc[-1, midx_e] - bsupumnc[-2, midx_e]
        bsupvmnc[-1, midx_e] = 2 * bsupvmnc[-1, midx_e] - bsupvmnc[-2, midx_e]
        bsupumnc[-1, midx_o] = w3 * bsupumnc[-1, midx_o] - bsupumnc[-2, midx_o]
        bsupvmnc[-1, midx_o] = w3 * bsupvmnc[-1, midx_o] - bsupvmnc[-2, midx_o]

        # inverse Fourier transform to (s,u,v) coordinates
        r_grid = ImportData.costransform(theta, phi, rmnc, xm, xn)  # R (m)
        z_grid = ImportData.sintransform(theta, phi, zmns, xm, xn)  # Z (m)
        ru_grid = ImportData.sintransform(theta, phi, rumns, xm, xn)  # dR/du (m/rad)
        zu_grid = ImportData.costransform(theta, phi, zumnc, xm, xn)  # dZ/du (m/rad)
        rv_grid = ImportData.sintransform(theta, phi, rvmns, xm, xn)  # dR/dv (m/rad)
        zv_grid = ImportData.costransform(theta, phi, zvmnc, xm, xn)  # dZ/dv (m/rad)
        bu_grid = ImportData.costransform(theta, phi, bsupumnc, xm_nyq, xn_nyq)  # B^u (T)
        bv_grid = ImportData.costransform(theta, phi, bsupvmnc, xm_nyq, xn_nyq)  # B^v (T)

        # magnetic axis
        axis_r = r_grid[0, 0, :]  # m
        axis_z = z_grid[0, 0, :]  # m

        #get lcfs
        lcfs_r = r_grid[-1, :, :] # m
        lcfs_z = z_grid[-1, :, :] # m

        # calculate B_R, B_phi, B_z
        br_grid = bu_grid * ru_grid + bv_grid * rv_grid
        bphi_grid = bv_grid * r_grid
        bz_grid = bu_grid * zu_grid + bv_grid * zv_grid

        # range in R and Z, make interpolation arrays
        rmin = np.amin(r_grid)
        rmax = np.amax(r_grid)
        zmin = np.amin(z_grid)
        zmax = np.amax(z_grid)
        r_1d = np.linspace(rmin, rmax, nr)  # m
        z_1d = np.linspace(zmin, zmax, nz)  # m
        z_2d, r_2d = np.meshgrid(z_1d, r_1d)
    
        # toroidal magnetic flux
        psi0 = psi_1d[0]  # axis
        psi1 = psi_1d[-1]  # LCFS
        psi_2d = np.tile(psi_1d, (ntheta, 1)).T  # Wb

        # interpolate psi, B_R, B_phi, B_Z to cylindircal coordinates
        psi = np.zeros([nr, nz, nphi])
        br = np.zeros([nr, nz, nphi])
        bphi = np.zeros([nr, nz, nphi])
        bz = np.zeros([nr, nz, nphi])

        # interpolate to cylindrical grid, iterate through toroidal angle
        for i in range(nphi):
            # interpolate data inside VMEC domain
            print(i)
            psi[:, :, i] = griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                psi_2d.flatten(),
                (r_2d, z_2d),
                fill_value=psi1,
            )
            br[:, :, i] = griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                br_grid[:, :, i].flatten(),
                (r_2d, z_2d),
            )
            bphi[:, :, i] = griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                bphi_grid[:, :, i].flatten(),
                (r_2d, z_2d),
            )
            bz[:, :, i] = griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                bz_grid[:, :, i].flatten(),
                (r_2d, z_2d),
            )

            #Replace br, bphi, bz NaN values outside LCFS with closest values
            data = br[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            br[:,:,i] = filled_data
       
            data = bz[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bz[:,:,i] = filled_data

            data = bphi[:,:,i]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bphi[:,:,i] = filled_data

        # change order from [R,Z,phi] to [R,phi,Z]
        psi = np.transpose(psi, (0, 2, 1))
        br = np.transpose(br, (0, 2, 1))
        bphi = np.transpose(bphi, (0, 2, 1))
        bz = np.transpose(bz, (0, 2, 1))

        #pad psi0 if needed
        if psipad != 0.0:
            print('Warning: Padding psi0 with',psipad)
            psi0 += psipad

        out = {
            "axis_phimin": phimin,  # deg
            "axis_phimax": np.rad2deg(phi[-1]),  # deg
            "axis_nphi": nphi,
            "axisr": axis_r,  # m
            "axisz": axis_z,  # m
            #"rlcfs": lcfs_r,  # m
            #"zlcfs": lcfs_z,  # m
            "b_rmin": rmin,  # m
            "b_rmax": rmax,  # m
            "b_nr": nr,
            "b_zmin": zmin,  # m
            "b_zmax": zmax,  # m
            "b_nz": nz,
            "b_phimin": phimin,  # deg
            "b_phimax": np.rad2deg(phi[-1]),  # deg
            "b_nphi": nphi,
            "br": br,  # T
            "bphi": bphi,  # T
            "bz": bz,  # T
            "psi": psi,  # Wb
            "psi0": psi0,  # Wb
            "psi1": psi1,  # Wb
            "psi_rmin": rmin,  # m
            "psi_rmax": rmax,  # m
            "psi_nr": nr,
            "psi_zmin": zmin,  # m
            "psi_zmax": zmax,  # m
            "psi_nz": nz,
            "psi_phimin": phimin,  # deg
            "psi_phimax": np.rad2deg(phi[-1]),  # deg
            "psi_nphi": nphi,
        }

        return (gtype, out)

    @staticmethod
    def desc_field(h5file,phimin=0,phimax=361,nphi=361,ntheta=120,
                   nr=100,nz=100,psipad=0.0):
        """Load magnetic field data from a DESC equilibrium.

        Notes
        -----
        The toroidal magnetic flux is saved in the DESC output as the variable `Psi`
    
        In B_STS B_STS_eval_rho defines psi as:
        `rho[0] = sqrt( (psi - Bdata->psi0) / delta );`
        So psi is the toroidal magnetic flux.

        Values outside the LCFS are interpolated to the nearest values. Simulations
        should not extend past rho>1.

        Parameters
        ----------
        h5file : str
           File path to DESC HDF5 output.
        phimin : float, optional
           Minimum toroidal angle phi (deg). Default = 0.
        phimax : float, optional
           Maximum toroidal angle phi (deg). Default = 360.
        nphi : int, optional
           Number of toroidal angle phi grid points. Default = 360.
        ntheta : int, optional
           Number of poloidal angle theta grid points. Default = 120.
        nr : int, optional
           Number of radial coordinate R grid points. Default = 100.
        nz : int, optional
           Number of vertical coordinate Z grid points. Default = 100.
        psipad : float, optional
           Padding to slightly alter flux on axis.

        Returns
        -------
        out : dict
        Dictionary with the following items:
           - `'axis_nphi'`, `'b_nphi'`, `'psi_nphi'`: nphi phi bins
           - `'b_nr'`, `'psi_nr'`: nr radial bins
           - `'b_nz'`, `'psi_nz'`: nz z-bins
           - `'axis_phimin'`, `'b_phimin'`, `'psi_phimin'`: phimin (deg)
           - `'axis_phimax'`, `'b_phimax'`, `'psi_phimax'`: (phimax-phimin)*(nphi-1)/nphi (deg)
           - `'b_rmin'`, `'psi_rmin'`: minimum radial coordinate R of output grids (m)
           - `'b_rmax'`, `'psi_rmax'`: maximum radial coordinate R of output grids (m)
           - `'b_zmin'`, `'psi_zmin'`: minimum vertical coordinate Z of output grids (m)
           - `'b_zmax'`, `'psi_zmax'`: maximum vertical coordinate Z of output grids (m)
           - `'axis_r'`: R(phi) on the magnetic axis (m)
           - `'axis_z'`: Z(phi) on the magnetic axis (m)
           - `'rlcfs'`: R(phi,len(u)) for the LCFS (m)
           - `'zlcfs'`: Z(phi,len(u)) for the LCFS (m)
           - `'psi0'`: toroidal magnetic flux on the magnetic axis (Wb)
           - `'psi1'`: toroidal magnetic flux through the last closed flux surface (Wb)
           - `'psi'`: toroidal magnetic flux psi(R,phi,Z) (Wb)
           - `'br'`: radial magnetic field B_R(R,phi,Z) (T)
           - `'bphi'`: toroidal magnetic field B_phi(R,phi,Z) (T)
           - `'bz'`: vertical magnetic field B_Z(R,phi,Z) (T)
        """

        #load DESC equilibrium
        if not dscio: raise ImportError("Package desc-opt not found")
        fam = dscio.load(h5file)
        try:  # if file is an EquilibriaFamily, use final Equilibrium
            eq = fam[-1]
        except:  # file is already an Equilibrium
            eq = fam
        eq.resolution_summary()

        # poloidal angle array
        theta = np.linspace(0, 2.0 * np.pi, ntheta)  # rad

        # toroidal angle array
        # note: phi should start at 0 and end on 360, inclusive
        phi = np.deg2rad(np.linspace(phimin, phimax, nphi, endpoint=False))  # rad

        # magnetic axis
        grid_axis = dscg.LinearGrid(rho=0.0, zeta=nphi, NFP=1)
        data_axis = eq.compute(["R", "Z"], grid=grid_axis)
        axis_r = data_axis["R"]  # m
        axis_z = data_axis["Z"]  # m
        psi0 = 0  # Wb

        # boundary
        grid = dscg.LinearGrid(rho=1.0, theta=ntheta, zeta=nphi,
                               NFP=1, sym=False,endpoint=True)
        data = eq.compute(["R", "Z"], grid=grid)
        theta_2D = grid.nodes[:, 1].reshape((grid.num_zeta, grid.num_theta),
                                            order="C").T
        phi_2D = grid.nodes[:, 2].reshape((grid.num_zeta, grid.num_theta),
                                          order="C").T
        bdry_r = data["R"].reshape((grid.num_zeta, grid.num_theta),
                                   order="C").T
        bdry_z = data["Z"].reshape((grid.num_zeta, grid.num_theta),
                                   order="C").T
        
        # boundaries
        rmin = np.min(bdry_r)  # m
        rmax = np.max(bdry_r)  # m
        zmin = np.min(bdry_z)  # m
        zmax = np.max(bdry_z)  # m
        psi1 = eq.Psi  # Wb

        # output domain
        R_1d = np.linspace(rmin, rmax, nr)  # m
        Z_1d = np.linspace(zmin, zmax, nz)  # m
        Z_2d, R_2d = np.meshgrid(Z_1d, R_1d)

        # interpolate psi, B_R, B_phi, B_Z to cylindircal coordinates
        psi = np.zeros([nr, nz, nphi])
        br = np.zeros([nr, nz, nphi])
        bphi = np.zeros([nr, nz, nphi])
        bz = np.zeros([nr, nz, nphi])

        # interpolate to cylindrical grid, iterate through toroidal angle
        for k in range(nphi):
            print(k)
            grid = dscg.ConcentricGrid(
                L=eq.L_grid, M=eq.M_grid, N=0, NFP=eq.NFP, node_pattern="linear")
            grid._nodes[:, 2] = phi[k]
            data = eq.compute(["R", "Z", "psi", "B_R", "B_phi", "B_Z"], grid=grid)
            
            # interpolate data inside DESC domain
            psi[:, :, k] = griddata(
                (data["R"], data["Z"]),
                data["psi"] * 2 * np.pi,  # DESC `psi` is normalized by 2 pi
                (R_2d, Z_2d),
                fill_value=psi1,
            )
            br[:,:,k] = griddata((data["R"],data["Z"]),data["B_R"],(R_2d,Z_2d))
            bphi[:,:,k]=griddata((data["R"],data["Z"]),data["B_phi"],(R_2d,Z_2d))
            bz[:,:,k] = griddata((data["R"],data["Z"]),data["B_Z"], (R_2d,Z_2d))

            #Replace br, bphi, bz NaN values outside LCFS with closest values
            data = br[:,:,k]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            br[:,:,k] = filled_data

            data = bz[:,:,k]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bz[:,:,k] = filled_data

            data = bphi[:,:,k]
            mask = np.where(~np.isnan(data))
            interp = NearestNDInterpolator(np.transpose(mask),data[mask])
            filled_data = interp(*np.indices(data.shape))
            bphi[:,:,k] = filled_data

        # change order from [R,Z,phiang] to [R,phiang,Z]
        psi = np.transpose(psi, (0, 2, 1))
        br = np.transpose(br, (0, 2, 1))
        bphi = np.transpose(bphi, (0, 2, 1))
        bz = np.transpose(bz, (0, 2, 1))

        #pad psi0 if needed
        if psipad != 0.0:
            print('Warning: Padding psi0 with',psipad)
            psi0 += psipad

        out = {
            "axis_phimin": phimin,  # deg
            "axis_phimax": np.rad2deg(phi[-1]),  # deg
            "axis_nphi": nphi,
            "axisr": axis_r,  # m
            "axisz": axis_z,  # m
            "rlcfs": bdry_r,  # m
            "zlcfs": bdry_z,  # m
            "b_rmin": rmin,  # m
            "b_rmax": rmax,  # m
            "b_nr": nr,
            "b_zmin": zmin,  # m
            "b_zmax": zmax,  # m
            "b_nz": nz,
            "b_phimin": phimin,  # deg
            "b_phimax": np.rad2deg(phi[-1]),  # deg
            "b_nphi": nphi,
            "br": br,  # T
            "bphi": bphi,  # T
            "bz": bz,  # T
            "psi": psi,  # Wb
            "psi0": psi0,  # Wb
            "psi1": psi1,  # Wb
            "psi_rmin": rmin,  # m
            "psi_rmax": rmax,  # m
            "psi_nr": nr,
            "psi_zmin": zmin,  # m
            "psi_zmax": zmax,  # m
            "psi_nz": nz,
            "psi_phimin": phimin,  # deg
            "psi_phimax": np.rad2deg(phi[-1]),  # deg
            "psi_nphi": nphi,
        }

        return out

    @staticmethod
    def extender_field(ncfile,extfile,ntheta=120,psipad=0.0):
        """Load magnetic field data from VMEC and EXTENDER outputs into ASCOT input format.

        Notes
        -----
        The toroidal magnetic flux is saved in the VMEC output as the variable `phi`
        
        In B_STS B_STS_eval_rho defines psi as:
        `rho[0] = sqrt( (psi - Bdata->psi0) / delta );`
        So psi is the toroidal magnetic flux.

        Psi values outside the LCFS are defined as psi1, so while markers may exist in the
        vacuum region (beyond LCFS) the mapping to the rho variable will not be valid
        as the rho=1 everywhere in the vacuum region. 

        The magnetic field (Br,Bphi,Bz) is defined via the EXTENDER input file while the
        toroidal flux (psi) is supplied via the VMEC input. Care should be taken to 
        make sure that the VMEC psi calculations agree with EXTENDER for rho<1. 

        EXTENDER input is assumed to cover one field period!!

        Parameters
        ----------
        ncfile : str
           File path to VMEC NetCDF output (free-boundary solution).
        extfile : str
           File path to EXTENDER NetCDF output that corresponds to the VMEC solution.
        ntheta : int, optional
           Number of poloidal angle theta grid points. Default = 120.
        psipad : float, optional
           Padding to slightly alter flux on axis.

        Returns
        -------
        out : dict
        Dictionary with the following items:
           - `'axis_nphi'`, `'b_nphi'`, `'psi_nphi'`: nphi phi bins
           - `'b_nr'`, `'psi_nr'`: nr radial bins
           - `'b_nz'`, `'psi_nz'`: nz z-bins
           - `'axis_phimin'`, `'b_phimin'`, `'psi_phimin'`: phimin (deg)
           - `'axis_phimax'`, `'b_phimax'`, `'psi_phimax'`: (phimax-phimin)*(nphi-1)/nphi (deg)
           - `'b_rmin'`, `'psi_rmin'`: minimum radial coordinate R of output grids (m)
           - `'b_rmax'`, `'psi_rmax'`: maximum radial coordinate R of output grids (m)
           - `'b_zmin'`, `'psi_zmin'`: minimum vertical coordinate Z of output grids (m)
           - `'b_zmax'`, `'psi_zmax'`: maximum vertical coordinate Z of output grids (m)
           - `'axis_r'`: R(phi) on the magnetic axis (m)
           - `'axis_z'`: Z(phi) on the magnetic axis (m)
           - `'rlcfs'`: R(phi,len(u)) for the LCFS (m)
           - `'zlcfs'`: Z(phi,len(u)) for the LCFS (m)
           - `'psi0'`: toroidal magnetic flux on the magnetic axis (Wb)
           - `'psi1'`: toroidal magnetic flux through the last closed flux surface (Wb)
           - `'psi'`: toroidal magnetic flux psi(R,phi,Z) (Wb)
           - `'br'`: radial magnetic field B_R(R,phi,Z) (T)
           - `'bphi'`: toroidal magnetic field B_phi(R,phi,Z) (T)
           - `'bz'`: vertical magnetic field B_Z(R,phi,Z) (T)
        """
        
        # load NetCDF data
        if not nc: raise ImportError("Package netCDF4 not found")
        vmec = nc.Dataset(ncfile)
        extender = nc.Dataset(extfile)
    
        # array dimenstions
        nr = int(extender.variables["ir"].getValue())
        nz = int(extender.variables["jz"].getValue())
        nphi = int(extender.variables["kp"].getValue())
        nfp = int(extender.variables["nfp"].getValue())
   
        # coordinate bounds (m)
        rmin = float(extender.variables["rmin"].getValue())
        rmax = float(extender.variables["rmax"].getValue())
        zmin = float(extender.variables["zmin"].getValue())
        zmax = float(extender.variables["zmax"].getValue())
    
        # magnetic field (T)
        br = np.array(extender.variables["br_001"])
        bphi = np.array(extender.variables["bp_001"])
        bz = np.array(extender.variables["bz_001"])

        # poloidal and toroidal angles (rad)
        theta = np.linspace(0, 2 * np.pi, ntheta, endpoint=True)
        phi = np.linspace(0, 2 * np.pi / nfp, nphi, endpoint=False)
    
        # VMEC spectral coefficients
        xm = np.array(vmec.variables["xm"])  # poloidal mode numbers
        xn = np.array(vmec.variables["xn"])  # toroidal mode numbers
        rmnc = np.array(vmec.variables["rmnc"])  # cos(mn) comp of cyl R, full mesh
        zmns = np.array(vmec.variables["zmns"])  # sin(mn) comp of cyl Z, full mesh
        psi_1d = np.array(vmec.variables["phi"])  # toroidal flux, full mesh

        # inverse Fourier transform to (psi,theta,phi) coordinates (m)
        r_grid = self.costransform(theta, phi, rmnc, xm, xn)
        z_grid = self.sintransform(theta, phi, zmns, xm, xn)

        # magnetic axis (m)
        axisr = r_grid[0, 0, :]
        axisz = z_grid[0, 0, :]

        #get lcfs (m)
        lcfs_r = r_grid[-1, :, :] 
        lcfs_z = z_grid[-1, :, :]

        # toroidal magnetic flux (Wb)
        psi0 = psi_1d[0] #axis
        psi1 = psi_1d[-1] #LCFS
        psi_2d = np.tile(psi_1d, (ntheta, 1)).T
    
        # R,Z interpolation points (m)
        r_1d = np.linspace(rmin, rmax, nr)
        z_1d = np.linspace(zmin, zmax, nz)
        z_2d, r_2d = np.meshgrid(z_1d, r_1d)

        # interpolate psi inside VMEC domain to cylindircal coordinates
        psi = np.zeros([nr, nz, nphi])
        for i in range(nphi):
            print(i)
            psi[:, :, i] = griddata(
                (r_grid[:, :, i].flatten(), z_grid[:, :, i].flatten()),
                psi_2d.flatten(),
                (r_2d, z_2d),
                fill_value=psi1,
            )
    
        # close NetCDF data
        vmec.close()
        extender.close()
    
        # change order from [R,Z,phi] to [phi,Z,R] for below
        psi = np.transpose(psi, (2, 1, 0))
    
        # repeat for each field period
        phidum = phi
        for i in range(1,nfp):
            phi = np.append(phi,phidum+i*2*np.pi/nfp)
        axisr = np.tile(axisr, nfp)
        axisz = np.tile(axisz, nfp)
        lcfs_r = np.tile(lcfs_r,(1,nfp))
        lcfs_z = np.tile(lcfs_z,(1,nfp))
        br = np.tile(br, (nfp, 1, 1))
        bphi = np.tile(bphi, (nfp, 1, 1))
        bz = np.tile(bz, (nfp, 1, 1))
        psi = np.tile(psi, (nfp, 1, 1))
       
        # repeat endpoint phi=0 == phi=360
        phi = np.append(phi,2*np.pi)
        nphi = len(phi)
        axisr = np.append(axisr, axisr[0])
        axisz = np.append(axisz, axisz[0])
        br = np.append(br, br[0, :, :][None, :], axis=0)
        bphi = np.append(bphi, bphi[0, :, :][None, :], axis=0)
        bz = np.append(bz, bz[0, :, :][None, :], axis=0)
        psi = np.append(psi, psi[0, :, :][None, :], axis=0)
    
        #dumb way to append to lcfs endpoint
        dumx,dumy = lcfs_r.shape
        new_lcfsr = np.zeros([dumx,dumy+1])
        new_lcfsr[:,0:dumy] = lcfs_r
        new_lcfsr[:,-1] = lcfs_r[:,0]
        new_lcfsz = np.zeros([dumx,dumy+1])
        new_lcfsz[:,0:dumy] = lcfs_z
        new_lcfsz[:,-1] = lcfs_z[:,0]
        lcfs_r = new_lcfsr
        lcfs_z = new_lcfsz

        # change order from [phi,Z,R] to [R,phi,Z]
        br = np.transpose(br, (2, 0, 1))
        bphi = np.transpose(bphi, (2, 0, 1))
        bz = np.transpose(bz, (2, 0, 1))
        psi = np.transpose(psi, (2, 0, 1))

        #pad psi0 if needed
        if psipad != 0.0:
            print('Warning: Padding psi0 with',psipad)
            psi0 += psipad
    
        out = {
            "axis_phimin": 0,  # deg
            "axis_phimax": np.rad2deg(phi[-1]),  # deg
            "axis_nphi": nphi,
            "axisr": axisr,  # m
            "axisz": axisz,  # m
            "rlcfs": lcfs_r, # m
            "zlcfs": lcfs_z, #m
            "b_rmin": rmin,  # m
            "b_rmax": rmax,  # m
            "b_nr": nr,
            "b_zmin": zmin,  # m
            "b_zmax": zmax,  # m
            "b_nz": nz,
            "b_phimin": 0,  # deg
            "b_phimax": np.rad2deg(phi[-1]),  # deg
            "b_nphi": nphi,
            "br": br,  # T
            "bphi": bphi,  # T
            "bz": bz,  # T
            "psi": psi,  # Wb
            "psi0": psi0,  # Wb
            "psi1": psi1,  # Wb
            "psi_rmin": rmin,  # m
            "psi_rmax": rmax,  # m
            "psi_nr": nr,
            "psi_zmin": zmin,  # m
            "psi_zmax": zmax,  # m
            "psi_nz": nz,
            "psi_phimin": 0,  # deg
            "psi_phimax": np.rad2deg(phi[-1]),  # deg
            "psi_nphi": nphi,
        }

        return out
