import os
import numpy as np
import warnings
from freeqdsk import geqdsk
import re

from a5py.physlib import cocos as cocosmod
import unyt

try:
    import adas
except ImportError:
   adas = None

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
            Mass in amu of fast ion in CX reaction.
        z2cx : int, optional
            Atomic number of background neutral (donor) in CX reaction.
        a2cx : int, optional
            Atomic mass number of background neutral in CX reaction.
        m2cx : float, optional
            Mass in amu of background neutral in CX reaction.
        ekinmincx : float, optional
            Minimum of energy abscissa for CX rate coefficients.
        ekinmaxcx : float, optional
            Maximum of energy abscissa for CX rate coefficients.
        nekincx : float, optional
            Number of points in energy abscissa for CX rate coefficients.
        tempmincx : float, optional
            Minimum of temperature abscissa for CX rate coefficients.
        tempmaxcx : float, optional
            Maximum of temperature abscissa for CX rate coefficients.
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
        with open(fn, "r") as f:
            eqd = geqdsk.read(f)

        cocos = cocosmod.assign(
            eqd["qpsi"][0], eqd["cpasma"], eqd["bcentr"], eqd["simagx"],
            eqd["sibdry"], phiclockwise, weberperrad)

        verbose and print("Eqdsk cocos: "+str(cocos))
        verbose and print("ASCOT cocos: "+str(cocosmod.COCOS_ASCOT))

        if cocos != cocosmod.COCOS_ASCOT:
            warnings.warn(
                "EQDSK COCOS is %d while ASCOT5 expects 3. Transforming COCOS"
                % cocos)
            eqd = cocosmod.fromCocosNtoCocosM(eqd, cocosmod.COCOS_ASCOT)

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
