import numpy as np
import scipy
import h5py
import warnings

class Ascot4Templates():

    def ascot4_particles(self, fn="input.particles", force=None):
        """Convert marker input from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 marker input.
        force : {"prt", "gc"}, optional
            Sometimes ASCOT4 input.particles has both particle and guiding
            center data, and this toggle can be used to choose between them.

            If not given, the markers will be converted to particles.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        headerLength = 0
        with open(fn,'r') as f:
            line = f.readline()
            headerLength += 1
            if line[:-1] != " PARTICLE INITIAL DATA FOR ASCOT" \
               and line[:-1] != " PARTICLE OUTPUT DATA FROM ASCOT":
                raise NameError("Unrecognized first line in \""+fn+"\".")

            line = f.readline()
            headerLength += 1
            if line.split()[0] != "4":
                raise NameError("Bad file version in \""+fn+"\".")
            line = f.readline()
            line = f.readline()
            headerLength += 2
            nComments = int(float(line.split()[0]))
            for i in range(nComments):
                f.readline()
                headerLength += 1

            line = f.readline()
            line = f.readline()
            headerLength += 2
            nParticles = int(float(line.split()[0]))

            line = f.readline()
            line = f.readline()
            headerLength += 2
            nFields = int(float(line.split()[0]))

            fieldNames = []
            for i in range(0,nFields):
                fieldNames.append(f.readline()[:10].strip())
                headerLength += 1

            line = f.readline()
            headerLength += 1

        # Read the actual data table
        # --------------------------
        columns = np.loadtxt(fn, skiprows=headerLength)

        # -1 means unknown number of lines.
        if nParticles == -1:
            nParticles = columns.shape[0]

        fields = {}
        for j in range(nFields):
            fields[fieldNames[j]] = -999.0*np.ones((nParticles,))
            for i in range(nParticles):
                fields[fieldNames[j]][i] = columns[i,j]

        if "charge" not in fieldNames and "Znum" in fieldNames:
            warnings.warn("Converting Znum to charge.")
            fields['charge'] = fields['Znum'].astype('float')
        if "mass" not in fieldNames and "Anum" in fieldNames:
            warnings.warn("Converting Anum to mass as mass = Anum * amu.")
            fields["mass"] = fields['Anum'].astype('float')
        if "id" not in fieldNames:
            warnings.warn("No IDs in Data. Generating unique ids.")
            fields["id"] = np.arange(1, fields["charge"].size + 1)
        if (min(fields["id"]) <= 0):
            zero_ind = np.where(fields["id"] == 0)[0]
            fields["id"][zero_ind] = max(fields["id"] ) + 1
            warnings.warn("Converting ID 0 to new unique ID.")

        if ("vphi" in fieldNames and force is None) or force == "prt":
            # We have particles
            warnings.warn("Set time to zero for all markers.")
            out = {     "n":fields["id"].size,    "ids":fields["id"],
                     "mass":fields["mass"],    "charge":fields["charge"],
                        "r":fields["Rprt"],       "phi":fields["phiprt"],
                        "z":fields["zprt"],        "vr":fields["vR"],
                     "vphi":fields["vphi"],        "vz":fields["vz"],
                     "anum":fields['Anum'],      "znum":fields['Znum'],
                   "weight":fields["weight"],    "time":fields["weight"]*0 }
            return ("prt", out)

        if ("energy" in fieldNames and force is None) or force == "gc":
            # We have guiding centers (theta is random)
            warnings.warn(
                "Set time to zero and randomizing zeta for all markers.")
            zeta = 2*np.pi*np.random.rand(fields["id"].size)
            out = {     "n":fields["id"].size,    "ids":fields["id"],
                     "mass":fields["mass"],    "charge":fields["charge"],
                        "r":fields["R"],          "phi":fields["phi"],
                        "z":fields["z"],       "energy":fields["energy"],
                    "pitch":fields["pitch"],     "zeta":zeta,
                     "anum":fields['Anum'],      "znum":fields['Znum'],
                   "weight":fields["weight"],    "time":fields["weight"]*0 }
            return ("gc", out)

    def ascot4_wall2d(self, fn="input.wall_2d"):
        """Convert 2D wall data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 3D wall data.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        with open(fn, "r") as f:
            txt = np.loadtxt(f,skiprows=1)
            # "divflag" : txt[:,2] Not needed but just FYI
            out = {"nelements" : txt.shape[0], "r" : txt[:,0], "z" : txt[:,1]}
        return ("wall_2D", out)

    def ascot4_wall3d(self, fn="input.wall_3d", hdf5=False):
        """Convert 3D wall data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 3D wall data.
        hdf5 : bool, optional
            If True, the data is read from HDF5 file instead of
            the input.wall_3d text format.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        if hdf5:
            data = {}
            with h5py.File(fn, 'r') as f: # Open for reading
                data["x1x2x3"] = f["/wall/3d/triangles_x1x2x3"][:]
                data["y1y2y3"] = f["/wall/3d/triangles_y1y2y3"][:]
                data["z1z2z3"] = f["/wall/3d/triangles_z1z2z3"][:]
                data["flag"]   = f["/wall/3d/triangles_flag"][:]

            data["nelements"] = data["flag"].size
            data["flag"] = np.reshape(data["flag"], (data["flag"].size, 1))
            return ("wall_3D", data)

        # Crudely check number of lines to get maximum size for data
        num_lines = sum(1 for line in open(fn))
        copysector = False

        with open(fn, "r") as f:
            data = {}
            n_sectors = int(f.readline().split()[0])

            # Read sector ids
            sector_id = {}
            f.readline()
            for sector in range(n_sectors):
                line = f.readline().split()
                sector_id[int(float(line[0]))] = int(float(line[1]))

            # Read data for each sector
            n_read = 0

            data["x1x2x3"] = np.zeros((num_lines, 3))
            data["y1y2y3"] = np.zeros((num_lines, 3))
            data["z1z2z3"] = np.zeros((num_lines, 3))
            ids = np.zeros((num_lines, 1))
            # This is not filled (for lack of test-files)
            #data["flag"] = np.zeros((num_lines, 1))
            for sector in range(1, n_sectors + 1):
                f.readline() # Skip empty line
                try:
                    n_elements = int(float(f.readline().split()[0]))
                except Exception:
                    warnings.warn(
                        "There are multiple sectors but data only for one "\
                        + "sector. The sector is copied and rotated." )
                    copysector = True
                    break
                for i in range(n_elements):
                    f.readline() # Skip header
                    f.readline() # Skip empty line
                    f.readline() # Skip element info
                    for j in range(3):
                        line = f.readline().split()
                        data["x1x2x3"][n_read,j] = line[0]
                        data["y1y2y3"][n_read,j] = line[1]
                        data["z1z2z3"][n_read,j] = line[2]
                    ids[n_read] = sector_id[sector]
                    n_read = n_read + 1

                    # Some data may contain rectangles instead of triangles.
                    # This just means we read the "extra" vertice and make
                    # another triangle from it.
                    pos  = f.tell()
                    line = f.readline()
                    if len(line) == 0 or len(line) == 1:
                        # This was an empty line (no vertice data). Return to
                        # previous line.
                        f.seek(pos)
                    else:
                        line = line.split()
                        data["x1x2x3"][n_read,0] = line[0]
                        data["y1y2y3"][n_read,0] = line[1]
                        data["z1z2z3"][n_read,0] = line[2]
                        data["x1x2x3"][n_read,1] = data["x1x2x3"][n_read-1,0]
                        data["y1y2y3"][n_read,1] = data["y1y2y3"][n_read-1,0]
                        data["z1z2z3"][n_read,1] = data["z1z2z3"][n_read-1,0]
                        data["x1x2x3"][n_read,2] = data["x1x2x3"][n_read-1,2]
                        data["y1y2y3"][n_read,2] = data["y1y2y3"][n_read-1,2]
                        data["z1z2z3"][n_read,2] = data["z1z2z3"][n_read-1,2]
                        ids[n_read] = sector_id[sector]
                        n_read = n_read + 1

        # Remove extra zeros
        data["x1x2x3"] = data["x1x2x3"][0:n_read,:]
        data["y1y2y3"] = data["y1y2y3"][0:n_read,:]
        data["z1z2z3"] = data["z1z2z3"][0:n_read,:]
        ids = ids[0:n_read,:]

        if copysector:
            r1r2r3 = np.sqrt(data["x1x2x3"]**2 + data["y1y2y3"]**2)
            t1t2t3 = np.arctan2(data["y1y2y3"], data["x1x2x3"])
            z1z2z3 = data["z1z2z3"]
            data["x1x2x3"] = np.zeros( (n_read*n_sectors,3) )
            data["y1y2y3"] = np.zeros( (n_read*n_sectors,3) )
            data["z1z2z3"] = np.zeros( (n_read*n_sectors,3) )

            for i in range(n_sectors):
                x1x2x3 = r1r2r3 * np.cos(t1t2t3 + i*2*np.pi/n_sectors)
                y1y2y3 = r1r2r3 * np.sin(t1t2t3 + i*2*np.pi/n_sectors)
                data["x1x2x3"][n_read*i:n_read*(i+1),:] = x1x2x3
                data["y1y2y3"][n_read*i:n_read*(i+1),:] = y1y2y3
                data["z1z2z3"][n_read*i:n_read*(i+1),:] = z1z2z3

            ids = np.arange(1, n_read*n_sectors+1)

        data["nelements"] = ids.size
        #data["flag"] = np.reshape(data["flag"], (data["flag"].size, 1))
        return ("wall_3D", data)

    def ascot4_tokamak(self, fn="input.magn_bkg", fnheader="input.magn_header",
                       interpolate_psi0=False):
        """Convert tokamak magnetic field data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 magnetic field data.
        fnheader : str, optional
            Filename of the ASCOT4 magnetic field header.
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
        data = {}
        with open(fn) as fh:
            tmp = list(map(float,fh.readline().split()))
            phi0       = tmp[0]
            nSector    = int(tmp[1])
            nPhi       = int(tmp[2])
            nCoil      = int(tmp[3])
            zeroAtCoil = int(tmp[4])

            r1,r2,nr = [float(number) for number in fh.readline().split()]
            nr = int(nr)
            z1,z2,nz = [float(number) for number in fh.readline().split()]
            nz = int(nz)

            r = np.linspace(r1, r2, nr)
            z = np.linspace(z1, z2, nz)
            if(nPhi > 1):
                dphi = 360.0 / ( nSector * nPhi )
                phi = np.linspace(phi0 + dphi * 0.5,
                                  phi0 + 360.0 / nSector - dphi * 0.5,
                                  nPhi)

            phimap_tor = np.array(
                [int(float(number)) for number in fh.readline().split()])
            phimap_pol = np.array(
                [int(float(number)) for number in fh.readline().split()])

            h5data = np.array(fh.read().split(), dtype=float).flatten()

            sz2d = nr*nz
            sz3d = nr*nz*nPhi
            data["psi"]  = h5data[:sz2d].reshape(nz,nr).T / (2*np.pi)
            data["br"]   = h5data[sz2d+0*sz3d:sz2d+1*sz3d].\
                reshape(nz,nPhi,nr).squeeze().T
            data["bphi"] = h5data[sz2d+1*sz3d:sz2d+2*sz3d].\
                reshape(nz,nPhi,nr).squeeze().T
            data["bz"]   = h5data[sz2d+2*sz3d:sz2d+3*sz3d].\
                reshape(nz,nPhi,nr).squeeze().T

        data["axisr"], data["axisz"], data["psi0"], data["psi1"] = \
            self._ascot4_magn_header(fnheader)
        if nPhi == 1:
            data["nr"]   = nr
            data["rmin"] = r[0]
            data["rmax"] = r[-1]
            data["nz"]   = nz
            data["zmin"] = z[0]
            data["zmax"] = z[-1]
            return ("B_2DS", data)

        data["b_nr"]     = nr
        data["b_rmin"]   = r[0]
        data["b_rmax"]   = r[-1]
        data["b_nz"]     = nz
        data["b_zmin"]   = z[0]
        data["b_zmax"]   = z[-1]
        data["b_nphi"]   = nPhi
        data["b_phimin"] = phi[0]
        data["b_phimax"] = phi[-1]

        if interpolate_psi0:
            from scipy.interpolate import RegularGridInterpolator
            psi0 = RegularGridInterpolator((r, z), data['psi'])(
                (data['axisr'], data['axisz']))
            self._ascot.input_init(bfield=data)
            data["axisr"], data["axisz"], data["psi0"] = \
                self._ascot.input_findpsi0(data["psi0"])
            self._ascot.input_free()
        return ("B_3DS", data)

    def ascot4_stellarator(self, fn="ascot4.h5", psi0=None, psi1=None):
        """Convert stellarator magnetic field data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            HDF5 file with the ASCOT4 stellarator magnetic field data.
        psi0 : float, optional
            Poloidal flux on axis.
        psi1 : float, optional
            Poloidal flux at separatrix.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        if psi0 is None:
            psi0 = 0
            warnings.warn(
                "Using dummy value for psi0. Find proper value with ascotpy")
        if psi1 is None:
            psi1 = 1
            warnings.warn(
                "Using dummy value for psi1. Find proper value with ascotpy")
        data = {}
        with h5py.File(fn, "r") as f:
            f = f["bfield/stellarator"]

            r    = f["r"][:]
            phi  = f["phi"][:] * 180 / np.pi
            z    = f["z"][:]

            data["br"]   = f["br"][:]
            data["bphi"] = f["bphi"][:]
            data["bz"]   = f["bz"][:]
            data["psi"]  = f["s"][:]

            axis_r   = f["axis_R"][:]
            axis_phi = f["axis_phi"][:] * 180 / np.pi
            axis_z   = f["axis_z"][:]

            nperiods = f["toroidalPeriods"][:]
            if "symmetrymode" in f:
                symmetrymode = f["symmetrymode"][:]
            else:
                warnings.warn(
                    "No symmetry mode specified. Defaulting to stellarator "
                    + "symmetry")
                symmetrymode = 0

        philim = 360 / nperiods / (2 if symmetrymode == 0 else 1)
        if(symmetrymode == 0):
            # Data is in the format f(z, phi, r)
            phi = np.concatenate( (phi[:-1], phi[:] + phi[-1]) )

            t            = data["br"][:, :-1, :]
            t_sym        = -data["br"][-1::-1, -1:0:-1, :]
            data["br"]   = np.concatenate((t, t_sym),1)

            t            = data["bphi"][:, :-1, :]
            t_sym        = data["bphi"][-1::-1, -1:0:-1, :]
            data["bphi"] = np.concatenate((t, t_sym),1)

            t            = data["bz"][:, :-1, :]
            t_sym        = data["bz"][-1::-1, -1:0:-1, :]
            data["bz"]   = np.concatenate((t, t_sym),1)

            t            = data["psi"][:, :-1, :]
            t_sym        = data["psi"][-1::-1, -1:0:-1, :]
            data["psi"]  = np.concatenate((t, t_sym),1)
        elif (phi[0] == np.mod(phi[-1], philim)):
            warnings.warn("Removing duplicate field data point")
            data["br"]   = data["br"][:, :-1, :]
            data["bphi"] = data["bphi"][:, :-1, :]
            data["bz"]   = data["bz"][:, :-1, :]
            data["psi"]  = data["psi"][:, :-1, :]
        if (axis_phi[0] == np.mod(axis_phi[-1], philim)):
            warnings.warn("Removing duplicate axis data point")
            axis_r = axis_r[:-1]
            axis_z = axis_z[:-1]

        # Transpose to Ascot5 format
        data["br"]   = np.transpose(data["br"],   (2,1,0))
        data["bphi"] = np.transpose(data["bphi"], (2,1,0))
        data["bz"]   = np.transpose(data["bz"],   (2,1,0))
        data["psi"]  = np.transpose(data["psi"],  (2,1,0))

        data["b_rmin"]   = r[0]
        data["b_rmax"]   = r[-1]
        data["b_nr"]     = r.size
        data["b_phimin"] = phi[0]
        data["b_phimax"] = phi[-1]
        data["b_nphi"]   = phi.size - 1
        data["b_zmin"]   = z[0]
        data["b_zmax"]   = z[-1]
        data["b_nz"]     = z.size

        data["psi0"] = psi0
        data["psi1"] = psi1
        data["axis_phimin"] = axis_phi[0]
        data["axis_phimax"] = axis_phi[-1]
        data["axis_nphi"]   = axis_phi.size-1
        data["axisr"]       = axis_r
        data["axisz"]       = axis_z

        return ("B_STS", data)

    def ascot4_erad(self, fn="input.erad", reff=1.0):
        """Convert 1D electric field data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 electric field data.
        reff : float, optional
            Effective minor radius of the plasma used to convert ``dv/drho`` to
            SI units as ``drho/dr=1/reff`` [m].

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        data = {}
        with open(fn, "r") as f:
            f.readline() # Skip comment line
            data["nrho"] = int(float(f.readline().split()[0]))

            h5data = np.loadtxt(f)

            rho = h5data[:,0]
            data["rhomin"]  = h5data[0,0]
            data["rhomax"]  = h5data[-1,0]
            data["dvdrho"] = h5data[:,1]
            # For data in format dV/rho, we can ignore effective minor radius

        # Make sure the input is linearly spaced. If not, interpolate
        tol = 1.0001
        diff = np.diff(rho)
        if ( max(diff)/min(diff) > tol):
            warnings.warn("Interpolating dV/drho to uniform grid")
            new_rho = np.linspace(np.amin(rho), np.amax(rho), data["nrho"])
            data["dvdrho"] = np.interp(new_rho, rho, data["dvdrho"])

        data["reff"] = reff
        return ("E_1DS", data)

    def ascot4_neutral1d(self, fn="input.neutral1d", anum=1, znum=1):
        """Convert 1D neutral data from ASCOT4 to ASCOT5.

        Note that ASCOT4 input is always for a single species.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 neutral data.
        anum : int, optional
            Atomic mass number.
        znum : int, optional
            Atomic mass number.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        data = {}
        data["nspecies"] = 1 #ASCOT4 doesn't support several neutral species
        data["anum"]     = np.ones((1,1) * anum)
        data["znum"]     = np.ones((1,1) * znum)
        with open(fn, "r") as fh:
            fh.readline()
            fh.readline()
            fh.readline()
            data["nrho"] = int(float(fh.readline()))
            data["density"]     = np.zeros((data["nrho"],1))
            data["temperature"] = np.zeros((data["nrho"],1))

            fh.readline() # ignore headers
            h5data = np.loadtxt(fh)
            rho    = np.array(h5data[:,0])
            data["density"][:,0]     = np.array(h5data[:,1])
            data["temperature"][:,0] = np.array(h5data[:,2])

            # Make sure the input is linearly spaced. If not, interpolate
            tol = 1.0001
            diff = np.diff(rho)
            if ( max(diff)/min(diff) > tol ):
                warnings.warn("Interpolating neutral data to uniform grid")
                new_rho = np.linspace(rho[0], rho[-1], data["nrho"])
                data["density"][:,0] = np.interp(
                    new_rho, rho, data["density"][:,0])
                data["temperature"][:,0] = np.interp(
                    new_rho, rho, data["temperature"][:,0])

            # Add extra data point outside rho=1 to avoid out of data range
            # errors
            if ( np.amax(rho) <= 1.0 ):
                warnings.warn("Adding small datapoint outside rho=1.0")
                data["nrho"] = data["nrho"] + 1
                rho = np.append(rho, 2*rho[-1]-rho[-2])
                data["density"] = np.append(
                    data["density"], data["density"][-1,:], axis=0)
                data["temperature"] = np.append(
                    data["temperature"], data["temperature"][-1,:], axis=0)

            data["rhomin"] = rho[0]
            data["rhomax"] = rho[-1]
        return ("N0_1D", data)

    def ascot4_alfven(self, fn="input.alfven"):
        """Convert Alfvén MHD data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 Alfvén eigenmode data.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        data = {}
        with open(fn) as fh:

            # Skip first three lines
            fh.readline()
            fh.readline()
            fh.readline()

            # AEs on (1) or off (0)
            ae = int(fh.readline().split()[0]) == 1

            # Total number of different modes
            data["nmode"] = int(fh.readline().split()[0])

            # Number of radial points for profiles
            data["nrho"] = int(fh.readline().split()[0])

            # Skip empty line
            fh.readline()

            # Poloidal mode numbers
            data["mmodes"] = np.array(fh.readline().split()[:data["nmode"]])
            data["mmodes"] = -data["mmodes"].astype("i4")

            # Toroidal mode numbers
            data["nmodes"] = np.array(fh.readline().split()[:data["nmode"]])
            data["nmodes"] = data["nmodes"].astype("i4")

            # Amplitudes
            data["amplitude"] = np.array(fh.readline().split()[:data["nmode"]])
            data["amplitude"] = data["amplitude"].astype("f8")

            # Angular frequencies (omega) [rad/s]
            data["omega"] = np.array(fh.readline().split()[:data["nmode"]])
            data["omega"] = data["omega"].astype("f8")

            # Phase not given, fix it at zero
            warnings.warn("Phase not given, setting it to zero for all modes")
            data["phase"] = data["omega"]*0

            # psin, alpha profile, phi profile, each line corresponds to one psi
            line = fh.readline()
            psi   = np.zeros( (data["nrho"],1) )
            alpha = np.zeros( (data["nrho"],data["nmode"]) )
            phi   = np.zeros( (data["nrho"],data["nmode"]) )
            line = fh.readline()

            # File ends in #EOF
            for i in range(data["nrho"]):
                line = line.split()
                psi[i]     = line[0]
                alpha[i,:] = line[1:data["nmode"]+1]
                phi[i,:]   = line[data["nmode"]+1:2*data["nmode"]+1]
                line = fh.readline()

            # Set data on uniform grid
            data["rhomin"] = psi[0]
            data["rhomax"] = psi[-1]

            psiq = np.linspace(data["rhomin"], data["rhomax"], data["nrho"])

            data["alpha"] = np.zeros( (data["nrho"],data["nmode"]) )
            data["phi"]   = np.zeros( (data["nrho"],data["nmode"]) )
            for i in range(data["nmode"]):
                f = scipy.interpolate.interp1d(psi.ravel(), alpha[:,i])
                data["alpha"][:,i] = f(psiq).ravel()

                f = scipy.interpolate.interp1d(psi.ravel(), phi[:,i])
                data["phi"][:,i] = f(psiq).ravel()

            for i in range(data["nmode"]):
                f = scipy.interpolate.interp1d(psiq.ravel(), data["alpha"][:,i])
                data["alpha"][:,i] = f(psiq*psiq).ravel()

                f = scipy.interpolate.interp1d(psiq.ravel(), data["phi"][:,i])
                data["phi"][:,i] = f(psiq*psiq).ravel()

        return ("MHD_STAT", data)

    def ascot4_plasma1d(self, fn="input.plasma_1d"):
        """Convert 1D plasma data from ASCOT4 to ASCOT5.

        Parameters
        ----------
        fn : str, optional
            Filename of the ASCOT4 plasma data.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        data = {}
        with open(fn) as fh:
            fh.readline()
            fh.readline()
            fh.readline()
            line = fh.readline().split()[:2]
            nrho = int(float(line[0]))
            nion = int(float(line[1]))
            line = fh.readline().split()
            data["znum"] = np.array(
                [int(float(znum)) for znum in line[:nion]])
            line = fh.readline().split()
            data["anum"] = np.array(
                [int(float(anum)) for anum in line[:nion]])
            fh.readline() # Ignore collisions on-off line
            fh.readline() # ignore headers
            h5data = np.loadtxt(fh)
            rho = np.array(h5data[:,0])
            te  = np.array(h5data[:,1])
            ne  = np.array(h5data[:,2])
            #data["vtor"] = np.array(h5data[:,3]) Rotation not yet implemented
            ti  = np.array(h5data[:,4])
            ni  = np.zeros((nion,nrho))
            for i in range(nion):
                ni[i,:] = np.array(h5data[:,5+i])
            # Make sure the input is linearly spaced. If not, interpolate
            tol = 1.0001
            diff = np.diff(rho)
            if ( max(diff)/min(diff) > tol):
                warnings.warn("Interpolating plasma data to uniform grid")
                new_rho = np.linspace(np.amin(rho), np.amax(rho), nrho)
                ne = np.interp(new_rho, rho, ne)
                te = np.interp(new_rho, rho, te)
                for i in range(nion):
                    ni[i,:] = np.interp(new_rho, rho, ni[i,:])
                ti  = np.interp(new_rho, rho, ti)
                rho = new_rho

            ni = np.transpose(ni)
            # Add extra data point outside rho=1 to avoid out of data range
            # errors
            if ( np.amax(rho) <= 1.0 ):
                warnings.warn("Adding small datapoint outside rho=1.0")
                nrho = nrho + 1
                rho  = np.append(rho, 2*rho[-1]-rho[-2])
                ne   = np.append(ne, ne[-1]*1e-10)
                te   = np.append(te, te[-1])
                ni   = np.append(ni, np.expand_dims(ni[-1,:]*1e-10, 1).T,
                                 axis=0)
                ti   = np.append(ti, ti[-1])

        warnings.warn(
            "Ascot4 data does not contain masses and charges explicitly. " + \
            "Assuming fully-ionized ions and that mass = Anum * amu")
        data.update({"nion":nion, "nrho":nrho, "rho":rho, "etemperature":te,
                     "itemperature":ti, "edensity":ne, "idensity":ni,
                     "mass":data["anum"], "charge":data["znum"]})
        return("plasma_1D", data)

    def _ascot4_magn_header(self, fn):
        """Read data from input.magn_header files.

        Parameters
        ----------
        fn : str
            Filename for input.magn_header.

        Returns
        -------
        axisr : float
            Magnetic axis R coordinate.
        axisz : float
            Magnetic axis z coordinate.
        psi0 : float
            Poloidal flux at axis.
        psi1 : float
            Poloidal flux at separatrix.
        """
        with open(fn) as fh:
            # first four lines contain nothing interesting
            fh.readline()
            fh.readline()
            fh.readline()
            fh.readline()

            # Next three lines contain axis psi, R, and z values
            tmp  = [float(number) for number in fh.readline().split()]
            psi0 = tmp[0] / (2*np.pi)
            psi1 = tmp[1] / (2*np.pi)

            tmp   = [float(number) for number in fh.readline().split()]
            axisr = tmp[0]

            tmp   = [float(number) for number in fh.readline().split()]
            axisz = tmp[0]

        return axisr, axisz, psi0, psi1
