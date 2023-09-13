import os
import numpy as np

class ImportAdas():

    def import_adas(self,
                    fn_CX_DD0   = "input.sigmavCX_DD0",
                    fn_BMS_H0H  = "input.BMSsigmav_H0H",
                    fn_BMS_H0He = "input.BMSsigmav_H0He",
                    fn_BMS_H0C  = "input.BMSsigmav_H0C"):
        """Import data from ADAS files.

        Note: this method should be redone completely as right now it is just
        the script asigma_loc_read2hdf5.py. The problem with that script
        (and this method) is that the input files are expected to have in-house
        format. Instead, what we would like to have is that this script would
        take the ADAS files directly.

        Parameters
        ----------
        files : str or list [str]
            Folder where input files are located or a list of input files.

        Returns
        -------
        gtype : str
            Type of the generated input data.
        data : dict
            Input data that can be passed to ``write_hdf5`` method of
            a corresponding type.
        """
        # Count the number of reactions for which data files can be found,
        # and account for possible multiples due to isotope combinations
        N_reac = 0
        CX_DD0_found = 0
        # Name HH0 instead of DD0 since we may, in the future, use this variable
        # for the total number of hydrogen-hydrogen CX reactions when different
        # isotope combinations are included.
        N_CX_HH0 = 0
        # NOTE: This list of BMS data filenames must be up to date with the
        # function parameter list.
        fn_BMSs = [fn_BMS_H0H,fn_BMS_H0He,fn_BMS_H0C]
        N_fn_BMS = len(fn_BMSs)
        BMS_founds = [0]*N_fn_BMS
        N_BMS = 0
        if(os.path.isfile(fn_CX_DD0)):
            CX_DD0_found = 1
            N_CX_HH0 = 1
        for i_BMS in range(0,N_fn_BMS):
            if(os.path.isfile(fn_BMSs[i_BMS])):
                BMS_founds[i_BMS] = 1
                N_BMS += 1
        N_reac = N_CX_HH0 + N_BMS

        # Initialize arrays for reaction identifiers and abscissae
        z_1       = np.zeros(N_reac, dtype=int)
        a_1       = np.zeros(N_reac, dtype=int)
        z_2       = np.zeros(N_reac, dtype=int)
        a_2       = np.zeros(N_reac, dtype=int)
        reac_type = np.zeros(N_reac, dtype=int)
        N_E       = np.zeros(N_reac, dtype=int)
        E_min     = np.zeros(N_reac, dtype=float)
        E_max     = np.zeros(N_reac, dtype=float)
        N_n       = np.zeros(N_reac, dtype=int)
        n_min     = np.zeros(N_reac, dtype=float)
        n_max     = np.zeros(N_reac, dtype=float)
        N_T       = np.zeros(N_reac, dtype=int)
        T_min     = np.zeros(N_reac, dtype=float)
        T_max     = np.zeros(N_reac, dtype=float)
        N         = np.zeros(N_reac, dtype=int)

        # Read atomic reaction data from local files
        #print("Reading atomic reaction data from local files")
        i_reac = 0
        # Read CX reaction data of the form sigmav(E,T)
        #print("  Trying to read D-D0 CX rate coefficient data")
        if(CX_DD0_found):
            f = open(fn_CX_DD0, "r")
            if f.mode == 'r':
                print("    Reading file " + fn_CX_DD0)
                lines = f.readlines()
                metadata = lines[0].split(' ')
                data_format = int(metadata[0]) # Not used
                z_1[i_reac]       = int(metadata[1])
                a_1[i_reac]       = int(metadata[2])
                z_2[i_reac]       = int(metadata[3])
                a_2[i_reac]       = int(metadata[4])
                reac_type[i_reac] = int(metadata[5])
                N_E[i_reac]       = int(metadata[6])
                N_n[i_reac]       = 1
                N_T[i_reac]       = int(metadata[7])
                N[i_reac]         = N_E[i_reac]*N_n[i_reac]*N_T[i_reac]
                abscissa1 = lines[1].split(' ')
                abscissa2 = lines[2].split(' ')
                ordinate  = lines[3].split(' ')
                # Currently, the ASCOT5 spline interpolation module only
                # supports uniform grid spacing. Hence, only min, max and N
                # is needed to reconstruct the abscissa on the ASCOT5 side.
                E_min[i_reac] = float(abscissa1[0])
                E_max[i_reac] = float(abscissa1[N_E[i_reac]-1])
                T_min[i_reac] = float(abscissa2[0])
                T_max[i_reac] = float(abscissa2[N_T[i_reac]-1])
                # Read the rate coefficient data
                sigmavCX_DD0 = np.zeros(N[i_reac])
                for i in range(0, N[i_reac]):
                    sigmavCX_DD0[i] = float(ordinate[i])
                # Increment the reaction index
                i_reac += 1
        else:
            print("    File not found: " + fn_CX_DD0)
        # Read BMS reaction data of the form BMSsigmav(E, n, T)
        #print("  Trying to read BMS data")
        # Initialize a list that can be filled with lists of unknown length below
        BMS_sigmavs = []
        i_BMS = 0
        for i_fn_BMS in range(0,N_fn_BMS):
            if(BMS_founds[i_fn_BMS]):
                f = open(fn_BMSs[i_fn_BMS], "r")
                if f.mode == 'r':
                    #print("    Reading file " + fn_BMSs[i_fn_BMS])
                    lines = f.readlines()
                    metadata = lines[0].split(' ')
                    data_format = int(metadata[0]) # Overwritten, because not used
                    z_1[i_reac]       = int(metadata[1])
                    a_1[i_reac]       = int(metadata[2])
                    z_2[i_reac]       = int(metadata[3])
                    a_2[i_reac]       = int(metadata[4])
                    reac_type[i_reac] = int(metadata[5])
                    N_E[i_reac]       = int(metadata[6])
                    N_n[i_reac]       = int(metadata[7])
                    N_T[i_reac]       = int(metadata[8])
                    N[i_reac]         = N_E[i_reac]*N_n[i_reac]*N_T[i_reac]
                    abscissa1 = lines[1].split(' ')
                    abscissa2 = lines[2].split(' ')
                    abscissa3 = lines[3].split(' ')
                    ordinate  = lines[4].split(' ')
                    # Currently, the ASCOT5 spline interpolation module only
                    # supports uniform grid spacing. Hence, only min, max and N is
                    # needed to reconstruct the abscissa on the ASCOT5 side.
                    # The density unit is converted (cm^-3 --> m^-3).
                    E_min[i_reac] =     float(abscissa1[0])
                    E_max[i_reac] =     float(abscissa1[N_E[i_reac]-1])
                    n_min[i_reac] = 1e6*float(abscissa2[0])
                    n_max[i_reac] = 1e6*float(abscissa2[N_n[i_reac]-1])
                    T_min[i_reac] =     float(abscissa3[0])
                    T_max[i_reac] =     float(abscissa3[N_T[i_reac]-1])
                    # Read the BMS coefficients and convert
                    # the unit (cm^3s^-1 --> m^3s^-1)
                    BMS_sigmavs.append(ordinate[0:N[i_reac]])
                    for i in range(0, N[i_reac]):
                        BMS_sigmavs[i_BMS][i] = 1e-6*float(BMS_sigmavs[i_BMS][i])
                    # Increment the BMS reaction index
                    i_BMS += 1
                    # Increment the reaction index
                    i_reac += 1
            else:
                print("    File not found: " + fn_BMSs[i_fn_BMS])

        # Move reaction data to one long array
        N_tot = 0
        for i_reac in range(0, N_reac):
            N_tot += N[i_reac]
        sigma = np.zeros((1,N_tot))
        i_reac  = 0
        i_sigma = 0
        if(CX_DD0_found):
            for i in range(0, N[i_reac]):
                sigma[0,i_sigma+i] = sigmavCX_DD0[i]
            i_sigma += N[i_reac]
            i_reac  += 1
        # We write no 'else' for this 'if' here because we already did above
        i_BMS = 0
        for i_fn_BMS in range(0,N_fn_BMS):
            if(BMS_founds[i_fn_BMS]):
                for i in range(0, N[i_reac]):
                    sigma[0,i_sigma+i] = BMS_sigmavs[i_BMS][i]
                i_sigma += N[i_reac]
                i_BMS   += 1
                i_reac  += 1

        out = {
            "nreac" : N_reac, "z1" : z_1, "a1" : a_1, "z2" : z_2, "a2" : a_2,
            "reactype" : reac_type, "nenergy" : N_E, "energymin" : E_min,
            "energymax" : E_max, "ndensity" : N_n, "densitymin" : n_min,
            "densitymax" : n_max, "ntemperature" : N_T,
            "temperaturemin" : T_min, "temperaturemax" : T_max, "sigma" : sigma
        }
        return ("asigma_loc", out)
