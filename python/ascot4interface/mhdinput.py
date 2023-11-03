import numpy as np
from scipy import interpolate

def read_alfven(fn):
    data = {}
    with open(fn) as fh:

        # Skip first three lines
        fh.readline()
        fh.readline()
        fh.readline()

        # AEs on (1) or off (0)
        data["ae"] = int(fh.readline().split()[0]) == 1

        # Total number of different modes
        data["nmode"] = int(fh.readline().split()[0])

        # Number of radial points for profiles
        data["nrho"] = int(fh.readline().split()[0])

        # Skip empty line
        fh.readline()

        # Poloidal mode numbers
        data["mmodes"] = np.array(fh.readline().split()[:data["nmode"]])
        data["mmodes"] = data["mmodes"].astype(int)

        # Toroidal mode numbers
        data["nmodes"] = np.array(fh.readline().split()[:data["nmode"]])
        data["nmodes"] = data["nmodes"].astype(int)

        # Amplitudes
        data["amplitude"] = np.array(fh.readline().split()[:data["nmode"]])
        data["amplitude"] = data["amplitude"].astype(np.float)

        # Angular frequencies (omega) [rad/s]
        data["omega"] = np.array(fh.readline().split()[:data["nmode"]])
        data["omega"] = data["omega"].astype(np.float)

        # Phase not given, fix it at zero
        data["phase"] = data["omega"]*0

        # psin, alpha profile, phi profile, each line corresponds to one psislot
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
            f = interpolate.interp1d(psi.ravel(), alpha[:,i])
            data["alpha"][:,i] = f(psiq).ravel()

            f = interpolate.interp1d(psi.ravel(), phi[:,i])
            data["phi"][:,i] = f(psiq).ravel()

        for i in range(data["nmode"]):
            f = interpolate.interp1d(psiq.ravel(), data["alpha"][:,i])
            data["alpha"][:,i] = f(psiq*psiq).ravel()

            f = interpolate.interp1d(psiq.ravel(), data["phi"][:,i])
            data["phi"][:,i] = f(psiq*psiq).ravel()

    return data
