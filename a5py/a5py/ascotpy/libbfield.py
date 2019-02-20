import numpy as np

from a5py.ascotpy.libascot import LibAscot

class LibBfield(LibAscot):

    def evaluate(self, R, phi, z, t, quantity, grid=False,
                 squeeze=[None, None, None, None]):
        R   = np.asarray(R).ravel()
        phi = np.asarray(phi).ravel()
        z   = np.asarray(z).ravel()
        t   = np.asarray(t).ravel()

        if grid:
            arrsize = (R.size, phi.size, z.size, t.size)
            R, phi, z, t = np.meshgrid(R, phi, z, t)
        else:
            # Not a grid so check that dimensions are correct (and make
            # single-valued vectors correct size)
            arrsize = np.amax(np.array([R.size, phi.size, z.size, t.size]))
            errmsg = "Input arrays have inconsistent sizes ({}, {}, {}, {})"
            assert (R.size == 1 or R.size == arrsize) and  \
                (phi.size == 1 or phi.size == arrsize) and \
                (z.size == 1 or z.size == arrsize) and     \
                (t.size == 1 or t.size == arrsize),        \
                errmsg.format(R.size, phi.size, z.size, t.size)

            if R.size == 1:
                R = R[0]*np.ones((arrsize,))
            if phi.size == 1:
                phi = phi[0]*np.ones((arrsize,))
            if z.size == 1:
                z = z[0]*np.ones((arrsize,))
            if t.size == 1:
                t = t[0]*np.ones((arrsize,))

        out = None
        if quantity in ["rho", "psi"]:
            out = self.eval_bfield(R, phi, z, t, evalrho=True, evalpsi=True)[quantity]

        if quantity in ["br", "bphi", "bz", "brdr", "brdphi", "brdz", "bphidr",
                        "bphidphi", "bphidz", "bzdr", "bzdphi", "bzdz"]:
            out = self.eval_bfield(R, phi, z, t, evalb=True)[quantity]

        if quantity == "divergence":
            out = self.eval_bfield(R, phi, z, t, evalb=True)
            out = out["br"]/R + out["brdr"] + out["bphidphi"]/R + out["bzdz"]

        assert out is not None, "Unknown quantity"

        if grid:
            out = np.reshape(out, arrsize)
            idx = [slice(1), slice(1), slice(1), slice(1)]
            for i in range(len(squeeze)):
                if squeeze[i] is not None:
                    # This does not work yet
                    out = np.apply_along_axis(squeeze[i], i, out)
                else:
                    idx[i] = slice(None)

            out = out[idx[0], idx[1], idx[2], idx[3]]

        return out


    def evaluates(self, R, phi, z, t, quantity, ):

        if Rvec is None:
            Rvec = np.linspace(0.01, 10, 1000)
        if zvec is None:
            zvec = np.linspace()

        self.evaluate(xgrid, ygrid)
        pass


    def plot(self):
        pass
