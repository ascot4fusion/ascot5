"""
Methods to evaluate quantities from boozer data.

File: libbbozer.py
"""

class LibBoozer():

    quantities = ["psi (bzr)", "theta", "zeta",
                  "dpsidr (bzr)", "dpsidphi (bzr)", "dpsidz (bzr)",
                  "dthetadr", "dthetadphi", "dthetadz",
                  "dzetadr", "dzetadphi", "dzetadz", "rho (bzr)", "qprof",
                  "jacobian", "jacobianb2"]

    def evaluate(self, R, phi, z, t, quantity, br=None, bphi=None, bz=None):

        out = None
        if quantity == "psi (bzr)":
            out = self.eval_boozer(R, phi, z, t)["psi"]
        elif quantity == "rho (bzr)":
            out = self.eval_boozer(R, phi, z, t)["rho"]
        elif quantity == "dpsidr (bzr)":
            out = self.eval_boozer(R, phi, z, t)["dpsidr"]
        elif quantity == "dpsidphi (bzr)":
            out = self.eval_boozer(R, phi, z, t)["dpsidphi"]
        elif quantity == "dpsidz (bzr)":
            out = self.eval_boozer(R, phi, z, t)["dpsidz"]
        elif quantity in ["theta", "zeta",
                          "dthetadr", "dthetadphi", "dthetadz", "dzetadr",
                          "dzetadphi", "dzetadz"]:
             out = self.eval_boozer(R, phi, z, t)[quantity]
        elif quantity in ["qprof", "jacobian", "jacobianb2"]:
            out = self.eval_boozer(R, phi, z, t, evalfun=True)[quantity]

        assert out is not None, "Unknown quantity"

        return out
