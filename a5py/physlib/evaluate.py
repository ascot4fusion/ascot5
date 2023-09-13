import numpy as np


def evaluate(quantity, m=None, q=None, r=None, phi=None, z=None, pvec=None,
             mu=None, ppar=None, zeta=None, b=None):
    """
    Evaluate derived quantities from the given quantities.

    Note that it is enough to pass enough parameters that are needed to evaluate
    the queried quantity.
    """

    if quantity in ["xprt", "xgc"]:
        return coords.xcoord(r=r, phi=phi)

    if quantity in ["yprt", "ygc"]:
        return coords.ycoord(r=r, phi=phi)

    if quantity in ["energyprt"]:
        return energy_momentum(m=m, p=pvec)

    if quantity in ["energygc"]:
        return energy_muppar(m=m, mu=mu, ppar=ppar, b=b)

    if quantity in ["pitchprt"]:
        return pitch_momentum(p=pvec, b=b)

    if quantity in ["pitchgc"]:
        return pitch_muppar(m=m, mu=mu, ppar=ppar, b=b)

    if quantity in ["muprt"]:
        return mu_momentum(m=m, p=pvec, b=b)

    if quantity in ["pparprt"]:
        return mu_momentum(m=m, p=pvec, b=b)

    if quantity in ["vparprt"]:
        return mu_momentum(m=m, p=pvec, b=b)

    if quantity in ["pnormprt"]:
        return np.sqrt( np.sum(pvec**2, axis=1) )

    if quantity in ["vvecprt", "vrprt", "vphiprt", "vzprt", "vnormprt"]:
        vvec = velocity_momentum(pvec)
        if quantity == "vvecprt":
            return vvec
        if quantity == "vrprt":
            return vvec[0]
        if quantity == "vphiprt":
            return vvec[1]
        if quantity == "vzprt":
            return vvec[2]
        if quantity == "vnormprt":
            return np.sqrt( np.sum(vvec**2) )

    if quantity in ["vrgc", "vphigc", "vzgc", "vnormgc", "vvecgc"]:
        if quantity == "vvecgc":
            return vvec
        if quantity == "vrgc":
            return vvec[0]
        if quantity == "vphigc":
            return vvec[1]
        if quantity == "vzgc":
            return vvec[2]
        if quantity == "vnormgc":
            return np.sqrt( np.sum(vvec**2) )
