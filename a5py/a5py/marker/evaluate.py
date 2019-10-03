"""
Routines from evaluating derived quantities from marker fields.

File: evaluate.py
"""
import numpy           as np
import scipy.constants as const

from .alias import get as alias

def eval_particle(quantity, mass=None, charge=None,
                  R=None, phi=None, z=None,
                  vR=None, vphi=None, vz=None,
                  BR=None, Bphi=None, Bz=None, psi=None):
    """
    Evaluated particle quantities.

    All inputs must be in SI units.

    Args:
        quantity : str <br>
            Quantity to be evaluated.

    Returns:
        Evaluated quantity in SI units.

    Raises:
        AssertionError if not all arguments necessary to evaluate the quantity
        are given.
        ValueError if the quantity is unknown.
    """

    quantity = alias(quantity)

    if quantity in ["x", "y", "xprt", "yprt"]:
        assert(R is not None and phi is not None)

        if quantity in ["x", "xprt"]:
            return R * np.cos(phi)
        if quantity in ["y", "yprt"]:
            return R * np.sin(phi)

    if quantity in ["phimod", "phimodprt"]:
        return np.mod(phi + 2*np.pi, 2*np.pi)

    if quantity in ["energy", "energyprt"]:
        assert(mass is not None and vR is not None and vphi is not None and
               vz is not None)

        return 0.5 * mass * (vR*vR + vphi*vphi + vz*vz)

    if quantity in ["gamma", "gammaprt"]:
        assert(vR is not None and vphi is not None and vz is not None)

        return 1/ np.sqrt( 1 - (vR*vR + vphi*vphi + vz*vz)
                           / (const.c * const.c) )

    if quantity in ["pr", "pphi", "pz"]:
        assert(mass is not None and vR is not None and vphi is not None and
               vz is not None)

        g = eval_particle("gamma", vR=vR, vphi=vphi, vz=vz)

        if quantity == "pr":
            return g*mass*vR
        if quantity == "pphi":
            return g*mass*vphi
        if quantity == "pz":
            return g*mass*vz

    if quantity in ["pitch", "pitchprt"]:
        assert(vR is not None and vphi is not None and vz is not None and
               BR is not None and Bphi is not None and Bz is not None)

        vnorm = eval_particle("vnorm", vR=vR, vphi=vphi, vz=vz)
        Bnorm = eval_particle("bnorm", BR=BR, Bphi=Bphi, Bz=Bz)
        return (vR*BR + vphi*Bphi + vz*Bz) / (vnorm*Bnorm)

    if quantity in ["mu", "muprt"]:
        assert(vR is not None and vphi is not None and vz is not None and
               BR is not None and Bphi is not None and Bz is not None and
               mass is not None)

        pitch = eval_particle("pitch", vR=vR, vphi=vphi, vz=vz,
                              BR=BR, Bphi=Bphi, Bz=Bz)
        vnorm = eval_particle("vnorm", vR=vR, vphi=vphi, vz=vz)
        g     = eval_particle("gamma", vR=vR, vphi=vphi, vz=vz)
        Bnorm = eval_particle("bnorm", BR=BR, Bphi=Bphi, Bz=Bz)
        return mass * np.power(g*vnorm, 2) * (1 - pitch*pitch) / (2*Bnorm)

    if quantity in ["vpar", "vparprt"]:
        assert(vR is not None and vphi is not None and vz is not None and
               BR is not None and Bphi is not None and Bz is not None)

        pitch = eval_particle("pitch", vR=vR, vphi=vphi, vz=vz,
                              BR=BR, Bphi=Bphi, Bz=Bz)
        vnorm = eval_particle("vnorm", vR=vR, vphi=vphi, vz=vz)
        return pitch * vnorm

    if quantity in ["ppar", "pparprt"]:
        assert(vR is not None and vphi is not None and vz is not None and
               BR is not None and Bphi is not None and Bz is not None and
               mass is not None)

        vpar = eval_particle("vpar", vR=vR, vphi=vphi, vz=vz,
                             BR=BR, Bphi=Bphi, Bz=Bz)
        g = eval_particle("gamma", vR=vR, vphi=vphi, vz=vz)
        return g*mass*vpar

    if quantity in ["vnorm", "vnormprt"]:
        assert(vR is not None and vphi is not None and vz is not None)

        return np.sqrt(vR*vR + vphi*vphi + vz*vz)

    if quantity in ["pnorm", "pnormprt"]:
        assert(mass is not None and vR is not None and vphi is not None and
               vz is not None)

        pR   = eval_particle("pR", mass=mass, vR=vR, vphi=vphi, vz=vz)
        pphi = eval_particle("pphi", mass=mass, vR=vR, vphi=vphi, vz=vz)
        pz   = eval_particle("pz", mass=mass, vR=vR, vphi=vphi, vz=vz)
        return np.sqrt(pR*pR + pphi*pphi + pz*pz)

    if quantity in ["bnorm", "bnormprt"]:
        assert(BR is not None and Bphi is not None and Bz is not None)

        return np.sqrt(BR*BR + Bphi*Bphi + Bz*Bz)

    if quantity in ["ptor", "ptorprt"]:
        assert(R is not None and charge is not None and mass is not None and
               vR is not None and vphi is not None and vz is not None and
               psi is not None)

        g = eval_particle("gamma", vR=vR, vphi=vphi, vz=vz)
        return g * mass * R * vphi + charge * psi

    raise ValueError(quantity + " is not a valid quantity.")


def eval_guidingcenter(quantity, mass=None, charge=None,
                       R=None, phi=None, z=None,
                       vpar=None, mu=None, theta=None,
                       BR=None, Bphi=None, Bz=None, psi=None):
    """
    Evaluated guiding center quantities.

    All inputs must be in SI units.

    Args:
        quantity : str <br>
            Quantity to be evaluated.

    Returns:
        Evaluated quantity in SI units.

    Raises:
        AssertionError if not all arguments necessary to evaluate the quantity
        are given.
        ValueError if the quantity is unknown.
    """
    quantity = alias(quantity)

    if quantity in ["x", "y", "xgc", "ygc"]:
        assert(R is not None and phi is not None)

        if quantity in ["x", "xgc"]:
            return R * np.cos(phi)
        if quantity in ["y", "ygc"]:
            return R * np.sin(phi)

    if quantity in ["phimod", "phimodgc"]:
        return np.mod(phi + 2*np.pi, 2*np.pi)

    if quantity in ["energy", "energygc"]:
        assert(mass is not None and vpar is not None and mu is not None and
               BR is not None and Bphi is not None and Bz is not None)

        Bnorm = eval_guidingcenter("bnorm", BR=BR, Bphi=Bphi, Bz=Bz)
        return mu * Bnorm + mass * vpar*vpar / 2

    if quantity in ["gamma", "gammagc"]:
        assert(mass is not None and mu is not None and vpar is not None and
               BR is not None and Bphi is not None and Bz is not None)

        Bnorm = eval_guidingcenter("bnorm", BR=BR, Bphi=Bphi, Bz=Bz)
        return np.sqrt(
            ( 1 + (2 * mu * Bnorm) / (mass * const.c * const.c) )
            / (1 - vpar * vpar / (const.c * const.c) ) )

    if quantity in ["vnorm", "vnormgc"]:
        assert(mass is not None and vpar is not None and mu is not None and
               BR is not None and Bphi is not None and Bz is not None)

        g     = eval_guidingcenter("gamma", mass=mass, mu=mu, vpar=vpar,
                                   BR=BR, Bphi=Bphi, Bz=Bz)
        return np.sqrt(1-1/g**2) * const.c

    if quantity in ["pitch", "pitchgc"]:
        assert(mass is not None and vpar is not None and mu is not None and
               BR is not None and Bphi is not None and Bz is not None)

        vnorm = eval_guidingcenter("vnorm", mass=mass, mu=mu, vpar=vpar,
                                   BR=BR, Bphi=Bphi, Bz=Bz)
        return vpar/vnorm

    if quantity in ["ppar", "ppargc"]:
        assert(mass is not None and vpar is not None and mu is not None and
               BR is not None and Bphi is not None and Bz is not None)

        g    = eval_guidingcenter("gamma", mass=mass, mu=mu, vpar=vpar,
                                  BR=BR, Bphi=Bphi, Bz=Bz)
        return g * mass * vpar

    if quantity in ["pnorm", "pnormgc"]:
        assert(mass is not None and vpar is not None and mu is not None and
               BR is not None and Bphi is not None and Bz is not None)

        vnorm = eval_guidingcenter("vnorm", mass=mass, mu=mu, vpar=vpar,
                                   BR=BR, Bphi=Bphi, Bz=Bz)
        g     = eval_guidingcenter("gamma", mass=mass, mu=mu, vpar=vpar,
                                   BR=BR, Bphi=Bphi, Bz=Bz)
        return g * mass * vnorm

    if quantity in ["bnorm", "bnormgc"]:
        assert(BR is not None and Bphi is not None and Bz is not None)

        return np.sqrt(BR*BR + Bphi*Bphi + Bz*Bz)

    if quantity in ["ptor", "ptorgc"]:
        assert(R is not None and charge is not None and mass is not None and
               vpar is not None and mu is not None and psi is not None and
               BR is not None and Bphi is not None and Bz is not None)

        g = eval_guidingcenter("gamma", mass=mass, mu=mu, vpar=vpar,
                               BR=BR, Bphi=Bphi, Bz=Bz)
        return g * mass * R * vpar + charge * psi

    raise ValueError(quantity + " is not a valid quantity.")
