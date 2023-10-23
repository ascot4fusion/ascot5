"""Functions for checking and transforming equilibrium COCOS.

https://www.sciencedirect.com/science/article/pii/S0010465512002962
"""
import numpy as np
import copy

COCOS_ASCOT = 3

def assign(q, ip, b0, psiaxis, psibndr, phiclockwise, weberperrad):
    """Assign correct COCOS for the given equilibrium.

    Parameters
    ----------
    q : float
        Safety factor (at any point, sign included).
    ip : float
        Plasma current (sign included).
    b0 : float
        Toroidal field (sign included).
    psiaxis : float
        Poloidal flux at the magnetic axis.
    psinbd : float
        Poloidal flux at the boundary.
    phiclockwise : bool
        If True, the toroidal angle increases clockwise when viewed from above.
    weberperrad : bool
        If True, the poloidal flux is not in Webers but Wb/rad (divided by 2pi).

    Returns
    -------
    cocos : int
        The COCOS number corresponding to this equilibrium.

    Raises
    ------
    ValueError
        If correct COCOS could not be determined.
    """
    sign_q  = np.sign(q)
    sign_ip = np.sign(ip)
    sign_b0 = np.sign(b0)
    cocos = set([1, 2, 3, 4, 5, 6, 7, 8])

    sigma_bp = np.sign(psibndr-psiaxis)/sign_ip
    if sigma_bp > 0:
        for i in [3,4,7,8]:
            cocos.discard(i)
    else:
        for i in [1,2,5,6]:
            cocos.discard(i)

    sigma_rhothetaphi = sign_q/(sign_ip*sign_b0)
    if sigma_rhothetaphi > 0:
        for i in [1,2,7,8]:
            cocos.discard(i)
    else:
        for i in [3,4,5,6]:
            cocos.discard(i)

    if phiclockwise:
        for i in [1,3,5,7]:
            cocos.discard(i)
    else:
        for i in [2,4,6,8]:
            cocos.discard(i)

    if len(cocos) > 1:
        raise ValueError("Could not determine COCOS")
    cocos = cocos.pop()
    if weberperrad:
        cocos += 10
    return cocos

def tococos3(eqd, cocos):
    """Transform equilibrium to COCOS3.

    Parameters
    ----------
    eqd : dict
        Dictionary from reading the EQDSK file.
    cocos : int
        Target COCOS.

    Returns
    -------
    eqdout : dict
        Equilibrium data converted to COCOS3.
    """

    def getparameters(cocos):
        """Retrieve sign conventions for a given COCOS.

        Parameters
        ----------
        cocos : int
            COCOS label.

        Returns
        -------
        sigma_Bp : int
            Is the sign of Bpol same as grad phi x grad psi (+1)
            or opposite (-1).
        sigma_Rphiz : int
            Are cylindrical coordinates right (+1) or left (-1) handed.
        sigma_rhothetaphi : int
            Are flux coordinates right (+1) or left (-1) handed.
        sigma_q_pos : int
            Is the sign of q same as dPsi_tor/dpsi (+1) or opposite (-1).
        sign_dpdpsi_pos : int
            The sign of dp/dpsi.
        exp_Bp : int
            Is the psi divided by 2pi when calculating Bpol or not.
        """
        # Replace the strings in this list with the corresponding values
        signs = ["sigma_Bp", "sigma_RphiZ", "sigma_rhothetaphi",
                 "sign_q_pos", "sign_pprime_pos", "exp_Bp"]
        if cocos > 10:
            signs[5] = 1
            cocos -= 10
        else:
            signs[5] = 0

        if cocos == 1:
            signs[:5] = [+1,+1,+1,+1,-1]
        elif cocos == 2:
            signs[:5] = [+1,-1,+1,+1,-1]
        elif cocos == 3:
            signs[:5] = [-1,+1,-1,-1,+1]
        elif cocos == 4:
            signs[:5] = [-1,-1,-1,-1,+1]
        elif cocos == 5:
            signs[:5] = [+1,+1,-1,-1,-1]
        elif cocos == 6:
            signs[:5] = [+1,-1,-1,-1,-1]
        elif cocos == 7:
            signs[:5] = [-1,+1,+1,+1,+1]
        elif cocos == 8:
            signs[:5] = [-1,-1,+1,+1,+1]

        return signs

    sigma_Bp, sigma_RphiZ, sigma_rhothetaphi, sign_q_pos, sign_pprime_pos, \
        exp_Bp = getparameters(cocos)

    #cocosin["sigma_ip"] = np.sign(eqd.Ip)
    #cocosin["sigma_b0"] = np.sign(eqd.B0EXP)

    #These cocos are for COCOS 3
    ascot_Bp, ascot_RphiZ, ascot_rhothetaphi, ascot_q_pos, ascot_pprime_pos, \
        ascot_Bp = getparameters(COCOS_ASCOT)

    #Checking the signs of the current and field desired as output
    #cocosout["sigma_ip"] = np.sign(eqd.Ip)    if sigma_ip == 0 else sigma_ip
    #cocosout["sigma_b0"] = np.sign(eqd.B0EXP) if sigma_b0 == 0 else sigma_b0

    # Define effective variables: sigma_Ip_eff, sigma_B0_eff, sigma_Bp_eff,
    # exp_Bp_eff as in Appendix C
    sigma_Bp_eff = sigma_Bp * ascot_Bp
    exp_Bp_eff   = ascot_Bp - exp_Bp
    sigma_rhothetaphi_eff = sigma_rhothetaphi * sigma_rhothetaphi
    #if 'sigma_ip' in cocosout.keys() and 'sigma_b0' in cocosout.keys():
    #    sigma_Ip_eff = cocosin['sigma_ip']*cocosout['sigma_ip']
    #    sigma_B0_eff = cocosin['sigma_b0']*cocosout['sigma_b0']
    #else:
    #    # No sign of Ip nor B0 requested
    sigma_Ip_eff = sigma_RphiZ * ascot_RphiZ
    sigma_B0_eff = sigma_RphiZ * ascot_RphiZ

    # Define input
    F_in       = eqd["fpol"]
    FFprime_in = eqd["ffprime"]
    pprime_in  = eqd["pprime"]

    psirz_in   = eqd["psi"]
    psiaxis_in = eqd["simagx"]
    psiedge_in = eqd["sibdry"]

    q_in  = eqd["qpsi"]
    b0_in = eqd["bcentr"]
    ip_in = eqd["cpasma"]

    # Transform
    F         = F_in       * sigma_B0_eff
    FFprime   = FFprime_in * sigma_Ip_eff * sigma_Bp_eff / (2*np.pi)**exp_Bp_eff
    pprime    = pprime_in  * sigma_Ip_eff * sigma_Bp_eff / (2*np.pi)**exp_Bp_eff

    _fact_psi = sigma_Ip_eff * sigma_Bp_eff * (2*np.pi)**exp_Bp_eff
    psirz     = psirz_in   * _fact_psi
    psiaxis   = psiaxis_in * _fact_psi
    psiedge   = psiedge_in * _fact_psi

    q  = q_in  * sigma_Ip_eff * sigma_B0_eff * sigma_rhothetaphi_eff
    b0 = b0_in * sigma_B0_eff
    ip = ip_in * sigma_Ip_eff

    # Define output
    eqdout = copy.deepcopy(eqd)
    eqdout["fpol"]    = F
    eqdout["ffprime"] = FFprime
    eqdout["pprime"]  = pprime

    eqdout["psi"]     = psirz
    eqdout["psiaxis"] = psiaxis
    eqdout["psiedge"] = psiedge

    eqdout["qpsi"]   = q
    eqdout["bcentr"] = b0
    eqdout["cpasma"] = ip
    return eqdout
