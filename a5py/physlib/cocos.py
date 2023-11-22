"""Functions for checking and transforming equilibrium COCOS.

https://www.sciencedirect.com/science/article/pii/S0010465512002962
"""
import numpy as np
import copy

COCOS_ASCOT = 3

class COCOS:
    """
    COCOS class

    `cocos::Int`           = COCOS ID number\\
    `exp_Bp::Int`          = 0 or 1, depending if psi is already divided by 2pi or not, respectively\\
    `sigma_Bp::Int`        = +1 or -1, depending if psi is increasing or decreasing with Ip and B0 positive\\
    `sigma_RpZ::Int`       = +1 or -1, depending if (R,phi,Z) is right-handed or (R,Z,phi), respectively\\
    `sigma_rhotp::Int`     = +1 or -1, depending if (rho, theta, phi) is right-handed or (rho,phi,theta), repectively\\
    `sign_q_pos::Int`      = +1 or -1, depending if q is positive or negative with Ip and B0 positive\\
    `sign_pprime_pos::Int` = +1 or -1, depending if dp/dpsi is positive or negative with Ip and B0 positive
    """
    def __init__(self, cocos, exp_Bp, sigma_Bp, sigma_RpZ, sigma_rhotp, sign_q_pos, sign_pprime_pos):
        self.cocos = cocos
        self.exp_Bp = exp_Bp
        self.sigma_Bp = sigma_Bp
        self.sigma_RpZ = sigma_RpZ
        self.sigma_rhotp = sigma_rhotp
        self.sign_q_pos = sign_q_pos
        self.sign_pprime_pos = sign_pprime_pos

def cocos(cocos_in):
    exp_Bp = 1 if cocos_in >= 11 else 0

    if cocos_in in (1, 11):
        # ITER, Boozer are COCOS=11
        return COCOS(cocos_in, exp_Bp, 1, 1, 1, 1, -1)
    elif cocos_in in (2, 12):
        # CHEASE, ONETWO, Hinton-Hazeltine, LION is COCOS=2
        return COCOS(cocos_in, exp_Bp, 1, -1, 1, 1, -1)
    elif cocos_in in (3, 13):
        # Freidberg, CAXE, KINX, EFIT are COCOS=3
        # EU-ITM up to end of 2011 is COCOS=13
        return COCOS(cocos_in, exp_Bp, -1, 1, -1, -1, 1)
    elif cocos_in in (4, 14):
        return COCOS(cocos_in, exp_Bp, -1, -1, -1, -1, 1)
    elif cocos_in in (5, 15):
        return COCOS(cocos_in, exp_Bp, 1, 1, -1, -1, -1)
    elif cocos_in in (6, 16):
        return COCOS(cocos_in, exp_Bp, 1, -1, -1, -1, -1)
    elif cocos_in in (7, 17):
        return COCOS(cocos_in, exp_Bp, -1, 1, 1, 1, 1)
    elif cocos_in in (8, 18):
        return COCOS(cocos_in, exp_Bp, -1, -1, 1, 1, 1)
    else:
        raise ValueError(f"COCOS = {cocos_in} does not exist")

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
    if sigma_rhothetaphi < 0:
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

def transform_cocos(cc_in: COCOS, cc_out: COCOS,
                    sigma_Ip = None,
                    sigma_B0 = None,
                    ld = (1, 1),
                    lB = (1, 1),
                    exp_mu0 = (0, 0)):
    """
    Returns a dictionary of the multiplicative factors to transform COCOS from `cc_in` to `cc_out`

    Parameters
    ----------
    sigma_Ip : Union[Tuple[int, int], None]
        A tuple of the (Input, Output) current sign or nothing
    sigma_B0 : Union{NTuple{2,Int},Nothing}` - A tuple of the (Input, Output) toroidal field sign or nothing
    ld : NTuple{2,<:Real}
        A tuple of the (Input, Output) length scale factor. Default = (1,1)
    lB : NTuple{2,<:Real}
        A tuple of the (Input, Output) magnetic field scale factor. Default = (1,1)
    exp_mu0 : NTuple{2,<:Real}
        A tuple of the (Input, Output) mu0 exponent (0, 1). Default = (0,0)

    Returns
    -------
    transforms : dict
        Transform multiplicative factors to be able to convert from `cc_in to `cc_out`
    """

    ld_eff = ld[1] / ld[0]
    lB_eff = lB[1] / lB[0]
    exp_mu0_eff = exp_mu0[1] - exp_mu0[0]

    sigma_RpZ_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ

    if sigma_Ip is None:
        sigma_Ip_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ
    else:
        sigma_Ip_eff = sigma_Ip[0] * sigma_Ip[1]

    if sigma_B0 is None:
        sigma_B0_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ
    else:
        sigma_B0_eff = sigma_B0[0] * sigma_B0[1]

    sigma_Bp_eff = cc_in.sigma_Bp * cc_out.sigma_Bp
    exp_Bp_eff = cc_out.exp_Bp - cc_in.exp_Bp
    sigma_rhotp_eff = cc_in.sigma_rhotp * cc_out.sigma_rhotp

    mu0 = 4 * 3.14159265358979323846 * 1e-7  # pi is used directly for more precision

    transforms = {}
    transforms["R"] = ld_eff
    transforms["Z"] = ld_eff
    transforms["PRES"] = (lB_eff ** 2) / (mu0 ** exp_mu0_eff)
    transforms["PSI"] = lB_eff * (ld_eff ** 2) * sigma_Ip_eff * sigma_Bp_eff * ((2 * 3.14159265358979323846) ** exp_Bp_eff) * (ld_eff ** 2) * lB_eff
    transforms["TOR"] = lB_eff * (ld_eff ** 2) * sigma_B0_eff
    transforms["PPRIME"] = (lB_eff / ((ld_eff ** 2) * (mu0 ** exp_mu0_eff))) * sigma_Ip_eff * sigma_Bp_eff / ((2 * 3.14159265358979323846) ** exp_Bp_eff)
    transforms["FFPRIME"] = lB_eff * sigma_Ip_eff * sigma_Bp_eff / ((2 * 3.14159265358979323846) ** exp_Bp_eff)
    transforms["B"] = lB_eff * sigma_B0_eff
    transforms["F"] = sigma_B0_eff * ld_eff * lB_eff
    transforms["I"] = sigma_Ip_eff * ld_eff * lB_eff / (mu0 ** exp_mu0_eff)
    transforms["J"] = sigma_Ip_eff * lB_eff / ((mu0 ** exp_mu0_eff) * ld_eff)
    transforms["Q"] = sigma_Ip_eff * sigma_B0_eff * sigma_rhotp_eff

    return transforms

def fromCocosNtoCocosM(eqd, cocos_m, phiclockwise=None, weberperrad=None):
    """Transform equilibrium from cocos_n (determined from eqd, see below) to cocos_m.

    Parameters
    ----------
    eqd : dict
        Dictionary from reading the EQDSK file.
    cocosm : int
        Target COCOS.

    Returns
    -------
    eqdout : dict
        Equilibrium data converted to cocos_m.
    """

    cocos_n = assign(eqd["qpsi"][0], eqd["cpasma"], eqd["bcentr"], eqd["simagx"], eqd["sibdry"], phiclockwise, weberperrad)

    transform_dict = transform_cocos(cocos(cocos_n), cocos(cocos_m))

    # Define output
    eqdout = copy.deepcopy(eqd)
    # eqdout["nx"] = eqd["nx"] # For clarity
    # eqdout["ny"] = eqd["ny"] # -||-
    eqdout["rdim"] = eqd["rdim"]*transform_dict["R"]
    eqdout["zdim"] = eqd["zdim"]*transform_dict["Z"]
    eqdout["rcentr"] = eqd["rcentr"]*transform_dict["R"]
    eqdout["rleft"] = eqd["rleft"]*transform_dict["R"]
    eqdout["zmid"] = eqd["zmid"]*transform_dict["Z"]
    eqdout["rmagx"] = eqd["rmagx"]*transform_dict["R"]
    eqdout["zmagx"] = eqd["zmagx"]*transform_dict["Z"]
    eqdout["simagx"] = eqd["simagx"]*transform_dict["PSI"]
    eqdout["sibdry"] = eqd["sibdry"]*transform_dict["PSI"]
    eqdout["bcentr"] = eqd["bcentr"]*transform_dict["B"]
    eqdout["cpasma"] = eqd["cpasma"]*transform_dict["I"]
    eqdout["fpol"] = eqd["fpol"]*transform_dict["F"]
    eqdout["pres"] = eqd["pres"]*transform_dict["PRES"]
    eqdout["ffprime"] = eqd["ffprime"]*transform_dict["FFPRIME"]
    eqdout["pprime"] = eqd["pprime"]*transform_dict["PPRIME"]
    eqdout["psi"] = eqd["psi"]*transform_dict["PSI"]
    eqdout["qpsi"] = eqd["qpsi"]*transform_dict["Q"]
    #eqdout["psiaxis"] = psiaxis # Wrong key? psiaxis is called simagx according to freeqdsk docs
    #eqdout["psiedge"] = psiedge # Wrong key? psiedge is called sibdry according to freeqdsk docs
    #eqdout["nbdry"] = eqd["nbdry"]
    #eqdout["nlim"] = eqd["nlim"]
    eqdout["rbdry"] = eqd["rbdry"]*transform_dict["R"]
    eqdout["zbdry"] = eqd["zbdry"]*transform_dict["Z"]
    eqdout["rlim"] = eqd["rlim"]*transform_dict["R"]
    eqdout["zlim"] = eqd["zlim"]*transform_dict["Z"]

    return eqdout