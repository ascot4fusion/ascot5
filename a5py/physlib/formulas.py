import numpy as np

import unyt
from unyt.physical_constants import speed_of_light as c

class Quantity:
    registry = {}
    """Registry of all known quantities."""

    def __init__(self, name, units):
        self.name = name
        self.units = units
        self._formulas = []
        Quantity.registry[name] = self

    def add_formula(self, func, deps):
        """Register a formula with explicit dependencies"""
        self._formulas.append((func, deps))

    @property
    def variables(self):
        """All possible variable combinations that can be used to compute
        this quantity."""
        return [tuple(deps) for _, deps in self._formulas]

    def compute(self, **variables):
        """Try to compute this quantity from the given variables.
        """
        needed, computethese = resolve_quantities(variables.keys(), [self.name])
        needed_qnt = {name: val for name, val in variables.items() if name in needed}
        return compute_quantities(needed_qnt, computethese)[self.name]


def resolve_quantities(available, requested):
    """Find all quantities needed to compute requested quantities and what are
    the minimum set needed to compute them from the available quantities.

    The available quantities are usually obtained from the ASCOT output or from
    interpolating the input data. Since this is both memory and time consuming,
    the motivation of this function is to minimize the number of quantities
    we need. This is also why the actual computation is done in separate
    function.

    Parameters
    ----------
    available : list[str]
        List of quantities available for computations.
    requested : list[str]
        List of requested quantities.

    Returns
    -------
    needed : list[str]
        A (sub)set of available quantities that are in fact needed to calculate
        the required quantities.
    computethese : list[str]
        List of all quantities necessary to be evaluated in order to calculate
        the requested quantities.

        This list also contains the requested quantities.
    """
    registry = Quantity.registry
    available = set(available)
    requested = set(requested)

    computable_quantities = set(available)
    quantity_vs_what_is_needed = {}

    def can_be_computed(qnt):
        """Check if `qnt` can be computed from the already computable
        quantities.

        If yes, the quantity is also included in the list of computable
        quantities.
        """
        for dependencies in qnt.variables:
            if all(dep in computable_quantities for dep in dependencies):
                computable_quantities.add(qnt.name)
                quantity_vs_what_is_needed[qnt.name] = dependencies
                return True
        return False

    def scan_for_computable_quantities():
        """Find if there are any new quantities that can be computed from
        the already computable quantities.

        If yes, the quantity is also included in the list of computable
        quantities.
        """
        found_new_quantities = False
        for qname, q in registry.items():
            if qname in computable_quantities:
                continue
            if can_be_computed(q):
                found_new_quantities = True
        return found_new_quantities

    found_new_quantities = True
    while found_new_quantities:
        found_new_quantities = scan_for_computable_quantities()

    if not requested.issubset(computable_quantities):
        raise RuntimeError(
            f"Cannot resolve {requested - computable_quantities}"
            )

    needed = set()
    computethese = set()

    def backtrack(qname):
        if qname in available:
            needed.add(qname)
        elif qname in quantity_vs_what_is_needed:
            for dep in quantity_vs_what_is_needed[qname]:
                backtrack(dep)
        computethese.add(qname)

    for qname in requested:
        backtrack(qname)

    return needed, computethese


def compute_quantities(available, requested):
    """Compute requested quantities from those that are available.

    Parameters
    ----------
    available : dict[str, float]
        Dictionary of available quantities and their values.
    requested : list[str]
        List of requested quantities.

    Returns
    -------
    computed :dict[str, float]
        Dictionary of requested quantities and their values.
    """
    computed = {}
    known = dict(available)  # working copy of available values

    # try until no new quantities can be computed
    progress = True
    while progress:
        progress = False
        for qname in list(requested):
            if qname in known:
                computed[qname] = known[qname]
                requested.remove(qname)
                progress = True
                continue

            q = Quantity.registry.get(qname)
            if q is None:
                continue  # unknown quantity

            for func, deps in q._formulas:
                if all(dep in known for dep in deps):
                    # compute the value
                    val = func(*(known[d] for d in deps))
                    known[qname] = val
                    computed[qname] = val
                    requested.remove(qname)
                    progress = True
                    break

    return computed


bnorm = Quantity("bnorm", "T")
"""Magnetic field norm."""

bnorm.add_formula(
    lambda br, bphi, bz: np.sqrt( br**2 + bphi**2 + bz**2 ),
    ["br", "bphi", "bz"]
    )


divb = Quantity("divb", "T/m")
"""Divergence of the magnetic field."""

divb.add_formula(
    lambda br, brdr, bphidphi, bzdz, r: br/r + brdr + bphidphi/r + bzdz,
    ["br", "brdr", "bphidphi", "bzdz", "r"]
    )


gradbr = Quantity("gradbr", "T/m")
r""":math:`R` component of the magnetic field gradient."""

gradbr.add_formula(
    lambda br, brdr, bphi, bphidr, bz, bzdr:
    (br*brdr + bphi*bphidr + bz*bzdr) / np.sqrt(br**2 + bphi**2 + bz**2),
    ["br", "brdr", "bphi", "bphidr", "bz", "bzdr"]
    )


gradbphi = Quantity("gradbphi", "T/m")
r""":math:`\phi` component  of the magnetic field gradient."""

gradbphi.add_formula(
    lambda br, brdphi, bphi, bphidphi, bz, bzdphi, r:
    (br*brdphi + bphi*bphidphi + bz*bzdphi) / (np.sqrt(br**2 + bphi**2 + bz**2)*r),
    ["br", "brdphi", "bphi", "bphidphi", "bz", "bzdphi", "r"]
    )


gradbz = Quantity("gradbz", "T/m")
r""":math:`z` component of the magnetic field gradient."""

gradbz.add_formula(
    lambda br, brdz, bphi, bphidz, bz, bzdz:
    (br*brdz + bphi*bphidz + bz*bzdz) / np.sqrt(br**2 + bphi**2 + bz**2),
    ["br", "brdz", "bphi", "bphidz", "bz", "bzdz"]
    )


curlbr = Quantity("curlbr", "T/m")
r""":math:`R` component of the magnetic field curl."""

curlbr.add_formula(
    lambda bzdphi, bphidz, r: bzdphi/r - bphidz,
    ["bzdphi", "bphidz", "r"]
    )


curlbphi = Quantity("curlbphi", "T/m")
r""":math:`\phi` component of the magnetic field curl."""

curlbphi.add_formula(
    lambda brdz, bzdr: brdz - bzdr,
    ["brdz", "bzdr"]
    )


curlbz = Quantity("curlbz", "T/m")
r""":math:`z` component of the magnetic field curl."""

curlbz.add_formula(
    lambda bphi, bphidr, brdphi, r: (r*bphidr + bphi - brdphi) / r,
    ["bphi", "bphidr", "brdphi", "r"]
    )


jr = Quantity("jr", "A/m**2")
r"""Plasma current flux :math:`R` component."""

jr.add_formula(
    lambda curlbr: curlbr / unyt.mu_0,
    ["curlbr"]
    )


jphi = Quantity("jphi", "A/m**2")
r"""Plasma current flux :math:`\phi` component."""

jphi.add_formula(
    lambda curlbphi: curlbphi / unyt.mu_0,
    ["curlbphi"]
    )


jz = Quantity("jz", "A/m**2")
r"""Plasma current flux :math:`z` component."""

jz.add_formula(
    lambda curlbz: curlbz / unyt.mu_0,
    ["curlbz"]
    )


jnorm = Quantity("jnorm", "A/m**2")
"""Plasma current density."""

jnorm.add_formula(
    lambda jr, jphi, jz: np.sqrt(jr**2 + jphi**2 + jz**2),
    ["jr", "jphi", "jz"]
    )


pnorm = Quantity("pnorm", "kg*m/s")
"""Momentum norm."""

pnorm.add_formula(
    lambda pr, pphi, pz: np.sqrt( pr**2 + pphi**2 + pz**2 ),
    ["pr", "pphi", "pz"]
    )
pnorm.add_formula(
    lambda ppar, pperp: np.sqrt( ppar**2 + pperp**2 ),
    ["ppar", "pperp"]
    )
pnorm.add_formula(
    lambda gamma, m: np.sqrt(gamma**2 - 1) * m * c,
    ["gamma", "mass"]
    )
pnorm.add_formula(
    lambda m, mu, ppar, b: np.sqrt( 2 * m * mu * b + ppar**2 ),
    ["mass", "mu", "ppar", "bnorm"]
    )


vnorm = Quantity("vnorm", "m/s")
"""Velocity norm."""

vnorm.add_formula(
    lambda p, m: p / (m * c),
    ["pnorm", "mass"]
    )

vnorm.add_formula(
    lambda gamma: np.sqrt(1 - 1.0 / gamma**2) * c,
    ["gamma"]
    )


energy = Quantity("ekin", "eV")
"""Kinetic energy."""

energy.add_formula(
    lambda p, m: p**2 / (2*m),
    ["pnorm", "mass"]
    )
energy.add_formula(
    lambda v, m: 0.5 * m * v**2,
    ["vnorm", "mass"]
    )
energy.add_formula(
    lambda gamma, m: (gamma - 1.0) * m * c**2,
    ["gamma", "mass"]
    )

pitch = Quantity("pitch", "1")
r"""Pitch :math:`v_\parallel/v`."""

pitch.add_formula(
    lambda pr, pphi, pz, br, bphi, bz: (pr*br + pphi*bphi + pz*bz) / (
        np.sqrt( pr**2 + pphi**2 + pz**2 ) * np.sqrt( br**2 + bphi**2 + bz**2 ) ),
    ["pr", "pphi", "pz", "br", "bphi", "bz"]
    )
pitch.add_formula(
    lambda ppar, b, m, mu: ppar / np.sqrt( 2 * m * mu * b + ppar**2 ),
    ["ppar", "bnorm", "mass", "mu"]
    )

gamma = Quantity("gamma", "1")
"""Lorentz factor."""

gamma.add_formula(
    lambda m, p: np.sqrt( 1 + p / ( m * c )**2 ),
    ["mass", "pnorm"]
    )
gamma.add_formula(
    lambda v: np.sqrt( 1.0 / ( ( 1 - v / c ) * ( 1 + v / c ) ) ),
    ["vnorm"]
    )
gamma.add_formula(
    lambda m, mu, ppar, b: np.sqrt( 1 + 2*mu*b / (m*c**2) + (ppar / (m*c))**2),
    ["mass", "mu", "ppar", "bnorm"]
    )
gamma.add_formula(
    lambda m, mu, vpar, b: np.sqrt( (1 + 2*mu*b / (m*c**2)) / ((1 - vpar/c) * (1 + vpar/c)) ),
    ["mass", "mu", "vpar", "bnorm"]
    )
gamma.add_formula(
    lambda m, energy: 1 + energy / (m*c**2),
    ["mass", "ekin"]
    )

vpar = Quantity("vpar", "m/s")
"""Parallel velocity."""

vpar.add_formula(
    lambda ppar, gamma, m: ppar / ( m*gamma ),
    ["ppar", "gamma", "mass"]
    )

ppar = Quantity("ppar", "kg*m/s")
"""Parallel momentum."""

ppar.add_formula(
    lambda pr, pphi, pz, br, bphi, bz: (pr*br + pphi*bphi + pz*bz) / (
        np.sqrt( br**2 + bphi**2 + bz**2 ) ),
    ["pr", "pphi", "pz", "br", "bphi", "bz"]
    )


mu = Quantity("mu", "T/eV")
"""Magnetic moment."""

mu.add_formula(
    lambda m, xi, p, b: ( 1 - xi**2 ) * p**2 / ( 2 * b * m ),
    ["mass", "pitch", "pnorm", "bnorm"]
    )


ptor = Quantity("ptor", "kg*m^2/s")
"""Toroidal canonical angular momentum."""

ptor.add_formula(
    lambda pphi, r, q, psi: pphi * r + q * psi,
    ["pphi", "r", "charge", "psi"]
    )
ptor.add_formula(
    lambda ppar, r, q, psi, br, bphi, bz: ppar * r * bphi / np.sqrt(br**2 + bphi**2 + bz**2) + q * psi,
    ["ppar", "r", "charge", "psi", "br", "bphi", "bz"]
    )


gyroradius = Quantity("gyroradius", "m")
"""Radius of the gyromotion in uniform magnetic field."""

gyroradius.add_formula(
    lambda q, p, xi, b: ( np.sqrt(1 - xi**2) * p / ( b * np.abs(q) ) ),
    ["charge", "pnorm", "pitch", "bnorm"]
    )


gyrofrequency = Quantity("gyrofrequency", "rad/s")
"""Frequency of the gyromotion in uniform magnetic field."""

gyrofrequency.add_formula(
    lambda m, q, b, gamma: np.abs(q) * b / ( gamma*m ),
    ["mass", "charge", "bnorm", "gamma"]
    )


funlib = {
    "cart2pol":lambda x, y, z=None:(
        np.sqrt(x**2 + y**2), np.arctan2(y, x), z
        ),
    "pol2cart":lambda r, phi, z=None: (
        r * np.cos(phi), r * np.sin(phi), z
        ),
    "cart2pol_vec":lambda x, y, vx, vy, vz=None: (
         vx * np.cos(np.arctan2( y, x )) + vy * np.sin(np.arctan2( y, x )),
        -vx * np.sin(np.arctan2( y, x )) + vy * np.cos(np.arctan2( y, x )),
        vz,
        ),
    "pol2cart_vec": lambda phi, vr, vphi, vz=None: (
        vr * np.cos(phi) - vphi * np.sin(phi),
        vr * np.sin(phi) + vphi * np.cos(phi),
        vz,
        ),
}

def cart2pol_vec(vx,x,vy,y,vz=None,z=None):
    phi = np.arctan2( y, x )
    vr  =  vx * np.cos(phi) + vy * np.sin(phi)
    vphi= -vx * np.sin(phi) + vy * np.cos(phi)
    mps = unyt.m / unyt.s
    if vz is not None:
        return vr*mps,vphi*mps,vz*mps
    return vr*mps,vphi*mps,vz

def mod0to360(angle):
    """Transfer arbitrary angle to between interval [0, 2pi].
    """
    return np.mod(angle + 2*np.pi, 2*np.pi)





def bouncefrequency(m, ekin, minorradius, majorradius, safetyfactor):
    """Estimate bounce frequency from energy and aspect ratio.
    """
    eps = minorradius / majorradius
    gamma = gamma_energy(m, ekin)
    vnorm = vnorm_gamma(gamma)
    return np.sqrt(0.5 * eps) * vnorm / ( safetyfactor * majorradius )

def collfreq_ei(mi, qi, ne, Te, clog):
    """Evaluate electron-ion collision frequency.
    """
    return ( ( np.sqrt(2 / np.pi) / 3 )
             * ( unyt.e * qi / ( 4 * np.pi * unyt.eps_0 ) )**2 * ne * clog
             * 4 * np.pi / np.sqrt( unyt.me * Te**3 ) ).to("1/s")

def collfreq_ie(mi, qi, ne, Te, clog):
    """Evaluate ion-electron collision frequency.
    """
    return ( (unyt.me / mi) * ( np.sqrt(2 / np.pi) / 3 )
             * ( unyt.e * qi / ( 4 * np.pi * unyt.eps_0 ) )**2 * ne * clog
             * 4 * np.pi / np.sqrt( unyt.me * Te**3 ) ).to("1/s")
