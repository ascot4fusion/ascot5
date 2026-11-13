"""Templates for importing atomic data (cross sections and rate cofficients)."""

from typing import Optional

import unyt
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import quad

from a5py.data.atomic import CXdata, BMSdata, IIdata, EIdata
from a5py.physlib import Species
from .template import InputTemplate


class ImportAtomicData(InputTemplate):
    """Import atomic data from different formats.

    Initializing this object will not import anything yet. Use the methods to
    import the data which is then stored to this object. Each reaction has to be
    imported individually.

    Attributes
    ----------
    ascot : :class:`Ascot`
        The Ascot object where the input will be created.
    EI : :class:`EIdata`
        Electron impact ionization rate coefficient.
    II : list of :class:`IIdata`
        Ion impact ionization rate coefficients.
    CX : list of :class:`CXdata`
        Charge exchange rate coefficients.
    BMS : list of :class:`BMSdata`
        Beam stopping rate coefficients.
    """

    def __init__(self, ascot):
        """Initialize empty template.

        Parameters
        ----------
        ascot : :class:`Ascot`
            The Ascot object where the input will be created.
        """
        self.EI = None
        self.II = []
        self.CX = []
        self.BMS = []
        super().__init__(ascot, "Atomic", {})

    def import_aladdin_ei(
        self,
        filename: str,
        beam: tuple[str, int],
        Tgrid: Optional[unyt.unyt_array] = None,
    ) -> None:
        """Import electron impact ionization rate coefficient from IAEA ALADDIN.

        Parameters
        ----------
        filename : str
            Name of the ALADDIN file.
        beam : tuple[str, int]
            Beam species and charge state.
        Tgrid : :class:`unyt.unyt_array`, *optional*
            Temperature grid for electrons to interpolate to.
        """
        try:
            temperature, sigmav_ref, _, _ = np.loadtxt(
                filename, skiprows=0, unpack=True
            )
        except ValueError:
            try:
                temperature, sigmav_ref, _, _ = np.loadtxt(
                    filename, skiprows=1, unpack=True
                )
            except Exception as e:
                raise e from None

        temperature *= unyt.eV
        sigmav_ref *= unyt.cm**3 / unyt.s
        beam_species, beam_charge = beam
        beam_species = Species.from_string(beam_species)

        if Tgrid is None:
            Tgrid = np.logspace(np.log10(1e2), np.log10(5e4), 20) * unyt.eV

        sigmav = np.interp(Tgrid, temperature, sigmav_ref).to("m**3/s")
        self.EI = EIdata(
            beam_species=beam_species,
            beam_charge=beam_charge,
            temperature=Tgrid,
            sigmav=sigmav,
        )

    def import_aladdin_ii(
        self,
        filename: str,
        beam: tuple[str, int],
        target: tuple[str, int],
        frame: str,
        egrid: Optional[unyt.unyt_array] = None,
        Tgrid: Optional[unyt.unyt_array] = None,
    ) -> None:
        """Import ion impact ionization cross section from IAEA ALADDIN.

        Parameters
        ----------
        filename : str
            Name of the ALADDIN file.
        beam : tuple[str, int]
            Beam species and charge state.
        target : tuple[str, int]
            Target species and charge state.
        frame : {"beam", "target"}
            Frame in which the cross section is given.
        egrid : :class:`unyt.unyt_array`, *optional*
            Energy grid for beam ions to interpolate to.
        Tgrid : :class:`unyt.unyt_array`, *optional*
            Temperature grid for target ions to interpolate to.
        """
        try:
            ekin, sigma, _, _ = np.loadtxt(filename, skiprows=0, unpack=True)
        except ValueError:
            try:
                ekin, sigma, _, _ = np.loadtxt(filename, skiprows=1, unpack=True)
            except Exception as e:
                raise e from None

        ekin = ekin * unyt.eV
        sigma = sigma * unyt.cm**2
        beam_species, beam_charge = beam
        target_species, target_charge = target
        beam_species = Species.from_string(beam_species)
        target_species = Species.from_string(target_species)

        if frame == "target":
            ekin = ekin * (beam_species.mass / target_species.mass)

        if egrid is None:
            egrid = np.logspace(np.log10(ekin[0]), np.log10(ekin[-1]), 20) * unyt.eV
        if Tgrid is None:
            Tgrid = np.logspace(np.log10(1e2), np.log10(5e4), 20) * unyt.eV

        sigmav = np.zeros((egrid.size, Tgrid.size)) * unyt.m**3 / unyt.s
        for i, T in enumerate(Tgrid):
            for j, E in enumerate(egrid):
                rate = calculate_maxwellian_average(
                    E.to("J").v,
                    T.to("J").v,
                    beam_species.mass.to("kg").v,
                    target_species.mass.to("kg").v,
                    ekin.to("J").v,
                    sigma.to("m**2").v,
                )
                sigmav[j, i] = rate
        self.II.append(
            IIdata(
                beam_species=beam_species,
                beam_charge=beam_charge,
                target_species=target_species,
                target_charge=target_charge,
                temperature=Tgrid,
                energy=egrid,
                sigmav=sigmav,
            )
        )

    def import_aladdin_cx(
        self,
        filename: str,
        beam: tuple[str, int],
        target_species: str,
        egrid: Optional[unyt.unyt_array] = None,
        Tgrid: Optional[unyt.unyt_array] = None,
    ) -> None:
        """Import charge exchange cross section from IAEA ALADDIN.

        Parameters
        ----------
        filename : str
            Name of the ALADDIN file.
        beam : tuple[str, int]
            Beam species and charge state.
        target_species : str
            Neutral species.
        egrid : :class:`unyt.unyt_array`, *optional*
            Energy grid for beam ions to interpolate to.
        Tgrid : :class:`unyt.unyt_array`, *optional*
            Temperature grid for target neutrals to interpolate to.
        """
        try:
            ekin, sigma, _, _ = np.loadtxt(filename, skiprows=0, unpack=True)
        except ValueError:
            try:
                ekin, sigma, _, _ = np.loadtxt(filename, skiprows=1, unpack=True)
            except Exception as e:
                raise e from None

        ekin = ekin * unyt.eV
        sigma = sigma * unyt.cm**2
        beam_species, beam_charge = beam
        beam_species = Species.from_string(beam_species)
        target_species = Species.from_string(target_species)

        if egrid is None:
            egrid = np.logspace(np.log10(ekin[0]), np.log10(ekin[-1]), 20) * unyt.eV
        if Tgrid is None:
            Tgrid = np.logspace(np.log10(1e2), np.log10(5e4), 20) * unyt.eV

        sigmav = np.zeros((egrid.size, Tgrid.size)) * unyt.m**3 / unyt.s
        for i, T in enumerate(Tgrid):
            for j, E in enumerate(egrid):
                rate = calculate_maxwellian_average(
                    E.to("J").v,
                    T.to("J").v,
                    beam_species.mass.to("kg").v,
                    target_species.mass.to("kg").v,
                    ekin.to("J").v,
                    sigma.to("m**2").v,
                )
                sigmav[j, i] = rate

        self.CX.append(
            CXdata(
                beam_species=beam_species,
                beam_charge=beam_charge,
                target_species=target_species,
                temperature=Tgrid,
                energy=egrid,
                sigmav=sigmav,
            )
        )

    def import_openadas_bms(
        self,
        filename: str,
        beam_species: str,
        target: tuple[str, int],
        egrid: Optional[unyt.unyt_array] = None,
        Tgrid: Optional[unyt.unyt_array] = None,
        ngrid: Optional[unyt.unyt_array] = None,
    ) -> None:
        """Import beam stopping rate coefficients from Open-ADAS format.

        Parameters
        ----------
        filename : str
            Name of the Open-ADAS file.
        beam_species : str
            Beam species.
        target : tuple[str, int]
            Target species and charge state.
        egrid : :class:`unyt.unyt_array`, *optional*
            Energy grid for beam neutrals to interpolate to.
        Tgrid : :class:`unyt.unyt_array`, *optional*
            Temperature grid for target ions to interpolate to.
        ngrid : :class:`unyt.unyt_array`, *optional*
            Density grid for target ions to interpolate to.
        """

        data = read_open_adas_bms(filename)

        egrid_ref = data["EB"] * unyt.eV
        ngrid_ref = data["DT"] * unyt.cm**-3
        Tgrid_ref = data["TT"] * unyt.eV
        Tref = data["TREF"] * unyt.eV
        temperature_scaling = data["SVT"]
        sigmav_ref = data["SV"] * unyt.cm**3 / unyt.s

        if egrid is None:
            egrid = (
                np.logspace(
                    np.log10(egrid_ref[0]), np.log10(egrid_ref[-1]), egrid_ref.size
                )
                * unyt.eV
            )
            egrid[0], egrid[-1] = egrid_ref[0], egrid_ref[-1]
        if ngrid is None:
            ngrid = (
                np.logspace(
                    np.log10(ngrid_ref[0]), np.log10(ngrid_ref[-1]), ngrid_ref.size
                )
                * unyt.cm**-3
            )
            ngrid[0], ngrid[-1] = ngrid_ref[0], ngrid_ref[-1]
        if Tgrid is None:
            Tgrid = (
                np.logspace(
                    np.log10(Tgrid_ref[0]), np.log10(Tgrid_ref[-1]), Tgrid_ref.size
                )
                * unyt.eV
            )
            Tgrid[0], Tgrid[-1] = Tgrid_ref[0], Tgrid_ref[-1]

        scaling = np.interp(Tgrid, Tgrid_ref, temperature_scaling) / np.interp(
            Tref, Tgrid_ref, temperature_scaling
        )

        interp = RegularGridInterpolator(
            (egrid_ref, ngrid_ref), sigmav_ref, method="linear", bounds_error=True
        )

        X, Y = np.meshgrid(egrid, ngrid, indexing="ij")
        sigmav = (
            interp((X, Y))[:, None, :] * scaling[None, :, None] * unyt.cm**3 / unyt.s
        )

        self.BMS.append(
            BMSdata(
                beam_species=beam_species,
                target_species=target[0],
                target_charge=target[1],
                density=ngrid,
                temperature=Tgrid,
                energy=egrid,
                sigmav=sigmav,
            )
        )

    def create_input(self) -> object:
        """Uses the Master to create the correct input object."""
        return self._ascot.data.create_atomic(self.CX, self.BMS, self.II, self.EI)


def read_open_adas_bms(filename: str) -> dict:
    """Read beam stopping coefficients, energy/density grids, and temperature
    dependent scaling factors from an Open-ADAS data file.

    Parameters
    ----------
    filename : str
        Path to the input data file.

    Returns
    -------
    dict
        Dictionary containing:
            ZT      : target ion charge
            SVREF   : reference stopping coefficient
            SPEC    : target specification
            DATE    : calculation date
            CODE    : source code name
            NEB     : number of beam energies
            NDT     : number of target densities
            TREF    : reference target temperature (eV)
            EB      : beam energy grid (eV/amu), shape (NEB,)
            DT      : target density grid (cm^-3), shape (NDT,)
            SV      : stopping coefficients (cm^3/s),
                      shape (NDT, NEB)
            NTT     : number of target temperatures
            EREF    : reference beam energy (eV/amu)
            DREF    : reference target density (cm^-3)
            TT      : target temperature grid (eV), shape (NTT,)
            SVT     : temperature scaling factors, shape (NTT,)
    """

    def read_values(f, n, dtype=float):
        """Read exactly n whitespace-separated values across any number of
        lines.
        """
        values = []
        while len(values) < n:
            line = f.readline()
            if not line:
                raise ValueError(
                    f"Unexpected end of file: expected {n} values, "
                    f"got {len(values)}"
                )
            values.extend(line.split())

        if len(values) > n:
            raise ValueError(
                f"Expected {n} values, but read {len(values)}. "
                "Check the file layout."
            )

        if dtype is float:
            return np.array(
                [float(v.replace("D", "E").replace("d", "e")) for v in values]
            )
        return np.array([dtype(v) for v in values])

    with open(filename, "r") as f:
        # First line: metadata
        ZT, SVREF, SPEC, DATE, CODE = f.readline().split()

        ZT = int(ZT)
        SVREF = float(SVREF.split("=")[1])

        # Second line: grid sizes and reference temperature
        f.readline()
        NEB, NDT, TREF = f.readline().split()
        NEB = int(NEB)
        NDT = int(NDT)
        TREF = float(TREF.split("=")[1])
        f.readline()

        # Beam energy grid
        EB = read_values(f, NEB)

        # Target density grid
        DT = read_values(f, NDT)

        f.readline()

        # Stopping coefficient matrix.
        # Each density has NEB values, one per beam energy.
        SV = np.empty((NDT, NEB))

        for idt in range(NDT):
            SV[idt, :] = read_values(f, NEB)

        # Temperature-dependent scaling data
        f.readline()
        NTT, EREF, DREF = f.readline().split()
        NTT = int(NTT)
        EREF = float(EREF.split("=")[1])
        DREF = float(DREF.split("=")[1])
        f.readline()
        TT = read_values(f, NTT)
        f.readline()
        SVT = read_values(f, NTT)

    return {
        "ZT": ZT,
        "SVREF": SVREF,
        "SPEC": SPEC,
        "DATE": DATE,
        "CODE": CODE,
        "NEB": NEB,
        "NDT": NDT,
        "TREF": TREF,
        "EB": EB,
        "DT": DT,
        "SV": SV,
        "NTT": NTT,
        "EREF": EREF,
        "DREF": DREF,
        "TT": TT,
        "SVT": SVT,
    }


def calculate_maxwellian_average(
    E: float,
    T: float,
    m_beam: float,
    m_target: float,
    energy: np.ndarray,
    sigma: np.ndarray,
) -> float:
    """
    Parameters
    ----------
    E : float
        beam energy [J]
    T : float
        Temperature of the Maxwellian species [J].
    m_beam : float
        Mass of the beam particle [kg].
    m_target : float
        Mass of the Maxwellian particle [kg].
    energy : np.ndarray
        Collision energy in the target frame [J].
    sigma : np.ndarray
        Cross section with respect to collision energy [m^2].

    Returns
    -------
    sigmav : float
        Rate coefficient [m^3/s].
    """

    mu = m_beam * m_target / (m_beam + m_target)
    v_beam = np.sqrt(2 * E / m_beam)
    v_th = np.sqrt(2 * T / m_target)

    def integrand(u):
        E_beam = (0.5 * mu * u**2) * m_beam / mu
        S = np.interp(E_beam, energy, np.log(sigma), right=np.nan)
        S = 0.0 if np.isnan(S) else np.exp(S)
        A = (
            S
            * u**2
            * (
                np.exp(-(((u - v_beam) / v_th) ** 2))
                - np.exp(-(((u + v_beam) / v_th) ** 2))
            )
            * np.sqrt(1 / np.pi)
            / (v_th * v_beam)
        )
        return A

    return quad(
        integrand,
        0,
        200 * v_th,
        epsabs=0,
        epsrel=1e-3,
    )[0]
