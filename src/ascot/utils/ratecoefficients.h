/**
 * Analytical models for various cross-sections and rate coefficients.
 *
 * Electron impact ionization
 * --------------------------
 *
 * Electron impact ionization is calculated using Lotz formula: Lotz W., "An
 * empirical formula for the electron impact ionization cross-section",
 * Zeitschrift für Physik, 1967, https://doi.org/10.1007/BF01325928. The rate
 * coefficient is calculated by integrating Maxwellian electron distribution
 * using Gauss–Laguerre quadrature.
 *
 * Only reaction considered is He^+1 + e- -> He^+2 + e- + e-.
 *
 * Ion impact ionization
 * ---------------------
 *
 * Ion impact ionization is calculated using Barnett's fit for the rate
 * coefficient, see Barnett C. and others, "Atomic data for fusion. Volume 1:
 * Collisions of H, H2, He and Li atoms and ions with atoms and molecules",
 * 1990, https://doi.org/10.2172/6570226.
 *
 * Only one reaction (D-52) is supported: He^+1 + H^+1 -> He^+2 + H^+1 + e-.
 *
 * Charge-exchange
 * ---------------
 *
 * Charge-exchange is calculated using Lindsay-Stebbings model that gives cross
 * section for H^+1 + H^0 -> H^0 + H^+1, see Lindsay B. G. and Stebbngs R. F.,
 * "Charge transfer cross sections for energetic neutral atom data analysis",
 * Journal of Geophysical Research: Space Physics, 2005,
 * https://www.doi.org/10.1029/2005JA011298. The rate coefficient is calculated
 * by integrating Maxwellian neutral distribution using Gauss–Legendre
 * quadrature.
 *
 * Fusion cross-section and reactivity
 * -----------------------------------
 *
 * The model is adapted from H.-S. Bosch and G.M. Hale, "Improved formulas for
 * fusion cross-sections and thermal reactivities", Nuclear Fusion, 1993,
 * https://www.doi.org/10.1088/0029-5515/32/4/I07.
 *
 * Beam-stopping rate coefficient
 * ------------------------------
 *
 * The data and model is from S. Suzuki and others, "Attenuation of high-energy
 * neutral hydrogen beams in high-density plasmas", Plasma Physics and
 * Controlled Fusion, 1998, https://www.doi.org/10.1088/0741-3335/40/12/009.
 *
 * The model calculates beam-stopping coefficients for hydrogen beams in
 * hydrogen plasma with impurities present.
 *
 * Available impurities are: He, B, C, N, O, Fe.
 */
#ifndef RATECOEFFICIENTS_H
#define RATECOEFFICIENTS_H

#include "defines.h"
#include "datatypes.h"

/**
 * Amount of energy released in a fusion reaction [eV].
 *
 * The order of the reactions is:
 * 1. H2 + H3 => He4 + n
 * 2. H2 + He3 => He4 + H1
 * 3. H2 + H2 => H3 + H1
 * 4. H2 + H2 => He3 + n
 */
extern const double Q_value[4];

/**
 * Evaluate a fusion reaction cross section.
 *
 * @param reaction Specifies the reaction.
 * @param ekin Ion energy [J].
 * @return Cross section [m^2].
 */
double boschhale_sigma(Reaction reaction, double ekin);

/**
 * Evaluate reactivity for two Maxwellian distributions with the same
 * temperature for a given fusion reaction.
 *
 * @param reaction Specifies the reaction.
 * @param Ti Temperature of the ion species' [J].
 * @return Reactivity [m^3/s].
 */
double boschhale_sigmav(Reaction reaction, double Ti);

/**
 * Calculate the electron impact ionization cross section using Lotz formula.
 *
 * @param nsubshell Number of subshells.
 * @param ekin Electron kinetic energy [eV].
 * @param occupancy Occupancy of subshells.
 * @param binding_energy Binding energies of subshells [eV].
 * @return Electron impact ionization cross section [m^2].
 */
real lotzformula_sigma(
    size_t nsubshell, real ekin, const int *occupancy,
    const real *binding_energy);

/**
 * Integrate the electron impact ionization rate coefficient assuming Maxwellian
 * electron distribution.
 *
 * The integral is calculated using Gauss–Laguerre quadrature.
 *
 * @param anum Atomic mass number of the target ion.
 * @param znum Atomic number of the target ion.
 * @param Te Electron temperature [J].
 * @param sigmav Evaluated rate coefficient [m^3/s].
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t lotzformula_sigmav(int anum, int znum, real Te, real *sigmav);

/**
 * Calculate the cross-section for charge-exchange reaction H^+1 + H^0 -> H^0 +
 * H^+1 using Lindsay-Steppings model.
 *
 * @param ekin Electron kinetic energy [keV].
 */
real lindsaystebbings_sigma(real ekin);

/**
 * Integrate the charge-exchange rate coefficient for H^+1 + H^0 -> H^0 + H^+1
 * assuming Maxwellian neutral  distribution.
 *
 * The integral is calculated using Gauss–Laguerre quadrature.
 *
 * @param ekin Electron kinetic energy [J].
 * @param Ti Electron temperature [J].
 * @param sigmav Evaluated rate coefficient [m^3/s].
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t lindsaystebbings_sigmav(real ekin, real Ti, real *sigmav);

/**
 * Calculate the rate coefficient for ion impact ionization He^+1 + H^+1 ->
 * He^+2 + H^+1 + e- using Barnett's fit.
 *
 * @param ekin He^+1 kinetic energy [J].
 * @param Ti Ion temperature [J].
 * @param sigmav Evaluated rate coefficient [m^3/s].
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t barnett_sigmav(real ekin, real Ti, real *sigmav);

/**
 * Calculate beam-stopping cross-section according to Suzuki model.
 *
 * @param sigmav Evaluated beam stopping cross section [m^2].
 * @param eperamu Test particle energy divided by its atomic mass number [J].
 * @param ne Electron particle density [m^-3].
 * @param te Electron temperature [J].
 * @param nion Number of ion species present in the plasma.
 * @param ni Densities of the ion species with at least single hydrogen
 *        species.
 * @param anum Background ion species' atomic mass number.
 * @param znum Background ion species' charge number.
 *
 * @return Zero if evaluation was successful.
 */
err_t suzuki_sigma(
    real *sigmav, real eperamu, real ne, real te, size_t nion,
    real *ni, const int *anum, const int *znum);

#endif