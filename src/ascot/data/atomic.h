/**
 * @file atomic.h
 * Atomic reaction data interface
 *
 * This is an interface through which atomic reaction data is initialized
 * and accessed. Reading, for example from local files, is done elsewhere.
 *
 * The name asigma is short for atomicsigma. The word sigma refers to
 * cross-section, a fundamental type of reaction probability data. Note
 * that the data is not necessarily in the form of pure cross-sections.
 * It might be in some derivative form, such as rate coefficients.
 *
 * To add a new atomic reaction data instance, make sure the functions
 * are implemented and called from this interface, and that asigma.h
 * contains enum type for the new instance.
 *
 * The interface checks which instance given data corresponds to from the
 * "type"-field in asigma_offload_data or asigma_data that is given
 * as an argument, and calls the relevant function for that instance.
 */
#ifndef ATOMIC_H
#define ATOMIC_H

#include "defines.h"
#include "utils/interp.h"
#include "atomic_data.h"

/**
 * Initialize local file atomic data and check inputs.
 *
 * @param data pointer to the data struct
 *
 * @return zero if initialization success
 */
int asigma_init(Atomic* atomic, size_t nreac, int *z1, int *a1, int *z2, int *a2,
    int *reactype, int *ne, real *emin, real *emax, int *nn, real *nmin,
    real *nmax, int *nT, real *Tmin, real *Tmax, real *sigma);

/**
 * Free allocated resources.
 *
 * @param data Pointer to the data struct.
 */
void asigma_free(Atomic* atomic);

/**
 * Offload data to the accelerator.
 *
 * @param atomic pointer to the data struct.
 */
void asigma_offload(Atomic* atomic);

/**
 * Toggle extrapolation when evaluating cross sections.
 *
 * In this context the extrapolation means values outside the abscissae
 * are set to zero instead of raising an error.
 *
 * @param extrapolate Flag whether to extrapolate.
 */
void asigma_extrapolate(int extrapolate);

/**
 * Evaluate atomic reaction cross-section.
 *
 * This function evaluates the cross-section (sigma) for the atomic reaction
 * corresponding to the reaction identifiers given as parameters at the
 * given mass-normalized collision energy.
 *
 * This is a SIMD function.
 *
 * @param sigma Pointer to evaluated cross-section.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param z_2 Atomic number of bulk particle.
 * @param a_2 Atomic mass number of bulk particle.
 * @param E_coll_per_amu Energy per amu corresponding to collision speed.
 * @param reac_type Reaction type.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
err_t asigma_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    asigma_reac_type reac_type, Atomic* atomic);

/**
 * Evaluate atomic reaction rate coefficient.
 *
 * This function evaluates the rate coefficient (<sigma*v>) for the atomic
 * reaction corresponding to the reaction identifiers given as parameters
 * at the given fast particle energy and bulk plasma conditions.
 *
 * This is a SIMD function.
 *
 * @param sigmav Pointer to evaluated rate coefficient.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param m_1 Mass of fast particle.
 * @param z_2 Atomic number of bulk particle.
 * @param a_2 Atomic mass number of bulk particle.
 * @param E Energy of fast particle.
 * @param T_e Electron temperature of bulk plasma.
 * @param T_0 Temperature of bulk neutrals.
 * @param n_i Ion density of bulk plasma.
 * @param reac_type Reaction type.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
err_t asigma_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, asigma_reac_type reac_type,
    Atomic* atomic);

/**
 * Evaluate charge exchange rate coefficient.
 *
 * This is a SIMD function.
 *
 * @param ratecoeff Pointer to evaluated rate coefficient.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param E Energy of fast particle.
 * @param mass Mass of fast particle.
 * @param nspec Number of bulk neutral species.
 * @param znum Atomic numbers of bulk particles.
 * @param anum Atomic mass numbers of bulk particles.
 * @param T_0 Temperature of bulk neutrals.
 * @param n_0 Densities of bulk neutrals.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
err_t asigma_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, size_t nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    Atomic* atomic);

/**
 * Evaluate beam stopping rate coefficient.
 *
 * This function first tries to evaluate BMS with ADAS data. If not present,
 * the Suzuki model is used instead.
 *
 * @param ratecoeff Pointer to evaluated rate coefficient.
 * @param z_1 Atomic number of fast particle.
 * @param a_1 Atomic mass number of fast particle.
 * @param E Energy of fast particle.
 * @param mass Mass of fast particle.
 * @param nion Number of bulk ion species.
 * @param znum Atomic numbers of bulk particles.
 * @param anum Atomic mass numbers of bulk particles.
 * @param T_e Electron temperature of bulk plasma.
 * @param n_i Densities of bulk ions.
 * @param asigma_data Pointer to atomic data struct.
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
DECLARE_TARGET_SIMD_UNIFORM(atomic)
err_t asigma_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, size_t nion,
    const int* znum, const int* anum, real T_e, real* n_i,
    Atomic* atomic);

#endif
