/**
 * @file atomic.h
 * Atomic reaction data interface
 */
#ifndef ATOMIC_H
#define ATOMIC_H

#include "atomic_data.h"
#include "defines.h"
#include "utils/interp.h"

/**
 * Offload data to the accelerator.
 *
 * @param atomic Data to be offloaded.
 */
void Atomic_offload(Atomic *atomic);

GPU_DECLARE_TARGET_SIMD_UNIFORM(z_1, a_1, extrapolate, atomic)
/**
 * Evaluate electron impact ionization rate coefficient
 *
 * If rate coefficient data is not available, Lotz model is used instead.
 *
 * @param ratecoeff Evaluated rate coefficient [m^3/s].
 * @param z_1 Test particle atomic number.
 * @param a_1 Test particle mass number.
 * @param q_1 Test particle charge state.
 * @param Te Electron temperature [J].
 * @param ne Electron density [m^3].
 * @param extrapolate Flag to extrapolate data with zero values.
 * @param atomic Atomic data.
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t Atomic_eval_ei(
    real *ratecoeff, int z_1, int a_1, int q_1, real Te,
    real ne, int extrapolate, Atomic *atomic);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(nspecies, znum, anum, extrapolate, atomic)
/**
 * Evaluate ion impact ionization rate coefficient
 *
 * If rate coefficient data is not available, Barnett model is used instead.
 *
 * @param ratecoeff Evaluated rate coefficient [m^3/s].
 * @param charge Test particle charge state.
 * @param ekin Test particle energy [J].
 * @param nspecies Number of background ion species.
 * @param znum Atomic numbers of background ions.
 * @param anum Mass numbers of background ions.
 * @param Ti Ion temperature [J].
 * @param ni Ion densities [m^-3].
 * @param extrapolate Flag to extrapolate data with zero values.
 * @param atomic Atomic data.
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t Atomic_eval_ii(
    real *ratecoeff, int charge, real ekin, size_t nspecies,
    const int *znum, const int *anum, real Ti, real *ni, int extrapolate,
    Atomic *atomic);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(nspecies, znum, anum, extrapolate, atomic)
/**
 * Evaluate charge exchange rate coefficient.
 *
 * If rate coefficient data is not available, Lindsay-Stebbings model is used
 * instead.
 *
 * @param ratecoeff Evaluated rate coefficient [m^3/s].
 * @param charge Test particle charge state.
 * @param ekin Test particle energy [J].
 * @param nspecies Number of background neutral species.
 * @param znum Atomic numbers of background neutrals.
 * @param anum Mass numbers of background neutrals.
 * @param T0 Neutral temperature [J].
 * @param n0 Neutral densities [m^-3].
 * @param extrapolate Flag to extrapolate data with zero values.
 * @param atomic Atomic data.
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t Atomic_eval_cx(
    real *ratecoeff, int charge, real ekin, size_t nspecies,
    const int *znum, const int *anum, real T0, real *n0, int extrapolate,
    Atomic *atomic);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(mass, nion, znum, anum, extrapolate, atomic)
/**
 * Evaluate beam stopping rate coefficient.
 *
 * If rate coefficient data is not available, Suzuki model is used.
 *
 * @param ratecoeff Evaluated rate coefficient [m^3/s].
 * @param ekin Test particle energy [J].
 * @param mass Test particle mass [kg].
 * @param nion Number of background ion species.
 * @param znum Atomic numbers of background ions.
 * @param anum Mass numbers of background ions.
 * @param Te Electron temperature [J].
 * @param ne Electron density [m^-3].
 * @param ni Background ion densities [m^-3].
 * @param extrapolate Flag to extrapolate data with zero values.
 * @param atomic Atomic data.
 * @return Non-zero err_t value if evaluation failed, zero otherwise.
 */
err_t Atomic_eval_bms(
    real ratecoeff[1], real ekin, real mass, size_t nion,
    const int *znum, const int *anum, real Te, real ne, real *ni,
    int extrapolate, Atomic *atomic);
GPU_DECLARE_TARGET_SIMD_UNIFORM_END

#endif
