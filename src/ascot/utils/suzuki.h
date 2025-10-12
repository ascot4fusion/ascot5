/**
 * @file suzuki.h
 * Suzuki beam-stopping coefficients
 *
 * The data and model is from
 * S Suzuki et al 1998 Plasma Phys. Control. Fusion 40 2097
 * doi.org/10.1088/0741-3335/40/12/009
 *
 * The model calculates beam-stopping coefficients for hydrogen beams in
 * hydrogen plasma with impurities present.
 *
 * Available impurities are: He, Li, Be, B, C, N, O, Fe.
 */
#ifndef SUZUKI_H
#define SUZUKI_H

#include "defines.h"

#define NIMPURITIES 9

/**
 * @brief Fitting parameters for Aijk in equation 28.
 *
 * Table 2a from Suzuki's paper for H, D, and T. Note the last column in the
 * paper is for different magnetic field magnitude, whose effect was negligible
 * and so that is not included here.
 *
 * These values are valid in range 100 < E (keV/amu) < 10000.
 */
extern const real A_highE[3][10];

/**
 * @brief Fitting parameters for Aijk in equation 28.
 *
 * Table 3a from Suzuki's paper for H, D, and T. Note the last column in the
 * paper is for different magnetic field magnitude, whose effect was negligible
 * and so that is not included here.
 *
 * These values are valid in range 10 < E (keV/amu) < 100.
 */
extern const real A_lowE[3][10];

/**
 * @brief Charge number corresponding to the tabulated values of Bijk
 */
extern const int Z_imp[NIMPURITIES];

/**
 * @brief Minimum valid Zeff corresponding to the tabulated values of Bijk
 */
extern const real Zeffmin_imp[NIMPURITIES];

/**
 * @brief Maximum valid Zeff corresponding to the tabulated values of Bijk
 */
extern const real Zeffmax_imp[NIMPURITIES];

/**
 * @brief Fitting parameters for Bijk in equation 28
 *
 * Table 2bc (partially) from Suzuki's paper for H, D, and T.
 *
 * These values are valid in range 100 < E (keV/amu) < 10000.
 */
extern const real B_highE[NIMPURITIES][12];

/**
 * @brief Fitting parameters for Bijk in equation 28.
 *
 * Table 3bc (partially) from Suzuki's paper for H, D, and T.
 *
 * These values are valid in range 10 < E (keV/amu) < 100.
 */
extern const real B_lowE[NIMPURITIES][12];

/**
 * @brief Calculate beam-stopping cross-section according to Suzuki model
 *
 * @param sigmav evaluated beam stopping cross section [m^2]
 * @param EperAmu test particle energy divided by its atomic mass number [J]
 * @param vnorm test particle velocity [m/s]
 * @param ne electron particle density [m^-3]
 * @param te electron temperature [J]
 * @param nion number of ion species present in the plasma
 * @param ni densities of the ion species with at least single hydrogen
 *           species
 * @param anum background ion species' atomic mass number
 * @param znum background ion species' charge number
 *
 * @return zero if evaluation was succesfull
 */
err_t suzuki_sigmav(
    real *sigmav, real EperAmu, real vnorm, real ne, real te, size_t nion,
    real *ni, const int *Anum, const int *Znum);

#endif
