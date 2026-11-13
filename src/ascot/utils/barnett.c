/**
 * Implements Barnett ion impact ionization rate coefficient in
 * "ratecoefficients.h".
 */
#include "consts.h"
#include "defines.h"
#include "ratecoefficients.h"
#include <math.h>

/**
 * The energy abscissa where coefficients are tabulated in terms of He energy
 * [eV/amu].
 */
static const real energy_grid[7] = {
    20e3, 40e3, 50e3, 80e3, 100e3, 200e3, 500e3,
};

/**
 * The coefficients A0–A7 for He^+1 + H^+1 -> He^+2 + H^+1 + e- in terms of
 * energy.
 */
static const real coefficients[7][7] = {
    {-45.6943, 1.57907, 0.710583, 0.175768, -0.00705451, -0.0294110,
     -0.0139117},
    {-41.8177, 0.512347, 0.311944, 0.126340, 0.0238011, -0.00849678,
     -0.00968659},
    {-40.2324, 0.239690, 0.153320, 0.0661005, 0.0168871, -0.00195984,
     -0.00376827},
    {-39.3054, 0.111744, 0.0654773, 0.0331280, 0.00699774, 0.00254954,
     -0.00222171},
    {-38.6923, 0.0227031, 0.0195819, 0.0149078, 0.00271427, 0.00119878,
     -0.00147599},
    {-37.8257, -0.0418823, 0.0248902, -0.0121876, -0.00555012, -0.00208690,
     -5.276582E-04},
    {-39.5936, -0.0146717, -0.00917239, -0.00524336, -0.00267631, -0.00119112,
     -4.070685E-04}};

err_t barnett_sigmav(real ekin, real Ti, real *sigmav)
{
    err_t err = 0;
    real ekin_eV_per_amu = ekin / (4.003 * CONST_E);
    size_t idx = 0;
    for (size_t i = 1; i < 7; i++)
        idx = ekin_eV_per_amu < energy_grid[i] ? i : idx;

    real k = (ekin_eV_per_amu - energy_grid[idx - 1]) /
             (energy_grid[idx] - energy_grid[idx - 1]);
    const real *Ai = coefficients[idx == 0 ? idx : idx - 1];
    const real *Af = coefficients[idx];
    const real A[7] = {
        k * Ai[0] + (1. - k) * Af[0], k * Ai[1] + (1. - k) * Af[1],
        k * Ai[2] + (1. - k) * Af[2], k * Ai[3] + (1. - k) * Af[3],
        k * Ai[4] + (1. - k) * Af[4], k * Ai[5] + (1. - k) * Af[5],
        k * Ai[6] + (1. - k) * Af[6]};

    real Tmin = 1e3, Tmax = 2e4;
    real T_eV_per_amu = Ti / CONST_E;
    real X = 2.0 * log(T_eV_per_amu / sqrt(Tmin * Tmax)) / log(Tmax / Tmin);
    real X2 = X * X;
    real X3 = X2 * X;
    real X4 = X3 * X;
    real X5 = X4 * X;
    real X6 = X5 * X;

    *sigmav = exp(A[0] / 2 + A[1] * X + A[2] * (2 * X2 - 1.0) +
                  A[3] * (4 * X3 - 3 * X) + A[4] * (8 * X4 - 8 * X2 + 1.0) +
                  A[5] * (16 * X5 - 20 * X3 + 5 * X) +
                  A[6] * (32 * X6 - 48 * X4 + 18 * X2 - 1.0)) *
              1e-6;
    return err;
}
