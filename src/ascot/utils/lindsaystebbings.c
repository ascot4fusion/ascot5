/**
 * Implements Lindsay-Stebbings rate coefficients for hydrogen charge-exchange
 * in "ratecoefficients.h".
 */
#include "consts.h"
#include "defines.h"
#include "mathlib.h"
#include "ratecoefficients.h"
#include <math.h>

real lindsaystebbings_sigma(real ekin)
{
    real A = 4.15 - 0.531 * log(ekin);
    return A * A * pow(-expm1(-67.3 / ekin), 4.5) * 1e-20;
}

err_t lindsaystebbings_sigmav(real ekin, real Ti, real *sigmav)
{
    err_t err = 0;

    real ub = sqrt(2.0 * ekin / CONST_M_P);
    real uT = sqrt(2.0 * Ti / CONST_M_P);
    real a = ub / uT;
    real xmax = a + 2;

    real integral = 0;
    for (size_t i = 0; i < GAUSS_LEGENDRE_NQUAD; i++)
    {
        real x = 0.5 * xmax * (mathlib_gauss_legendre_node[i] + 1.0);
        real Erel_keV = Ti * x * x / (1.0e3 * CONST_E);
        real d1 = x - a;
        real d2 = x + a;
        integral += mathlib_gauss_legendre_weight[i] * x * x *
                    lindsaystebbings_sigma(Erel_keV) *
                    (exp(-d1 * d1) - exp(-d2 * d2));
    }

    *sigmav = 0.5 * xmax * uT * integral / (CONST_SQRTPI * a);

    return err;
}
