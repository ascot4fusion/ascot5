
/**
 * Implements "lotzformula.h".
 */
#include "ratecoefficients.h"

#include "consts.h"
#include "defines.h"
#include "mathlib.h"
#include <math.h>
#include <stddef.h>

/** Classical Lotz coefficient [m^2 eV^2]. */
static const double LOTZ_A = 4.5e-18;

/* Number of equivalent electrons in a given shell for Helium. */
const int OCCUPANCY_HE[1] = {1};

/* Electron binding energies for Helium [eV]. */
const real EBINDING_HE[1] = {54.418};

/* Number of subshells for: Helium. */
const int nsubshells[1] = {1};

real lotzformula_sigma(
    size_t nsubshell, real ekin, const int *occupancy,
    const real *binding_energy)
{
    real sigma = 0.0;
    for (size_t i = 0; i < nsubshell; ++i)
    {
        real q = occupancy[i];
        real I = binding_energy[i];
        sigma += (ekin >= I) * LOTZ_A * q * log(ekin / I) / (ekin * I);
    }
    return sigma;
}

err_t lotzformula_sigmav(int anum, int znum, real Te, real *sigmav)
{
    err_t err = 0;
    const int *occupancy = NULL;
    const real *binding_energy = NULL;
    size_t nsubshell = 0;
    if (anum == 4 && znum == 2)
    {
        occupancy = OCCUPANCY_HE;
        binding_energy = EBINDING_HE;
        nsubshell = nsubshells[0];
    }

    Te = Te / CONST_E;
    for (size_t i = 0; i < GAUSS_LAGUERRE_NQUAD; i++)
    {
        real ekin = mathlib_gauss_laguerre_node[i] * Te;
        real sigma =
            lotzformula_sigma(nsubshell, ekin, occupancy, binding_energy);
        sigmav[0] += mathlib_gauss_laguerre_weight[i] * sigma *
                     mathlib_gauss_laguerre_node[i];
    }

    sigmav[0] *= sqrt(8.0 * CONST_E * Te / (CONST_PI * CONST_M_E));

    return err;
}
