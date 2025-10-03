/**
 * Implements suzuki.h.
 */
#include "suzuki.h"
#include "ascot5.h"
#include "consts.h"
#include "error.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

a5err suzuki_sigmav(
    real *sigmav, real EperAmu, real vnorm, real ne, real te, integer nion,
    real *ni, const int *anum, const int *znum)
{
    a5err err = 0;
    /* Convert eperamu to keV and te to eV */
    EperAmu /= (1e3 * CONST_E);
    te /= CONST_E;

    int ind_H[3]; /* Indices of possible hydrogen species H, D, and T */
    int n_H = 0;  /* Total number of hydrogen species (max 3)         */
    int ind_Z[MAX_SPECIES]; /* Indices of impurity species            */
    int n_Z = 0; /* Total number of impurity species                 */

    /* Separate ions into hydrogen species and impurities and calculate
     * their total densities and Zeff */
    real dens_H = 0.0, dens_Z = 0.0;
    real Zeff_sum1 = 0.0, Zeff_sum2 = 0.0;
    for (int i = 0; i < nion; i++)
    {
        if (znum[i] == 1)
        {
            ind_H[n_H] = i;
            dens_H += ni[i];
            n_H++;
        }
        else
        {
            ind_Z[n_Z] = i;
            dens_Z += ni[i];
            n_Z++;
        }

        Zeff_sum1 += ni[i] * znum[i] * znum[i];
        Zeff_sum2 += ni[i] * znum[i];
    }
    real Zeff = Zeff_sum1 / Zeff_sum2;

    if (n_H == 0)
    {
        /* No hydrogen species present */
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_SUZUKI);
    }

    /* Select low- or high-energy coefficient tables */
    real(*A)[10];
    real(*B)[12];
    if (EperAmu >= 9.0 && EperAmu < 100.0)
    {
        A = A_lowE;
        B = B_lowE;
    }
    else if (EperAmu < 10000.0)
    {
        A = A_highE;
        B = B_highE;
    }
    else
    {
        /* Outside energy range, just put some data to avoid segfaults */
        A = A_lowE;
        B = B_lowE;
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_SUZUKI);
    }

    real logE = log(EperAmu);
    real N = ne * 1.0e-19;
    real logN = log(N);
    real U = log(te * 1.0e-3);

    /* Equation 28 for sigma_H */
    real sigma_H = 0.0;
    for (int i = 0; i < n_H; i++)
    {
        real *Aijk = A[anum[ind_H[i]] - 1];
        /* Weight with density in case we have multiple hydrogen species */
        sigma_H +=
            ni[ind_H[i]] * (Aijk[0] * 1.e-16 / EperAmu) *
            (1.0 + Aijk[1] * logE + Aijk[2] * logE * logE) *
            (1.0 + pow(1.0 - exp(-Aijk[3] * N), Aijk[4]) *
                       (Aijk[5] + Aijk[6] * logE + Aijk[7] * logE * logE)) *
            (1.0 + Aijk[8] * U + Aijk[9] * U * U);
    }
    sigma_H /= dens_H;

    /* Equations 26 & 27 for Sz */
    real sigma_Z = 0.0;
    for (int i = 0; i < n_Z; i++)
    {
        int ind_B = -1;
        for (int j = 0; j < 9; j++)
        {
            if (Z_imp[j] == znum[ind_Z[i]] && Zeff > Zeffmin_imp[j] &&
                Zeff < Zeffmax_imp[j])
            {
                ind_B = j;
            }
        }
        if (ind_B < 0 && !err)
        {
            /* Data for the input species not tabulated */
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_SUZUKI);
        }
        sigma_Z +=
            ni[ind_Z[i]] / ne * znum[ind_Z[i]] *
            (B[ind_B][0] + B[ind_B][1] * U + B[ind_B][2] * logN +
             B[ind_B][3] * logN * U + B[ind_B][4] * logE +
             B[ind_B][5] * logE * U + B[ind_B][6] * logE * logN +
             B[ind_B][7] * logE * logN * U + B[ind_B][8] * logE * logE +
             B[ind_B][9] * logE * logE * U + B[ind_B][10] * logE * logE * logN +
             B[ind_B][11] * logE * logE * logN * U);
    }

    /* Equation 24 and convert cm^2 to m^2*/
    *sigmav = sigma_H * (1 + (Zeff - 1) * sigma_Z) * 1e-4;
    /* Multiply with velocity to get correct units */
    *sigmav *= vnorm;

    if (err)
    {
        *sigmav = 0.0;
    }
    return err;
}
