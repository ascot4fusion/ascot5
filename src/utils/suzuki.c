/**
 * Implements suzuki.h.
 */
#include "suzuki.h"
#include "defines.h"
#include "consts.h"
#include "errors.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

const real A_highE[3][10] = {
    {12.7, 1.25, 0.452, 0.0105, 0.547, -0.102, 0.360, -0.0298, -0.0959,
     4.21e-3},
    {14.1, 1.11, 0.408, 0.0105, 0.547, -0.0403, 0.345, -0.0288, -0.0971,
     4.74e-3},
    {12.7, 1.26, 0.449, 0.0105, 0.547, -0.00577, 0.336, -0.0282, -0.0974,
     4.87e-3}};

const real A_lowE[3][10] = {
    {-52.9, -1.36, 0.0719, 0.0137, 0.454, 0.403, -0.220, 0.0666, -0.0677,
     -1.48e-3},
    {-67.9, -1.22, 0.0814, 0.0139, 0.454, 0.465, -0.273, 0.0751, -0.0630,
     -5.08e-4},
    {-74.2, -1.18, 0.0843, 0.0139, 0.453, 0.491, -0.294, 0.0788, -0.0612,
     -1.85e-4}};

const integer Z_imp[NIMPURITIES] = {2, 6, 6, 4, 8, 7, 3, 5, 26};

const real Zeffmin_imp[NIMPURITIES] = {
    1.0, 1.0, 5.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

const real Zeffmax_imp[NIMPURITIES] = {
    2.1, 5.0, 6.0, 4.0, 5.0, 5.0, 3.0, 5.0, 5.0};

const real B_highE[NIMPURITIES][12] = {
    {0.231, 0.343, -0.185, -0.162e-1, 0.105, -0.703e-1, 0.531e-1, 0.342e-2,
     -0.838e-2, 0.415e-2, -0.335e-2, -0.221e-3},
    {-0.101e1, -0.865e-2, -0.124, -0.145e-1, 0.391, 0.161e-1, 0.298e-1,
     0.332e-2, -0.248e-1, -0.104e-2, -0.152e-2, -0.189e-3},
    {-0.100e1, -0.255e-1, -0.125, -0.142e-1, 0.388, 0.206e-1, 0.297e-1,
     0.326e-2, -0.246e-1, -0.131e-2, -0.148e-2, -0.180e-3},
    {-0.613, 0.552e-1, -0.167, -0.159e-1, 0.304, 0.154e-2, 0.436e-1, 0.378e-2,
     -0.201e-1, -0.216e-3, -0.251e-2, -0.227e-3},
    {-0.102e1, -0.148e-1, -0.674e-1, -0.917e-2, 0.359, 0.143e-1, 0.139e-1,
     0.184e-2, -0.209e-1, -0.732e-3, -0.502e-3, -0.949e-4},
    {-0.102e1, -0.139e-1, -0.979e-1, -0.117e-1, 0.375, 0.156e-1, 0.224e-1,
     0.254e-2, -0.226e-1, -0.889e-3, -0.104e-2, -0.139e-3},
    {-0.441, 0.129, -0.170, -0.162e-1, 0.277, -0.156e-1, 0.466e-1, 0.379e-2,
     -0.193e-1, 0.753e-3, -0.286e-2, -0.239e-3},
    {-0.732, 0.183e-1, -0.155, -0.172e-1, 0.321, 0.946e-2, 0.397e-1, 0.420e-2,
     -0.204e-1, -0.619e-3, -0.224e-2, -0.254e-3},
    {-0.820, -0.636e-2, 0.542e-1, 0.395e-2, 0.202, 0.806e-3, -0.200e-2,
     -0.178e-2, -0.610e-2, 0.651e-3, 0.175e-2, 0.146e-3}};

const real B_lowE[NIMPURITIES][12] = {
    {-0.792, 0.420e-1, 0.530e-1, -0.139e-1, 0.301, -0.264e-1, -0.299e-1,
     0.607e-2, 0.272e-3, 0.611e-2, 0.347e-2, -0.919e-3},
    {0.161, 0.598e-1, -0.336e-2, -0.426e-2, -0.157, -0.396e-1, 0.460e-2,
     0.219e-2, 0.391e-1, 0.711e-2, -0.144e-2, -0.385e-3},
    {0.158, 0.554e-1, -0.431e-2, -0.335e-2, -0.155, -0.374e-1, 0.537e-2,
     0.174e-2, 0.388e-1, 0.683e-2, -0.160e-2, -0.322e-3},
    {0.112, 0.495e-1, 0.116e-1, -0.286e-2, -0.149, -0.331e-1, -0.426e-2,
     0.980e-3, 0.447e-1, 0.652e-2, -0.356e-3, -0.203e-3},
    {0.111, 0.541e-1, -0.346e-3, -0.368e-2, -0.108, -0.347e-1, 0.193e-2,
     0.181e-2, 0.280e-1, 0.604e-2, -0.841e-3, -0.317e-3},
    {0.139, 0.606e-1, -0.306e-2, -0.455e-2, -0.133, -0.394e-1, 0.399e-2,
     0.236e-2, 0.335e-1, 0.690e-2, -0.124e-2, -0.405e-3},
    {0.112, 0.495e-1, 0.116e-1, -0.286e-2, -0.149, -0.331e-1, -0.426e-2,
     0.980e-3, 0.447e-1, 0.652e-2, -0.356e-3, -0.203e-3},
    {0.122, 0.527e-1, -0.430e-3, -0.318e-2, -0.151, -0.364e-1, 0.343e-2,
     0.151e-2, 0.420e-1, 0.692e-2, -0.141e-2, -0.290e-3},
    {-0.110e-1, 0.202e-1, 0.946e-3, -0.409e-2, -0.666e-2, -0.117e-1, -0.236e-3,
     0.202e-2, 0.408e-2, 0.185e-2, -0.648e-4, -0.313e-3}};

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
    const real(*A)[10];
    const real(*B)[12];
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
        const real *Aijk = A[anum[ind_H[i]] - 1];
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
