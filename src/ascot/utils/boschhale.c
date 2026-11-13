/**
 * Implements Bosch-Hale fusion cross-section and rate coefficients in
 * "ratecoefficients.h".
 */
#include "ratecoefficients.h"
#include "consts.h"
#include "defines.h"
#include <math.h>
#include <stdio.h>

const double Q_value[4] = {17.6e6, 18.3e6, 4.03e6, 3.27e6};

/**
 * Parameters for calculating fusion reactivities.
 *
 * The order of the reactions is:
 * 1. H2 + H3 => He4 + n
 * 2. H2 + He3 => He4 + H1
 * 3. H2 + H2 => H3 + H1
 * 4. H2 + H2 => He3 + n
 *
 * The parameters themselves are in the order of B_G, m_r*c^2, C_1, C_2, C_3,
 * C_4, C_5, C_6, C_7.
 */
const double sigmav_parameters[4][9] = {
    {34.3827, 1124656.0, 1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3,
     1.35000e-2, -1.06750e-4, 1.36600e-5},
    {68.7508, 1124572, 5.51036e-10, 6.41918e-3, -2.02896e-3, -1.91080e-5,
     1.35776e-4, 0., 0.},
    {31.3970, 937814, 5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0,
     0.0},
    {31.3970, 937814, 5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6,
     0.0, 0.0}};

/**
 * Parameters for calculating fusion reactivities.
 *
 * The order of the reactions is:
 * 1. H2 + H3 => He4 + n
 * 2. H2 + He3 => He4 + H1
 * 3. H2 + H2 => H3 + H1
 * 4. H2 + H2 => He3 + n
 *
 * The parameters themselves are in the order of B_G, A_1, A_2, A_3, A_4, A_5,
 * B_1, B_2, B_3, B_4.
 */
const double sigma_parameters[6][10] = {
    {34.3827, 6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0., 6.38e1, -9.95e-1,
     6.981e-5, 1.728e-4},
    {34.3827, -1.4714e6, 0., 0., 0., 0., -8.4127e-3, 4.7983e-6, -1.0748e-9,
     8.5184e-14},
    {68.7508, 5.7501e6, 2.5226e3, 4.5566e1, 0., 0., -3.1995e-3, -8.5530e-6,
     5.9014e-8, 0.},
    {68.7508, -8.3993e5, 0., 0., 0., 0., -2.6830e-3, 1.1633e-6, -2.1332e-10,
     1.4250e-14},
    {31.3970, 5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10, 0., 0., 0.,
     0.},
    {31.3970, 5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9, 0., 0., 0.,
     0.}};

const double ekin_limits[4][2] = {
    {0.5, 4700}, {0.3, 4800}, {0.5, 5000}, {0.5, 4900}};

double boschhale_sigma(Reaction reaction, double ekin)
{
    double J_to_keV = 1. / (1.e3 * CONST_E);
    ekin *= J_to_keV;
    double E_min = ekin_limits[reaction][0], E_max = ekin_limits[reaction][1];

    size_t index = 0;
    switch (reaction)
    {
    case DT_He4n:
        index = ekin <= 550 ? 0 : 1;
        break;
    case DHe3_He4p:
        index = ekin <= 900 ? 2 : 3;
        break;
    case DD_Tp:
        index = 4;
        break;
    case DD_He3n:
        index = 5;
    }
    double BG = sigma_parameters[index][0];
    const double *A = &sigma_parameters[index][1];
    const double *B = &sigma_parameters[index][6];

    /* Cap energy for astrophysical S-factor */
    double E2 = ekin;
    if (E2 > E_max)
        E2 = E_max;

    double S = (A[0] + E2 * (A[1] + E2 * (A[2] + E2 * (A[3] + E2 * A[4])))) /
               (1. + E2 * (B[0] + E2 * (B[1] + E2 * (B[2] + E2 * B[3]))));

    /* "With ekin in keV, the sigma is given in millibarns", hence 1e-31 */
    double sigma = S / (ekin * exp(BG / sqrt(ekin))) * 1e-31;

    /* Check for underflow and minimum energy limit */
    if (BG / sqrt(E2) > 700 || ekin <= E_min)
        sigma = 0;

    return sigma;
}

double boschhale_sigmav(Reaction reaction, double Ti)
{
    Ti /= 1e3 * CONST_E;
    double BG = sigmav_parameters[reaction][0];
    double MRC2 = sigmav_parameters[reaction][1];
    const double *C = &sigmav_parameters[reaction][2];

    double cm3_to_m3 = 1.e-6;
    double theta = Ti / (1 - Ti * (C[1] + Ti * (C[3] + Ti * C[5])) /
                                 (1 + Ti * (C[2] + Ti * (C[4] + Ti * C[6]))));
    double xi = pow((BG * BG / (4 * theta)), 1.0 / 3.0);
    double sigmav = C[0] * theta * sqrt(xi / (MRC2 * Ti * Ti * Ti)) *
                    exp(-3 * xi) * cm3_to_m3;

    return sigmav;
}
