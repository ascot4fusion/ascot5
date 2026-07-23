/**
 * @file boschhale.c
 * @brief Formulas for fusion cross-sections and thermal reactivities.
 *
 * The model is adapted from here:
 * https://www.doi.org/10.1088/0029-5515/32/4/I07
 */
#include "ascot5.h"
#include <math.h>
#include "consts.h"
#include "boschhale.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Get masses and charges of particles participating in the reaction and
 * the released energy
 *
 * @param reaction reaction enum
 * @param m1 mass of the first reactant [kg]
 * @param q1 charge of the first reactant [C]
 * @param m2 mass of the second reactant [kg]
 * @param q2 charge of the second reactant [C]
 * @param mprod1 mass of the first product [kg]
 * @param qprod1 charge of the first product [C]
 * @param mprod2 mass of the second product [kg]
 * @param qprod2 charge of the second product [C]
 * @param Q energy released [J]
 */
void boschhale_reaction(
    Reaction reaction, real* m1, real* q1, real* m2, real* q2,
    real* mprod1, real* qprod1, real* mprod2, real* qprod2, real* Q) {
    switch(reaction) {
        case DT_He4n:
            *m1     = 3.344e-27; // D
            *q1     = CONST_E;
            *m2     = 5.008e-27; // T
            *q2     = CONST_E;
            *mprod1 = 6.645e-27; // He4
            *qprod1 = 2*CONST_E;
            *mprod2 = 1.675e-27; // n
            *qprod2 = 0.0;
            *Q      = 17.6e6*CONST_E;
            break;
        case DHe3_He4p:
            *m1     = 3.344e-27; // D
            *q1     = CONST_E;
            *m2     = 5.008e-27; // He3
            *q2     = 2*CONST_E;
            *mprod1 = 6.645e-27; // He4
            *qprod1 = 2*CONST_E;
            *mprod2 = 1.673e-27; // p
            *qprod2 = CONST_E;
            *Q      = 18.3e6*CONST_E;
            break;
        case DD_Tp:
            *m1     = 3.344e-27; // D
            *q1     = CONST_E;
            *m2     = 3.344e-27; // D
            *q2     = CONST_E;
            *mprod1 = 5.008e-27; // T
            *qprod1 = CONST_E;
            *mprod2 = 1.673e-27; // p
            *qprod2 = CONST_E;
            *Q      = 4.03e6*CONST_E;
            break;
        case DD_He3n:
            *m1     = 3.344e-27; // D
            *q1     = CONST_E;
            *m2     = 3.344e-27; // D
            *q2     = CONST_E;
            *mprod1 = 5.008e-27; // He3
            *qprod1 = 2*CONST_E;
            *mprod2 = 1.675e-27; // n
            *qprod2 = 0.0;
            *Q      = 3.27e6*CONST_E;
            break;
    }
}

/**
 * @brief Estimate cross-section for a given fusion reaction.
 *
 * See: Bosch and Hale, 1992, Nuclear Fusion. Vol. 32, No.4. Section 4.2
 *
 * @param reaction reaction for which the cross-section is estimated.
 * @param E Sum of ion kinetic energies in the CM frame [J].
 *
 * @return cross-section [m^2].
 */
real boschhale_sigma(Reaction reaction, real E) {

    real BG, A[5], B[4];
    real E_min, E_max;
    E = E / (1.e3 * CONST_E); // Convert to keV

    switch(reaction) {

    case DT_He4n:
        if(E <= 550) {
            BG = 34.3827;
            A[0] = 6.927e4;
            A[1] = 7.454e8;
            A[2] = 2.050e6;
            A[3] = 5.2002e4;
            A[4] = 0.0;
            B[0] = 6.38e1;
            B[1] = -9.95e-1;
            B[2] = 6.981e-5;
            B[3] = 1.728e-4;
        }
        else {
            BG = 34.3827;
            A[0] = -1.4714e6;
            A[1] = 0.0;
            A[2] = 0.0;
            A[3] = 0.0;
            A[4] = 0.0;
            B[0] = -8.4127e-3;
            B[1] = 4.7983e-6;
            B[2] = -1.0748e-9;
            B[3] = 8.5184e-14;
        }
        E_min = 0.5;
        E_max = 4700;
        break;

    case DHe3_He4p:
        if(E <= 900) {
            BG = 68.7508;
            A[0] = 5.7501e6;
            A[1] = 2.5226e3;
            A[2] = 4.5566e1;
            A[3] = 0.0;
            A[4] = 0.0;
            B[0] = -3.1995e-3;
            B[1] = -8.5530e-6;
            B[2] = 5.9014e-8;
            B[3] = 0.0;
        }
        else {
            BG = 68.7508;
            A[0] = -8.3993e5;
            A[1] = 0.0;
            A[2] = 0.0;
            A[3] = 0.0;
            A[4] = 0.0;
            B[0] = -2.6830e-3;
            B[1] = 1.1633e-6;
            B[2] = -2.1332e-10;
            B[3] = 1.4250e-14;
        }
        E_min = 0.3;
        E_max = 4800;
        break;

    case DD_Tp:
        BG = 31.3970;
        A[0] = 5.5576e4;
        A[1] = 2.1054e2;
        A[2] = -3.2638e-2;
        A[3] = 1.4987e-6;
        A[4] = 1.8181e-10;
        B[0] = 0.0;
        B[1] = 0.0;
        B[2] = 0.0;
        B[3] = 0.0;
        E_min = 0.5;
        E_max = 5000;
        break;

    case DD_He3n:
        BG = 31.3970;
        A[0] = 5.3701e4;
        A[1] = 3.3027e2;
        A[2] = -1.2706e-1;
        A[3] = 2.9327e-5;
        A[4] = -2.5151e-9;
        B[0] = 0.0;
        B[1] = 0.0;
        B[2] = 0.0;
        B[3] = 0.0;
        E_min = 0.5;
        E_max = 4900;
        break;

    default:
        return -1;
    }

    if(E <= E_min) {
        return 0;
    }

    /* Cap energy for astrophysical S-factor */
    real E2 = E;
    if(E2 > E_max) {
        E2 = E_max;
    }

    real S = (A[0] + E2*(A[1] + E2*(A[2] + E2*(A[3] + E2*A[4]))))
             / (1 + E2*(B[0] + E2*(B[1] + E2*(B[2]+E2*B[3]))));

    /* Check for underflow */
    if(BG / sqrt(E2) > 700) {
        return 0;
    }

    /* "With E in keV, the sigma is given in millibarns", hence 1e-31 */
    real sigma = S / (E * exp(BG / sqrt(E))) * 1e-31;

    return sigma;
}

/**
 * @brief Estimate reactivity for a given fusion reaction.
 *
 * For two Maxwellian distributions with temperature Ti.
 *
 * See: Bosch and Hale, 1992, Nuclear Fusion. Vol. 32, No.4. Section 5.2
 *
 * @param reaction reaction for which the reactivity is estimated.
 * @param Ti ion temperature [keV].
 *
 * @return reactivity [m^3/s].
 */
real boschhale_sigmav(Reaction reaction, real Ti) {

    real BG, MRC2, C1, C2, C3, C4, C5, C6, C7;

    switch(reaction) {

    case DT_He4n:
        BG = 34.3827;
        MRC2 = 1124656;
        C1 = 1.17302E-9;
        C2 = 1.51361E-2;
        C3 = 7.51886E-2;
        C4 = 4.60643E-3;
        C5 = 1.35000E-2;
        C6 = -1.06750E-4;
        C7 = 1.36600E-5;
        break;

    case DHe3_He4p:
        BG = 68.7508;
        MRC2 = 1124572;
        C1 = 5.51036E-10;
        C2 = 6.41918E-3;
        C3 = -2.02896E-3;
        C4 = -1.91080E-5;
        C5 = 1.35776E-4;
        C6 = 0.0;
        C7 = 0.0;
        break;

    case DD_Tp:
        BG = 31.3970;
        MRC2 = 937814;
        C1 = 5.65718E-12;
        C2 = 3.41267E-3;
        C3 = 1.99167E-3;
        C4 = 0.0;
        C5 = 1.05060E-5;
        C6 = 0.0;
        C7 = 0.0;
        break;

    case DD_He3n:
        BG = 31.3970;
        MRC2 = 937814;
        C1 = 5.43360E-12;
        C2 = 5.85778E-3;
        C3 = 7.68222E-3;
        C4 = 0.0;
        C5 = -2.96400E-6;
        C6 = 0.0;
        C7 = 0.0;
        break;

    default:
        return -1;
    }

    real theta = Ti / (1 - Ti*(C2 + Ti*(C4 + Ti*C6))
                           / (1 + Ti*(C3 + Ti*(C5 + Ti*C7))));

    real xi = pow((BG*BG / (4*theta)), 1.0/3.0);

    real sigmav = C1 * theta * sqrt(xi / (MRC2 * Ti*Ti*Ti))
                  * exp(-3*xi) * 1.e-6;

    return sigmav;
}

/**
 * @brief Estimate reactivity for a given fusion reaction for single ion-bulk
 *
 * The bulk plasma is assumed to have a Maxwellian distribution with thermal
 * speed vt. For evaluating the reactivity integral, Bosch-Hale cross section is
 * used.
 *
 * See: Bosch and Hale, 1992, Nuclear Fusion. Vol. 32, No.4. Section 4.2
 *
 * @param reaction reaction for which the reactivity is estimated.
 * @param vt Bulk plasma thermal speed sqrt(2kT/m) [m/s]
 * @param vf Fast ion speed in the rest frame of the background plasma [m/s]
 * @param N Number of intervals for trapezoidal integration
 *
 * @return reactivity <sigma*u> [m^3/s].
 */
real boschhale_sigmav_beam_bulk(
    Reaction reaction,
    real vt,
    real vf,
    int N)
{
    if (N < 2) {
        return 0.0;
    }

    // Based on reaction, get the constants
    real mu = -1.0;   // Reduced mass of the fusion reactants [kg]
    real Emin = -1.0; // Minimum energy for which Bosch-Hale is valid [J]
    real Emax = -1.0; // Maximum energy for which Bosch-Hale is valid [J]

    switch(reaction) {
    case DT_He4n:
        mu = 2.0048659575986528e-27;
        Emin = 0.5e3 * CONST_E;
        Emax = 4700e3 * CONST_E;
        break;

    case DHe3_He4p:
        mu = 2.004714615900972e-27;
        Emin = 0.3e3 * CONST_E;
        Emax = 4800e3 * CONST_E;
        break;

    case DD_Tp:
        mu = 1.671791818550634e-27;
        Emin = 0.5e3 * CONST_E;
        Emax = 5000e3 * CONST_E;
        break;

    case DD_He3n:
        mu = 1.671791818550634e-27;
        Emin = 0.5e3 * CONST_E;
        Emax = 4900e3 * CONST_E;
        break;

    default:
        return -1;
    }

    const real dE = (Emax - Emin) / (N - 1);

    const real prefactor =
        sqrt(2.0) / (sqrt(M_PI) * vt * vf * mu * sqrt(mu));

    real integral = 0.0;

    const real inv_vt2 = 1.0 / (vt * vt);

    for (int i = 0; i < N; ++i) {

        const real E = Emin + i * dE;

        const real sigma = boschhale_sigma(reaction, E);

        /*
        const real a = (vf / vt) * (vf / vt);
        const real b = 2.0 * E / (mu * vt * vt);
        const real c = 2.0 * vf / (vt * vt)
                       * sqrt(2.0 * E / mu);

        const real f =
            sqrt(E) * sigma *
            (exp(c - a - b) - exp(-c - a - b))
            * prefactor;
        */

        const real u = sqrt(2.0 * E / mu);
        const real f =
            sqrt(E) * sigma *
            (exp(-inv_vt2 * (vf - u)*(vf - u)) - exp(-inv_vt2 * (vf + u)*(vf + u)));

        /* trapezoidal weights */
        if (i == 0 || i == N - 1) {
            integral += 0.5 * f;
        } else {
            integral += f;
        }
    }

    return prefactor * integral * dE;
}
