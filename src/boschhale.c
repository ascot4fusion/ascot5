/**
 * @file boschhale.c
 * @brief Formulas for fusion cross-sections and thermal reactivities.
 *
 * The model is adapted from here:
 * https://www.doi.org/10.1088/0029-5515/32/4/I07
 */
#include "ascot5.h"
#include <math.h>
#include "boschhale.h"

/**
 * @brief Estimate cross-section for a given fusion reaction.
 *
 * Reaction indices:
 *
 * 1 = T(d,n)^4He
 * 2 = ^3He(d,p)^4He
 * 3 = D(d,p)T
 * 4 = D(d,n)^3He
 *
 * @param reaction reaction for which the cross-section is estimated.
 * @param E ion energy [keV].
 * @return cross-section [m^2].
 */
real boschhale_sigma(int reaction, real E) {

    real BG, A[5], B[4];
    real E_min, E_max;

    switch(reaction) {

    case 1:
        if(E <= 530) {
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

    case 2:
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

    case 3:
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

    case 4:
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

    real sigma = S / (E * exp(BG / sqrt(E))) * 1e-31;

    return sigma;

}

/**
 * @brief Estimate reactivity for a given fusion reaction.
 *
 * Reaction indices:
 *
 * 1 = T(d,n)^4He
 * 2 = ^3He(d,p)^4He
 * 3 = D(d,p)T
 * 4 = D(d,n)^3He
 *
 * @param reaction reaction for which the reactivity is estimated.
 * @param Ti ion temperature [keV].
 * @return reactivity.
 */
real boschhale_sigmav(int reaction, real Ti) {

    real BG, MRC2, C1, C2, C3, C4, C5, C6, C7;

    switch(reaction) {

    case 1:
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

    case 2:
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

    case 3:
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

    case 4:
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
