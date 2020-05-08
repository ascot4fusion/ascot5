/**
 * @file step_gceom_mhd.h
 * @brief Guiding center equations of motion with MHD activity
 */
#ifndef STEP_GCEOM_MHD_H
#define STEP_GCEOM_MHD_H

#include <math.h>
#include "../../ascot5.h"
#include "../../math.h"
#include "../../physlib.h"
#include "../../consts.h"

#pragma omp declare target


/**
 * @brief Calculate guiding center equations of motion for a single particle
 *
 * @param i particle index that is calculated
 * @param ydot output right hand side of the equations of motion in a
 *             6-length array (rdot, phidot, zdot, ppardot, mudot, chidot)
 * @param y input coordinates in a 5-length array (r, phi, z, rhopar, mu)
 * @param mass mass
 * @param charge charge
 * @param B_dB magnetic field and derivatives at the guiding center location
 * @param E electric field at the guiding center location
 * @param mhd_dmhd mhd perturbation information evaluated by mhd.c
 */
#pragma omp declare simd
static void step_gceom_mhd(real* ydot, real* y, real mass, real charge,
                           real* B_dB, real* E, real* mhd_dmhd) {

    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[4];
    B[2] = B_dB[8];

    real normB = sqrt(math_dot(B, B));
    real gamma = physlib_gamma_ppar(mass, y[4], y[3], normB);

    real gradB[3];
    gradB[0] = (B[0]*B_dB[1] + B[1]*B_dB[5] + B[2]*B_dB[9]) / normB;
    gradB[1] = (B[0]*B_dB[2] + B[1]*B_dB[6] + B[2]*B_dB[10])
               / (normB * y[0]);
    gradB[2] = (B[0]*B_dB[3] + B[1]*B_dB[7] + B[2]*B_dB[11]) / normB;

    real gradBcrossB[3];
    math_cross(gradB, B, gradBcrossB);

    real curlB[3];
    curlB[0] = B_dB[10] / y[0] - B_dB[7];
    curlB[1] = B_dB[3] - B_dB[9];
    curlB[2] = (B[1] - B_dB[2]) / y[0] + B_dB[5];

    real gradalpha[3];
    real alpha = mhd_dmhd[0];
    gradalpha[0] = mhd_dmhd[2];
    gradalpha[1] = mhd_dmhd[3];
    gradalpha[2] = mhd_dmhd[4];

    real alphadot = mhd_dmhd[1]
        +        ydot[0] * gradalpha[0]
        + y[0] * ydot[1] * gradalpha[1]
        +        ydot[2] * gradalpha[2];

    real Bstar[3];
    Bstar[0] = B[0] + (y[3] / ( charge * normB ) + alpha)
                      * curlB[0];
    Bstar[1] = B[1] + (y[3] / ( charge * normB ) + alpha)
                      * curlB[1];
    Bstar[2] = B[2] + (y[3] / ( charge * normB ) + alpha)
                      * curlB[2];

    real Estar[3];
    Estar[0] = E[0] - mhd_dmhd[7]
        + normB * y[3] * gradalpha[0] / ( gamma * mass )
        - ( y[4] + y[3]*y[3] / ( mass * normB) ) * gradB[0]
        / ( charge * gamma );
    Estar[1] = E[0] - mhd_dmhd[8]
        + normB * y[3] * gradalpha[1] / ( gamma * mass )
        - ( y[4] + y[3]*y[3] / ( mass * normB) ) * gradB[1]
        / ( charge * gamma );
    Estar[2] = E[0] - mhd_dmhd[9]
        + normB * y[3] * gradalpha[2] / ( gamma * mass )
        - ( y[4] + y[3]*y[3] / ( mass * normB) ) * gradB[2]
        / ( charge * gamma );

    real Bhat[3];
    Bhat[0] = B[0] / normB;
    Bhat[1] = B[1] / normB;
    Bhat[2] = B[2] / normB;

    real BhatDotBstar = math_dot(Bhat, Bstar);

    real EstarcrossBhat[3];
    math_cross(Estar, Bhat, EstarcrossBhat);

    ydot[0] = (y[3] * Bstar[0] / ( gamma * mass )
               + EstarcrossBhat[0]) / BhatDotBstar;
    ydot[1] = (y[3] * Bstar[1] / ( gamma * mass )
               + EstarcrossBhat[1]) / (y[0]*BhatDotBstar);
    ydot[2] = (y[3] * Bstar[2] / ( gamma * mass )
               + EstarcrossBhat[2]) / BhatDotBstar;
    ydot[3] = charge * math_dot(Bstar,Estar) / BhatDotBstar
        - alphadot * charge * normB;
    ydot[4] = 0;
    ydot[5] = charge * normB / (gamma*mass);

}

#pragma omp end declare target

#endif
