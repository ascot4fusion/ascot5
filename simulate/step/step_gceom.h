/**
 * @file step_gceom.h
 * @brief Guiding center equations of motion
 */
#ifndef STEP_GCEOM_H
#define STEP_GCEOM_H

#include <math.h>
#include "../../ascot5.h"
#include "../../math.h"
#include "../../physlib.h"

#pragma omp declare target

/**
 * @brief Calculate guiding center equations of motion for a single particle
 *
 * @param ydot output right hand side of the equations of motion in a
 *             5-length array (rdot, phidot, zdot, vpardot, mudot)
 * @param y input coordinates in a 5-length array (r, phi, z, vpar, mu)
 * @param mass mass [kg]
 * @param charge charge [C]
 * @param B_dB magnetic field and derivatives at the guiding center location
 * @param E electric field at the guiding center location
 */
#pragma omp declare simd
static void step_gceom(real* ydot, real* y, real mass, real charge,
                       real* B_dB, real* E) {

    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[4];
    B[2] = B_dB[8];

    real normB = sqrt(math_dot(B, B));
    real gamma = physlib_gamma_vpar(mass, y[4], y[3], normB);

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

    real Bstar[3];
    Bstar[0] = B[0] + (mass * y[3] * gamma / charge)
                      * (curlB[0] / normB - gradBcrossB[0] / (normB*normB));
    Bstar[1] = B[1] + (mass * y[3] * gamma / charge)
                      * (curlB[1] / normB - gradBcrossB[1] / (normB*normB));
    Bstar[2] = B[2] + (mass * y[3] * gamma / charge)
                      * (curlB[2] / normB - gradBcrossB[2] / (normB*normB));

    real Estar[3];
    Estar[0] = E[0] - y[4] * gradB[0] / (charge * gamma);
    Estar[1] = E[1] - y[4] * gradB[1] / (charge * gamma);
    Estar[2] = E[2] - y[4] * gradB[2] / (charge * gamma);

    real Bhat[3];
    Bhat[0] = B[0] / normB;
    Bhat[1] = B[1] / normB;
    Bhat[2] = B[2] / normB;

    real BhatDotBstar = math_dot(Bhat, Bstar);

    real EstarcrossBhat[3];
    math_cross(Estar, Bhat, EstarcrossBhat);

    ydot[0] = (y[3]*Bstar[0]+EstarcrossBhat[0])/BhatDotBstar;
    ydot[1] = (y[3]*Bstar[1]+EstarcrossBhat[1])/(y[0]*BhatDotBstar);
    ydot[2] = (y[3]*Bstar[2]+EstarcrossBhat[2])/BhatDotBstar;
    ydot[3] = (charge/(mass*gamma)) * math_dot(Bstar,Estar)/BhatDotBstar;
    ydot[4] = 0;
    ydot[5] = charge * normB/(gamma*mass);

}

#pragma omp end declare target

#endif
