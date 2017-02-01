/**
 * @file step_gc_rk4.c
 * @brief Calculate a guiding center step for a struct of particles with RK4
 **/
#include <math.h>
#include "ascot5.h"
#include "step_gc_rk4.h"
#include "B_field.h"
#include "math.h"
#include "particle.h"

/**
 * @brief Calculate guiding center equations of motion for a single particle
 *
 * @param i particle index that is calculated
 * @param ydot output right hand side of the equations of motion in a
 *             5-length array (rdot, phidot, zdot, vpardot, mudot)
 * @param t time
 * @param yprev input coordinates in a 5-length array (r, phi, z, vpar, mu)
 * @param mass mass 
 * @param charge charge
 * @param B_dB magnetic field and derivatives at the particle location
 *             (Br = B_dB[0][i], dBr/dr = B_dB[1][i], dBr/dphi = B_dB[2][i], 
 *             dBr/dz = B_dB[3][i], Bphi = B_dB[4][i], dBphi/dr = B_dB[5][i],
 *             dBphi/dphi = B_dB[6][i],dBphi/dz = B[7][i], Bz = B[8][i],
 *             dBz/dr = B[9][i], dBz/dphi = B[10][i], dBz/dz = B[11][i])
 */
inline void ydot_gc(real ydot[], real t, real y[], real mass[], real charge[], real B_dB[]) {

    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[4];
    B[2] = B_dB[8];

    real normB = sqrt(math_dot(B, B));

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
    Bstar[0] = B[0] + (mass[0] * y[3] / charge[0])
                      * (curlB[0] / normB - gradBcrossB[0] / (normB*normB));
    Bstar[1] = B[1] + (mass[0] * y[3] / charge[0])
                      * (curlB[1] / normB - gradBcrossB[1] / (normB*normB));
    Bstar[2] = B[2] + (mass[0] * y[3] / charge[0])
                      * (curlB[2] / normB - gradBcrossB[2] / (normB*normB));

    real Estar[3];
    Estar[0] = -y[4] * gradB[0] / charge[0];
    Estar[1] = -y[4] * gradB[1] / charge[0];
    Estar[2] = -y[4] * gradB[2] / charge[0];

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
    ydot[3] = (charge[0]/mass[0]) * math_dot(Bstar,Estar)/BhatDotBstar;
    ydot[4] = 0;
}

/**
 * @brief Integrate a guiding center step for a struct of particles with RK4
 *
 * This function calculates a guiding center step for a struct of NSIMD 
 * particles with RK4 simultaneously using SIMD instructions. All arrays in the 
 * function are of NSIMD length so vectorization can be performed directly 
 * without gather and scatter operations.
 *
 * @param p particle struct that will be updated
 * @param t time
 * @param h length of time step
 * @param Bdata pointer to magnetic field data
 */
void step_gc_rk4(particle_simd_gc* p, real t, real h, B_field_data* Bdata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd 
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real k1[5];
            real k2[5];
            real k3[5];
            real k4[5];
            real tempy[5];
            real yprev[5];
            real y[5];

           /* Mass and charge need to have the same structure as in the previous
               arrays so that the compiler can vectorize the function call, but
               we only use the first row for actual data */
            real mass[5];
            real charge[5];

            real B[3];
            real B_dB[12];

            /* Coordinates are copied from the struct into an array to make 
             * passing parameters easier */
            yprev[0] = p->r[i];
            yprev[1] = p->phi[i];
            yprev[2] = p->z[i];
            yprev[3] = p->vpar[i];
            yprev[4] = p->mu[i];
            mass[0] = p->mass[i];
            charge[0] = p->charge[i];

            B_field_eval_B_dB(B_dB, yprev[0], yprev[1], yprev[2], Bdata);
            ydot_gc(k1, t, yprev, mass, charge, B_dB);
            int j;
            /* particle coordinates for the subsequent ydot evaluations are
             * stored in tempy */
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h/2.0*k1[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            ydot_gc(k2, t+h/2.0, tempy, mass, charge, B_dB);
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h/2.0*k2[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            ydot_gc(k3, t+h/2.0, tempy, mass, charge, B_dB);
            for(j = 0; j < 5; j++) {
                tempy[j] = yprev[j] + h*k3[j];
            }

            B_field_eval_B_dB(B_dB, tempy[0], tempy[1], tempy[2], Bdata);
            ydot_gc(k4, t+h, tempy, mass, charge, B_dB);
            for(j = 0; j < 5; j++) {
                y[j] = yprev[j]
                    + h/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
            } 

            p->r[i] = y[0];
            p->phi[i] = y[1];
            p->z[i] = y[2];
            p->vpar[i] = y[3];

            p->time[i] = p->time[i] + h;

            /* Update other particle parameters to be consistent */
            B_field_eval_B(B, y[0], y[1], y[2], Bdata);
            p->B_r[i] = B[0];
            p->B_phi[i] = B[1];
            p->B_z[i] = B[2]; 

            p->prev_r[i] = yprev[0];
            p->prev_phi[i] = yprev[1];
            p->prev_z[i] = yprev[2];

        }
    }
}
