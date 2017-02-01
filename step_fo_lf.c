/**
 * @file step_fo_lf.c
 * @brief Calculate a full orbit step for a struct of particles with leap-frog
 **/
#include <math.h>
#include "math.h"
#include "ascot5.h"
#include "step_fo_lf.h"
#include "B_field.h"
#include "particle.h"

/**
 * @brief Integrate a full orbit step for a struct of particles with leap-frog
 *
 * The integration is performed for a struct of NSIMD particles using the 
 * leap frog algorithm (ASCOT paper (11)).
 *
 * @param p particle struct that will be updated
 * @param t time
 * @param h length of time step
 * @param Bdata pointer to magnetic field data
 */
void step_fo_lf(particle_simd_fo* p, real t, real h, B_field_data* Bdata) {

    int i;
    /* Following loop will be executed simultaneously for all i */
    #pragma omp simd 
    for(i = 0; i < NSIMD; i++) {
        if(!p->running[i]) {
            /* Convert velocity to cartesian coordinates */
            real vprevxyz[3];
            vprevxyz[0] = p->rdot[i] * cos(p->phi[i]) - p->phidot[i] * sin(p->phi[i]);
            vprevxyz[1] = p->rdot[i] * sin(p->phi[i]) + p->phidot[i] * cos(p->phi[i]);
            vprevxyz[2] = p->zdot[i];

            real Brpz[3];
            Brpz[0] = p->B_r[i];
            Brpz[1] = p->B_phi[i];
            Brpz[2] = p->B_z[i];

            /* Magnetic field to cartesian coordinates */
            real Bxyz[3];
            Bxyz[0] = Brpz[0] * cos(p->phi[i]) - Brpz[1] * sin(p->phi[i]);
            Bxyz[1] = Brpz[0] * sin(p->phi[i]) + Brpz[1] * cos(p->phi[i]);
            Bxyz[2] = Brpz[2];

            /* Precompute some values that will be used repeatedly */
            real a = h * p->charge[i] / p->mass[i];
            real BdotB = Bxyz[0]*Bxyz[0] + Bxyz[1]*Bxyz[1] + Bxyz[2]*Bxyz[2];
            real vprevdotB = vprevxyz[0]*Bxyz[0] + vprevxyz[1]*Bxyz[1]
                         + vprevxyz[2]*Bxyz[2];
            real vprevcrossB[3];
            vprevcrossB[0] = vprevxyz[1]*Bxyz[2] - vprevxyz[2]*Bxyz[1];
            vprevcrossB[1] = vprevxyz[2]*Bxyz[0] - vprevxyz[0]*Bxyz[2];
            vprevcrossB[2] = vprevxyz[0]*Bxyz[1] - vprevxyz[1]*Bxyz[0];

            /* Update velocities */
            real vxyz[3];
            vxyz[0] = (vprevxyz[0] + a*vprevcrossB[0]
                       + a*a/2*vprevdotB*Bxyz[0] - a*a/4*BdotB*vprevxyz[0])
                      / (1 + a*a/4*BdotB);
            vxyz[1] = (vprevxyz[1] + a*vprevcrossB[1]
                      + a*a/2*vprevdotB*Bxyz[1] - a*a/4*BdotB*vprevxyz[1])
                      / (1 + a*a/4*BdotB);
            vxyz[2] = (vprevxyz[2] + a*vprevcrossB[2]
                       + a*a/2*vprevdotB*Bxyz[2] - a*a/4*BdotB*vprevxyz[2])
                      / (1 + a*a/4*BdotB);

            /* Positions to cartesian coordinates */
            real prevxyz[3];
            prevxyz[0] = p->r[i] * cos(p->phi[i]);
            prevxyz[1] = p->r[i] * sin(p->phi[i]);
            prevxyz[2] = p->z[i];

            /* Update positions */
            real xyz[3];
            xyz[0] = prevxyz[0] + h*vxyz[0];
            xyz[1] = prevxyz[1] + h*vxyz[1];
            xyz[2] = prevxyz[2] + h*vxyz[2];

            /* Back to cylindrical coordinates */
            p->r[i] = sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]);
            p->phi[i] = atan2(xyz[1], xyz[0]);
            p->z[i] = xyz[2];
            p->rdot[i] = vxyz[0] * cos(p->phi[i]) + vxyz[1] * sin(p->phi[i]);
            p->phidot[i] = -vxyz[0] * sin(p->phi[i]) + vxyz[1] * cos(p->phi[i]);
            p->zdot[i] = vxyz[2];

            /* Evaluate magnetic field at new position */
            B_field_eval_B(Brpz, p->r[i], p->phi[i], p->z[i], Bdata);
            p->B_r[i] = Brpz[0];
            p->B_phi[i] = Brpz[1];
            p->B_z[i] = Brpz[2];
        }
    }
}
