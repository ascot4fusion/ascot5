/**
 * @file flow.c
 * @brief Rotating plasma
 */
#include "ascot5.h"
#include "B_field.h"
#include "particle.h"
#include "physlib.h"
#include "plasma.h"

void flow_gc_to_plasma(particle_simd_gc* p, B_field_data* Bdata,
                       plasma_data* pdata) {

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real vr, vphi, vz;

            plasma_eval_rotation(&vr, &vphi, &vz, p->rho[i], p->r[i], p->phi[i],
                                 p->z[i], p->time[i], pdata);

            real vpar = physlib_vnorm_pnorm(p->mass[i], p->ppar[i]);

            vpar -= vphi * p->B_phi[i] / math_normc(p->B_r[i], p->B_phi[i],
                                                    p->B_z[i]);

            p->ppar[i] = physlib_pnorm_vnorm(p->mass[i], vpar);
        }
    }
}

void flow_gc_to_lab(particle_simd_gc* p, B_field_data* Bdata,
                    plasma_data* pdata) {

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real vr, vphi, vz;

            plasma_eval_rotation(&vr, &vphi, &vz, p->rho[i], p->r[i], p->phi[i],
                                 p->z[i], p->time[i], pdata);

            real vpar = physlib_vnorm_pnorm(p->mass[i], p->ppar[i]);

            vpar += vphi * p->B_phi[i] / math_normc(p->B_r[i], p->B_phi[i],
                                                    p->B_z[i]);

            p->ppar[i] = physlib_pnorm_vnorm(p->mass[i], vpar);
        }
    }
}

void flow_fo_to_plasma(particle_simd_gc* p, B_field_data* Bdata,
                       plasma_data* pdata) {

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real vr, vphi, vz;

            plasma_eval_rotation(&vr, &vphi, &vz, p->rho[i], p->r[i], p->phi[i],
                                 p->z[i], p->time[i], pdata);

            real vpar = physlib_vnorm_pnorm(p->mass[i], p->ppar[i]);

            vpar -= vphi * p->B_phi[i] / math_normc(p->B_r[i], p->B_phi[i],
                                                    p->B_z[i]);

            p->ppar[i] = physlib_pnorm_vnorm(p->mass[i], vpar);
        }
    }
}

void flow_fo_to_lab(particle_simd_gc* p, B_field_data* Bdata,
                    plasma_data* pdata) {

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real vr, vphi, vz;

            plasma_eval_rotation(&vr, &vphi, &vz, p->rho[i], p->r[i], p->phi[i],
                                 p->z[i], p->time[i], pdata);

            real vpar = physlib_vnorm_pnorm(p->mass[i], p->ppar[i]);

            vpar += vphi * p->B_phi[i] / math_normc(p->B_r[i], p->B_phi[i],
                                                    p->B_z[i]);

            p->ppar[i] = physlib_pnorm_vnorm(p->mass[i], vpar);
        }
    }
}
