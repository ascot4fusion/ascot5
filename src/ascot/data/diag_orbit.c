/**
 * Implements diag_orbit.h.
 */
#include "utils/mathlib.h"
#include "diag_orbit.h"
#include "consts.h"
#include "defines.h"
#include "diag.h"
#include "bfield.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

void DiagOrbit_offload(DiagOrbit *orbit)
{
    (void) orbit;
}

void DiagOrbit_update_go(
    DiagOrbit *orbit, Bfield *bfield, MarkerGyroOrbit *mrk_f,
    MarkerGyroOrbit *mrk_i)
{
    size_t nloop = 1;
    if(orbit->interval < 0)
        nloop = orbit->ntoroidal + orbit->npoloidal + orbit->nradial;

    for(size_t j = 0; j < nloop; j++) {
        real beta;
        if(orbit->interval >= 0)
            beta = 0;
        else if(j < orbit->ntoroidal)
            beta = orbit->toroidal[j];
        else
            beta = orbit->poloidal[j - orbit->ntoroidal];

        #pragma omp simd
        for (size_t i = 0; i < mrk_f->n_mrk; i++)
        {
            if(mrk_f->id[i] > 0) {
                size_t imrk = mrk_f->index[i];
                size_t ipoint = orbit->idx[imrk];
                size_t idx = imrk * orbit->npoint + ipoint;
                double k = -1.0;
                if(orbit->interval < 0 && mrk_f->mileage[i] != mrk_i->mileage[i])
                {
                    int istoroidal = j < orbit->ntoroidal ? 1 : 0;
                    real alpha2 = istoroidal ? mrk_f->phi[i] : mrk_f->theta[i];
                    real alpha1 = istoroidal ? mrk_i->phi[i] : mrk_i->theta[i];
                    k = math_crossed_plane(alpha1, alpha2, beta);
                }
                else if(orbit->interval >= 0 &&
                    (( fabs(mrk_f->mileage[i] - orbit->stamp[imrk]) > orbit->interval) ||
                    ( orbit->id[imrk * orbit->npoint] == 0 )))
                    k = 1;
                real d = 1 - k;

                orbit->charge[idx] = k >= 0 ? rint(mrk_f->charge[i]/CONST_E) : orbit->charge[idx];
                orbit->r[idx] = k >= 0 ? k * mrk_f->r[i] + d * mrk_i->r[i] : orbit->r[idx];
                orbit->z[idx] = k >= 0 ? k * mrk_f->z[i] + d * mrk_i->z[i] : orbit->z[idx];
                orbit->p1[idx] = k * mrk_f->p_r[i] + d * mrk_i->p_r[i];
                orbit->p3[idx] = k * mrk_f->p_z[i] + d * mrk_i->p_z[i];
                orbit->p2[idx] = k * mrk_f->p_phi[i] + d * mrk_i->p_phi[i];
                orbit->phi[idx] = k >= 0 ? k * mrk_f->phi[i] + d * mrk_i->phi[i] : orbit->phi[idx];
                orbit->mileage[idx] = k >= 0 ? k * mrk_f->mileage[i] + d * mrk_i->mileage[i] : orbit->mileage[idx];

                orbit->id[idx] = k >= 0 ? mrk_f->id[i] : orbit->id[idx];
                orbit->poincare[idx] = k >= 0 ? j : orbit->poincare[idx];
                orbit->simmode[idx] = k >= 0 ? DIAG_ORB_GC : orbit->simmode[idx];

                ipoint++;
                ipoint = ipoint == orbit->npoint ? 0 : ipoint;
                orbit->idx[imrk] = k >= 0 ? ipoint : orbit->idx[imrk];
                orbit->stamp[imrk] = k >= 0 ? mrk_f->mileage[i] : orbit->stamp[imrk];
            }
        }
    }
}

void DiagOrbit_update_gc(
    DiagOrbit *orbit, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i)
{
    size_t nloop = 1;
    if(orbit->interval < 0)
        nloop = orbit->ntoroidal + orbit->npoloidal + orbit->nradial;

    for(size_t j = 0; j < nloop; j++) {
        real beta;
        if(orbit->interval >= 0)
            beta = 0;
        else if(j < orbit->ntoroidal)
            beta = orbit->toroidal[j];
        else
            beta = orbit->poloidal[j - orbit->ntoroidal];

        #pragma omp simd
        for (size_t i = 0; i < mrk_f->n_mrk; i++)
        {
            if(mrk_f->id[i] > 0) {
                size_t imrk = mrk_f->index[i];
                size_t ipoint = orbit->idx[imrk];
                size_t idx = imrk * orbit->npoint + ipoint;
                double k = -1.0;
                if(orbit->interval < 0 && mrk_f->mileage[i] != mrk_i->mileage[i]) {
                    int istoroidal = j < orbit->ntoroidal ? 1 : 0;
                    real alpha2 = istoroidal ? mrk_f->phi[i] : mrk_f->theta[i];
                    real alpha1 = istoroidal ? mrk_i->phi[i] : mrk_i->theta[i];
                    k = math_crossed_plane(alpha1, alpha2, beta);
                }
                else if(orbit->interval >= 0 &&
                    (( fabs(mrk_f->mileage[i] - orbit->stamp[imrk]) > orbit->interval) ||
                    ( orbit->id[imrk * orbit->npoint] == 0 ))) {
                    k = 1;
                }
                real d = 1 - k;

                orbit->charge[idx] = k >= 0 ? rint(mrk_f->charge[i]/CONST_E) : orbit->charge[idx];
                orbit->r[idx] = k >= 0 ? k * mrk_f->r[i] + d * mrk_i->r[i] : orbit->r[idx];
                orbit->z[idx] = k >= 0 ? k * mrk_f->z[i] + d * mrk_i->z[i] : orbit->z[idx];
                orbit->p1[idx] = k * mrk_f->ppar[i] + d * mrk_i->ppar[i];
                orbit->p3[idx] = k * mrk_f->mu[i] + d * mrk_i->mu[i];
                orbit->p2[idx] = k * mrk_f->zeta[i] + d * mrk_i->zeta[i];
                orbit->phi[idx] = k >= 0 ? k * mrk_f->phi[i] + d * mrk_i->phi[i] : orbit->phi[idx];
                orbit->mileage[idx] = k >= 0 ? k * mrk_f->mileage[i] + d * mrk_i->mileage[i] : orbit->mileage[idx];

                orbit->id[idx] = k >= 0 ? mrk_f->id[i] : orbit->id[idx];
                orbit->poincare[idx] = k >= 0 ? j : orbit->poincare[idx];
                orbit->simmode[idx] = k >= 0 ? DIAG_ORB_GC : orbit->simmode[idx];

                ipoint++;
                ipoint = ipoint == orbit->npoint ? 0 : ipoint;
                orbit->idx[imrk] = k >= 0 ? ipoint : orbit->idx[imrk];
                orbit->stamp[imrk] = k >= 0 ? mrk_f->mileage[i] : orbit->stamp[imrk];
            }
        }
    }
}

void DiagOrbit_update_fl(
    DiagOrbit *orbit, Bfield *bfield, MarkerFieldLine *mrk_f,
    MarkerFieldLine *mrk_i)
{
    size_t nloop = 1;
    if(orbit->interval < 0)
        nloop = orbit->ntoroidal + orbit->npoloidal + orbit->nradial;

    for(size_t j = 0; j < nloop; j++) {
        real beta;
        if(orbit->interval >= 0)
            beta = 0;
        else if(j < orbit->ntoroidal)
            beta = orbit->toroidal[j];
        else
            beta = orbit->poloidal[j - orbit->ntoroidal];

        #pragma omp simd
        for (size_t i = 0; i < mrk_f->n_mrk; i++)
        {
            if(mrk_f->id[i] > 0) {
                size_t imrk = mrk_f->index[i];
                size_t ipoint = orbit->idx[imrk];
                size_t idx = imrk * orbit->npoint + ipoint;
                double k = -1.0;
                if(orbit->interval < 0 && mrk_f->mileage[i] != mrk_i->mileage[i]) {
                    int istoroidal = j < orbit->ntoroidal ? 1 : 0;
                    real alpha2 = istoroidal ? mrk_f->phi[i] : mrk_f->theta[i];
                    real alpha1 = istoroidal ? mrk_i->phi[i] : mrk_i->theta[i];
                    k = math_crossed_plane(alpha1, alpha2, beta);
                }
                else if(orbit->interval >= 0 &&
                    (( fabs(mrk_f->mileage[i] - orbit->stamp[imrk]) > orbit->interval) ||
                    ( orbit->id[imrk * orbit->npoint] == 0 ))) {
                    k = 1;
                }
                real d = 1 - k;

                orbit->r[idx] = k >= 0 ? k * mrk_f->r[i] + d * mrk_i->r[i] : orbit->r[idx];
                orbit->z[idx] = k >= 0 ? k * mrk_f->z[i] + d * mrk_i->z[i] : orbit->z[idx];
                orbit->p1[idx] = k >= 0 ? k * mrk_f->pitch[i] + d * mrk_i->pitch[i] : orbit->p1[idx];
                orbit->phi[idx] = k >= 0 ? k * mrk_f->phi[i] + d * mrk_i->phi[i] : orbit->phi[idx];
                orbit->mileage[idx] = k >= 0 ? k * mrk_f->mileage[i] + d * mrk_i->mileage[i] : orbit->mileage[idx];

                orbit->id[idx] = k >= 0 ? mrk_f->id[i] : orbit->id[idx];
                orbit->poincare[idx] = k >= 0 ? j : orbit->poincare[idx];
                orbit->simmode[idx] = k >= 0 ? DIAG_ORB_ML : orbit->simmode[idx];

                ipoint++;
                ipoint = ipoint == orbit->npoint ? 0 : ipoint;
                orbit->idx[imrk] = k >= 0 ? ipoint : orbit->idx[imrk];
                orbit->stamp[imrk] = k >= 0 ? mrk_f->mileage[i] : orbit->stamp[imrk];
            }
        }
    }
}
