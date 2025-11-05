/**
 * Implements diag_orbit.h.
 */
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
        #pragma omp simd
        for (size_t i = 0; i < mrk_f->n_mrk; i++)
        {
            if(mrk_f->id[i] > 0) {
                size_t imrk = mrk_f->index[i];
                size_t ipoint = orbit->idx[imrk];
                size_t idx = imrk * orbit->npoint + ipoint;
                real k = -1;
                if(orbit->interval >= 0 &&
                    ( fabs(mrk_f->mileage[i] - orbit->stamp[imrk]) > orbit->interval) ||
                    ( orbit->id[imrk * orbit->npoint] == 0 )) {
                    k = 1;
                }
                real d = 1 - k;

                orbit->r[idx] = k >= 0 ? k * mrk_f->r[i] + d * mrk_i->r[i] : orbit->r[idx];
                orbit->z[idx] = k * mrk_f->z[i] + d * mrk_i->z[i];
                orbit->p1[idx] = k * mrk_f->p_r[i] + d * mrk_i->p_r[i];
                orbit->p3[idx] = k * mrk_f->p_z[i] + d * mrk_i->p_z[i];
                orbit->p2[idx] = k * mrk_f->p_phi[i] + d * mrk_i->p_phi[i];
                orbit->phi[idx] = k * mrk_f->phi[i] + d * mrk_i->phi[i];
                orbit->mileage[idx] = k * mrk_f->mileage[i] + d * mrk_i->mileage[i];

                orbit->id[idx] = mrk_f->id[i];
                orbit->charge[idx] = mrk_i->charge[i] / CONST_E;
                orbit->poincare[idx] = j;
                orbit->simmode[idx] = DIAG_ORB_FO;

                ipoint++;
                ipoint = ipoint == orbit->npoint ? 0 : ipoint;
                orbit->idx[imrk] = ipoint;
                orbit->stamp[imrk] = mrk_f->mileage[i];
            }
        }
    }
}

void DiagOrbit_update_gc(
    DiagOrbit *orbit, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i)
{

}

void DiagOrbit_update_fl(
    DiagOrbit *orbit, Bfield *bfield, MarkerFieldLine *mrk_f,
    MarkerFieldLine *mrk_i)
{
}

int crossed_poincare_plane(double alpha0, double alpha1, double beta)
{
    const double TWO_PI = 2.0 * CONST_PI;

    /* Normalize to [0, 2π) */
    alpha0 = fmod(fmod(alpha0, TWO_PI) + TWO_PI, TWO_PI);
    alpha1 = fmod(fmod(alpha1, TWO_PI) + TWO_PI, TWO_PI);

    /* Handle 2pi -> 0 transition */
    double dphi = alpha1 - alpha0;
    dphi -= dphi > CONST_PI ? TWO_PI : 0;
    dphi += dphi < -CONST_PI ? TWO_PI : 0;

    double rel_plane = beta - alpha0;
    rel_plane += rel_plane < 0 ? TWO_PI : 0;

    return (dphi > 0 && rel_plane <= dphi) || (dphi < 0 && rel_plane + dphi >= 0);
}
