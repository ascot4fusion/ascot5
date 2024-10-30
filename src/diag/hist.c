/**
 * @file hist.c
 * @brief Diagnostic for collecting histograms representing the distribution
 *        function.
 *
 * The histogram can have an arbitrary number of dimensions consisting of any
 * of the predefined coordinate axes.
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "../math.h"
#include "../consts.h"
#include "../offload.h"
#include "../particle.h"
#include "hist.h"

/**
 * @brief Initialize the histogram object.
 */
int hist_init(histogram* data, int dimensions, hist_coordinate* coordinates,
              real* binmin, real* binmax, size_t* nbin) {
    data->axes[0].name = R;
    data->axes[1].name = PHI;
    data->axes[2].name = Z;
    data->axes[3].name = RHO;
    data->axes[4].name = THETA;
    data->axes[5].name = PPAR;
    data->axes[6].name = PPERP;
    data->axes[7].name = PR;
    data->axes[8].name = PPHI;
    data->axes[9].name = PZ;
    data->axes[10].name = EKIN;
    data->axes[11].name = XI;
    data->axes[12].name = MU;
    data->axes[13].name = PTOR;
    data->axes[14].name = TIME;
    data->axes[15].name = CHARGE;

    data->nbin = 1;
    for(int i = HIST_ALLDIM-1; i >= 0; i--) {
        data->axes[i].n = 0;
        data->axes[i].min = 0;
        data->axes[i].max = 1;
        for(int k = 0; k < dimensions; k++) {
            if(coordinates[k] == data->axes[i].name) {
                data->axes[i].min = binmin[k];
                data->axes[i].max = binmax[k];
                data->axes[i].n = nbin[k];
                data->nbin *= nbin[k];
            }
        }
        if(i < HIST_ALLDIM-1) {
            data->strides[i] = data->axes[i+1].n ? data->axes[i+1].n : 1;
            if( i < HIST_ALLDIM-2) {
                data->strides[i] *= data->strides[i+1];
            }
        }
    }
    data->bins = (real*) calloc(data->nbin, sizeof(real));
    return 0;
}

/**
 * @brief Free allocated resources
 */
void hist_free(histogram* data) {
    free(data->bins);
    data->nbin = 0;
}

/**
 * @brief Offload the data to the accelerator
 */
void hist_offload(histogram* data) {
    GPU_MAP_TO_DEVICE(
        data->axes[0:data->dimensions], data->bins[0:data->ns],
    )
}

/**
 * @brief Update the histogram in the particle picture
 */
void hist_update_fo(histogram* hist, particle_simd_fo* p_f,
                    particle_simd_fo* p_i) {
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for(int i = 0; i < p_f->n_mrk; i++) {
        hist_axis* axis;

        real phi = fmod(p_f->phi[i], 2*CONST_PI);
        phi += (phi<0) * 2*CONST_PI;

        real theta = fmod(p_f->theta[i], 2*CONST_PI);
        theta += (theta<0) * 2*CONST_PI;

        real ppar = (  p_f->p_r[i]   * p_f->B_r[i]
                     + p_f->p_phi[i] * p_f->B_phi[i]
                     + p_f->p_z[i]   * p_f->B_z[i])
                      / sqrt(  p_f->B_r[i]  * p_f->B_r[i]
                             + p_f->B_phi[i]* p_f->B_phi[i]
                             + p_f->B_z[i]  * p_f->B_z[i]);

        real pperp = sqrt(
                    p_f->p_r[i]   * p_f->p_r[i]
                  + p_f->p_phi[i] * p_f->p_phi[i]
                  + p_f->p_z[i]   * p_f->p_z[i]
                  - ppar * ppar);

        axis = &hist->axes[15];
        size_t nq = axis->n;
        size_t i15 = 0;

        axis = &hist->axes[14];
        size_t i14 = math_bin_index(
            p_f->time[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[13];
        size_t nptor = axis->n;
        size_t i13 = 0;

        axis = &hist->axes[12];
        size_t nmu = axis->n;
        size_t i12 = 0;

        axis = &hist->axes[11];
        size_t nxi = axis->n;
        size_t i11 = 0;

        axis = &hist->axes[10];
        size_t nekin = axis->n;
        size_t i10 = 0;

        axis = &hist->axes[9];
        size_t i9 = math_bin_index(p_f->p_z[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[8];
        size_t i8 = math_bin_index(
            p_f->p_phi[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[7];
        size_t i7 = math_bin_index(p_f->p_r[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[6];
        size_t i6 = math_bin_index(pperp, axis->n, axis->min, axis->max);

        axis = &hist->axes[5];
        size_t i5 = math_bin_index(ppar, axis->n, axis->min, axis->max);

        axis = &hist->axes[4];
        size_t i4 = math_bin_index(theta, axis->n, axis->min, axis->max);

        axis = &hist->axes[3];
        size_t i3 = math_bin_index(p_f->rho[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[2];
        size_t i2 = math_bin_index(p_f->z[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[1];
        size_t i1 = math_bin_index(phi, axis->n, axis->min, axis->max);

        axis = &hist->axes[0];
        size_t i0 = math_bin_index(p_f->r[i], axis->n, axis->min, axis->max);

        size_t* n = hist->strides;
        size_t index =  i0*n[0]  + i1*n[1]   + i2*n[2]   + i3*n[3]   + i4*n[4]
                     +  i5*n[5]  + i6*n[6]   + i7*n[7]   + i8*n[8]   + i9*n[9]
                     + i10*n[10] + i11*n[11] + i12*n[12] + i13*n[13] + i14*n[14]
                     + i15;
        GPU_ATOMIC
        hist->bins[index] += p_f->weight[i];
    }
}

/**
 * @brief Update the histogram in the GC picture
 */
void hist_update_gc(histogram* hist, particle_simd_gc* p_f,
                    particle_simd_gc* p_i) {
    //TODO
}