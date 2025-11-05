/**
 * Implements hist.h.
 */
#include "diag_hist.h"
#include "consts.h"
#include "defines.h"
#include "marker.h"
#include "parallel.h"
#include "utils/mathlib.h"
#include <stdlib.h>


void DiagHist_offload(DiagHist *hist)
{
    SUPPRESS_UNUSED_WARNING(hist);
    GPU_MAP_TO_DEVICE(
        hist->axes[0:hist->dimensions], hist->bins[0:hist->ns], )
}

void DiagHist_update_go(
    DiagHist *hist, MarkerGyroOrbit *mrk_f, MarkerGyroOrbit *mrk_i)
{
    (void)mrk_i;
#ifdef GPU
    size_t index;
    real weight;
#else
    size_t index[NSIMD];
    real weight[NSIMD];
#endif
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (int i = 0; i < mrk_f->n_mrk; i++)
    {
        HistAxis *axis;

        real phi = fmod(mrk_f->phi[i], 2 * CONST_PI);
        phi += (phi < 0) * 2 * CONST_PI;

        real theta = fmod(mrk_f->theta[i], 2 * CONST_PI);
        theta += (theta < 0) * 2 * CONST_PI;

        real ppar =
            (mrk_f->p_r[i] * mrk_f->B_r[i] + mrk_f->p_phi[i] * mrk_f->B_phi[i] +
             mrk_f->p_z[i] * mrk_f->B_z[i]) /
            sqrt(
                mrk_f->B_r[i] * mrk_f->B_r[i] + mrk_f->B_phi[i] * mrk_f->B_phi[i] +
                mrk_f->B_z[i] * mrk_f->B_z[i]);

        real pperp = sqrt(
            mrk_f->p_r[i] * mrk_f->p_r[i] + mrk_f->p_phi[i] * mrk_f->p_phi[i] +
            mrk_f->p_z[i] * mrk_f->p_z[i] - ppar * ppar);

        axis = &hist->axes[15];
        // size_t nq = axis->n;
        size_t i15 = 0;

        axis = &hist->axes[14];
        size_t i14 =
            math_bin_index(mrk_f->time[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[13];
        // size_t nptor = axis->n;
        size_t i13 = 0;

        axis = &hist->axes[12];
        // size_t nmu = axis->n;
        size_t i12 = 0;

        axis = &hist->axes[11];
        // size_t nxi = axis->n;
        size_t i11 = 0;

        axis = &hist->axes[10];
        // size_t nekin = axis->n;
        size_t i10 = 0;

        axis = &hist->axes[9];
        size_t i9 = math_bin_index(mrk_f->p_z[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[8];
        size_t i8 =
            math_bin_index(mrk_f->p_phi[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[7];
        size_t i7 = math_bin_index(mrk_f->p_r[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[6];
        size_t i6 = math_bin_index(pperp, axis->n, axis->min, axis->max);

        axis = &hist->axes[5];
        size_t i5 = math_bin_index(ppar, axis->n, axis->min, axis->max);

        axis = &hist->axes[4];
        size_t i4 = math_bin_index(theta, axis->n, axis->min, axis->max);

        axis = &hist->axes[3];
        size_t i3 = math_bin_index(mrk_f->rho[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[2];
        size_t i2 = math_bin_index(mrk_f->z[i], axis->n, axis->min, axis->max);

        axis = &hist->axes[1];
        size_t i1 = math_bin_index(phi, axis->n, axis->min, axis->max);

        axis = &hist->axes[0];
        size_t i0 = math_bin_index(mrk_f->r[i], axis->n, axis->min, axis->max);

#ifdef GPU
        size_t *n = hist->strides;
        index = i0 * n[0] + i1 * n[1] + i2 * n[2] + i3 * n[3] + i4 * n[4] +
                i5 * n[5] + i6 * n[6] + i7 * n[7] + i8 * n[8] + i9 * n[9] +
                i10 * n[10] + i11 * n[11] + i12 * n[12] + i13 * n[13] +
                i14 * n[14] + i15;
        weight = mrk_f->weight[i];
        GPU_ATOMIC
        hist->bins[index] += mrk_f->weight[i];
#else
        size_t *n = hist->strides;
        index[i] = i0 * n[0] + i1 * n[1] + i2 * n[2] + i3 * n[3] + i4 * n[4] +
                   i5 * n[5] + i6 * n[6] + i7 * n[7] + i8 * n[8] + i9 * n[9] +
                   i10 * n[10] + i11 * n[11] + i12 * n[12] + i13 * n[13] +
                   i14 * n[14] + i15;
        weight[i] = mrk_f->weight[i];
#endif
    }
#ifndef GPU
    for (size_t i = 0; i < mrk_f->n_mrk; i++)
    {
        if (mrk_f->running[i] && index[i] < hist->nbin)
        {
            GPU_ATOMIC
            hist->bins[index[i]] += weight[i];
        }
    }
#endif
}

void DiagHist_update_gc(
    DiagHist *hist, MarkerGuidingCenter *mrk_f, MarkerGuidingCenter *mrk_i)
{
    // TODO
    (void)hist;
    (void)mrk_f;
    (void)mrk_i;
}