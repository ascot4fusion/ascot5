/**
 * Implements diag_hist.h.
 */
#include "diag_hist.h"
#include "consts.h"
#include "data/bfield.h"
#include "defines.h"
#include "marker.h"
#include "parallel.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <stdlib.h>

void DiagHist_offload(DiagHist *hist)
{
    SUPPRESS_UNUSED_WARNING(hist);
    GPU_MAP_TO_DEVICE(hist->axes [0:HIST_NDIM], hist->bins [0:hist->ns], )
}

void DiagHist_onload(DiagHist *hist)
{
    SUPPRESS_UNUSED_WARNING(hist);
    GPU_MAP_FROM_DEVICE(hist->axes [0:HIST_NDIM], hist->bins [0:hist->ns], )
}

void DiagHist_update_go(
    DiagHist *hist, Bfield *bfield, MarkerGyroOrbit *mrk_f,
    MarkerGyroOrbit *mrk_i)
{
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk_f->size; i++)
    {
        HistAxis *axis;

        real psi;
        Bfield_eval_psi(
            &psi, mrk_f->r[i], mrk_f->phi[i], mrk_f->z[i], mrk_f->time[i],
            bfield);
        real bnorm = math_normc(mrk_f->br[i], mrk_f->bphi[i], mrk_f->bz[i]);

        real phi = fmod(mrk_f->phi[i], 2 * CONST_PI);
        phi += (phi < 0) * 2 * CONST_PI;

        real theta = fmod(mrk_f->theta[i], 2 * CONST_PI);
        theta += (theta < 0) * 2 * CONST_PI;

        real ppar =
            (mrk_f->p_r[i] * mrk_f->br[i] + mrk_f->p_phi[i] * mrk_f->bphi[i] +
             mrk_f->p_z[i] * mrk_f->bz[i]) /
            sqrt(
                mrk_f->br[i] * mrk_f->br[i] +
                mrk_f->bphi[i] * mrk_f->bphi[i] +
                mrk_f->bz[i] * mrk_f->bz[i]);

        real pperp = sqrt(
            mrk_f->p_r[i] * mrk_f->p_r[i] + mrk_f->p_phi[i] * mrk_f->p_phi[i] +
            mrk_f->p_z[i] * mrk_f->p_z[i] - ppar * ppar);

        real charge = mrk_f->charge[i] * CONST_E;
        real pnorm = sqrt(ppar * ppar + pperp * pperp);
        real gamma = physlib_gamma_pnorm(mrk_f->mass, pnorm);
        real ekin = physlib_Ekin_gamma(mrk_f->mass, gamma);
        real pitch = ppar / pnorm;
        real mu = physlib_gc_mu(mrk_f->mass, pnorm, pitch, bnorm);
        real ptor = phys_ptoroid_fo(
            charge, mrk_f->r[i], mrk_f->p_phi[i], psi);

        int valid = 1;
        axis = &hist->axes[15];
        size_t i15 = ((mrk_f->charge[i] - axis->min) /
                      (axis->max - axis->min)) *
                     axis->n;
        valid *= axis->n ? 1
                         : (mrk_f->charge[i] >= axis->min &&
                            mrk_f->charge[i] <= axis->max);

        axis = &hist->axes[14];
        size_t i14 =
            math_bin_index(mrk_f->time[i], axis->n, axis->min, axis->max);
        valid *= axis->n ? 1
                         : (mrk_f->time[i] > axis->min &&
                            mrk_f->time[i] < axis->max);

        axis = &hist->axes[13];
        size_t i13 = math_bin_index(ptor, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (ptor > axis->min && ptor < axis->max);

        axis = &hist->axes[12];
        size_t i12 = math_bin_index(mu, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (mu > axis->min && mu < axis->max);

        axis = &hist->axes[11];
        size_t i11 = math_bin_index(pitch, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (pitch > axis->min && pitch < axis->max);

        axis = &hist->axes[10];
        size_t i10 = math_bin_index(ekin, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (ekin > axis->min && ekin < axis->max);

        axis = &hist->axes[9];
        size_t i9 =
            math_bin_index(mrk_f->p_z[i], axis->n, axis->min, axis->max);
        valid *= axis->n
                     ? 1
                     : (mrk_f->p_z[i] > axis->min && mrk_f->p_z[i] < axis->max);

        axis = &hist->axes[8];
        size_t i8 =
            math_bin_index(mrk_f->p_phi[i], axis->n, axis->min, axis->max);
        valid *= axis->n ? 1
                         : (mrk_f->p_phi[i] > axis->min &&
                            mrk_f->p_phi[i] < axis->max);

        axis = &hist->axes[7];
        size_t i7 =
            math_bin_index(mrk_f->p_r[i], axis->n, axis->min, axis->max);
        valid *= axis->n
                     ? 1
                     : (mrk_f->p_r[i] > axis->min && mrk_f->p_r[i] < axis->max);

        axis = &hist->axes[6];
        size_t i6 = math_bin_index(pperp, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (pperp > axis->min && pperp < axis->max);

        axis = &hist->axes[5];
        size_t i5 = math_bin_index(ppar, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (ppar > axis->min && ppar < axis->max);

        axis = &hist->axes[4];
        size_t i4 = math_bin_index(theta, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (theta > axis->min && theta < axis->max);

        axis = &hist->axes[3];
        size_t i3 =
            math_bin_index(mrk_f->rho[i], axis->n, axis->min, axis->max);
        valid *= axis->n
                     ? 1
                     : (mrk_f->rho[i] > axis->min && mrk_f->rho[i] < axis->max);

        axis = &hist->axes[2];
        size_t i2 = math_bin_index(mrk_f->z[i], axis->n, axis->min, axis->max);
        valid *=
            axis->n ? 1 : (mrk_f->z[i] > axis->min && mrk_f->z[i] < axis->max);

        axis = &hist->axes[1];
        size_t i1 = math_bin_index(phi, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (phi > axis->min && phi < axis->max);

        axis = &hist->axes[0];
        size_t i0 = math_bin_index(mrk_f->r[i], axis->n, axis->min, axis->max);
        valid *=
            axis->n ? 1 : (mrk_f->r[i] > axis->min && mrk_f->r[i] < axis->max);

        size_t *n = hist->strides;
        size_t index = i0 * n[0] + i1 * n[1] + i2 * n[2] + i3 * n[3] +
                       i4 * n[4] + i5 * n[5] + i6 * n[6] + i7 * n[7] +
                       i8 * n[8] + i9 * n[9] + i10 * n[10] + i11 * n[11] +
                       i12 * n[12] + i13 * n[13] + i14 * n[14] + i15;
        real weight = mrk_f->weight[i] * (mrk_f->time[i] - mrk_i->time[i]);
        index = valid ? index : 0;

        GPU_ATOMIC
        hist->values[index] += valid * weight;
    }
}

void DiagHist_update_gc(
    DiagHist *hist, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i)
{
    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk_f->size; i++)
    {
        HistAxis *axis;

        real psi;
        Bfield_eval_psi(
            &psi, mrk_f->r[i], mrk_f->phi[i], mrk_f->z[i], mrk_f->time[i],
            bfield);
        real bnorm = math_normc(mrk_f->br[i], mrk_f->bphi[i], mrk_f->bz[i]);

        real phi = fmod(mrk_f->phi[i], 2 * CONST_PI);
        phi += (phi < 0) * 2 * CONST_PI;

        real theta = fmod(mrk_f->theta[i], 2 * CONST_PI);
        theta += (theta < 0) * 2 * CONST_PI;

        real charge = mrk_f->charge[i] * CONST_E;
        real pnorm =
            physlib_gc_p(mrk_f->mass, mrk_f->mu[i], mrk_f->ppar[i], bnorm);
        real pperp = sqrt(pnorm * pnorm - mrk_f->ppar[i] * mrk_f->ppar[i]);
        real gamma = physlib_gamma_pnorm(mrk_f->mass, pnorm);
        real ekin = physlib_Ekin_gamma(mrk_f->mass, gamma);
        real pitch = mrk_f->ppar[i] / pnorm;
        real ptor = phys_ptoroid_gc(
            charge, mrk_f->r[i], mrk_f->ppar[i], psi, bnorm,
            mrk_f->bphi[i]);

        int valid = 1;
        axis = &hist->axes[15];
        size_t i15 = ((mrk_f->charge[i] - axis->min) /
                      (axis->max - axis->min)) *
                     axis->n;
        valid *= axis->n ? 1
                         : (mrk_f->charge[i] >= axis->min &&
                            mrk_f->charge[i] <= axis->max);

        axis = &hist->axes[14];
        size_t i14 =
            math_bin_index(mrk_f->time[i], axis->n, axis->min, axis->max);
        valid *= axis->n ? 1
                         : (mrk_f->time[i] > axis->min &&
                            mrk_f->time[i] < axis->max);

        axis = &hist->axes[13];
        size_t i13 = math_bin_index(ptor, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (ptor > axis->min && ptor < axis->max);

        axis = &hist->axes[12];
        size_t i12 = math_bin_index(mrk_f->mu[i], axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (mrk_f->mu[i] > axis->min && mrk_f->mu[i] < axis->max);

        axis = &hist->axes[11];
        size_t i11 = math_bin_index(pitch, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (pitch > axis->min && pitch < axis->max);

        axis = &hist->axes[10];
        size_t i10 = math_bin_index(ekin, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (ekin > axis->min && ekin < axis->max);

        // In GC mode pr, pphi, and pz are not defined
        size_t i9 = 0, i8 = 0, i7 = 0;

        axis = &hist->axes[6];
        size_t i6 = math_bin_index(pperp, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (pperp > axis->min && pperp < axis->max);

        axis = &hist->axes[5];
        size_t i5 =
            math_bin_index(mrk_f->ppar[i], axis->n, axis->min, axis->max);
        valid *= axis->n ? 1
                         : (mrk_f->ppar[i] > axis->min &&
                            mrk_f->ppar[i] < axis->max);

        axis = &hist->axes[4];
        size_t i4 = math_bin_index(theta, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (theta > axis->min && theta < axis->max);

        axis = &hist->axes[3];
        size_t i3 =
            math_bin_index(mrk_f->rho[i], axis->n, axis->min, axis->max);
        valid *= axis->n
                     ? 1
                     : (mrk_f->rho[i] > axis->min && mrk_f->rho[i] < axis->max);

        axis = &hist->axes[2];
        size_t i2 = math_bin_index(mrk_f->z[i], axis->n, axis->min, axis->max);
        valid *=
            axis->n ? 1 : (mrk_f->z[i] > axis->min && mrk_f->z[i] < axis->max);

        axis = &hist->axes[1];
        size_t i1 = math_bin_index(phi, axis->n, axis->min, axis->max);
        valid *= axis->n ? 1 : (phi > axis->min && phi < axis->max);

        axis = &hist->axes[0];
        size_t i0 = math_bin_index(mrk_f->r[i], axis->n, axis->min, axis->max);
        valid *=
            axis->n ? 1 : (mrk_f->r[i] > axis->min && mrk_f->r[i] < axis->max);

        size_t *n = hist->strides;
        size_t index = i0 * n[0] + i1 * n[1] + i2 * n[2] + i3 * n[3] +
                       i4 * n[4] + i5 * n[5] + i6 * n[6] + i7 * n[7] +
                       i8 * n[8] + i9 * n[9] + i10 * n[10] + i11 * n[11] +
                       i12 * n[12] + i13 * n[13] + i14 * n[14] + i15;
        real weight = mrk_f->weight[i] * (mrk_f->time[i] - mrk_i->time[i]);
        index = valid ? index : 0;

        GPU_ATOMIC
        hist->values[index] += valid * weight;
    }
}
