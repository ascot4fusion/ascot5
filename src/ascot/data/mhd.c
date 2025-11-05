/**
 * Implements mhd.h.
 */
#include "mhd.h"
#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd_dynamic.h"
#include "mhd_stationary.h"
#include "utils/mathlib.h"
#include <stdlib.h>

void Mhd_free(Mhd *mhd)
{
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        MhdStationary_free(mhd->stationary);
        break;
    case MHD_DYNAMIC:
        MhdDynamic_free(mhd->dynamic);
        break;
    }
}

void Mhd_offload(Mhd *mhd)
{
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        MhdStationary_offload(mhd->stationary);
        break;
    case MHD_DYNAMIC:
        MhdDynamic_offload(mhd->dynamic);
        break;
    }
}

err_t Mhd_eval_alpha_Phi(
    real alpha[5], real Phi[5], real r, real phi, real z, real t,
    size_t include_mode, Mhd *mhd, Bfield *bfield, Boozer *boozer)
{
    err_t err = 0;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        err = MhdStationary_eval_alpha_Phi(
            alpha, Phi, r, phi, z, t, include_mode, mhd->stationary, bfield,
            boozer);
        break;
    case MHD_DYNAMIC:
        err = MhdDynamic_eval_alpha_Phi(
            alpha, Phi, r, phi, z, t, include_mode, mhd->dynamic, bfield,
            boozer);
        break;
    }

    Phi[0] = err ? 0.0 : Phi[0];
    Phi[1] = err ? 0.0 : Phi[1];
    Phi[2] = err ? 0.0 : Phi[2];
    Phi[3] = err ? 0.0 : Phi[3];
    Phi[4] = err ? 0.0 : Phi[4];
    alpha[0] = err ? 0.0 : alpha[0];
    alpha[1] = err ? 0.0 : alpha[1];
    alpha[2] = err ? 0.0 : alpha[2];
    alpha[3] = err ? 0.0 : alpha[3];
    alpha[4] = err ? 0.0 : alpha[4];

    return err;
}

err_t Mhd_eval_perturbation(
    real b[3], real e[3], real Phi[1], real r, real phi, real z, real t,
    int include_background, size_t include_mode, Mhd *mhd, Bfield *bfield,
    Boozer *boozer)
{
    err_t err = 0;
    real alpha_dalpha[5], Phi_dPhi[5], B_dB[15];
    err = err ? err : Bfield_eval_b_db(B_dB, r, phi, z, t, bfield);

    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[1];
    B[2] = B_dB[2];

    real curlB[3];
    curlB[0] = B_dB[10] / r - B_dB[8];
    curlB[1] = B_dB[5] - B_dB[9];
    curlB[2] = (B_dB[1] - B_dB[4]) / r + B_dB[6];

    err = err ? err
              : Mhd_eval_alpha_Phi(
                    alpha_dalpha, Phi_dPhi, r, phi, z, t, include_mode,
                    mhd, bfield, boozer);

    real gradalphacrossB[3], *gradalpha = &alpha_dalpha[2];
    math_cross(gradalpha, B, gradalphacrossB);

    b[0] = alpha_dalpha[0] * curlB[0] + gradalphacrossB[0];
    b[1] = alpha_dalpha[0] * curlB[1] + gradalphacrossB[1];
    b[2] = alpha_dalpha[0] * curlB[2] + gradalphacrossB[2];
    b[0] += include_background ? B[0] : 0;
    b[1] += include_background ? B[1] : 0;
    b[2] += include_background ? B[2] : 0;

    e[0] = -Phi_dPhi[2] - B[0] * alpha_dalpha[1];
    e[1] = -Phi_dPhi[3] - B[1] * alpha_dalpha[1];
    e[2] = -Phi_dPhi[4] - B[2] * alpha_dalpha[1];
    Phi[0] = Phi_dPhi[0];

    b[0] = err ? 1.0 : b[0];
    b[1] = err ? 0.0 : b[1];
    b[2] = err ? 0.0 : b[2];
    e[0] = err ? 0.0 : e[0];
    e[1] = err ? 0.0 : e[1];
    e[2] = err ? 0.0 : e[2];
    Phi[0] = err ? 0.0 : Phi[0];

    return err;
}

size_t Mhd_get_n_modes(Mhd *mhd)
{
    int val = 0;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        val = mhd->stationary->n;
        break;
    case MHD_DYNAMIC:
        val = mhd->dynamic->n;
        break;
    }
    return val;
}

const int *Mhd_get_nmode(Mhd *mhd)
{
    const int *val = NULL;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        val = mhd->stationary->nmode;
        break;
    case MHD_DYNAMIC:
        val = mhd->dynamic->nmode;
        break;
    }
    return val;
}

const int *Mhd_get_mmode(Mhd *mhd)
{
    const int *val = NULL;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        val = mhd->stationary->mmode;
        break;
    case MHD_DYNAMIC:
        val = mhd->dynamic->mmode;
        break;
    }
    return val;
}

const real *Mhd_get_amplitude(Mhd *mhd)
{
    const real *val = NULL;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        val = mhd->stationary->amplitude;
        break;
    case MHD_DYNAMIC:
        val = mhd->dynamic->amplitude;
        break;
    }
    return val;
}

const real *Mhd_get_frequency(Mhd *mhd)
{
    const real *val = NULL;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        val = mhd->stationary->omega;
        break;
    case MHD_DYNAMIC:
        val = mhd->dynamic->omega;
        break;
    }
    return val;
}

const real *Mhd_get_phase(Mhd *mhd)
{
    const real *val = NULL;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        val = mhd->stationary->phase;
        break;
    case MHD_DYNAMIC:
        val = mhd->dynamic->phase;
        break;
    }
    return val;
}
