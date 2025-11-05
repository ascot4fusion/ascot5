/**
 * Implements marker.h.
 */
#include "marker.h"
#include "bfield.h"
#include "consts.h"
#include "datatypes.h"
#include "defines.h"
#include "efield.h"
#include "utils/gctransform.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * Allocate marker vector field or goto fail if malloc failed.
 *
 * Assumes variables ``err`` and ``vector_size`` are defined.
 *
 * @param a Pointer to the field to allocate.
 */
#define allocate_field(a)                                                      \
    do                                                                         \
    {                                                                          \
        a = xmalloc(&err, vector_size * sizeof(a));                            \
        if (err != 0)                                                          \
            goto fail;                                                         \
    } while (0)

int MarkerGyroOrbit_allocate(MarkerGyroOrbit *mrk, size_t vector_size)
{
    int err = 0;
    allocate_field(mrk->r);
    allocate_field(mrk->phi);
    allocate_field(mrk->z);
    allocate_field(mrk->p_r);
    allocate_field(mrk->p_phi);
    allocate_field(mrk->p_z);
    allocate_field(mrk->mass);
    allocate_field(mrk->charge);
    allocate_field(mrk->time);
    allocate_field(mrk->znum);
    allocate_field(mrk->anum);
    allocate_field(mrk->B_r);
    allocate_field(mrk->B_phi);
    allocate_field(mrk->B_z);
    allocate_field(mrk->B_r_dr);
    allocate_field(mrk->B_phi_dr);
    allocate_field(mrk->B_z_dr);
    allocate_field(mrk->B_r_dphi);
    allocate_field(mrk->B_phi_dphi);
    allocate_field(mrk->B_z_dphi);
    allocate_field(mrk->B_r_dz);
    allocate_field(mrk->B_phi_dz);
    allocate_field(mrk->B_z_dz);
    allocate_field(mrk->cputime);
    allocate_field(mrk->rho);
    allocate_field(mrk->theta);
    allocate_field(mrk->id);
    allocate_field(mrk->bounces);
    allocate_field(mrk->endcond);
    allocate_field(mrk->walltile);
    allocate_field(mrk->mileage);
    allocate_field(mrk->running);
    allocate_field(mrk->err);
    allocate_field(mrk->index);
    mrk->n_mrk = vector_size;
    return 0;

fail:
    MarkerGyroOrbit_deallocate(mrk);
    return 1;
}

void MarkerGyroOrbit_deallocate(MarkerGyroOrbit *mrk)
{
    free(mrk->r);
    free(mrk->phi);
    free(mrk->z);
    free(mrk->p_r);
    free(mrk->p_phi);
    free(mrk->p_z);
    free(mrk->mass);
    free(mrk->charge);
    free(mrk->time);
    free(mrk->znum);
    free(mrk->anum);
    free(mrk->B_r);
    free(mrk->B_phi);
    free(mrk->B_z);
    free(mrk->B_r_dr);
    free(mrk->B_phi_dr);
    free(mrk->B_z_dr);
    free(mrk->B_r_dphi);
    free(mrk->B_phi_dphi);
    free(mrk->B_z_dphi);
    free(mrk->B_r_dz);
    free(mrk->B_phi_dz);
    free(mrk->B_z_dz);
    free(mrk->cputime);
    free(mrk->rho);
    free(mrk->theta);
    free(mrk->id);
    free(mrk->bounces);
    free(mrk->endcond);
    free(mrk->walltile);
    free(mrk->mileage);
    free(mrk->running);
    free(mrk->err);
    free(mrk->index);
    mrk->n_mrk = 0;
}

void MarkerGyroOrbit_offload(MarkerGyroOrbit *p)
{
    SUPPRESS_UNUSED_WARNING(p);
    GPU_MAP_TO_DEVICE(
        p [0:1], p->running [0:p->n_mrk], p->r [0:p->n_mrk],
        p->phi [0:p->n_mrk], p->p_r [0:p->n_mrk], p->p_phi [0:p->n_mrk],
        p->p_z [0:p->n_mrk], p->mileage [0:p->n_mrk], p->z [0:p->n_mrk],
        p->charge [0:p->n_mrk], p->mass [0:p->n_mrk], p->B_r [0:p->n_mrk],
        p->B_r_dr [0:p->n_mrk], p->B_r_dphi [0:p->n_mrk],
        p->B_r_dz [0:p->n_mrk], p->B_phi [0:p->n_mrk], p->B_phi_dr [0:p->n_mrk],
        p->B_phi_dphi [0:p->n_mrk], p->B_phi_dz [0:p->n_mrk],
        p->B_z [0:p->n_mrk], p->B_z_dr [0:p->n_mrk], p->B_z_dphi [0:p->n_mrk],
        p->B_z_dz [0:p->n_mrk], p->rho [0:p->n_mrk], p->theta [0:p->n_mrk],
        p->err [0:p->n_mrk], p->time [0:p->n_mrk], p->weight [0:p->n_mrk],
        p->cputime [0:p->n_mrk], p->id [0:p->n_mrk], p->endcond [0:p->n_mrk],
        p->walltile [0:p->n_mrk], p->index [0:p->n_mrk], p->znum [0:p->n_mrk],
        p->anum [0:p->n_mrk], p->bounces [0:p->n_mrk])
}

void MarkerGyroOrbit_onload(MarkerGyroOrbit *p)
{
    SUPPRESS_UNUSED_WARNING(p);
    GPU_UPDATE_FROM_DEVICE(
        p->running [0:p->n_mrk], p->r [0:p->n_mrk], p->phi [0:p->n_mrk],
        p->p_r [0:p->n_mrk], p->p_phi [0:p->n_mrk], p->p_z [0:p->n_mrk],
        p->mileage [0:p->n_mrk], p->z [0:p->n_mrk], p->charge [0:p->n_mrk],
        p->mass [0:p->n_mrk], p->B_r [0:p->n_mrk], p->B_r_dr [0:p->n_mrk],
        p->B_r_dphi [0:p->n_mrk], p->B_r_dz [0:p->n_mrk], p->B_phi [0:p->n_mrk],
        p->B_phi_dr [0:p->n_mrk], p->B_phi_dphi [0:p->n_mrk],
        p->B_phi_dz [0:p->n_mrk], p->B_z [0:p->n_mrk], p->B_z_dr [0:p->n_mrk],
        p->B_z_dphi [0:p->n_mrk], p->B_z_dz [0:p->n_mrk], p->rho [0:p->n_mrk],
        p->theta [0:p->n_mrk], p->err [0:p->n_mrk], p->time [0:p->n_mrk],
        p->weight [0:p->n_mrk], p->cputime [0:p->n_mrk], p->id [0:p->n_mrk],
        p->endcond [0:p->n_mrk], p->walltile [0:p->n_mrk],
        p->index [0:p->n_mrk], p->znum [0:p->n_mrk], p->anum [0:p->n_mrk],
        p->bounces [0:p->n_mrk])
}

void MarkerGyroOrbit_copy(
    MarkerGyroOrbit *copy, MarkerGyroOrbit *original, size_t index)
{
    copy->r[index] = original->r[index];
    copy->phi[index] = original->phi[index];
    copy->z[index] = original->z[index];
    copy->p_r[index] = original->p_r[index];
    copy->p_phi[index] = original->p_phi[index];
    copy->p_z[index] = original->p_z[index];
    copy->time[index] = original->time[index];
    copy->mileage[index] = original->mileage[index];
    copy->cputime[index] = original->cputime[index];
    copy->rho[index] = original->rho[index];
    copy->weight[index] = original->weight[index];
    copy->cputime[index] = original->cputime[index];
    copy->rho[index] = original->rho[index];
    copy->theta[index] = original->theta[index];
    copy->mass[index] = original->mass[index];
    copy->charge[index] = original->charge[index];
    copy->znum[index] = original->znum[index];
    copy->anum[index] = original->anum[index];
    copy->id[index] = original->id[index];
    copy->bounces[index] = original->bounces[index];
    copy->running[index] = original->running[index];
    copy->endcond[index] = original->endcond[index];
    copy->walltile[index] = original->walltile[index];
    copy->B_r[index] = original->B_r[index];
    copy->B_phi[index] = original->B_phi[index];
    copy->B_z[index] = original->B_z[index];
    copy->B_r_dr[index] = original->B_r_dr[index];
    copy->B_r_dphi[index] = original->B_r_dphi[index];
    copy->B_r_dz[index] = original->B_r_dz[index];
    copy->B_phi_dr[index] = original->B_phi_dr[index];
    copy->B_phi_dphi[index] = original->B_phi_dphi[index];
    copy->B_phi_dz[index] = original->B_phi_dz[index];
    copy->B_z_dr[index] = original->B_z_dr[index];
    copy->B_z_dphi[index] = original->B_z_dphi[index];
    copy->B_z_dz[index] = original->B_z_dz[index];
}

int MarkerGyroOrbit_from_queue(
    MarkerGyroOrbit *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield)
{
    State *p = queue->p[queue_index];
    err_t err = p->err;

    real B_dB[15], psi[1], rho[2];
    if (!err)
        err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);

    if (!err)
    {
        mrk->r[mrk_index] = p->rprt;
        mrk->phi[mrk_index] = p->phiprt;
        mrk->z[mrk_index] = p->zprt;
        mrk->p_r[mrk_index] = p->pr;
        mrk->p_phi[mrk_index] = p->pphi;
        mrk->p_z[mrk_index] = p->pz;
        mrk->mass[mrk_index] = p->mass;
        mrk->charge[mrk_index] = p->charge;
        mrk->znum[mrk_index] = p->znum;
        mrk->anum[mrk_index] = p->anum;
        mrk->bounces[mrk_index] = 0;
        mrk->weight[mrk_index] = p->weight;
        mrk->time[mrk_index] = p->time;
        mrk->theta[mrk_index] = p->theta;
        mrk->id[mrk_index] = p->id;
        mrk->endcond[mrk_index] = p->endcond;
        mrk->walltile[mrk_index] = p->walltile;
        mrk->mileage[mrk_index] = p->mileage;
        mrk->rho[mrk_index] = rho[0];
        mrk->B_r[mrk_index] = B_dB[0];
        mrk->B_r_dr[mrk_index] = B_dB[3];
        mrk->B_r_dphi[mrk_index] = B_dB[4];
        mrk->B_r_dz[mrk_index] = B_dB[5];
        mrk->B_phi[mrk_index] = B_dB[1];
        mrk->B_phi_dr[mrk_index] = B_dB[6];
        mrk->B_phi_dphi[mrk_index] = B_dB[7];
        mrk->B_phi_dz[mrk_index] = B_dB[8];
        mrk->B_z[mrk_index] = B_dB[2];
        mrk->B_z_dr[mrk_index] = B_dB[9];
        mrk->B_z_dphi[mrk_index] = B_dB[10];
        mrk->B_z_dz[mrk_index] = B_dB[11];
        mrk->running[mrk_index] = p->endcond == 0;
        mrk->cputime[mrk_index] = p->cputime;
        mrk->index[mrk_index] = queue_index;
        mrk->err[mrk_index] = 0;
    }
    if (err)
        p->err = err;

    return err > 0;
}

void MarkerGyroOrbit_to_queue(
    MarkerQueue *queue, MarkerGyroOrbit *mrk, size_t index, Bfield *bfield)
{
    err_t err = 0;
    State *p = queue->p[mrk->index[index]];
    p->rprt = mrk->r[index];
    p->phiprt = mrk->phi[index];
    p->zprt = mrk->z[index];
    p->pr = mrk->p_r[index];
    p->pphi = mrk->p_phi[index];
    p->pz = mrk->p_z[index];
    p->mass = mrk->mass[index];
    p->charge = mrk->charge[index];
    p->znum = mrk->znum[index];
    p->anum = mrk->anum[index];
    p->weight = mrk->weight[index];
    p->time = mrk->time[index];
    p->theta = mrk->theta[index];
    p->id = mrk->id[index];
    p->endcond = mrk->endcond[index];
    p->walltile = mrk->walltile[index];
    p->cputime = mrk->cputime[index];
    p->mileage = mrk->mileage[index];

    /* Particle to guiding center */
    real B_dB[15], psi[1], rho[2], ppar, mu;
    rho[0] = mrk->rho[index];
    B_dB[0] = mrk->B_r[index];
    B_dB[1] = mrk->B_r_dr[index];
    B_dB[2] = mrk->B_r_dphi[index];
    B_dB[3] = mrk->B_r_dz[index];
    B_dB[4] = mrk->B_phi[index];
    B_dB[5] = mrk->B_phi_dr[index];
    B_dB[6] = mrk->B_phi_dphi[index];
    B_dB[7] = mrk->B_phi_dz[index];
    B_dB[8] = mrk->B_z[index];
    B_dB[9] = mrk->B_z_dr[index];
    B_dB[10] = mrk->B_z_dphi[index];
    B_dB[11] = mrk->B_z_dz[index];

    gctransform_particle2guidingcenter(
        p->mass, p->charge, B_dB, p->rprt, p->phiprt, p->zprt, p->pr, p->pphi, p->pz,
        &p->r, &p->phi, &p->z, &ppar, &mu, &p->zeta);

    if (!err)
        err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);

    real Bnorm = math_normc(B_dB[0], B_dB[1], B_dB[2]);
    p->ekin = physlib_Ekin_ppar(p->mass, mu, ppar, Bnorm);
    p->pitch = physlib_gc_xi(p->mass, mu, ppar, Bnorm);

    /* If marker already has error flag, make sure it is not overwritten here */
    p->err = mrk->err[index] ? mrk->err[index] : err;
}

int MarkerGuidingCenter_allocate(MarkerGuidingCenter *mrk, size_t vector_size)
{
    int err = 0;
    allocate_field(mrk->r);
    allocate_field(mrk->phi);
    allocate_field(mrk->z);
    allocate_field(mrk->ppar);
    allocate_field(mrk->mu);
    allocate_field(mrk->zeta);
    allocate_field(mrk->B_r);
    allocate_field(mrk->B_phi);
    allocate_field(mrk->B_z);
    allocate_field(mrk->B_r_dr);
    allocate_field(mrk->B_phi_dr);
    allocate_field(mrk->B_z_dr);
    allocate_field(mrk->B_r_dphi);
    allocate_field(mrk->B_phi_dphi);
    allocate_field(mrk->B_z_dphi);
    allocate_field(mrk->B_r_dz);
    allocate_field(mrk->B_phi_dz);
    allocate_field(mrk->B_z_dz);
    allocate_field(mrk->cputime);
    allocate_field(mrk->rho);
    allocate_field(mrk->theta);
    allocate_field(mrk->id);
    allocate_field(mrk->bounces);
    allocate_field(mrk->endcond);
    allocate_field(mrk->walltile);
    allocate_field(mrk->mileage);
    allocate_field(mrk->running);
    allocate_field(mrk->err);
    allocate_field(mrk->index);
    mrk->n_mrk = vector_size;
    return 0;

fail:
    MarkerGuidingCenter_deallocate(mrk);
    return 1;
}

void MarkerGuidingCenter_deallocate(MarkerGuidingCenter *mrk)
{
    free(mrk->r);
    free(mrk->phi);
    free(mrk->z);
    free(mrk->ppar);
    free(mrk->mu);
    free(mrk->zeta);
    free(mrk->B_r);
    free(mrk->B_phi);
    free(mrk->B_z);
    free(mrk->B_r_dr);
    free(mrk->B_phi_dr);
    free(mrk->B_z_dr);
    free(mrk->B_r_dphi);
    free(mrk->B_phi_dphi);
    free(mrk->B_z_dphi);
    free(mrk->B_r_dz);
    free(mrk->B_phi_dz);
    free(mrk->B_z_dz);
    free(mrk->cputime);
    free(mrk->rho);
    free(mrk->theta);
    free(mrk->id);
    free(mrk->bounces);
    free(mrk->endcond);
    free(mrk->walltile);
    free(mrk->mileage);
    free(mrk->running);
    free(mrk->err);
    free(mrk->index);
    mrk->n_mrk = 0;
}

void MarkerGuidingCenter_offload(MarkerGuidingCenter *mrk) { (void)mrk; }

void MarkerGuidingCenter_onload(MarkerGuidingCenter *mrk) { (void)mrk; }

void MarkerGuidingCenter_copy(
    MarkerGuidingCenter *copy, MarkerGuidingCenter *original, size_t index)
{
    copy->r[index] = original->r[index];
    copy->phi[index] = original->phi[index];
    copy->z[index] = original->z[index];
    copy->ppar[index] = original->ppar[index];
    copy->mu[index] = original->mu[index];
    copy->zeta[index] = original->zeta[index];

    copy->time[index] = original->time[index];
    copy->mileage[index] = original->mileage[index];
    copy->weight[index] = original->weight[index];
    copy->cputime[index] = original->cputime[index];
    copy->rho[index] = original->rho[index];
    copy->theta[index] = original->theta[index];

    copy->mass[index] = original->mass[index];
    copy->charge[index] = original->charge[index];

    copy->id[index] = original->id[index];
    copy->bounces[index] = original->bounces[index];
    copy->running[index] = original->running[index];
    copy->endcond[index] = original->endcond[index];
    copy->walltile[index] = original->walltile[index];

    copy->B_r[index] = original->B_r[index];
    copy->B_phi[index] = original->B_phi[index];
    copy->B_z[index] = original->B_z[index];

    copy->B_r_dr[index] = original->B_r_dr[index];
    copy->B_r_dphi[index] = original->B_r_dphi[index];
    copy->B_r_dz[index] = original->B_r_dz[index];

    copy->B_phi_dr[index] = original->B_phi_dr[index];
    copy->B_phi_dphi[index] = original->B_phi_dphi[index];
    copy->B_phi_dz[index] = original->B_phi_dz[index];

    copy->B_z_dr[index] = original->B_z_dr[index];
    copy->B_z_dphi[index] = original->B_z_dphi[index];
    copy->B_z_dz[index] = original->B_z_dz[index];
}

int MarkerGuidingCenter_from_queue(
    MarkerGuidingCenter *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield)
{
    State *p = queue->p[queue_index];
    err_t err = p->err;

    real B_dB[15], psi[1], rho[2];
    if (!err)
        err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);

    if (!err)
    {
        real Bnorm = math_normc(B_dB[0], B_dB[1], B_dB[2]);
        real gamma = physlib_gamma_Ekin(p->mass, p->ekin);
        real pnorm = physlib_pnorm_gamma(p->mass, gamma);
        mrk->r[mrk_index] = p->r;
        mrk->phi[mrk_index] = p->phi;
        mrk->z[mrk_index] = p->z;
        mrk->zeta[mrk_index] = p->zeta;
        mrk->mass[mrk_index] = p->mass;
        mrk->charge[mrk_index] = p->charge;
        mrk->time[mrk_index] = p->time;
        mrk->bounces[mrk_index] = 0;
        mrk->weight[mrk_index] = p->weight;
        mrk->theta[mrk_index] = p->theta;
        mrk->id[mrk_index] = p->id;
        mrk->endcond[mrk_index] = p->endcond;
        mrk->walltile[mrk_index] = p->walltile;
        mrk->mileage[mrk_index] = p->mileage;
        mrk->rho[mrk_index] = rho[0];
        mrk->B_r[mrk_index] = B_dB[0];
        mrk->B_r_dr[mrk_index] = B_dB[3];
        mrk->B_r_dphi[mrk_index] = B_dB[4];
        mrk->B_r_dz[mrk_index] = B_dB[5];
        mrk->B_phi[mrk_index] = B_dB[1];
        mrk->B_phi_dr[mrk_index] = B_dB[6];
        mrk->B_phi_dphi[mrk_index] = B_dB[7];
        mrk->B_phi_dz[mrk_index] = B_dB[8];
        mrk->B_z[mrk_index] = B_dB[2];
        mrk->B_z_dr[mrk_index] = B_dB[9];
        mrk->B_z_dphi[mrk_index] = B_dB[10];
        mrk->B_z_dz[mrk_index] = B_dB[11];
        mrk->mu[mrk_index] = physlib_gc_mu(p->mass, pnorm, p->pitch, Bnorm);
        mrk->ppar[mrk_index] =
            phys_ppar_Ekin(p->mass, p->ekin, mrk->mu[mrk_index], Bnorm);
        mrk->running[mrk_index] = p->endcond == 0;
        mrk->cputime[mrk_index] = p->cputime;
        mrk->index[mrk_index] = queue_index;
        mrk->err[mrk_index] = 0;
    }
    if (err)
        p->err = err;

    return err > 0;
}

void MarkerGuidingCenter_to_queue(
    MarkerQueue *queue, MarkerGuidingCenter *mrk, size_t index, Bfield *bfield)
{
    err_t err = 0;
    State *p = queue->p[mrk->index[index]];
    p->r = mrk->r[index];
    p->phi = mrk->phi[index];
    p->z = mrk->z[index];

    p->mass = mrk->mass[index];
    p->charge = mrk->charge[index];
    p->time = mrk->time[index];
    p->weight = mrk->weight[index];
    p->id = mrk->id[index];
    p->cputime = mrk->cputime[index];
    p->theta = mrk->theta[index];
    p->endcond = mrk->endcond[index];
    p->walltile = mrk->walltile[index];
    p->mileage = mrk->mileage[index];

    /* Guiding center to particle transformation */
    real B_dB[15];
    B_dB[0] = mrk->B_r[index];
    B_dB[3] = mrk->B_r_dr[index];
    B_dB[4] = mrk->B_r_dphi[index];
    B_dB[5] = mrk->B_r_dz[index];
    B_dB[1] = mrk->B_phi[index];
    B_dB[6] = mrk->B_phi_dr[index];
    B_dB[7] = mrk->B_phi_dphi[index];
    B_dB[8] = mrk->B_phi_dz[index];
    B_dB[2] = mrk->B_z[index];
    B_dB[9] = mrk->B_z_dr[index];
    B_dB[10] = mrk->B_z_dphi[index];
    B_dB[11] = mrk->B_z_dz[index];

    real Bnorm =
        math_normc(mrk->B_r[index], mrk->B_phi[index], mrk->B_z[index]);
    p->ekin = physlib_Ekin_ppar(
        mrk->mass[index], mrk->mu[index], mrk->ppar[index], Bnorm);
    p->pitch = physlib_gc_xi(
        mrk->mass[index], mrk->mu[index], mrk->ppar[index], Bnorm);
    p->zeta = mrk->zeta[index];

    real pparprt, muprt, zetaprt;
    gctransform_guidingcenter2particle(
        p->mass, p->charge, B_dB, p->r, p->phi, p->z, mrk->ppar[index],
        mrk->mu[index], p->zeta, &p->rprt, &p->phiprt, &p->zprt, &pparprt,
        &muprt, &zetaprt);

    if (!err)
        err = Bfield_eval_b_db(B_dB, p->rprt, p->phiprt, p->zprt, p->time, bfield);

    gctransform_pparmuzeta2prpphipz(
        p->mass, p->charge, B_dB, p->phiprt, pparprt, muprt, zetaprt, &p->pr,
        &p->pphi, &p->pz);

    /* If marker already has error flag, make sure it is not overwritten here */
    p->err = mrk->err[index] ? mrk->err[index] : err;
}

int MarkerFieldLine_allocate(MarkerFieldLine *mrk, size_t vector_size)
{
    int err = 0;
    allocate_field(mrk->r);
    allocate_field(mrk->phi);
    allocate_field(mrk->z);
    allocate_field(mrk->pitch);
    allocate_field(mrk->B_r);
    allocate_field(mrk->B_phi);
    allocate_field(mrk->B_z);
    allocate_field(mrk->B_r_dr);
    allocate_field(mrk->B_phi_dr);
    allocate_field(mrk->B_z_dr);
    allocate_field(mrk->B_r_dphi);
    allocate_field(mrk->B_phi_dphi);
    allocate_field(mrk->B_z_dphi);
    allocate_field(mrk->B_r_dz);
    allocate_field(mrk->B_phi_dz);
    allocate_field(mrk->B_z_dz);
    allocate_field(mrk->cputime);
    allocate_field(mrk->rho);
    allocate_field(mrk->theta);
    allocate_field(mrk->id);
    allocate_field(mrk->time);
    allocate_field(mrk->endcond);
    allocate_field(mrk->walltile);
    allocate_field(mrk->mileage);
    allocate_field(mrk->running);
    allocate_field(mrk->err);
    allocate_field(mrk->index);
    mrk->n_mrk = vector_size;
    return 0;

fail:
    MarkerFieldLine_deallocate(mrk);
    return 1;
}

void MarkerFieldLine_deallocate(MarkerFieldLine *mrk)
{
    free(mrk->r);
    free(mrk->phi);
    free(mrk->z);
    free(mrk->pitch);
    free(mrk->B_r);
    free(mrk->B_phi);
    free(mrk->B_z);
    free(mrk->B_r_dr);
    free(mrk->B_phi_dr);
    free(mrk->B_z_dr);
    free(mrk->B_r_dphi);
    free(mrk->B_phi_dphi);
    free(mrk->B_z_dphi);
    free(mrk->B_r_dz);
    free(mrk->B_phi_dz);
    free(mrk->B_z_dz);
    free(mrk->cputime);
    free(mrk->rho);
    free(mrk->theta);
    free(mrk->id);
    free(mrk->time);
    free(mrk->endcond);
    free(mrk->walltile);
    free(mrk->mileage);
    free(mrk->running);
    free(mrk->err);
    free(mrk->index);
    mrk->n_mrk = 0;
}

void MarkerFieldLine_offload(MarkerFieldLine *mrk) { (void)mrk; }

void MarkerFieldLine_onload(MarkerFieldLine *mrk) { (void)mrk; }

void MarkerFieldLine_copy(
    MarkerFieldLine *copy, MarkerFieldLine *original, size_t index)
{
    copy->r[index] = original->r[index];
    copy->phi[index] = original->phi[index];
    copy->z[index] = original->z[index];
    copy->pitch[index] = original->pitch[index];
    copy->time[index] = original->time[index];
    copy->mileage[index] = original->mileage[index];
    copy->cputime[index] = original->cputime[index];
    copy->rho[index] = original->rho[index];
    copy->theta[index] = original->theta[index];
    copy->id[index] = original->id[index];
    copy->running[index] = original->running[index];
    copy->endcond[index] = original->endcond[index];
    copy->walltile[index] = original->walltile[index];
    copy->B_r[index] = original->B_r[index];
    copy->B_phi[index] = original->B_phi[index];
    copy->B_z[index] = original->B_z[index];
    copy->B_r_dr[index] = original->B_r_dr[index];
    copy->B_r_dphi[index] = original->B_r_dphi[index];
    copy->B_r_dz[index] = original->B_r_dz[index];
    copy->B_phi_dr[index] = original->B_phi_dr[index];
    copy->B_phi_dphi[index] = original->B_phi_dphi[index];
    copy->B_phi_dz[index] = original->B_phi_dz[index];
    copy->B_z_dr[index] = original->B_z_dr[index];
    copy->B_z_dphi[index] = original->B_z_dphi[index];
    copy->B_z_dz[index] = original->B_z_dz[index];
    copy->err[index] = original->err[index];
    copy->index[index] = original->index[index];
}

int MarkerFieldLine_from_queue(
    MarkerFieldLine *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield)
{
    State *p = queue->p[queue_index];
    err_t err = p->err;

    real B_dB[15], psi[1], rho[2];
    if (!err)
        err = Bfield_eval_b_db(B_dB, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err) {
        err = Bfield_eval_rho(rho, psi[0], bfield);
    }

    if (!err)
    {
        mrk->r[mrk_index] = p->r;
        mrk->phi[mrk_index] = p->phi;
        mrk->z[mrk_index] = p->z;
        mrk->pitch[mrk_index] = 2 * (p->pitch >= 0) - 1.0;
        mrk->time[mrk_index] = p->time;
        mrk->id[mrk_index] = p->id;
        mrk->cputime[mrk_index] = p->cputime;
        mrk->theta[mrk_index] = p->theta;
        mrk->endcond[mrk_index] = p->endcond;
        mrk->walltile[mrk_index] = p->walltile;
        mrk->mileage[mrk_index] = p->mileage;
        mrk->rho[mrk_index] = rho[0];
        mrk->B_r[mrk_index] = B_dB[0];
        mrk->B_r_dr[mrk_index] = B_dB[3];
        mrk->B_r_dphi[mrk_index] = B_dB[4];
        mrk->B_r_dz[mrk_index] = B_dB[5];
        mrk->B_phi[mrk_index] = B_dB[1];
        mrk->B_phi_dr[mrk_index] = B_dB[6];
        mrk->B_phi_dphi[mrk_index] = B_dB[7];
        mrk->B_phi_dz[mrk_index] = B_dB[8];
        mrk->B_z[mrk_index] = B_dB[2];
        mrk->B_z_dr[mrk_index] = B_dB[9];
        mrk->B_z_dphi[mrk_index] = B_dB[10];
        mrk->B_z_dz[mrk_index] = B_dB[11];
        mrk->running[mrk_index] = p->endcond == 0;
        mrk->index[mrk_index] = queue_index;
        mrk->err[mrk_index] = 0;
    }
    if (err)
        p->err = err;

    return err > 0;
}

void MarkerFieldLine_to_queue(
    MarkerQueue *queue, MarkerFieldLine *mrk, size_t index)
{
    State *p = queue->p[mrk->index[index]];
    p->rprt = mrk->r[index];
    p->phiprt = mrk->phi[index];
    p->zprt = mrk->z[index];
    p->pr = 0;
    p->pphi = 0;
    p->pz = 0;
    p->r = mrk->r[index];
    p->phi = mrk->phi[index];
    p->z = mrk->z[index];
    p->ekin = 0;
    p->pitch = mrk->pitch[index];
    p->zeta = 0;
    p->mass = 0;
    p->charge = 0;
    p->time = mrk->time[index];
    p->id = mrk->id[index];
    p->cputime = mrk->cputime[index];
    p->theta = mrk->theta[index];
    p->endcond = mrk->endcond[index];
    p->walltile = mrk->walltile[index];
    p->mileage = mrk->mileage[index];
    p->err = mrk->err[index];
}

int marker_go_to_gc(
    MarkerGyroOrbit *p_go, size_t j, MarkerGuidingCenter *p_gc, Bfield *bfield)
{
    err_t err = p_go->err[j];
    real axisrz[2];
    int simerr = 0; /* Error has already occurred */
    if (err)
    {
        simerr = 1;
    }
    p_gc->id[j] = p_go->id[j];
    p_gc->index[j] = p_go->index[j];

    real r, phi, z, ppar, mu, zeta, B_dB[15];
    if (!err)
    {
        real Rprt = p_go->r[j];
        real phiprt = p_go->phi[j];
        real zprt = p_go->z[j];
        real pr = p_go->p_r[j];
        real pphi = p_go->p_phi[j];
        real pz = p_go->p_z[j];
        real mass = p_go->mass[j];
        real charge = p_go->charge[j];

        p_gc->mass[j] = p_go->mass[j];
        p_gc->charge[j] = p_go->charge[j];
        p_gc->weight[j] = p_go->weight[j];
        p_gc->time[j] = p_go->time[j];
        p_gc->mileage[j] = p_go->mileage[j];
        p_gc->endcond[j] = p_go->endcond[j];
        p_gc->running[j] = p_go->running[j];
        p_gc->walltile[j] = p_go->walltile[j];
        p_gc->cputime[j] = p_go->cputime[j];

        B_dB[0] = p_go->B_r[j];
        B_dB[3] = p_go->B_r_dr[j];
        B_dB[4] = p_go->B_r_dphi[j];
        B_dB[5] = p_go->B_r_dz[j];
        B_dB[1] = p_go->B_phi[j];
        B_dB[6] = p_go->B_phi_dr[j];
        B_dB[7] = p_go->B_phi_dphi[j];
        B_dB[8] = p_go->B_phi_dz[j];
        B_dB[2] = p_go->B_z[j];
        B_dB[9] = p_go->B_z_dr[j];
        B_dB[10] = p_go->B_z_dphi[j];
        B_dB[11] = p_go->B_z_dz[j];

        /* Guiding center transformation */
        gctransform_particle2guidingcenter(
            mass, charge, B_dB, Rprt, phiprt, zprt, pr, pphi, pz, &r, &phi, &z,
            &ppar, &mu, &zeta);
    }
    if (!err && r <= 0)
    {
        err = ERROR_RAISE(ERR_UNPHYSICAL_MARKER, DATA_MARKER_C);
    }
    if (!err && mu < 0)
    {
        err = ERROR_RAISE(ERR_UNPHYSICAL_MARKER, DATA_MARKER_C);
    }

    real psi[1], rho[2];
    if (!err)
    {
        err = Bfield_eval_b_db(B_dB, r, phi, z, p_go->time[j], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_psi(psi, r, phi, z, p_go->time[j], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_rho(rho, psi[0], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_axis_rz(axisrz, bfield, p_gc->phi[j]);
    }

    if (!err)
    {
        p_gc->r[j] = r;
        p_gc->phi[j] = phi;
        p_gc->z[j] = z;
        p_gc->mu[j] = mu;
        p_gc->zeta[j] = zeta;
        p_gc->ppar[j] = ppar;
        p_gc->rho[j] = rho[0];

        /* Evaluate pol angle so that it is cumulative and at gc position */
        p_gc->theta[j] = p_go->theta[j];
        p_gc->theta[j] += atan2(
            (p_go->r[j] - axisrz[0]) * (p_gc->z[j] - axisrz[1]) -
                (p_go->z[j] - axisrz[1]) * (p_gc->r[j] - axisrz[0]),
            (p_go->r[j] - axisrz[0]) * (p_gc->r[j] - axisrz[0]) +
                (p_go->z[j] - axisrz[1]) * (p_gc->z[j] - axisrz[1]));

        p_gc->B_r[j] = B_dB[0];
        p_gc->B_r_dr[j] = B_dB[3];
        p_gc->B_r_dphi[j] = B_dB[4];
        p_gc->B_r_dz[j] = B_dB[5];

        p_gc->B_phi[j] = B_dB[1];
        p_gc->B_phi_dr[j] = B_dB[6];
        p_gc->B_phi_dphi[j] = B_dB[7];
        p_gc->B_phi_dz[j] = B_dB[8];

        p_gc->B_z[j] = B_dB[2];
        p_gc->B_z_dr[j] = B_dB[9];
        p_gc->B_z_dphi[j] = B_dB[10];
        p_gc->B_z_dz[j] = B_dB[11];
    }
    if (!simerr)
    {
        err = err;
    }
    p_gc->err[j] = err;
    if (p_gc->err[j])
    {
        p_gc->running[j] = 0;
        p_gc->endcond[j] = 0;
    }

    return err > 0;
}

#undef allocate_field
