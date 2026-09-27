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
    allocate_field(mrk->charge);
    allocate_field(mrk->time);
    allocate_field(mrk->br);
    allocate_field(mrk->bphi);
    allocate_field(mrk->bz);
    allocate_field(mrk->dbrdr);
    allocate_field(mrk->dbphidr);
    allocate_field(mrk->dbzdr);
    allocate_field(mrk->dbrdphi);
    allocate_field(mrk->dbphidphi);
    allocate_field(mrk->dbzdphi);
    allocate_field(mrk->dbrdz);
    allocate_field(mrk->dbphidz);
    allocate_field(mrk->dbzdz);
    allocate_field(mrk->cputime);
    allocate_field(mrk->rho);
    allocate_field(mrk->theta);
    allocate_field(mrk->weight);
    allocate_field(mrk->id);
    allocate_field(mrk->bounces);
    allocate_field(mrk->endcond);
    allocate_field(mrk->walltile);
    allocate_field(mrk->mileage);
    allocate_field(mrk->running);
    allocate_field(mrk->err);
    allocate_field(mrk->index);
    mrk->size = vector_size;
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
    free(mrk->charge);
    free(mrk->time);
    free(mrk->br);
    free(mrk->bphi);
    free(mrk->bz);
    free(mrk->dbrdr);
    free(mrk->dbphidr);
    free(mrk->dbzdr);
    free(mrk->dbrdphi);
    free(mrk->dbphidphi);
    free(mrk->dbzdphi);
    free(mrk->dbrdz);
    free(mrk->dbphidz);
    free(mrk->dbzdz);
    free(mrk->cputime);
    free(mrk->rho);
    free(mrk->theta);
    free(mrk->weight);
    free(mrk->id);
    free(mrk->bounces);
    free(mrk->endcond);
    free(mrk->walltile);
    free(mrk->mileage);
    free(mrk->running);
    free(mrk->err);
    free(mrk->index);
    mrk->size = 0;
}

void MarkerGyroOrbit_offload(MarkerGyroOrbit *p)
{
    SUPPRESS_UNUSED_WARNING(p);
    GPU_MAP_TO_DEVICE(
        p [0:1], p->running [0:p->size], p->r [0:p->size], p->phi [0:p->size],
        p->p_r [0:p->size], p->p_phi [0:p->size], p->p_z [0:p->size],
        p->mileage [0:p->size], p->z [0:p->size], p->charge [0:p->size],
        p->br [0:p->size], p->dbrdr [0:p->size], p->dbrdphi [0:p->size],
        p->dbrdz [0:p->size], p->bphi [0:p->size], p->dbphidr [0:p->size],
        p->dbphidphi [0:p->size], p->dbphidz [0:p->size], p->bz [0:p->size],
        p->dbzdr [0:p->size], p->dbzdphi [0:p->size], p->dbzdz [0:p->size],
        p->rho [0:p->size], p->theta [0:p->size], p->err [0:p->size],
        p->time [0:p->size], p->weight [0:p->size], p->cputime [0:p->size],
        p->id [0:p->size], p->endcond [0:p->size], p->walltile [0:p->size],
        p->index [0:p->size], p->bounces [0:p->size])
}

void MarkerGyroOrbit_onload(MarkerGyroOrbit *p)
{
    SUPPRESS_UNUSED_WARNING(p);
    GPU_UPDATE_FROM_DEVICE(
        p->running [0:p->size], p->r [0:p->size], p->phi [0:p->size],
        p->p_r [0:p->size], p->p_phi [0:p->size], p->p_z [0:p->size],
        p->mileage [0:p->size], p->z [0:p->size], p->charge [0:p->size],
        p->br [0:p->size], p->dbrdr [0:p->size], p->dbrdphi [0:p->size],
        p->dbrdz [0:p->size], p->bphi [0:p->size], p->dbphidr [0:p->size],
        p->dbphidphi [0:p->size], p->dbphidz [0:p->size], p->bz [0:p->size],
        p->dbzdr [0:p->size], p->dbzdphi [0:p->size], p->dbzdz [0:p->size],
        p->rho [0:p->size], p->theta [0:p->size], p->err [0:p->size],
        p->time [0:p->size], p->weight [0:p->size], p->cputime [0:p->size],
        p->id [0:p->size], p->endcond [0:p->size], p->walltile [0:p->size],
        p->index [0:p->size], p->bounces [0:p->size])
}

void MarkerGyroOrbit_copy(
    MarkerGyroOrbit *copy, MarkerGyroOrbit *original, size_t index)
{
    copy->mass = original->mass;
    copy->znum = original->znum;
    copy->anum = original->anum;
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
    copy->charge[index] = original->charge[index];
    copy->id[index] = original->id[index];
    copy->bounces[index] = original->bounces[index];
    copy->running[index] = original->running[index];
    copy->endcond[index] = original->endcond[index];
    copy->walltile[index] = original->walltile[index];
    copy->br[index] = original->br[index];
    copy->bphi[index] = original->bphi[index];
    copy->bz[index] = original->bz[index];
    copy->dbrdr[index] = original->dbrdr[index];
    copy->dbrdphi[index] = original->dbrdphi[index];
    copy->dbrdz[index] = original->dbrdz[index];
    copy->dbphidr[index] = original->dbphidr[index];
    copy->dbphidphi[index] = original->dbphidphi[index];
    copy->dbphidz[index] = original->dbphidz[index];
    copy->dbzdr[index] = original->dbzdr[index];
    copy->dbzdphi[index] = original->dbzdphi[index];
    copy->dbzdz[index] = original->dbzdz[index];
}

int MarkerGyroOrbit_from_queue(
    MarkerGyroOrbit *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield)
{
    State *p = queue->p[queue_index];
    err_t err = p->err;

    real b_db[15], psi[1], rho[2];
    if (!err)
        err = Bfield_eval_b_db(b_db, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);
    if (!err)
    {
        mrk->mass = p->mass;
        mrk->znum = p->znum;
        mrk->anum = p->anum;
        mrk->r[mrk_index] = p->rprt;
        mrk->phi[mrk_index] = p->phiprt;
        mrk->z[mrk_index] = p->zprt;
        mrk->p_r[mrk_index] = p->pr;
        mrk->p_phi[mrk_index] = p->pphi;
        mrk->p_z[mrk_index] = p->pz;
        mrk->charge[mrk_index] = p->charge;
        mrk->bounces[mrk_index] = 0;
        mrk->weight[mrk_index] = p->weight;
        mrk->time[mrk_index] = p->time;
        mrk->theta[mrk_index] = p->theta;
        mrk->id[mrk_index] = p->id;
        mrk->endcond[mrk_index] = p->endcond;
        mrk->walltile[mrk_index] = p->walltile;
        mrk->mileage[mrk_index] = p->mileage;
        mrk->rho[mrk_index] = rho[0];
        mrk->br[mrk_index] = b_db[0];
        mrk->dbrdr[mrk_index] = b_db[3];
        mrk->dbrdphi[mrk_index] = b_db[4];
        mrk->dbrdz[mrk_index] = b_db[5];
        mrk->bphi[mrk_index] = b_db[1];
        mrk->dbphidr[mrk_index] = b_db[6];
        mrk->dbphidphi[mrk_index] = b_db[7];
        mrk->dbphidz[mrk_index] = b_db[8];
        mrk->bz[mrk_index] = b_db[2];
        mrk->dbzdr[mrk_index] = b_db[9];
        mrk->dbzdphi[mrk_index] = b_db[10];
        mrk->dbzdz[mrk_index] = b_db[11];
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
    p->znum = mrk->znum;
    p->anum = mrk->anum;
    p->mass = mrk->mass;
    p->rprt = mrk->r[index];
    p->phiprt = mrk->phi[index];
    p->zprt = mrk->z[index];
    p->pr = mrk->p_r[index];
    p->pphi = mrk->p_phi[index];
    p->pz = mrk->p_z[index];
    p->charge = mrk->charge[index];
    p->weight = mrk->weight[index];
    p->time = mrk->time[index];
    p->theta = mrk->theta[index];
    p->id = mrk->id[index];
    p->endcond = mrk->endcond[index];
    p->walltile = mrk->walltile[index];
    p->cputime = mrk->cputime[index];
    p->mileage = mrk->mileage[index];

    /* Particle to guiding center */
    real b_db[15], psi[1], rho[2], ppar, mu;
    rho[0] = mrk->rho[index];
    b_db[0] = mrk->br[index];
    b_db[1] = mrk->dbrdr[index];
    b_db[2] = mrk->dbrdphi[index];
    b_db[3] = mrk->dbrdz[index];
    b_db[4] = mrk->bphi[index];
    b_db[5] = mrk->dbphidr[index];
    b_db[6] = mrk->dbphidphi[index];
    b_db[7] = mrk->dbphidz[index];
    b_db[8] = mrk->bz[index];
    b_db[9] = mrk->dbzdr[index];
    b_db[10] = mrk->dbzdphi[index];
    b_db[11] = mrk->dbzdz[index];

    gctransform_particle2guidingcenter(
        p->mass, p->charge * CONST_E, b_db, p->rprt, p->phiprt, p->zprt, p->pr, p->pphi,
        p->pz, &p->r, &p->phi, &p->z, &ppar, &mu, &p->zeta);

    if (!err)
        err = Bfield_eval_b_db(b_db, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);

    real Bnorm = math_normc(b_db[0], b_db[1], b_db[2]);
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
    allocate_field(mrk->charge);
    allocate_field(mrk->time);
    allocate_field(mrk->br);
    allocate_field(mrk->bphi);
    allocate_field(mrk->bz);
    allocate_field(mrk->dbrdr);
    allocate_field(mrk->dbphidr);
    allocate_field(mrk->dbzdr);
    allocate_field(mrk->dbrdphi);
    allocate_field(mrk->dbphidphi);
    allocate_field(mrk->dbzdphi);
    allocate_field(mrk->dbrdz);
    allocate_field(mrk->dbphidz);
    allocate_field(mrk->dbzdz);
    allocate_field(mrk->bounces);
    allocate_field(mrk->weight);
    allocate_field(mrk->cputime);
    allocate_field(mrk->rho);
    allocate_field(mrk->theta);
    allocate_field(mrk->id);
    allocate_field(mrk->endcond);
    allocate_field(mrk->walltile);
    allocate_field(mrk->mileage);
    allocate_field(mrk->running);
    allocate_field(mrk->err);
    allocate_field(mrk->index);
    mrk->size = vector_size;
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
    free(mrk->charge);
    free(mrk->time);
    free(mrk->br);
    free(mrk->bphi);
    free(mrk->bz);
    free(mrk->dbrdr);
    free(mrk->dbphidr);
    free(mrk->dbzdr);
    free(mrk->dbrdphi);
    free(mrk->dbphidphi);
    free(mrk->dbzdphi);
    free(mrk->dbrdz);
    free(mrk->dbphidz);
    free(mrk->dbzdz);
    free(mrk->bounces);
    free(mrk->weight);
    free(mrk->cputime);
    free(mrk->rho);
    free(mrk->theta);
    free(mrk->id);
    free(mrk->endcond);
    free(mrk->walltile);
    free(mrk->mileage);
    free(mrk->running);
    free(mrk->err);
    free(mrk->index);
    mrk->size = 0;
}

void MarkerGuidingCenter_offload(MarkerGuidingCenter *mrk)
{
    SUPPRESS_UNUSED_WARNING(mrk);
    GPU_MAP_TO_DEVICE(
        mrk [0:1], mrk->running [0:mrk->size], mrk->r [0:mrk->size],
        mrk->phi [0:mrk->size], mrk->ppar [0:mrk->size], mrk->mu [0:mrk->size],
        mrk->zeta [0:mrk->size], mrk->mileage [0:mrk->size],
        mrk->z [0:mrk->size], mrk->charge [0:mrk->size], mrk->br [0:mrk->size],
        mrk->dbrdr [0:mrk->size], mrk->dbrdphi [0:mrk->size],
        mrk->dbrdz [0:mrk->size], mrk->bphi [0:mrk->size],
        mrk->dbphidr [0:mrk->size], mrk->dbphidphi [0:mrk->size],
        mrk->dbphidz [0:mrk->size], mrk->bz [0:mrk->size],
        mrk->dbzdr [0:mrk->size], mrk->dbzdphi [0:mrk->size],
        mrk->dbzdz [0:mrk->size], mrk->rho [0:mrk->size],
        mrk->theta [0:mrk->size], mrk->err [0:mrk->size],
        mrk->time [0:mrk->size], mrk->weight [0:mrk->size],
        mrk->cputime [0:mrk->size], mrk->id [0:mrk->size],
        mrk->endcond [0:mrk->size], mrk->walltile [0:mrk->size],
        mrk->index [0:mrk->size], mrk->bounces [0:mrk->size])
}

void MarkerGuidingCenter_onload(MarkerGuidingCenter *mrk)
{
    SUPPRESS_UNUSED_WARNING(mrk);
    GPU_MAP_FROM_DEVICE(
        mrk [0:1], mrk->running [0:mrk->size], mrk->r [0:mrk->size],
        mrk->phi [0:mrk->size], mrk->ppar [0:mrk->size], mrk->mu [0:mrk->size],
        mrk->zeta [0:mrk->size], mrk->mileage [0:mrk->size],
        mrk->z [0:mrk->size], mrk->charge [0:mrk->size], mrk->br [0:mrk->size],
        mrk->dbrdr [0:mrk->size], mrk->dbrdphi [0:mrk->size],
        mrk->dbrdz [0:mrk->size], mrk->bphi [0:mrk->size],
        mrk->dbphidr [0:mrk->size], mrk->dbphidphi [0:mrk->size],
        mrk->dbphidz [0:mrk->size], mrk->bz [0:mrk->size],
        mrk->dbzdr [0:mrk->size], mrk->dbzdphi [0:mrk->size],
        mrk->dbzdz [0:mrk->size], mrk->rho [0:mrk->size],
        mrk->theta [0:mrk->size], mrk->err [0:mrk->size],
        mrk->time [0:mrk->size], mrk->weight [0:mrk->size],
        mrk->cputime [0:mrk->size], mrk->id [0:mrk->size],
        mrk->endcond [0:mrk->size], mrk->walltile [0:mrk->size],
        mrk->index [0:mrk->size], mrk->bounces [0:mrk->size])
}

void MarkerGuidingCenter_copy(
    MarkerGuidingCenter *copy, MarkerGuidingCenter *original, size_t index)
{
    copy->znum = original->znum;
    copy->anum = original->anum;
    copy->mass = original->mass;
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
    copy->charge[index] = original->charge[index];
    copy->id[index] = original->id[index];
    copy->err[index] = original->err[index];
    copy->index[index] = original->index[index];
    copy->bounces[index] = original->bounces[index];
    copy->running[index] = original->running[index];
    copy->endcond[index] = original->endcond[index];
    copy->walltile[index] = original->walltile[index];
    copy->br[index] = original->br[index];
    copy->bphi[index] = original->bphi[index];
    copy->bz[index] = original->bz[index];
    copy->dbrdr[index] = original->dbrdr[index];
    copy->dbrdphi[index] = original->dbrdphi[index];
    copy->dbrdz[index] = original->dbrdz[index];
    copy->dbphidr[index] = original->dbphidr[index];
    copy->dbphidphi[index] = original->dbphidphi[index];
    copy->dbphidz[index] = original->dbphidz[index];
    copy->dbzdr[index] = original->dbzdr[index];
    copy->dbzdphi[index] = original->dbzdphi[index];
    copy->dbzdz[index] = original->dbzdz[index];
}

int MarkerGuidingCenter_from_queue(
    MarkerGuidingCenter *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield)
{
    State *p = queue->p[queue_index];
    err_t err = p->err;

    real b_db[15], psi[1], rho[2];
    if (!err)
        err = Bfield_eval_b_db(b_db, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);

    if (!err)
    {
        real Bnorm = math_normc(b_db[0], b_db[1], b_db[2]);
        real gamma = physlib_gamma_Ekin(p->mass, p->ekin);
        real pnorm = physlib_pnorm_gamma(p->mass, gamma);
        mrk->mass = p->mass;
        mrk->anum = p->anum;
        mrk->znum = p->znum;
        mrk->r[mrk_index] = p->r;
        mrk->phi[mrk_index] = p->phi;
        mrk->z[mrk_index] = p->z;
        mrk->zeta[mrk_index] = p->zeta;
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
        mrk->br[mrk_index] = b_db[0];
        mrk->dbrdr[mrk_index] = b_db[3];
        mrk->dbrdphi[mrk_index] = b_db[4];
        mrk->dbrdz[mrk_index] = b_db[5];
        mrk->bphi[mrk_index] = b_db[1];
        mrk->dbphidr[mrk_index] = b_db[6];
        mrk->dbphidphi[mrk_index] = b_db[7];
        mrk->dbphidz[mrk_index] = b_db[8];
        mrk->bz[mrk_index] = b_db[2];
        mrk->dbzdr[mrk_index] = b_db[9];
        mrk->dbzdphi[mrk_index] = b_db[10];
        mrk->dbzdz[mrk_index] = b_db[11];
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

    p->mass = mrk->mass;
    p->anum = mrk->anum;
    p->znum = mrk->znum;
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
    real b_db[15];
    b_db[0] = mrk->br[index];
    b_db[3] = mrk->dbrdr[index];
    b_db[4] = mrk->dbrdphi[index];
    b_db[5] = mrk->dbrdz[index];
    b_db[1] = mrk->bphi[index];
    b_db[6] = mrk->dbphidr[index];
    b_db[7] = mrk->dbphidphi[index];
    b_db[8] = mrk->dbphidz[index];
    b_db[2] = mrk->bz[index];
    b_db[9] = mrk->dbzdr[index];
    b_db[10] = mrk->dbzdphi[index];
    b_db[11] = mrk->dbzdz[index];

    real Bnorm = math_normc(mrk->br[index], mrk->bphi[index], mrk->bz[index]);
    p->ekin =
        physlib_Ekin_ppar(mrk->mass, mrk->mu[index], mrk->ppar[index], Bnorm);
    p->pitch =
        physlib_gc_xi(mrk->mass, mrk->mu[index], mrk->ppar[index], Bnorm);
    p->zeta = mrk->zeta[index];

    real pparprt, muprt, zetaprt;
    gctransform_guidingcenter2particle(
        p->mass, p->charge * CONST_E, b_db, p->r, p->phi, p->z, mrk->ppar[index],
        mrk->mu[index], p->zeta, &p->rprt, &p->phiprt, &p->zprt, &pparprt,
        &muprt, &zetaprt);

    if (!err)
        err = Bfield_eval_b_db(
            b_db, p->rprt, p->phiprt, p->zprt, p->time, bfield);

    gctransform_pparmuzeta2prpphipz(
        p->mass, p->charge * CONST_E, b_db, p->phiprt, pparprt, muprt, zetaprt, &p->pr,
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
    allocate_field(mrk->br);
    allocate_field(mrk->bphi);
    allocate_field(mrk->bz);
    allocate_field(mrk->dbrdr);
    allocate_field(mrk->dbphidr);
    allocate_field(mrk->dbzdr);
    allocate_field(mrk->dbrdphi);
    allocate_field(mrk->dbphidphi);
    allocate_field(mrk->dbzdphi);
    allocate_field(mrk->dbrdz);
    allocate_field(mrk->dbphidz);
    allocate_field(mrk->dbzdz);
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
    mrk->size = vector_size;
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
    free(mrk->br);
    free(mrk->bphi);
    free(mrk->bz);
    free(mrk->dbrdr);
    free(mrk->dbphidr);
    free(mrk->dbzdr);
    free(mrk->dbrdphi);
    free(mrk->dbphidphi);
    free(mrk->dbzdphi);
    free(mrk->dbrdz);
    free(mrk->dbphidz);
    free(mrk->dbzdz);
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
    mrk->size = 0;
}

void MarkerFieldLine_offload(MarkerFieldLine *mrk)
{
    SUPPRESS_UNUSED_WARNING(mrk);
    GPU_MAP_TO_DEVICE(
        mrk [0:1], mrk->running [0:mrk->size], mrk->r [0:mrk->size],
        mrk->phi [0:mrk->size], mrk->z [0:mrk->size],
        mrk->mileage [0:mrk->size], mrk->br [0:mrk->size],
        mrk->dbrdr [0:mrk->size], mrk->dbrdphi [0:mrk->size],
        mrk->dbrdz [0:mrk->size], mrk->bphi [0:mrk->size],
        mrk->dbphidr [0:mrk->size], mrk->dbphidphi [0:mrk->size],
        mrk->dbphidz [0:mrk->size], mrk->bz [0:mrk->size],
        mrk->dbzdr [0:mrk->size], mrk->dbzdphi [0:mrk->size],
        mrk->dbzdz [0:mrk->size], mrk->rho [0:mrk->size],
        mrk->theta [0:mrk->size], mrk->err [0:mrk->size],
        mrk->time [0:mrk->size], mrk->cputime [0:mrk->size],
        mrk->id [0:mrk->size], mrk->endcond [0:mrk->size],
        mrk->walltile [0:mrk->size], mrk->index [0:mrk->size],
        mrk->bounces [0:mrk->size])
}

void MarkerFieldLine_onload(MarkerFieldLine *mrk)
{
    SUPPRESS_UNUSED_WARNING(mrk);
    GPU_MAP_FROM_DEVICE(
        mrk [0:1], mrk->running [0:mrk->size], mrk->r [0:mrk->size],
        mrk->phi [0:mrk->size], mrk->z [0:mrk->size],
        mrk->mileage [0:mrk->size], mrk->br [0:mrk->size],
        mrk->dbrdr [0:mrk->size], mrk->dbrdphi [0:mrk->size],
        mrk->dbrdz [0:mrk->size], mrk->bphi [0:mrk->size],
        mrk->dbphidr [0:mrk->size], mrk->dbphidphi [0:mrk->size],
        mrk->dbphidz [0:mrk->size], mrk->bz [0:mrk->size],
        mrk->dbzdr [0:mrk->size], mrk->dbzdphi [0:mrk->size],
        mrk->dbzdz [0:mrk->size], mrk->rho [0:mrk->size],
        mrk->theta [0:mrk->size], mrk->err [0:mrk->size],
        mrk->time [0:mrk->size], mrk->cputime [0:mrk->size],
        mrk->id [0:mrk->size], mrk->endcond [0:mrk->size],
        mrk->walltile [0:mrk->size], mrk->index [0:mrk->size],
        mrk->bounces [0:mrk->size])
}

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
    copy->br[index] = original->br[index];
    copy->bphi[index] = original->bphi[index];
    copy->bz[index] = original->bz[index];
    copy->dbrdr[index] = original->dbrdr[index];
    copy->dbrdphi[index] = original->dbrdphi[index];
    copy->dbrdz[index] = original->dbrdz[index];
    copy->dbphidr[index] = original->dbphidr[index];
    copy->dbphidphi[index] = original->dbphidphi[index];
    copy->dbphidz[index] = original->dbphidz[index];
    copy->dbzdr[index] = original->dbzdr[index];
    copy->dbzdphi[index] = original->dbzdphi[index];
    copy->dbzdz[index] = original->dbzdz[index];
    copy->err[index] = original->err[index];
    copy->index[index] = original->index[index];
}

int MarkerFieldLine_from_queue(
    MarkerFieldLine *mrk, MarkerQueue *queue, size_t mrk_index,
    size_t queue_index, Bfield *bfield)
{
    State *p = queue->p[queue_index];
    err_t err = p->err;

    real b_db[15], psi[1], rho[2];
    if (!err)
        err = Bfield_eval_b_db(b_db, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_psi(psi, p->r, p->phi, p->z, p->time, bfield);
    if (!err)
        err = Bfield_eval_rho(rho, psi[0], bfield);

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
        mrk->br[mrk_index] = b_db[0];
        mrk->dbrdr[mrk_index] = b_db[3];
        mrk->dbrdphi[mrk_index] = b_db[4];
        mrk->dbrdz[mrk_index] = b_db[5];
        mrk->bphi[mrk_index] = b_db[1];
        mrk->dbphidr[mrk_index] = b_db[6];
        mrk->dbphidphi[mrk_index] = b_db[7];
        mrk->dbphidz[mrk_index] = b_db[8];
        mrk->bz[mrk_index] = b_db[2];
        mrk->dbzdr[mrk_index] = b_db[9];
        mrk->dbzdphi[mrk_index] = b_db[10];
        mrk->dbzdz[mrk_index] = b_db[11];
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
    p->anum = 0;
    p->znum = 0;
    p->time = mrk->time[index];
    p->id = mrk->id[index];
    p->cputime = mrk->cputime[index];
    p->theta = mrk->theta[index];
    p->endcond = mrk->endcond[index];
    p->walltile = mrk->walltile[index];
    p->mileage = mrk->mileage[index];
    p->err = mrk->err[index];
}

int MarkerGyroOrbit_to_MarkerGuidingCenter(
    MarkerGuidingCenter *mrk_gc, const MarkerGyroOrbit *mrk_go, size_t index,
    Bfield *bfield)
{
    err_t err = mrk_go->err[index];
    real axisrz[2];
    int simerr = 0; /* Error has already occurred */
    if (err)
    {
        simerr = 1;
    }
    mrk_gc->id[index] = mrk_go->id[index];
    mrk_gc->index[index] = mrk_go->index[index];

    real r, phi, z, ppar, mu, zeta, b_db[15];
    if (!err)
    {
        real Rprt = mrk_go->r[index];
        real phiprt = mrk_go->phi[index];
        real zprt = mrk_go->z[index];
        real pr = mrk_go->p_r[index];
        real pphi = mrk_go->p_phi[index];
        real pz = mrk_go->p_z[index];
        real charge = mrk_go->charge[index];

        mrk_gc->mass = mrk_go->mass;
        mrk_gc->anum = mrk_go->anum;
        mrk_gc->znum = mrk_go->znum;
        mrk_gc->charge[index] = mrk_go->charge[index];
        mrk_gc->weight[index] = mrk_go->weight[index];
        mrk_gc->time[index] = mrk_go->time[index];
        mrk_gc->mileage[index] = mrk_go->mileage[index];
        mrk_gc->endcond[index] = mrk_go->endcond[index];
        mrk_gc->running[index] = mrk_go->running[index];
        mrk_gc->walltile[index] = mrk_go->walltile[index];
        mrk_gc->cputime[index] = mrk_go->cputime[index];

        b_db[0] = mrk_go->br[index];
        b_db[3] = mrk_go->dbrdr[index];
        b_db[4] = mrk_go->dbrdphi[index];
        b_db[5] = mrk_go->dbrdz[index];
        b_db[1] = mrk_go->bphi[index];
        b_db[6] = mrk_go->dbphidr[index];
        b_db[7] = mrk_go->dbphidphi[index];
        b_db[8] = mrk_go->dbphidz[index];
        b_db[2] = mrk_go->bz[index];
        b_db[9] = mrk_go->dbzdr[index];
        b_db[10] = mrk_go->dbzdphi[index];
        b_db[11] = mrk_go->dbzdz[index];

        /* Guiding center transformation */
        gctransform_particle2guidingcenter(
            mrk_go->mass, charge, b_db, Rprt, phiprt, zprt, pr, pphi, pz, &r,
            &phi, &z, &ppar, &mu, &zeta);
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
        err = Bfield_eval_b_db(b_db, r, phi, z, mrk_go->time[index], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_psi(psi, r, phi, z, mrk_go->time[index], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_rho(rho, psi[0], bfield);
    }
    if (!err)
    {
        err = Bfield_eval_axis_rz(axisrz, bfield, mrk_gc->phi[index]);
    }

    if (!err)
    {
        mrk_gc->r[index] = r;
        mrk_gc->phi[index] = phi;
        mrk_gc->z[index] = z;
        mrk_gc->mu[index] = mu;
        mrk_gc->zeta[index] = zeta;
        mrk_gc->ppar[index] = ppar;
        mrk_gc->rho[index] = rho[0];

        /* Evaluate pol angle so that it is cumulative and at gc position */
        mrk_gc->theta[index] = mrk_go->theta[index];
        mrk_gc->theta[index] += atan2(
            (mrk_go->r[index] - axisrz[0]) * (mrk_gc->z[index] - axisrz[1]) -
                (mrk_go->z[index] - axisrz[1]) * (mrk_gc->r[index] - axisrz[0]),
            (mrk_go->r[index] - axisrz[0]) * (mrk_gc->r[index] - axisrz[0]) +
                (mrk_go->z[index] - axisrz[1]) *
                    (mrk_gc->z[index] - axisrz[1]));

        mrk_gc->br[index] = b_db[0];
        mrk_gc->dbrdr[index] = b_db[3];
        mrk_gc->dbrdphi[index] = b_db[4];
        mrk_gc->dbrdz[index] = b_db[5];

        mrk_gc->bphi[index] = b_db[1];
        mrk_gc->dbphidr[index] = b_db[6];
        mrk_gc->dbphidphi[index] = b_db[7];
        mrk_gc->dbphidz[index] = b_db[8];

        mrk_gc->bz[index] = b_db[2];
        mrk_gc->dbzdr[index] = b_db[9];
        mrk_gc->dbzdphi[index] = b_db[10];
        mrk_gc->dbzdz[index] = b_db[11];
    }
    if (!simerr)
    {
        err = err;
    }
    mrk_gc->err[index] = err;
    if (mrk_gc->err[index])
    {
        mrk_gc->running[index] = 0;
        mrk_gc->endcond[index] = 0;
    }

    return err > 0;
}

size_t MarkerQueue_cycle(
    size_t *next_in_queue, MarkerQueue *q, size_t nmrk, size_t start,
    size_t ids[nmrk], int running[nmrk])
{
    size_t idx;
    for (idx = start; idx < nmrk; idx++)
    {
        int marker_finished = ids[idx] > 0 && !running[idx];
        int vector_initial_fill = ids[idx] == 0 && q->next < q->n;
        if (marker_finished)
        {
#pragma omp critical
            q->finished++;
            break;
        }
        if (vector_initial_fill)
            break;
    }

    *next_in_queue = q->n;
    int vector_has_empty_slot = idx < nmrk;
    if (vector_has_empty_slot)
    {
        size_t next;
#pragma omp critical
        next = q->next++;
        int empty_queue = next >= q->n;
        if (!empty_queue)
            *next_in_queue = next;
    }
    return idx;
}

#undef allocate_field
