/**
 * @file mhd.c
 * @brief MHD module interface
 *
 * This is an interface through which MHD data is initialized and accessed.
 * Reading e.g. from disk is done elsewhere. The MHD module produces helical
 * EM perturbations in to the EM field using the boozer module in making the
 * coordinate transformations between cylindrical and straight-field-line
 * coordinates.
 *
 * To add a new MHD instance, make sure these functions are implemented and
 * called from this interface, and that mhd.h contains enum type for the new
 * instance.
 *
 * The interface checks which instance given data corresponds to from
 * mhd_offload_data.type and Mhd.type from the struct that is given
 * as an argument, and calls the relevant function for that instance.
 */
#include "mhd.h"
#include "bfield.h"
#include "boozer.h"
#include "defines.h"
#include "mhd_dynamic.h"
#include "mhd_stationary.h"
#include <stdlib.h>

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
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

/**
 * @brief Offload mhd to the accelerator.
 *
 * @param mhd pointer to the mhd struct
 */
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

/**
 * @brief Evaluate the needed quantities from MHD mode for orbit following
 *
 * The quantities to be evaluated are alpha, phi, grad alpha, grad phi,
 * partial t alpha, partial t phi
 *
 * The values are stored in the given array as:
 * - mhd_dmhd[0] = alpha
 * - mhd_dmhd[1] = dalpha/dt
 * - mhd_dmhd[2] = grad alpha, r component
 * - mhd_dmhd[3] = grad alpha, phi component
 * - mhd_dmhd[4] = grad alpha, z component
 * - mhd_dmhd[5] = phi
 * - mhd_dmhd[6] = dphi/dt
 * - mhd_dmhd[7] = grad phi, r component
 * - mhd_dmhd[8] = grad phi, phi component
 * - mhd_dmhd[9] = grad phi, z component
 *
 * @param mhd_dmhd
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param includemode mode number to be included or MHD_INCLUDE_ALL
 * @param boozerdata pointer to boozer data
 * @param mhd pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
err_t Mhd_eval(
    real mhd_dmhd[10], real r, real phi, real z, real t, size_t includemode,
    Boozer *boozer, Mhd *mhd, Bfield *bfield)
{
    err_t err = 0;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        err = MhdStationary_eval(
            mhd_dmhd, r, phi, z, t, includemode, boozer, mhd->stationary,
            bfield);
        break;
    case MHD_DYNAMIC:
        err = MhdDynamic_eval(
            mhd_dmhd, r, phi, z, t, includemode, boozer, mhd->dynamic, bfield);
        break;
    }

    mhd_dmhd[0] = err ? 0.0 : mhd_dmhd[0];
    mhd_dmhd[1] = err ? 0.0 : mhd_dmhd[1];
    mhd_dmhd[2] = err ? 0.0 : mhd_dmhd[2];
    mhd_dmhd[3] = err ? 0.0 : mhd_dmhd[3];
    mhd_dmhd[4] = err ? 0.0 : mhd_dmhd[4];
    mhd_dmhd[5] = err ? 0.0 : mhd_dmhd[5];
    mhd_dmhd[6] = err ? 0.0 : mhd_dmhd[6];
    mhd_dmhd[7] = err ? 0.0 : mhd_dmhd[7];
    mhd_dmhd[8] = err ? 0.0 : mhd_dmhd[8];
    mhd_dmhd[9] = err ? 0.0 : mhd_dmhd[9];

    return err;
}

/**
 * @brief Evaluate perturbed fields Btilde, Etilde and potential Phi explicitly
 *
 * The values are stored in the given array as
 * - pert_field[0] = BtildeR
 * - pert_field[1] = BtildePhi
 * - pert_field[2] = BtildeZ
 * - pert_field[3] = EtildeR
 * - pert_field[4] = EtildePhi
 * - pert_field[5] = EtildeZ
 * - pert_field[6] = Phi
 *
 * Only the perturbation values for the magnetic field are returned if
 * pertonly=1, otherwise, the total perturbed field is returned. This is done to
 * avoid double evaluation of the magnetic field e.g. in field line tracing.
 * For electric field only the perturbation component is returned always.
 *
 * @param pert_field perturbation field components
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param t time coordinate [s]
 * @param pertonly flag whether to return the whole field or only perturbation
 * @param includemode mode number to be included or MHD_INCLUDE_ALL
 * @param boozerdata pointer to boozer data
 * @param mhd pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero err_t value if evaluation failed, zero otherwise
 */
err_t Mhd_perturbations(
    real pert_field[7], real r, real phi, real z, real t, int pertonly,
    size_t includemode, Boozer *boozer, Mhd *mhd, Bfield *bfield)
{
    err_t err = 0;
    switch (mhd->type)
    {
    case MHD_STATIONARY:
        err = MhdStationary_perturbations(
            pert_field, r, phi, z, t, pertonly, includemode, boozer,
            mhd->stationary, bfield);
        break;
    case MHD_DYNAMIC:
        err = MhdDynamic_perturbations(
            pert_field, r, phi, z, t, pertonly, includemode, boozer,
            mhd->dynamic, bfield);
        break;
    }

    pert_field[0] = err ? 0.0 : pert_field[0];
    pert_field[1] = err ? 0.0 : pert_field[1];
    pert_field[2] = err ? 0.0 : pert_field[2];
    pert_field[3] = err ? 0.0 : pert_field[3];
    pert_field[4] = err ? 0.0 : pert_field[4];
    pert_field[5] = err ? 0.0 : pert_field[5];
    pert_field[6] = err ? 0.0 : pert_field[6];

    return err;
}

/**
 * Get number of modes.
 *
 * @param mhd MHD data.
 *
 * @return Number of modes.
 */
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

/**
 * Get mode toroidal numbers.
 *
 * @param mhd MHD data.
 *
 * @return Mode toroidal numbers.
 */
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

/**
 * Get mode poloidal numbers.
 *
 * @param mhd MHD data.
 *
 * @return Mode poloidal numbers.
 */
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

/**
 * Get mode amplitudes.
 *
 * @param mhd MHD data.
 *
 * @return Mode amplitudes [1].
 */
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

/**
 * Get mode frequencies.
 *
 * @param mhd MHD data.
 *
 * @return Mode frequencies [rad/s].
 */
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

/**
 * Get mode phases.
 *
 * @param mhd MHD data.
 *
 * @return Mode phases [rad].
 */
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
