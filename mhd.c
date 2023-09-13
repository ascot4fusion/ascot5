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
 * mhd_offload_data.type and mhd_data.type from the struct that is given
 * as an argument, and calls the relevant function for that instance.
 */
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "mhd.h"
#include "B_field.h"
#include "boozer.h"
#include "mhd/mhd_stat.h"
#include "mhd/mhd_nonstat.h"

/**
 * @brief Load MHD data and prepare parameters
 *
 * This function fills the relevant MHD offload struct with
 * parameters and allocates and fills the offload array.
 *
 * The offload data has to have a type when this function is called as it should
 * be set when the offload data is constructed from inputs.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded
 */
int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array) {
    int err = 0;

    switch(offload_data->type) {

        case mhd_type_stat:
            err = mhd_stat_init_offload(&(offload_data->stat), offload_array);
            offload_data->offload_array_length =
                offload_data->stat.offload_array_length;
            break;

        case mhd_type_nonstat:
            err = mhd_nonstat_init_offload(&(offload_data->nonstat),
                                           offload_array);
            offload_data->offload_array_length =
                offload_data->nonstat.offload_array_length;
            break;

        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized MHD type.");
            err = 1;
            break;
    }

    if(!err) {
        print_out(VERBOSE_IO, "Estimated memory usage %.1f MB\n",
                  offload_data->offload_array_length
                  * sizeof(real) / (1024.0*1024.0) );
    }

    return err;
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * This function is host only.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void mhd_free_offload(mhd_offload_data* offload_data,
                      real** offload_array) {
    switch(offload_data->type) {

        case mhd_type_stat:
            mhd_stat_free_offload(&(offload_data->stat), offload_array);
            break;

        case mhd_type_nonstat:
            mhd_nonstat_free_offload(&(offload_data->nonstat), offload_array);
            break;
    }
}

/**
 * @brief Initialize MHD data struct on target
 *
 * This function copies the electric field parameters from the offload struct
 * to the struct on target and sets the MHD data pointers to correct offsets
 * in the offload array.
 *
 * This function returns error if the offload data has not been initialized.
 * The instances themselves should not return an error since all they do is
 * assign pointers and values.
 *
 * @param mhddata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array the offload array
 *
 * @return Non-zero integer if offload was not initialized beforehand
 */
int mhd_init(mhd_data* mhddata, mhd_offload_data* offload_data,
             real* offload_array) {
    int err = 0;

    switch(offload_data->type) {

        case mhd_type_stat:
            mhd_stat_init(&(mhddata->stat), &(offload_data->stat),
                          offload_array);
            break;

        case mhd_type_nonstat:
            mhd_nonstat_init(&(mhddata->nonstat), &(offload_data->nonstat),
                             offload_array);
            break;

        default:
            /* Unregonized input. Produce error. */
            print_err("Error: Unregonized electric field type.\n");
            err = 1;
            break;
    }
    mhddata->type = offload_data->type;

    return err;
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
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_eval(real mhd_dmhd[10], real r, real phi, real z, real t,
               int includemode, boozer_data* boozerdata, mhd_data* mhddata,
               B_field_data* Bdata) {
    a5err err = 0;

    switch(mhddata->type) {

        case mhd_type_stat:
            err = mhd_stat_eval(mhd_dmhd, r, phi, z, t, includemode,
                                boozerdata, &(mhddata->stat), Bdata);
            break;

        case mhd_type_nonstat:
            err = mhd_nonstat_eval(mhd_dmhd, r, phi, z, t, includemode,
                                   boozerdata, &(mhddata->nonstat), Bdata);
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_MHD );
            break;
    }

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
 * @param boozerdata pointer to boozer data
 * @param mhddata pointer to mhd data
 * @param Bdata pointer to magnetic field data
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err mhd_perturbations(real pert_field[7], real r, real phi, real z,
                        real t, int pertonly, int includemode,
                        boozer_data* boozerdata, mhd_data* mhddata,
                        B_field_data* Bdata) {
    a5err err = 0;

    switch(mhddata->type) {

        case mhd_type_stat:
            err =  mhd_stat_perturbations(pert_field, r, phi, z, t, pertonly,
                                          includemode, boozerdata,
                                          &(mhddata->stat), Bdata);
            break;

        case mhd_type_nonstat:
            err =  mhd_nonstat_perturbations(pert_field, r, phi, z, t, pertonly,
                                             includemode, boozerdata,
                                             &(mhddata->nonstat), Bdata);
            break;

        default:
            /* Unregonized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_MHD );
            break;
    }

    return err;
}

/**
 * @brief Return number of modes.
 *
 * @param mhddata pointer to mhd data
 *
 * @return number of modes
 */
int mhd_get_n_modes(mhd_data* mhddata) {
    int val = 0;
    switch(mhddata->type) {
        case mhd_type_stat:
            val = mhddata->stat.n_modes;
            break;
        case mhd_type_nonstat:
            val = mhddata->nonstat.n_modes;
            break;
    }
    return val;
}

/**
 * @brief Return mode toroidal numbers.
 *
 * @param mhddata pointer to mhd data
 *
 * @return mode n numbers
 */
const int* mhd_get_nmode(mhd_data* mhddata) {
    const int* val = NULL;
    switch(mhddata->type) {
        case mhd_type_stat:
            val = mhddata->stat.nmode;
            break;
        case mhd_type_nonstat:
            val = mhddata->nonstat.nmode;
            break;
    }
    return val;
}

/**
 * @brief Return mode poloidal numbers.
 *
 * @param mhddata pointer to mhd data
 *
 * @return mode m numbers
 */
const int* mhd_get_mmode(mhd_data* mhddata) {
    const int* val = NULL;
    switch(mhddata->type) {
        case mhd_type_stat:
            val = mhddata->stat.mmode;
            break;
        case mhd_type_nonstat:
            val = mhddata->nonstat.mmode;
            break;
    }
    return val;
}

/**
 * @brief Return mode amplitudes.
 *
 * @param mhddata pointer to mhd data
 *
 * @return mode amplitudes
 */
const real* mhd_get_amplitude(mhd_data* mhddata) {
    const real* val = NULL;
    switch(mhddata->type) {
        case mhd_type_stat:
            val = mhddata->stat.amplitude_nm;
            break;
        case mhd_type_nonstat:
            val = mhddata->nonstat.amplitude_nm;
            break;
    }
    return val;
}

/**
 * @brief Return mode frequencies.
 *
 * @param mhddata pointer to mhd data
 *
 * @return mode omega
 */
const real* mhd_get_frequency(mhd_data* mhddata) {
    const real* val = NULL;
    switch(mhddata->type) {
        case mhd_type_stat:
            val = mhddata->stat.omega_nm;
            break;
        case mhd_type_nonstat:
            val = mhddata->nonstat.omega_nm;
            break;
    }
    return val;
}

/**
 * @brief Return mode phases.
 *
 * @param mhddata pointer to mhd data
 *
 * @return mode phases
 */
const real* mhd_get_phase(mhd_data* mhddata) {
    const real* val = NULL;
    switch(mhddata->type) {
        case mhd_type_stat:
            val = mhddata->stat.phase_nm;
            break;
        case mhd_type_nonstat:
            val = mhddata->nonstat.phase_nm;
            break;
    }
    return val;
}
