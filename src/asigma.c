/**
 * @file asigma.c
 * @brief Atomic reaction data interface
 *
 * This is an interface through which atomic reaction data is initialized
 * and accessed. Reading, for example from local files, is done elsewhere.
 *
 * The name asigma is short for atomicsigma. The word sigma refers to
 * cross-section, a fundamental type of reaction probability data. Note
 * that the data is not necessarily in the form of pure cross-sections.
 * It might be in some derivative form, such as rate coefficients.
 *
 * To add a new atomic reaction data instance, make sure the functions
 * are implemented and called from this interface, and that asigma.h
 * contains enum type for the new instance.
 *
 * The interface checks which instance given data corresponds to from the
 * "type"-field in asigma_offload_data or asigma_data that is given
 * as an argument, and calls the relevant function for that instance.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "print.h"
#include "error.h"
#include "math.h"
#include "asigma.h"
#include "asigma/asigma_loc.h"
#include "consts.h"

#pragma omp declare target
/** Set values outside abscissae to zero instead of raising an error. */
static int ASIGMA_EXTRAPOLATE = 0;
#pragma omp end declare target

/**
 * @brief Toggle extrapolation when evaluating cross sections.
 *
 * In this context the extrapolation means values outside the abscissae
 * are set to zero instead of raising an error.
 *
 * @param extrapolate flag whether to extrapolate
 */
void asigma_extrapolate(int extrapolate) {
    ASIGMA_EXTRAPOLATE = extrapolate;
}

/**
 * @brief Load atomic reaction data and prepare parameters
 *
 * This function fills the relevant atomic sigma offload struct with parameters
 * and allocates and fills the offload array. Sets offload array length in the
 * offload struct.
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
int asigma_init_offload(asigma_offload_data* offload_data,
                        real** offload_array) {
    int err = 0;

    switch(offload_data->type) {

        case asigma_type_loc:
            err = asigma_loc_init_offload(&(offload_data->asigma_loc),
                                          offload_array);
            offload_data->offload_array_length =
                offload_data->asigma_loc.offload_array_length;
            break;

        default:
            /* Unrecognized input. Produce error. */
            print_err("Error: Unrecognized atomic reaction data type.");
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
void asigma_free_offload(asigma_offload_data* offload_data,
                         real** offload_array) {
    switch(offload_data->type) {
        case asigma_type_loc:
            asigma_loc_free_offload(&(offload_data->asigma_loc), offload_array);
            break;
    }
}

/**
 * @brief Initializes atomic reaction data struct on target
 *
 * This function copies some atomic reaction parameters from the offload
 * struct to the struct on target and sets the atomic reaction data pointers
 * to correct offsets in the offload array.
 *
 * @param asigma_data pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @return zero if initialization succeeded
 */
int asigma_init(asigma_data* asigma_data, asigma_offload_data* offload_data,
                real* offload_array) {
    int err = 0;
    switch(offload_data->type) {
        case asigma_type_loc:
            asigma_loc_init(&(asigma_data->asigma_loc),
                            &(offload_data->asigma_loc), offload_array);
            break;

        default:
            /* Unrecognized input. Produce error. */
            print_err("Error: Unrecognized atomic reaction data type.");
            err = 1;
            break;
    }
    asigma_data->type = offload_data->type;

    return err;
}

/**
 * @brief Evaluate atomic reaction cross-section
 *
 * This function evaluates the cross-section (sigma) for the atomic reaction
 * corresponding to the reaction identifiers given as parameters at the
 * given mass-normalized collision energy.
 *
 * This is a SIMD function.
 *
 * @param sigma pointer to evaluated cross-section
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param z_2 atomic number of bulk particle
 * @param a_2 atomic mass number of bulk particle
 * @param E_coll_per_amu energy per amu corresponding to collision speed
 * @param reac_type reaction type
 * @param asigma_data pointer to atomic data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err asigma_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    asigma_reac_type reac_type, asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_sigma(
                sigma, z_1, a_1, z_2, a_2, E_coll_per_amu, reac_type,
                ASIGMA_EXTRAPOLATE, &(asigma_data->asigma_loc));
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || sigma[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        sigma[0] = 0.0;
    }

    return err;
}

/**
 * @brief Evaluate atomic reaction rate coefficient
 *
 * This function evaluates the rate coefficient (<sigma*v>) for the atomic
 * reaction corresponding to the reaction identifiers given as parameters
 * at the given fast particle energy and bulk plasma conditions.
 *
 * This is a SIMD function.
 *
 * @param sigmav pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param m_1 mass of fast particle
 * @param z_2 atomic number of bulk particle
 * @param a_2 atomic mass number of bulk particle
 * @param E energy of fast particle
 * @param T_e electron temperature of bulk plasma
 * @param T_0 temperature of bulk neutrals
 * @param n_i ion density of bulk plasma
 * @param reac_type reaction type
 * @param asigma_data pointer to atomic data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err asigma_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, asigma_reac_type reac_type,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_sigmav(
                sigmav, z_1, a_1, m_1, z_2, a_2, E, T_e, T_0, n_i,
                reac_type, ASIGMA_EXTRAPOLATE, &(asigma_data->asigma_loc));
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || sigmav[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        sigmav[0] = 0.0;
    }

    return err;
}

/**
 * @brief Evaluate charge exchange rate coefficient
 *
 * This is a SIMD function.
 *
 * @param ratecoeff pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param E energy of fast particle
 * @param mass mass of fast particle
 * @param nspec number of bulk neutral species
 * @param znum atomic numbers of bulk particles
 * @param anum atomic mass numbers of bulk particles
 * @param T_0 temperature of bulk neutrals
 * @param n_0 densities of bulk neutrals
 * @param asigma_data pointer to atomic data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err asigma_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_cx(
                    ratecoeff, z_1, a_1, E, mass, nspec, znum, anum, T_0, n_0,
                    ASIGMA_EXTRAPOLATE, &(asigma_data->asigma_loc));
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || ratecoeff[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        ratecoeff[0] = 0.0;
    }

    return err;
}

/**
 * @brief Evaluate beam stopping rate coefficient
 *
 * This is a SIMD function.
 *
 * @param ratecoeff pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param E energy of fast particle
 * @param mass mass of fast particle
 * @param nion number of bulk ion species
 * @param znum atomic numbers of bulk particles
 * @param anum atomic mass numbers of bulk particles
 * @param T_e electron temperature of bulk plasma
 * @param n_i densities of bulk ions
 * @param asigma_data pointer to atomic data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err asigma_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_bms(
                    ratecoeff, z_1, a_1, E, mass, nion, znum, anum, T_e, n_i,
                    ASIGMA_EXTRAPOLATE, &(asigma_data->asigma_loc));
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || ratecoeff[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        ratecoeff[0] = 0.0;
    }

    return err;
}
