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
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void asigma_free(asigma_data* data) {
    switch(data->type) {
        case asigma_type_loc:
            asigma_loc_free(&data->asigma_loc);
            break;
    }
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void asigma_offload(asigma_data* data) {
    switch(data->type) {
        case asigma_type_loc:
            asigma_loc_offload(&data->asigma_loc);
            break;
    }
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
    real* ratecoeff, int z_1, int a_1, int q_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_cx(
                    ratecoeff, z_1, a_1, q_1, E, mass, nspec, znum, anum, T_0,
                    n_0, ASIGMA_EXTRAPOLATE, &(asigma_data->asigma_loc));
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


/**
 * @brief Evaluate electron impact ionization rate coefficient
 *
 * This is a SIMD function.
 *
 * @param ratecoeff pointer to evaluated rate coefficient
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param E energy of fast particle
 * @param mass mass of fast particle
 * @param Te electron temperature
 * @param Te electron density
 * @param asigma_data pointer to atomic data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err asigma_eval_eii(
    real* ratecoeff, int z_1, int a_1, int q_1, real E, real mass, real Te,
    real ne, asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_eii(
                    ratecoeff, z_1, a_1, q_1, E, mass, Te, ne,
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
