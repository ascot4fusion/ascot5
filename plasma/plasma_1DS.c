/**
 * @file plasma_1DS.c
 * @brief 1D spline plasma evaluation functions
 */
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "../ascot5.h"
#include "../print.h"
#include "../error.h"
#include "../consts.h"
#include "../spline/interp.h"
#include "plasma_1DS.h"

/**
 * @brief Flag to determine method to prevent negative values in interpolation
 *
 * It is possible that interpolating function f(x) with splines yield negative
 * values which would be unphysical if f is density or temperature. To prevent
 * this, one can either interpolate sqrt(f(x)) or log(f(x)) instead because
 * inverse of these functions are always positive.
 *
 * Set this parameter to
 *  - 0 : Do nothing
 *  - 1 : Take logarithm before interpolating
 *  - 2 : Take square root before interpolating
 */
#define PLASMA_1DS_NONEG 1

/** Logarithm flag */
#define PLASMA_1DS_LOG 1

/** Square root flag */
#define PLASMA_1DS_SQRT 2

/**
 * @brief Initialize 1DS plasma data and check inputs
 *
 * Before calling this function, the offload struct is expected to be fully
 * initialized.
 *
 * The offload array is expected to hold plasma data as
 *   &(*offload_array)[n_rho*0] = electron temperature
 *   &(*offload_array)[n_rho*1] = ion temperature
 *   &(*offload_array)[n_rho*2] = electron density
 *   &(*offload_array)[n_rho*3] = ion density
 *
 * This function initializes splines to plasma profiles and prints some values
 * as sanity checks.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded
 */
int plasma_1DS_init_offload(plasma_1DS_offload_data* offload_data,
                            real** offload_array) {

    /* Spline initialization */
    int err = 0;
    int n_rho     = offload_data->n_rho;
    int n_species = offload_data->n_species;

    /* Allocate enough space for two temperature and n_species density arrays */
    real* coeff_array = (real*) malloc((2 + n_species) * NSIZE_COMP1D
                                       * n_rho * sizeof(real));

    for(int i=0; i< ((2 + n_species)*n_rho); i++) {
        if(PLASMA_1DS_NONEG == PLASMA_1DS_LOG) {
            (*offload_array)[i] = log( (*offload_array)[i] + DBL_EPSILON);
        }
        if(PLASMA_1DS_NONEG == PLASMA_1DS_SQRT) {
            (*offload_array)[i] = sqrt( (*offload_array)[i] );
        }
    }

    /* Evaluate spline coefficients */

    /* Te */
    err += interp1Dcomp_init_coeff(
        coeff_array + 0*n_rho*NSIZE_COMP1D,
        *offload_array + 0*n_rho,
        offload_data->n_rho, NATURALBC,
        offload_data->rho_min, offload_data->rho_max);

    /* Ti */
    err += interp1Dcomp_init_coeff(
        coeff_array + 1*n_rho*NSIZE_COMP1D,
        *offload_array + 1*n_rho,
        offload_data->n_rho, NATURALBC,
        offload_data->rho_min, offload_data->rho_max);

    /* Densities */
    for(int i=0; i < n_species; i++) {
        err += interp1Dcomp_init_coeff(
            coeff_array + (2 +i)*n_rho*NSIZE_COMP1D,
            *offload_array + (2 + i)*n_rho,
            offload_data->n_rho, NATURALBC,
            offload_data->rho_min, offload_data->rho_max);
    }

    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = (2 + n_species) * NSIZE_COMP1D * n_rho;

    if(err) {
        return err;
    }

    print_out(VERBOSE_IO,
              "\n1D plasma profiles (P_1D)\n"
              " Min/Max rho               = %1.2le / %1.2le\n"
              " Number of rho grid points = %d\n"
              " Number of ion species     = %d\n",
              offload_data->rho_min, offload_data->rho_max,
              n_rho,
              n_species-1);

    return 0;
}

/**
 * @brief Free offload array and reset parameters
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void plasma_1DS_free_offload(plasma_1DS_offload_data* offload_data,
                             real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize magnetic field data struct on target
 *
 * This function copies the magnetic field parameters from the offload struct
 * to the struct on target and sets the plasma data pointers to
 * correct offsets in the offload array.
 *
 * @param plasma_data pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void plasma_1DS_init(plasma_1DS_data* plasma_data,
                     plasma_1DS_offload_data* offload_data,
                     real* offload_array) {
    plasma_data->n_species = offload_data->n_species;

    for(int i = 0; i < plasma_data->n_species; i++) {
        plasma_data->mass[i]   = offload_data->mass[i];
        plasma_data->charge[i] = offload_data->charge[i];
        plasma_data->znum[i]   = offload_data->znum[i];
        plasma_data->anum[i]   = offload_data->anum[i];
    }

    int n_rho = offload_data->n_rho;
    interp1Dcomp_init_spline(&(plasma_data->temp[0]),
                             &(offload_array[0*n_rho]),
                             offload_data->n_rho, NATURALBC,
                             offload_data->rho_min,
                             offload_data->rho_max);
    interp1Dcomp_init_spline(&(plasma_data->temp[1]),
                             &(offload_array[1*n_rho*NSIZE_COMP1D]),
                             offload_data->n_rho, NATURALBC,
                             offload_data->rho_min,
                             offload_data->rho_max);

    for(int i=0; i<offload_data->n_species; i++) {
        interp1Dcomp_init_spline(&(plasma_data->dens[i]),
                                 &(offload_array[(2+i)*n_rho*NSIZE_COMP1D]),
                                 offload_data->n_rho, NATURALBC,
                                 offload_data->rho_min,
                                 offload_data->rho_max);
    }
}

/**
 * @brief Evaluate plasma temperature
 *
 * @param temp temperature value will be stored in temp[0]
 * @param rho radial coordinate
 * @param species index of plasma species
 * @param plasma_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_1DS_eval_temp(real* temp, real rho, int species,
                           plasma_1DS_data* plasma_data) {
    int interperr = 0;
    interperr += interp1Dcomp_eval_f(temp, &plasma_data->temp[species>0], rho);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1DS );
    }

    if(PLASMA_1DS_NONEG == PLASMA_1DS_LOG) {
        *temp = exp(*temp);
    }
    else if(PLASMA_1DS_NONEG == PLASMA_1DS_SQRT) {
        *temp = (*temp) * (*temp);
    }
    else if(!err && *temp < 0){
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1DS );
    }

    return err;
}

/**
 * @brief Evaluate plasma density
 *
 * @param dens density value will be stored in dens[0]
 * @param rho radial coordinate
 * @param species index of plasma species
 * @param plasma_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_1DS_eval_dens(real* dens, real rho, int species,
                           plasma_1DS_data* plasma_data) {

    int interperr = 0;
    interperr += interp1Dcomp_eval_f(dens, &plasma_data->dens[species], rho);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1DS );
    }

    if(PLASMA_1DS_NONEG == PLASMA_1DS_LOG) {
        *dens = exp(*dens);
    }
    else if(PLASMA_1DS_NONEG == PLASMA_1DS_SQRT) {
        *dens = (*dens) * (*dens);
    }
    else if(!err && *dens < 0){
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1DS );
    }

    return err;
}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 * This function evaluates the density and temperature of all plasma species at
 * the given radial coordinate using spline interpolation.
 *
 * @param dens pointer to where interpolated densities [m^-3] are stored
 * @param temp pointer to where interpolated temperatures [J] are stored
 * @param rho radial coordinate
 * @param plasma_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_1DS_eval_densandtemp(real* dens, real* temp, real rho,
                                  plasma_1DS_data* plasma_data) {
    int interperr = 0;

    /* Evaluate electron temperature and density */
    interperr += interp1Dcomp_eval_f(&temp[0], &plasma_data->temp[0], rho);
    interperr += interp1Dcomp_eval_f(&dens[0], &plasma_data->dens[0], rho);

    /* Evaluate ion temperature (same for all ions) and densities */
    interperr += interp1Dcomp_eval_f(&temp[1], &plasma_data->temp[1], rho);
    for(int i=1; i<plasma_data->n_species; i++) {
        temp[i]    = temp[1];
        interperr += interp1Dcomp_eval_f(&dens[i], &plasma_data->dens[i], rho);
    }

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1DS );
    }

    if(PLASMA_1DS_NONEG == PLASMA_1DS_LOG) {
        for(int i=0; i<plasma_data->n_species; i++) {
            dens[i] = exp(dens[i]);
            temp[i] = exp(temp[i]);
        }
    }
    else if(PLASMA_1DS_NONEG == PLASMA_1DS_SQRT) {
        for(int i=0; i<plasma_data->n_species; i++) {
            dens[i] = dens[i]*dens[i];
            temp[i] = temp[i]*temp[i];
        }
    }
    else {
        for(int i=0; i<plasma_data->n_species; i++) {
            if(!err && (dens[i] < 0 || temp[i] < 0) ) {
                err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                   EF_PLASMA_1DS );
            }
        }
    }

    return err;
}

/**
 * @brief Return number of plasma species
 *
 * @param pls_data pointer to plasma data
 *
 * @return number of plasma species
 */
int plasma_1DS_get_n_species(plasma_1DS_data* pls_data) {
    return pls_data->n_species;
}

/**
 * @brief Return pointer to array storing species mass
 *
 * @param pls_data pointer to plasma data
 *
 * @return pointer to immutable MAX_SPECIES length array containing masses
 */
const real* plasma_1DS_get_species_mass(plasma_1DS_data* pls_data) {
    return pls_data->mass;
}

/**
 * @brief Return pointer to array storing species charge
 *
 * @param pls_data pointer to plasma data
 *
 * @return pointer to immutable MAX_SPECIES length array containing charges
 */
const real* plasma_1DS_get_species_charge(plasma_1DS_data* pls_data) {
    return pls_data->charge;
}

/**
 * @brief Return pointer to array storing species atomic number
 *
 * @param pls_data pointer to plasma data
 *
 * @return pointer to immutable MAX_SPECIES length array containing
 *         atomic numbers
 */
const int* plasma_1DS_get_species_znum(plasma_1DS_data* pls_data) {
    return pls_data->znum;
}

/**
 * @brief Return pointer to array storing species mass number
 *
 * @param pls_data pointer to plasma data
 *
 * @return pointer to immutable MAX_SPECIES length array containing
 *         mass numbers
 */
const int* plasma_1DS_get_species_anum(plasma_1DS_data* pls_data) {
    return pls_data->anum;
}
