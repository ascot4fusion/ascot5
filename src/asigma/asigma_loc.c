/**
 * @file asigma_loc.c
 * @brief Atomic reaction data from local files
 *
 * Atomic reaction data (sigmas) originating from local files and
 * interpolated using splines. If the data for a reaction is missing,
 * a corresponding analytical model implemented in ASCOT5 might be used.
 */
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../print.h"
#include "../error.h"
#include "asigma_loc.h"
#include "../asigma.h"
#include "../spline/interp.h"
#include "../consts.h"
#include "../math.h"

/**
 * @brief Initialize local file atomic data and check inputs
 *
 * Before calling this function, the offload struct is expected to be fully
 * initialized.
 *
 * The offload array is expected to hold atomic data as
 *
 *   [0*N_reac] = min value of energy abscissa [eV]
 *   [1*N_reac] = max value of energy abscissa [eV]
 *   [2*N_reac] = min value of density abscissa [m^-3]
 *   [3*N_reac] = max value of density abscissa [m^-3]
 *   [4*N_reac] = min value of temperature abscissa [eV]
 *   [5*N_reac] = max value of temperature abscissa [eV]
 *   [6*N_reac] = reaction probability data [(depends on reaction data type)]
 *
 * Each piece of data listed above is repeated for each reaction included.
 * Hence the N_reac interval between the 6 first pieces of data. The memory
 * space required for the last, reaction probability data, depends on the
 * dimensionality of the abscissae. The memory requirement for the reaction
 * data of the i_reac:th reaction is N_E[i_reac]*N_n[i_reac]*N_T[i_reac].
 *
 * This function initializes splines to atomic reaction data, thus increasing
 * the length of the offload array, and prints some values as sanity checks.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization success
 */
int asigma_loc_init_offload(asigma_loc_offload_data* offload_data,
                            real** offload_array) {
    int err = 0;
    int N_reac = offload_data->N_reac;

    /* Find how much space is needed for the output array */
    int temp_arr_length = 6*N_reac;
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        int N_E = offload_data->N_E[i_reac];
        int N_n = offload_data->N_n[i_reac];
        int N_T = offload_data->N_T[i_reac];
        int dim = (N_E > 1) + (N_n > 1) + (N_T > 1);
        switch(dim) {
            case 1:
                temp_arr_length += N_E * NSIZE_COMP1D;
                break;

            case 2:
                temp_arr_length += N_E * N_T * NSIZE_COMP2D;
                break;

            case 3:
                temp_arr_length += N_E * N_n * N_T * NSIZE_COMP3D;
                break;

            default:
                /* Unrecognized dimensionality. Produce error. */
                print_err("Error: Unrecognized abscissa dimensionality\n");
                err = 1;
                break;
        }
    }
    real* temp_array = (real*) malloc(temp_arr_length*sizeof(real));

    /* Helper pointers to keep track of the current memory positions in the
       reaction data parts of the arrays */
    real * temp_arr_pos    = temp_array + 5 * N_reac;
    real * offload_arr_pos = *offload_array + 5 * N_reac;

    /* Copy over reaction identifiers and abscissae parameters, and evaluate
       spline coefficients according to dimensionality of reaction data */
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        int N_E = offload_data->N_E[i_reac];
        int N_n = offload_data->N_n[i_reac];
        int N_T = offload_data->N_T[i_reac];
        int dim = (N_E > 1) + (N_n > 1) + (N_T > 1);

        real E_min = (*offload_array)[0*N_reac+i_reac];
        real E_max = (*offload_array)[1*N_reac+i_reac];
        real n_min = (*offload_array)[2*N_reac+i_reac];
        real n_max = (*offload_array)[3*N_reac+i_reac];
        real T_min = (*offload_array)[4*N_reac+i_reac];
        real T_max = (*offload_array)[5*N_reac+i_reac];
        temp_array[ 0*N_reac+i_reac] = E_min;
        temp_array[ 1*N_reac+i_reac] = E_max;
        temp_array[ 2*N_reac+i_reac] = n_min;
        temp_array[ 3*N_reac+i_reac] = n_max;
        temp_array[ 4*N_reac+i_reac] = T_min;
        temp_array[ 5*N_reac+i_reac] = T_max;
        switch(dim) {
            case 1:
                err += interp1Dcomp_init_coeff(
                    temp_arr_pos, offload_arr_pos,
                    N_E,
                    NATURALBC,
                    E_min, E_max);
                temp_arr_pos    += N_E * NSIZE_COMP1D;
                offload_arr_pos += N_E;
                break;

            case 2:
                err += interp2Dcomp_init_coeff(
                    temp_arr_pos, offload_arr_pos,
                    N_E, N_T,
                    NATURALBC, NATURALBC,
                    E_min, E_max,
                    T_min, T_max);
                temp_arr_pos    += N_E * N_T * NSIZE_COMP2D;
                offload_arr_pos += N_E * N_T;
                break;

            case 3:
                err += interp3Dcomp_init_coeff(
                    temp_arr_pos, offload_arr_pos,
                    N_E, N_n, N_T,
                    NATURALBC, NATURALBC, NATURALBC,
                    E_min, E_max,
                    n_min, n_max,
                    T_min, T_max);
                temp_arr_pos    += N_E * N_n * N_T * NSIZE_COMP3D;
                offload_arr_pos += N_E * N_n * N_T;
            break;

            default:
                /* Unrecognized dimensionality. Produce error. */
                print_err("Error: Unrecognized abscissa dimensionality\n");
                err = 1;
                break;
        }
    }

    free(*offload_array);
    *offload_array = temp_array;
    offload_data->offload_array_length = temp_arr_length;

    if(err) {
        return err;
    }

    print_out(VERBOSE_IO,"\nFound data for %d atomic reactions:\n", N_reac);
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        print_out(VERBOSE_IO,
              "Reaction number / Total number of reactions = %d / %d\n"
              "  Reactant species Z_1 / A_1, Z_2 / A_2     = %d / %d, %d / %d\n"
              "  Min/Max energy                            = %1.2le / %1.2le\n"
              "  Min/Max density                           = %1.2le / %1.2le\n"
              "  Min/Max temperature                       = %1.2le / %1.2le\n"
              "  Number of energy grid points              = %d\n"
              "  Number of density grid points             = %d\n"
              "  Number of temperature grid points         = %d\n",
              i_reac+1, N_reac,
              offload_data->z_1[i_reac], offload_data->a_1[i_reac],
              offload_data->z_2[i_reac], offload_data->a_2[i_reac],
              (*offload_array)[0 * N_reac + i_reac],
              (*offload_array)[1 * N_reac + i_reac],
              (*offload_array)[2 * N_reac + i_reac],
              (*offload_array)[3 * N_reac + i_reac],
              (*offload_array)[4 * N_reac + i_reac],
              (*offload_array)[5 * N_reac + i_reac],
              offload_data->N_E[i_reac],
              offload_data->N_n[i_reac],
              offload_data->N_T[i_reac]);
    }

    return 0;
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
*/
void asigma_loc_free_offload(asigma_loc_offload_data* offload_data,
                             real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize atomic reaction data struct on target
 *
 * This function copies atomic reaction data from the offload struct
 * and the offload array to the struct on the target, and uses
 * interp?Dcomp_init_spline() functions to initialize the precalculated
 * spline parameters of the reaction data in the spline structs within
 * the data struct on the target.
 *
 * @param asgm_loc_data pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void asigma_loc_init(
    asigma_loc_data* asgm_loc_data, asigma_loc_offload_data* offload_data,
    real* offload_array) {
    /* Copy over number of reactions and store it in a helper variable */
    int N_reac =  offload_data->N_reac;
    asgm_loc_data->N_reac = N_reac;

    /* Helper pointer to keep track of position in offload array */
    real* offload_arr_pos = offload_array+5*N_reac;

    /* Copy data from offload array to atomic sigma struct,
       initialize spline structs and determine reaction availability */
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {

        /* Reaction identifiers */
        asgm_loc_data->z_1[i_reac] = offload_data->z_1[i_reac];
        asgm_loc_data->a_1[i_reac] = offload_data->a_1[i_reac];
        asgm_loc_data->z_2[i_reac] = offload_data->z_2[i_reac];
        asgm_loc_data->a_2[i_reac] = offload_data->a_2[i_reac];
        asgm_loc_data->reac_type[i_reac] = offload_data->reac_type[i_reac];

        /* Initialize spline struct according to dimensionality of
           reaction data (and mark reaction availability) */
        int  N_E   = offload_data->N_E[i_reac];
        real E_min = offload_array[0*N_reac+i_reac];
        real E_max = offload_array[1*N_reac+i_reac];
        int  N_n   = offload_data->N_n[i_reac];
        real n_min = offload_array[2*N_reac+i_reac];
        real n_max = offload_array[3*N_reac+i_reac];
        int  N_T   = offload_data->N_T[i_reac];
        real T_min = offload_array[4*N_reac+i_reac];
        real T_max = offload_array[5*N_reac+i_reac];
        int  dim   = (N_E > 1) + (N_n > 1) + (N_T > 1);
        switch(dim) {
            case 1:
                interp1Dcomp_init_spline(
                    &(asgm_loc_data->sigma[i_reac]), offload_arr_pos,
                    N_E,
                    NATURALBC,
                    E_min, E_max);
                offload_arr_pos += N_E *NSIZE_COMP1D;
                break;

            case 2:
                interp2Dcomp_init_spline(
                    &(asgm_loc_data->sigmav[i_reac]), offload_arr_pos,
                    N_E, N_T,
                    NATURALBC, NATURALBC,
                    E_min, E_max,
                    T_min, T_max);
                offload_arr_pos += N_E * N_T * NSIZE_COMP2D;
                break;

            case 3:
                interp3Dcomp_init_spline(
                    &(asgm_loc_data->BMSsigmav[i_reac]), offload_arr_pos,
                    N_E, N_n, N_T,
                    NATURALBC, NATURALBC, NATURALBC,
                    E_min, E_max,
                    n_min, n_max,
                    T_min, T_max);
                offload_arr_pos += N_E * N_n * N_T * NSIZE_COMP3D;
                break;

            default:
                /* Unrecognized dimensionality. Produce error. */
                print_err("Error: Unrecognized abscissa dimensionality\n");
                break;
        }
    }
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
 * @param reac_type reaction type
 * @param asgm_loc_data pointer to atomic data struct
 * @param E_coll_per_amu energy per amu corresponding to collision speed
 * @param enable_atomic pointer to atomic enable and functionality flag
 *
 * @return zero if evaluation succeeded
 */
a5err asigma_loc_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, int reac_type,
    asigma_loc_data* asgm_loc_data, real E_coll_per_amu, int* enable_atomic) {
    a5err err = 0;

    /* We look for a match of the reaction identifiers in asgm_loc_data to
       determine if the reaction of interest has been initialized */
    int reac_found = -1, i_reac;
    for(i_reac = 0; i_reac < asgm_loc_data->N_reac; i_reac++) {
        if(z_1       == asgm_loc_data->z_1[i_reac] &&
           a_1       == asgm_loc_data->a_1[i_reac] &&
           z_2       == asgm_loc_data->z_2[i_reac] &&
           a_2       == asgm_loc_data->a_2[i_reac] &&
           reac_type == asgm_loc_data->reac_type[i_reac]) {
            reac_found = i_reac;
        }
    }
    i_reac = reac_found;

    /* The cross-section is evaluated if reaction data was found,
       is available, and its interpolation implemented. Otherwise,
       the cross-section is set to zero to avoid further problems. */
    if(reac_found < 0) {
        /* Reaction not found. Raise error. */
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
    } else {
        if(asgm_loc_data->reac_type[i_reac] == sigma_ioniz  ||
           asgm_loc_data->reac_type[i_reac] == sigma_recomb ||
           asgm_loc_data->reac_type[i_reac] == sigma_CX) {
            int interperr = 0;
            interperr += interp1Dcomp_eval_f(sigma,
                                             &asgm_loc_data->sigma[i_reac],
                                             E_coll_per_amu);
            if(interperr) {
                /* Energy is outside spline domain */
                if(*enable_atomic == 2) {
                    /* Set cross-section to 0 to avoid further problems */
                    *sigma = 0.0;
                } else {
                    /* Raise an error */
                    err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                       EF_ASIGMA_LOC );
                }
            }
        } else {
            /* Interpolation of cross-section not implemented. Raise error. */
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
        }
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
 * @param reac_type reaction type
 * @param asgm_loc_data pointer to atomic data struct
 * @param E energy of fast particle
 * @param T_e electron temperature of bulk plasma
 * @param T_0 temperature of bulk neutrals
 * @param n_i ion density of bulk plasma
 * @param enable_atomic pointer to atomic enable and functionality flag
 *
 * @return zero if evaluation succeeded
 */
a5err asigma_loc_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    int reac_type, asigma_loc_data* asgm_loc_data, real E,
    real T_e, real T_0, real n_i, int* enable_atomic) {
    a5err err = 0;

    /* Find the matching reaction. Note that BMS data is same for all
     * isotopes, so we don't compare anums */
    int reac_found = -1, i_reac;
    for(i_reac = 0; i_reac < asgm_loc_data->N_reac; i_reac++) {
        if(reac_type == sigmav_BMS &&
           z_1       == asgm_loc_data->z_1[i_reac] &&
           z_2       == asgm_loc_data->z_2[i_reac] &&
           reac_type == asgm_loc_data->reac_type[i_reac]) {
            reac_found = i_reac;
        } else if(z_1       == asgm_loc_data->z_1[i_reac] &&
                  a_1       == asgm_loc_data->a_1[i_reac] &&
                  z_2       == asgm_loc_data->z_2[i_reac] &&
                  a_2       == asgm_loc_data->a_2[i_reac] &&
                  reac_type == asgm_loc_data->reac_type[i_reac]) {
            reac_found = i_reac;
        }
    }
    i_reac = reac_found;

    /* The rate coefficient is evaluated if reaction data was found,
       is available, and its interpolation works. */
    if(reac_found < 0) {
        /* Reaction not found. Raise error. */
        err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
    } else {
        if(reac_type == sigmav_ioniz  ||
           reac_type == sigmav_recomb ||
           reac_type == sigmav_CX) {
            /* We remember unit conversions for parameters in the
               function call */
            int interperr = 0;
            interperr += interp2Dcomp_eval_f(sigmav,
                                             &asgm_loc_data->sigmav[i_reac],
                                             E/CONST_E, T_0/CONST_E);
            if(interperr) {
                /* An abscissa is outside spline domain */
                if(*enable_atomic == 2) {
                    /* Set rate coefficient to 0 to avoid further problems */
                    *sigmav = 0.0;
                } else {
                    /* Raise an error */
                    err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                       EF_ASIGMA_LOC );
                }
            }
        } else if(reac_type == sigmav_BMS) {
            /* Remember unit conversions for parameters as well as
               normalization of the fast particle energy with respect
               to mass (J --> eV/amu) in the function call. */
            int interperr = 0;
            interperr += interp3Dcomp_eval_f(
                sigmav, &asgm_loc_data->BMSsigmav[i_reac],
                E/CONST_E/(m_1/CONST_U), z_2*n_i, T_e/CONST_E);
            if(interperr) {
                /* An abscissa is outside spline domain */
                if(*enable_atomic == 2) {
                    /* Set rate coefficient to 0 to avoid further problems */
                    *sigmav = 0.0;
                } else {
                    /* Raise an error */
                    err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                       EF_ASIGMA_LOC );
                }
            }
        } else {
            /* Interpolation of rate coefficient not implemented.
               Raise error. */
            err = error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_ASIGMA_LOC);
        }
    }

    return err;
}
