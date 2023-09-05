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
 *           [0] = znum for fast particle species
 *      [N_reac] = anum for fast particle species
 *    [2*N_reac] = znum for bulk particle species
 *    [3*N_reac] = anum for bulk particle species
 *    [4*N_reac] = reaction type
 *    [5*N_reac] = dimension of energy abscissa (N_E)
 *    [6*N_reac] = min value of energy abscissa [eV]
 *    [7*N_reac] = max value of energy abscissa [eV]
 *    [8*N_reac] = dimension of density abscissa (N_n)
 *    [9*N_reac] = min value of density abscissa [m^-3]
 *   [10*N_reac] = max value of density abscissa [m^-3]
 *   [11*N_reac] = dimension of temperature abscissa (N_T)
 *   [12*N_reac] = min value of temperature abscissa [eV]
 *   [13*N_reac] = max value of temperature abscissa [eV]
 *   [14*N_reac] = reaction probability data [(depends on reaction data type)]
 * Each piece of data listed above is repeated for each reaction included.
 * Hence the N_reac interval between the 14 first pieces of data. The memory
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
 *
 * @todo In the part where the spline interpolations of the reaction data
 *   are initialized, it is assumed that, when there is one abscissa, it is
 *   energy, when there are two abscissae, they are, in order, energy and
 *   temperature, and, when there are three abscissae, they are, in order,
 *   energy, density and temperature. In the future, this may need to be
 *   made more general. For example, where we check the N:s, we can also
 *   extract and reorder the min:s and max:s depending on which N > 1.
 */
int asigma_loc_init_offload(asigma_loc_offload_data* offload_data,
                            real** offload_array) {
    int err = 0;
    int N_reac = offload_data->N_reac;

    /* Extract and store abscissa dimensions so we know how much memory
       needs to be allocated */
    int N_E_arr[N_reac];
    int N_n_arr[N_reac];
    int N_T_arr[N_reac];
    int N_arr[N_reac];
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        N_E_arr[i_reac] = (int) (*offload_array)[ 5*N_reac+i_reac];
        N_n_arr[i_reac] = (int) (*offload_array)[ 8*N_reac+i_reac];
        N_T_arr[i_reac] = (int) (*offload_array)[11*N_reac+i_reac];
        N_arr[i_reac] = N_E_arr[i_reac]*N_n_arr[i_reac]*N_T_arr[i_reac];
    }

    /* Allocate space for reaction identifiers, abscissa parameters and
       reaction data splines in a temporary array that will become the new
       offload array. The fixed number of needed memory space in the
       beginning of the temporary array (offload array) is 14*N_reac,
       because there are 5 reaction identifiers, namely z_1, a_1, z_2, a_2
       and reac_type, and 9 abscissa parameters, namely N_E, E_min, E_max,
       N_n, n_min, n_max, N_T, T_min and T_max. */
    int temp_arr_length = 14*N_reac;
    int dim[N_reac];
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        dim[i_reac] =
            (N_E_arr[i_reac]>1) + (N_n_arr[i_reac]>1) + (N_T_arr[i_reac] > 1);
        switch(dim[i_reac]) {
            case 1:
                temp_arr_length += N_arr[i_reac]*NSIZE_COMP1D;
                break;

            case 2:
                temp_arr_length += N_arr[i_reac]*NSIZE_COMP2D;
                break;

            case 3:
                temp_arr_length += N_arr[i_reac]*NSIZE_COMP3D;
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
    real * temp_arr_pos = temp_array+14*N_reac;
    real * offload_arr_pos = *offload_array+14*N_reac;

    /* Copy over reaction identifiers and abscissae parameters, and evaluate
       spline coefficients according to dimensionality of reaction data */
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        temp_array[ 0*N_reac+i_reac] = (*offload_array)[ 0*N_reac+i_reac];
        temp_array[ 1*N_reac+i_reac] = (*offload_array)[ 1*N_reac+i_reac];
        temp_array[ 2*N_reac+i_reac] = (*offload_array)[ 2*N_reac+i_reac];
        temp_array[ 3*N_reac+i_reac] = (*offload_array)[ 3*N_reac+i_reac];
        temp_array[ 4*N_reac+i_reac] = (*offload_array)[ 4*N_reac+i_reac];
        temp_array[ 5*N_reac+i_reac] = (*offload_array)[ 5*N_reac+i_reac];
        temp_array[ 6*N_reac+i_reac] = (*offload_array)[ 6*N_reac+i_reac];
        temp_array[ 7*N_reac+i_reac] = (*offload_array)[ 7*N_reac+i_reac];
        temp_array[ 8*N_reac+i_reac] = (*offload_array)[ 8*N_reac+i_reac];
        temp_array[ 9*N_reac+i_reac] = (*offload_array)[ 9*N_reac+i_reac];
        temp_array[10*N_reac+i_reac] = (*offload_array)[10*N_reac+i_reac];
        temp_array[11*N_reac+i_reac] = (*offload_array)[11*N_reac+i_reac];
        temp_array[12*N_reac+i_reac] = (*offload_array)[12*N_reac+i_reac];
        temp_array[13*N_reac+i_reac] = (*offload_array)[13*N_reac+i_reac];
        switch(dim[i_reac]) {
            case 1:
                err += interp1Dcomp_init_coeff(
                    temp_arr_pos, offload_arr_pos,
                    N_E_arr[i_reac],
                    NATURALBC,
                    temp_array[ 6*N_reac+i_reac], temp_array[ 7*N_reac+i_reac]);
                /* Update positions in arrays */
                temp_arr_pos    += N_arr[i_reac]*NSIZE_COMP1D;
                offload_arr_pos += N_arr[i_reac];
                break;

            case 2:
                err += interp2Dcomp_init_coeff(
                    temp_arr_pos, offload_arr_pos,
                    N_E_arr[i_reac], N_T_arr[i_reac],
                    NATURALBC, NATURALBC,
                    temp_array[ 6*N_reac+i_reac], temp_array[ 7*N_reac+i_reac],
                    temp_array[12*N_reac+i_reac], temp_array[13*N_reac+i_reac]);
                /* Update positions in arrays */
                temp_arr_pos    += N_arr[i_reac]*NSIZE_COMP2D;
                offload_arr_pos += N_arr[i_reac];
                break;

            case 3:
                err += interp3Dcomp_init_coeff(
                    temp_arr_pos, offload_arr_pos,
                    N_E_arr[i_reac], N_n_arr[i_reac], N_T_arr[i_reac],
                    NATURALBC, NATURALBC, NATURALBC,
                    temp_array[ 6*N_reac+i_reac], temp_array[ 7*N_reac+i_reac],
                    temp_array[ 9*N_reac+i_reac], temp_array[10*N_reac+i_reac],
                    temp_array[12*N_reac+i_reac], temp_array[13*N_reac+i_reac]);
                /* Update positions in arrays */
                temp_arr_pos    += N_arr[i_reac]*NSIZE_COMP3D;
                offload_arr_pos += N_arr[i_reac];
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
              (int) *(*offload_array+ 0*N_reac+i_reac),
              (int) *(*offload_array+ 1*N_reac+i_reac),
              (int) *(*offload_array+ 2*N_reac+i_reac),
              (int) *(*offload_array+ 3*N_reac+i_reac),
              *(*offload_array+ 6*N_reac+i_reac),
              *(*offload_array+ 7*N_reac+i_reac),
              *(*offload_array+ 9*N_reac+i_reac),
              *(*offload_array+10*N_reac+i_reac),
              *(*offload_array+12*N_reac+i_reac),
              *(*offload_array+13*N_reac+i_reac),
              N_E_arr[i_reac],
              N_n_arr[i_reac],
              N_T_arr[i_reac]);
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
 * interpXDcomp_init_spline() functions to initialize the precalculated
 * spline parameters of the reaction data in the spline structs within
 * the data struct on the target.
 *
 * @param asgm_loc_data pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @todo Should we change from using reac_avail flag to having in
 *   asgm_data only those reactions that are available and
 *   truly initialized? Actually, that is already done, reac_avail
 *   is a relic. However, we might want to use reac_avail for something.
 * @todo Can we avoid the excess memory usage currently resulting from
 *   allocating all three spline structs for each reaction?
 * @todo (Same as above in asigma_loc_init_offload()) In the part where
 *   the spline interpolations of the reaction data
 *   are initialized, it is assumed that, when there is one abscissa, it is
 *   energy, when there are two abscissae, they are, in order, energy and
 *   temperature, and, when there are three abscissae, they are, in order,
 *   energy, density and temperature. In the future, this may need to be
 *   made more general. For example, where we check the N:s, we can also
 *   extract and reorder the min:s and max:s depending on which N > 1.
 */
void asigma_loc_init(asigma_loc_data* asgm_loc_data,
                     asigma_loc_offload_data* offload_data,
                     real* offload_array) {
    /* Copy over number of reactions and store it in a helper variable */
    int N_reac =  offload_data->N_reac;
    asgm_loc_data->N_reac = N_reac;

    /* Allocate memory for atomic sigma struct entries */
    asgm_loc_data->z_1        = malloc(N_reac*sizeof(int));
    asgm_loc_data->a_1        = malloc(N_reac*sizeof(int));
    asgm_loc_data->z_2        = malloc(N_reac*sizeof(int));
    asgm_loc_data->a_2        = malloc(N_reac*sizeof(int));
    asgm_loc_data->reac_type  = malloc(N_reac*sizeof(int));
    asgm_loc_data->reac_avail = malloc(N_reac*sizeof(int));

    /* Extract abscissa parameters into helper arrays, and construct
       a helper variable for abscissa dimensionalities */
    int  N_E_arr[N_reac];
    real E_min_arr[N_reac];
    real E_max_arr[N_reac];
    int  N_n_arr[N_reac];
    real n_min_arr[N_reac];
    real n_max_arr[N_reac];
    int  N_T_arr[N_reac];
    real T_min_arr[N_reac];
    real T_max_arr[N_reac];
    int dim[N_reac];
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        N_E_arr[i_reac]   = (int) offload_array[ 5*N_reac+i_reac];
        E_min_arr[i_reac] =       offload_array[ 6*N_reac+i_reac];
        E_max_arr[i_reac] =       offload_array[ 7*N_reac+i_reac];
        N_n_arr[i_reac]   = (int) offload_array[ 8*N_reac+i_reac];
        n_min_arr[i_reac] =       offload_array[ 9*N_reac+i_reac];
        n_max_arr[i_reac] =       offload_array[10*N_reac+i_reac];
        N_T_arr[i_reac]   = (int) offload_array[11*N_reac+i_reac];
        T_min_arr[i_reac] =       offload_array[12*N_reac+i_reac];
        T_max_arr[i_reac] =       offload_array[13*N_reac+i_reac];
        dim[i_reac] =
            (N_E_arr[i_reac]>1) + (N_n_arr[i_reac]>1) + (N_T_arr[i_reac]>1);
    }
    /* Allocate memory for spline data for the different reactions.
       NOTE: Currently all three spline structs are allocated for
       each reaction, which results in excess memory usage. This
       full allocation is done to make sure that the i_reac counting
       here and elsewhere works. */
    asgm_loc_data->sigma     = malloc(N_reac*sizeof(interp1D_data));
    asgm_loc_data->sigmav    = malloc(N_reac*sizeof(interp2D_data));
    asgm_loc_data->BMSsigmav = malloc(N_reac*sizeof(interp3D_data));

    /* Helper pointer to keep track of position in offload array */
    real* offload_arr_pos = offload_array+14*N_reac;

    /* Copy data from offload array to atomic sigma struct,
       initialize spline structs and determine reaction availability */
    for(int i_reac = 0; i_reac < N_reac; i_reac++) {
        /* Reaction identifiers */
        asgm_loc_data->z_1[i_reac]       = (int) offload_array[0*N_reac+i_reac];
        asgm_loc_data->a_1[i_reac]       = (int) offload_array[1*N_reac+i_reac];
        asgm_loc_data->z_2[i_reac]       = (int) offload_array[2*N_reac+i_reac];
        asgm_loc_data->a_2[i_reac]       = (int) offload_array[3*N_reac+i_reac];
        asgm_loc_data->reac_type[i_reac] = (int) offload_array[4*N_reac+i_reac];
        /* Initialize spline struct according to dimensionality of
           reaction data (and mark reaction availability) */
        switch(dim[i_reac]) {
            case 1:
                interp1Dcomp_init_spline(
                    &(asgm_loc_data->sigma[i_reac]), offload_arr_pos,
                    N_E_arr[i_reac],
                    NATURALBC,
                    E_min_arr[i_reac], E_max_arr[i_reac]);
                asgm_loc_data->reac_avail[i_reac] = 1;
                /* Update position in array (2 of the N:s should be = 1) */
                offload_arr_pos += N_E_arr[i_reac]*
                                   N_n_arr[i_reac]*
                                   N_T_arr[i_reac]*NSIZE_COMP1D;
                break;

            case 2:
                interp2Dcomp_init_spline(
                    &(asgm_loc_data->sigmav[i_reac]), offload_arr_pos,
                    N_E_arr[i_reac], N_T_arr[i_reac],
                    NATURALBC, NATURALBC,
                    E_min_arr[i_reac], E_max_arr[i_reac],
                    T_min_arr[i_reac], T_max_arr[i_reac]);
                asgm_loc_data->reac_avail[i_reac] = 1;
                /* Update position in array (1 of the N:s should be = 1) */
                offload_arr_pos += N_E_arr[i_reac]*
                                   N_n_arr[i_reac]*
                                   N_T_arr[i_reac]*NSIZE_COMP2D;
                break;

            case 3:
                interp3Dcomp_init_spline(
                    &(asgm_loc_data->BMSsigmav[i_reac]), offload_arr_pos,
                    N_E_arr[i_reac], N_n_arr[i_reac], N_T_arr[i_reac],
                    NATURALBC, NATURALBC, NATURALBC,
                    E_min_arr[i_reac], E_max_arr[i_reac],
                    n_min_arr[i_reac], n_max_arr[i_reac],
                    T_min_arr[i_reac], T_max_arr[i_reac]);
                asgm_loc_data->reac_avail[i_reac] = 1;
                /* Update position in array */
                offload_arr_pos += N_E_arr[i_reac]*
                                   N_n_arr[i_reac]*
                                   N_T_arr[i_reac]*NSIZE_COMP3D;
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
a5err asigma_loc_eval_sigma(real* sigma,
                            int z_1, int a_1,
                            int z_2, int a_2,
                            int reac_type,
                            asigma_loc_data* asgm_loc_data,
                            real E_coll_per_amu,
                            int* enable_atomic) {
    a5err err = 0;

    /* We look for a match of the reaction identifiers in asgm_loc_data to
       determine if the reaction of interest has been initialized */
    int reac_found = 0;
    int i_reac;
    for(i_reac = 0; i_reac < asgm_loc_data->N_reac; i_reac++) {
        if(z_1       == asgm_loc_data->z_1[i_reac] &&
           a_1       == asgm_loc_data->a_1[i_reac] &&
           z_2       == asgm_loc_data->z_2[i_reac] &&
           a_2       == asgm_loc_data->a_2[i_reac] &&
           reac_type == asgm_loc_data->reac_type[i_reac]) {
            reac_found = 1;
            break;
        }
    }

    /* The cross-section is evaluated if reaction data was found,
       is available, and its interpolation implemented. Otherwise,
       the cross-section is set to zero to avoid further problems. */
    if(!reac_found) {
        /* Reaction not found. Raise error. */
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                           EF_ASIGMA_LOC );
    } else if(!asgm_loc_data->reac_avail[i_reac]) {
        /* Reaction not available. Raise error. */
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                           EF_ASIGMA_LOC );
    } else {
        if(asgm_loc_data->reac_type[i_reac] == reac_type_sigma_ion ||
           asgm_loc_data->reac_type[i_reac] == reac_type_sigma_rec ||
           asgm_loc_data->reac_type[i_reac] == reac_type_sigma_CX) {
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
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                               EF_ASIGMA_LOC );
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
 *
 * @todo Is there a cleaner and more general way to handle the
 *   isotopic invariance of beam-stopping data?
 */
a5err asigma_loc_eval_sigmav(real* sigmav,
                             int z_1, int a_1, real m_1,
                             int z_2, int a_2,
                             int reac_type,
                             asigma_loc_data* asgm_loc_data,
                             real E,
                             real T_e, real T_0, real n_i,
                             int* enable_atomic) {
    a5err err = 0;

    /* We look for a match of the reaction identifiers in asgm_loc_data to
       determine if the reaction of interest has been initialized.
       NOTE: Beam-stopping is a special case. The reaction data function
       for a certain element covers all of its different isotopes, as long
       as the energy parameter is given in units of energy/amu. Hence, for
       beam-stopping, the atomic mass numbers need not match below. */
    int reac_found = 0;
    int i_reac;
    for(i_reac = 0; i_reac < asgm_loc_data->N_reac; i_reac++) {
        if(reac_type == reac_type_BMS_sigmav &&
           z_1       == asgm_loc_data->z_1[i_reac] &&
           z_2       == asgm_loc_data->z_2[i_reac] &&
           reac_type == asgm_loc_data->reac_type[i_reac]) {
            reac_found = 1;
            break;
        } else if(z_1       == asgm_loc_data->z_1[i_reac] &&
                  a_1       == asgm_loc_data->a_1[i_reac] &&
                  z_2       == asgm_loc_data->z_2[i_reac] &&
                  a_2       == asgm_loc_data->a_2[i_reac] &&
                  reac_type == asgm_loc_data->reac_type[i_reac]) {
            reac_found = 1;
            break;
        }
    }

    /* The rate coefficient is evaluated if reaction data was found,
       is available, and its interpolation works. */
    if(!reac_found) {
        /* Reaction not found. Raise error. */
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                           EF_ASIGMA_LOC );
    } else if(!asgm_loc_data->reac_avail[i_reac]) {
        /* Reaction not available. Raise error. */
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                           EF_ASIGMA_LOC );
    } else {
        if(reac_type == reac_type_sigmav_ion ||
           reac_type == reac_type_sigmav_rec ||
           reac_type == reac_type_sigmav_CX) {
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
        } else if(reac_type == reac_type_BMS_sigmav) {
            /* Remember unit conversions for parameters as well as
               normalization of the fast particle energy with respect
               to mass (J --> eV/amu) in the function call. */
            int interperr = 0;
            interperr += interp3Dcomp_eval_f(sigmav,
					     &asgm_loc_data->BMSsigmav[i_reac],
					     E/CONST_E/(m_1/CONST_U),
					     z_2*n_i,
					     T_e/CONST_E);
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
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                               EF_ASIGMA_LOC );
        }
    }

    return err;
}
