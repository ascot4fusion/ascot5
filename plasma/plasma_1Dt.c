/**
 * @file plasma_1Dt.c
 * @brief 1D time-dependent plasma with linear interpolation
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../error.h"
#include "../consts.h"
#include "../print.h"
#include "plasma_1Dt.h"


/**
 * @brief Initialize 1D plasma data and check inputs
 *
 * Before calling this function, the offload struct is expected to be fully
 * initialized.
 *
 * The offload array is expected to hold plasma data as
 *   -                              [0] = rho grid
 *   -                          [n_rho] = electron temperature [J]
 *   -                        [n_rho*2] = ion temperature [J]
 *   -         [n_rho*2 + n_rho*n_ions] = electron density [m^-3]
 *   - [n_rho*2 + n_rho*n_ions + n_rho] = ion density [m^-3]
 *
 * Since this data requires no initialization, the only thing this function does
 * is that it prints some values as sanity check.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succes
 */
int plasma_1Dt_init_offload(plasma_1Dt_offload_data* offload_data,
                            real** offload_array) {

    int n_rho = offload_data->n_rho;
    int n_time = offload_data->n_time;
    int n_ions = offload_data->n_species -1;
    print_out(VERBOSE_IO, "\n1D plasma profiles (P_1Dt)\n");
    print_out(VERBOSE_IO,
              "Min rho = %1.2le, Max rho = %1.2le,"
              " Number of rho grid points = %d\n",
              (*offload_array)[0], (*offload_array)[n_rho-1], n_rho);
    print_out(VERBOSE_IO,
              "Min time = %1.2le, Max time = %1.2le,"
              " Number of time points = %d\n",
              (*offload_array)[n_rho], (*offload_array)[n_rho+n_time-1],n_time);
    print_out(VERBOSE_IO,
              "Number of ion species = %d\n", n_ions);
    return 0;
}

/**
 * @brief Free offload array and reset parameters
 *
 *This function deallocates the offload_array.

 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
*/
void plasma_1Dt_free_offload(plasma_1Dt_offload_data* offload_data,
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
 * @param pls_data pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
*/
void plasma_1Dt_init(plasma_1Dt_data* pls_data,
                     plasma_1Dt_offload_data* offload_data,
                     real* offload_array) {

    pls_data->n_rho = offload_data->n_rho;
    pls_data->n_time = offload_data->n_time;
    pls_data->n_species = offload_data->n_species;

    for(int i = 0; i < pls_data->n_species; i++) {
        pls_data->mass[i] = offload_data->mass[i];
        pls_data->charge[i] = offload_data->charge[i];
    }
    pls_data->rho  = &offload_array[0];
    pls_data->time = &offload_array[pls_data->n_rho];
    pls_data->temp = &offload_array[pls_data->n_rho+pls_data->n_time];
    pls_data->dens = &offload_array[pls_data->n_rho+pls_data->n_time
                                    +2*pls_data->n_rho*pls_data->n_time];
}

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param temp pointer to where evaluated temperature [J] is stored
 * @param rho radial coordinate
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_1Dt_eval_temp(real* temp, real rho, real t, int species,
                          plasma_1Dt_data* pls_data) {

    real temp_dens[MAX_SPECIES], temp_temp[MAX_SPECIES];

    a5err err = plasma_1Dt_eval_densandtemp(temp_dens, temp_temp, rho, t,
                                            pls_data);

    *temp = temp_temp[species];

    return err;

}

/**
 * @brief Evaluate plasma density
 *
 * This function evaluates the density of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param dens pointer to where evaluated density [m^-3] is stored
 * @param rho radial coordinate
 * @param species index of plasma species
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_1Dt_eval_dens(real* dens, real rho, real t, int species,
                           plasma_1Dt_data* pls_data) {
    real temp_dens[MAX_SPECIES], temp_temp[MAX_SPECIES];

    a5err err = plasma_1Dt_eval_densandtemp(temp_dens, temp_temp, rho, t,
                                            pls_data);

    *dens = temp_dens[species];

    return err;

}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 * This function evaluates the density and temperature of all plasma species at
 * the given radial coordinate using linear interpolation.
 *
 * @param dens pointer to where interpolated densities [m^-3] are stored
 * @param temp pointer to where interpolated temperatures [J] are stored
 * @param rho radial coordinate
 * @param pls_data pointer to plasma data struct
 *
 * @return zero if evaluation succeeded
 */
a5err plasma_1Dt_eval_densandtemp(real* dens, real* temp, real rho, real t,
                                  plasma_1Dt_data* pls_data) {

    a5err err = 0;
    if(rho < pls_data->rho[0]) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D );
    }
    else if(rho >= pls_data->rho[pls_data->n_rho-1]) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D );
    }
    else {
        int i_rho = 0;
        while(i_rho < pls_data->n_rho-1 && pls_data->rho[i_rho] <= rho) {
            i_rho++;
        }
        i_rho--;

        real t_rho = (rho - pls_data->rho[i_rho])
                 / (pls_data->rho[i_rho+1] - pls_data->rho[i_rho]);

        int i_time = 0;
        while(i_time < pls_data->n_time-1 && pls_data->time[i_time] <= t) {
            i_time++;
        }
        i_time--;

        real t_time = (t - pls_data->time[i_time])
                 / (pls_data->time[i_time+1] - pls_data->time[i_time]);

        if(i_time < 0) {
            /* time < t[0], use first profile */
            i_time = 0;
            t_time = 0;
        }
        else if(i_time >= pls_data->n_time-2) {
            /* time > t[n_time-1], use last profile */
            i_time = pls_data->n_time-2;
            t_time = 1;
        }

        for(int i = 0; i < pls_data->n_species; i++) {
            real p11, p12, p21, p22, p1, p2;

            p11 = pls_data->dens[i_time*pls_data->n_species*pls_data->n_rho
                                 + i*pls_data->n_rho
                                 + i_rho];
            p12 = pls_data->dens[i_time*pls_data->n_species*pls_data->n_rho
                                 + i*pls_data->n_rho
                                 + i_rho + 1];
            p21 = pls_data->dens[(i_time+1)*pls_data->n_species*pls_data->n_rho 
                                 + i*pls_data->n_rho
                                 + i_rho];
            p22 = pls_data->dens[(i_time+1)*pls_data->n_species*pls_data->n_rho 
                                 + i*pls_data->n_rho
                                 + i_rho + 1];
            
            p1 = p11 + t_rho * (p12 - p11);
            p2 = p21 + t_rho * (p22 - p21);

            dens[i] = p1 + t_time * (p2 - p1);

            if(i < 2) {
                /* Electron and ion temperature */
                p11 = pls_data->temp[i_time*2*pls_data->n_rho
                                     + i*pls_data->n_rho
                                     +i_rho];
                p12 = pls_data->temp[i_time*2*pls_data->n_rho
                                     + i*pls_data->n_rho
                                     + i_rho + 1];
                p21 = pls_data->temp[(i_time+1)*2*pls_data->n_rho
                                     + i*pls_data->n_rho
                                     + i_rho];
                p22 = pls_data->temp[(i_time+1)*2*pls_data->n_rho
                                     + i*pls_data->n_rho
                                     + i_rho + 1];

                p1 = p11 + t_rho * (p12 - p11);
                p2 = p21 + t_rho * (p22 - p21);

                temp[i] = p1 + t_time * (p2 - p1);
            }
            else {
                /* Temperature is same for all ion species */
                temp[i] = temp[1];
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
int plasma_1Dt_get_n_species(plasma_1Dt_data* pls_data) {
    return pls_data->n_species;
}

/**
 * @brief Return pointer to array storing species mass
 *
 * @param pls_data pointer to plasma data
 *
 * @return pointer to immutable MAX_SPECIES length array containing masses
 */
const real* plasma_1Dt_get_species_mass(plasma_1Dt_data* pls_data) {
    return pls_data->mass;
}

/**
 * @brief Return pointer to array storing species charge
 *
 * @param pls_data pointer to plasma data
 *
 * @return pointer to immutable MAX_SPECIES length array containing charges
 */
const real* plasma_1Dt_get_species_charge(plasma_1Dt_data* pls_data) {
    return pls_data->charge;
}
