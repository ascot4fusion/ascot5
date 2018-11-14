/**
 * @file plasma_1D.c
 * @brief 1D linearly interpolated plasma
 *
 * Plasma data which is defined in a 1D uniform grid from which the values are
 * interpolated linearly. The coordinate is the normalized poloidal flux.
 */
#include <stdio.h>
#include <stdlib.h>
#include "../ascot5.h"
#include "../error.h"
#include "../consts.h"
#include "../print.h"
#include "plasma_1D.h"


/**
 * @brief Initialize 1D plasma data and check inputs
 *
 * Before calling this function, the offload struct is expected to be fully
 * initialized.
 *
 * The offload array is expected to hold plasma data as
 *   &(*offload_array)[0] = rho grid
 *   &(*offload_array)[n_rho] = electron temperature
 *   &(*offload_array)[n_rho*2] = ion temperature
 *   &(*offload_array)[n_rho*2 + n_rho*n_ions] = electron density
 *   &(*offload_array)[n_rho*2 + n_rho*n_ions + n_rho] = ion density
 *
 * Since this data requires no initialization, the only thing this function does
 * is that it prints some values as sanity check.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
int plasma_1D_init_offload(plasma_1D_offload_data* offload_data,
                            real** offload_array) {
    // Do no initialization

    int n_rho = offload_data->n_rho;
    int n_ions = offload_data->n_species -1;
    print_out(VERBOSE_IO, "\n1D plasma profiles (P_1D)\n");
    print_out(VERBOSE_IO,
              "Min rho = %1.2le, Max rho = %1.2le,"
              " Number of rho grid points = %d,"
              " Number of ion species = %d\n",
              (*offload_array)[0], (*offload_array)[n_rho-1], n_rho, n_ions);
    print_out(VERBOSE_IO,
              "Species Z/A    Density [m^-3] at Min/Max rho"
              "    Temperature [eV] at Min/Max rho\n");
    print_out(VERBOSE_IO,
              "     Electrons              %1.2le/%1.2le          "
              "        %1.2le/%1.2le       \n",
              (*offload_array)[n_rho*2 + n_rho*n_ions],
              (*offload_array)[n_rho*3 + n_rho*n_ions - 1],
              (*offload_array)[n_rho] * CONST_KB / CONST_E,
              (*offload_array)[n_rho*2-1] * CONST_KB / CONST_E);
    for(int i=0; i < n_ions; i++) {
        print_out(VERBOSE_IO,
                  "      %3d/%3d               %1.2le/%1.2le     "
                  "             %1.2le/%1.2le       \n",
                  (int)(offload_data->charge[i+1]/CONST_E),
                  (int)(offload_data->mass[i+1]/CONST_U),
                  (*offload_array)[n_rho*(3+i) + n_rho*n_ions],
                  (*offload_array)[n_rho*(4+i) + n_rho*n_ions - 1],
                  (*offload_array)[n_rho*2] * CONST_KB / CONST_E,
                  (*offload_array)[n_rho*3-1] * CONST_KB / CONST_E);
    }
    real quasineutrality = 0;
    for(int k = 0; k <n_rho; k++) {
        real ele_qdens = (*offload_array)[n_rho*2 + n_rho*n_ions + k] * CONST_E;
        real ion_qdens = 0;
        for(int i=0; i < n_ions; i++) {
            ion_qdens +=
                (*offload_array)[n_rho*(3+i) + n_rho*n_ions + k]
                * offload_data->charge[i+1];
        }
        quasineutrality = fmax( quasineutrality,
                                fabs( 1 - ion_qdens / ele_qdens ) );
    }
    print_out(VERBOSE_IO, "Quasi-neutrality is (electron / ion charge density)"
              " %.2f\n", 1+quasineutrality);
    return 0;
}

/**
 * @brief Free offload array and reset parameters
 *
 *This function deallocates the offload_array.

 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
*/
void plasma_1D_free_offload(plasma_1D_offload_data* offload_data,
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
int plasma_1D_init(plasma_1D_data* pls_data,
                    plasma_1D_offload_data* offload_data,
                    real* offload_array) {

    int err = 0;

    pls_data->n_rho = offload_data->n_rho;
    pls_data->n_species = offload_data->n_species;
    int i;
    for(i = 0; i < pls_data->n_species; i++) {
        pls_data->mass[i] = offload_data->mass[i];
        pls_data->charge[i] = offload_data->charge[i];
    }
    pls_data->rho = &offload_array[0];
    pls_data->temp = &offload_array[pls_data->n_rho];
    pls_data->dens = &offload_array[pls_data->n_rho
                                  + pls_data->n_rho*pls_data->n_species];

    return err;
}

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param rho radial coordinate
 * @param species index of plasma species
 * @param plasma_data pointer to plasma data struct
 */
real plasma_1D_eval_temp(real rho, int species, plasma_1D_data* pls_data) {

    /* As the plasma data may be provided at irregular intervals, we must
     * search for the correct grid index */
    /** @todo Implement a more efficient search algorithm */
    real p = 0;
    if(rho < pls_data->rho[0]) {
        p = pls_data->temp[species*pls_data->n_rho];
    }
    else if(rho >= pls_data->rho[pls_data->n_rho-1]) {
        p = pls_data->temp[species*pls_data->n_rho + pls_data->n_rho - 1];
    }
    else {
        int i_rho = 0;
        while(i_rho < pls_data->n_rho - 1 && pls_data->rho[i_rho] <= rho) {
            i_rho++;
        }
        i_rho--;
        real t_rho = (rho - pls_data->rho[i_rho])
            / (pls_data->rho[i_rho+1] - pls_data->rho[i_rho]);

        real p1 = pls_data->temp[species*pls_data->n_rho + i_rho];
        real p2 = pls_data->temp[species*pls_data->n_rho + i_rho+1];
        p = p1 + t_rho * (p2 - p1);
    }

    return p;
}

/**
 * @brief Evaluate plasma density
 *
 * @see plasma_1d_eval_temp
 */
real plasma_1D_eval_dens(real rho, int species, plasma_1D_data* pls_data) {
    real p = 0;
    if(rho < pls_data->rho[0]) {
        p = pls_data->dens[species*pls_data->n_rho];
    }
    else if(rho >= pls_data->rho[pls_data->n_rho-1]) {
        p = pls_data->dens[species*pls_data->n_rho + pls_data->n_rho - 1];
    }
    else {
        int i_rho = 0;
        while(i_rho < pls_data->n_rho - 1 && pls_data->rho[i_rho] <= rho) {
            i_rho++;
        }
        i_rho--;
        real t_rho = (rho - pls_data->rho[i_rho])
                 / (pls_data->rho[i_rho+1] - pls_data->rho[i_rho]);

        real p1 = pls_data->dens[species*pls_data->n_rho + i_rho];
        real p2 = pls_data->dens[species*pls_data->n_rho + i_rho+1];
        p = p1 + t_rho * (p2 - p1);
    }

    return p;
}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 */
a5err plasma_1D_eval_densandtemp(real rho, plasma_1D_data* pls_data,
                                 real* dens, real* temp) {
    a5err err = 0;
    real p1, p2;
    if(rho < pls_data->rho[0]) {
        for(int i = 0; i < pls_data->n_species; i++) {
            dens[i] = pls_data->dens[i*pls_data->n_rho];

            if(i < 2) {
                /* Electron and ion temperature */
                temp[i] = pls_data->temp[i*pls_data->n_rho];
            }
            else {
                /* Temperature is same for all ion species */
                temp[i] = temp[1];
            }
        }
    }
    else if(rho >= pls_data->rho[pls_data->n_rho-1]) {
        for(int i = 0; i < pls_data->n_species; i++) {
            dens[i] = pls_data->dens[i*pls_data->n_rho + pls_data->n_rho - 1];
            if(i < 2) {
                /* Electron and ion temperature */
                temp[i] =
                    pls_data->temp[i*pls_data->n_rho + pls_data->n_rho - 1];
            }
            else {
                /* Temperature is same for all ion species */
                temp[i] = temp[1];
            }
        }
    }
    else {
        int i_rho = 0;
        while(i_rho < pls_data->n_rho-1 && pls_data->rho[i_rho] <= rho) {
            i_rho++;
        }
        i_rho--;

        real t_rho = (rho - pls_data->rho[i_rho])
                 / (pls_data->rho[i_rho+1] - pls_data->rho[i_rho]);

        for(int i = 0; i < pls_data->n_species; i++) {
            p1 = pls_data->dens[i*pls_data->n_rho + i_rho];
            p2 = pls_data->dens[i*pls_data->n_rho + i_rho+1];
            dens[i] = p1 + t_rho * (p2 - p1);

            if(i < 2) {
                /* Electron and ion temperature */
                p1 = pls_data->temp[i*pls_data->n_rho + i_rho];
                p2 = pls_data->temp[i*pls_data->n_rho + i_rho+1];
                temp[i] = p1 + t_rho * (p2 - p1);
            }
            else {
                /* Temperature is same for all ion species */
                temp[i] = temp[1];
            }
        }
    }

    return err;
}

int plasma_1D_get_n_species(plasma_1D_data* pls_data) {
    return pls_data->n_species;
}

real* plasma_1D_get_species_mass(plasma_1D_data* pls_data) {
    return pls_data->mass;
}

real* plasma_1D_get_species_charge(plasma_1D_data* pls_data) {
    return pls_data->charge;
}
