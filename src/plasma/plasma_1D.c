/**
 * @file plasma_1D.c
 * @brief 1D linearly interpolated plasma
 *
 * Plasma data which is defined in a 1D uniform grid from which the values are
 * interpolated linearly. The coordinate is the normalized poloidal flux.
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../ascot5.h"
#include "../error.h"
#include "../consts.h"
#include "../print.h"
#include "plasma_1D.h"


/**
 * @brief Initialize 1D plasma data and check inputs
 *
 * @param data pointer to the data struct
 *
 * @return zero if initialization succes
 */
int plasma_1D_init(plasma_1D_data* data, int nrho, int nion, real* rho,
                   int* anum, int* znum, real* mass, real* charge,
                   real* Te, real* Ti, real* ne, real* ni) {

    data->n_rho = nrho;
    data->n_species = nion + 1;

    data->anum = (int*) malloc( nion*sizeof(int) );
    data->znum = (int*) malloc( nion*sizeof(int) );
    data->mass = (real*) malloc( (nion+1)*sizeof(real) );
    data->charge = (real*) malloc( (nion+1)*sizeof(real) );
    for(int i = 0; i < data->n_species; i++) {
        if(i < nion) {
            data->znum[i] = znum[i];
            data->anum[i] = anum[i];
        }
        data->mass[i] = mass[i];
        data->charge[i] = charge[i];
    }
    data->rho = (real*) malloc( nrho*sizeof(real) );
    data->temp = (real*) malloc( 2*nrho*sizeof(real) );
    data->dens = (real*) malloc( (nion+1)*nrho*sizeof(real) );
    for(int i = 0; i < data->n_rho; i++) {
        data->rho[i] = rho[i];
        data->temp[i] = Te[i];
        data->temp[nrho + i] = Ti[i];
        data->dens[i] = ne[i];
        for(int j = 0; j < nion; j++) {
            data->dens[(j+1) * nrho + i] = ni[j*nrho + i];
        }
    }

    print_out(VERBOSE_IO, "\n1D plasma profiles (P_1D)\n");
    print_out(VERBOSE_IO,
              "Min rho = %1.2le, Max rho = %1.2le,"
              " Number of rho grid points = %d,"
              " Number of ion species = %d\n",
              data->rho[0], data->rho[data->n_rho-1], data->n_rho, nion);
    print_out(VERBOSE_IO,
              "Species Z/A  charge [e]/mass [amu] Density [m^-3] at Min/Max rho"
              "    Temperature [eV] at Min/Max rho\n");
    for(int i=0; i < nion; i++) {
        print_out(VERBOSE_IO,
                  " %3d  /%3d   %3d  /%7.3f             %1.2le/%1.2le     "
                  "           %1.2le/%1.2le       \n",
                  data->znum[i], data->anum[i],
                  (int)round(data->charge[i+1]/CONST_E),
                  data->mass[i+1]/CONST_U,
                  data->dens[(i+1)*nrho], data->dens[(i+1)*nrho - 1],
                  data->temp[nrho] / CONST_E, data->temp[2*nrho-1] / CONST_E);
    }
    print_out(VERBOSE_IO,
              "[electrons]  %3d  /%7.3f             %1.2le/%1.2le          "
              "      %1.2le/%1.2le       \n",
              -1, CONST_M_E/CONST_U,
              data->dens[0], data->dens[nrho - 1],
              data->temp[0] / CONST_E, data->temp[nrho-1] / CONST_E);
    real quasineutrality = 0;
    for(int k = 0; k < nrho; k++) {
        real ele_qdens = data->dens[k] * CONST_E;
        real ion_qdens = 0;
        for(int i = 0; i < nion; i++) {
            ion_qdens += data->dens[(i+1)*nrho + k] * data->charge[i+1];
        }
        quasineutrality = fmax( quasineutrality,
                                fabs( 1 - ion_qdens / ele_qdens ) );
    }
    print_out(VERBOSE_IO, "Quasi-neutrality is (electron / ion charge density)"
              " %.2f\n", 1+quasineutrality);
    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void plasma_1D_free(plasma_1D_data* data) {
    free(data->mass);
    free(data->charge);
    free(data->anum);
    free(data->znum);
    free(data->rho);
    free(data->temp);
    free(data->dens);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void plasma_1D_offload(plasma_1D_data* data) {
    GPU_MAP_TO_DEVICE(
        data->mass[0:data->n_species], data->charge[0:data->n_species], \
        data->anum[0:data->n_species-1], data->znum[0:data->n_species-1], \
        data->rho[0:data->n_rho], data->temp[0:2*data->n_rho], \
        data->dens[0:data->n_rho*data->n_species]
    )
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
a5err plasma_1D_eval_temp(real* temp, real rho, int species,
                          plasma_1D_data* pls_data) {

    a5err err = 0;
    if(rho < pls_data->rho[0]) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D );
    }
    else if(rho >= pls_data->rho[pls_data->n_rho-1]) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D );
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
        temp[0] = p1 + t_rho * (p2 - p1);
    }

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
a5err plasma_1D_eval_dens(real* dens, real rho, int species,
                          plasma_1D_data* pls_data) {

    a5err err = 0;
    if(rho < pls_data->rho[0]) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D );
    }
    else if(rho >= pls_data->rho[pls_data->n_rho-1]) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_PLASMA_1D );
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
        dens[0] = p1 + t_rho * (p2 - p1);
    }

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
a5err plasma_1D_eval_densandtemp(real* dens, real* temp, real rho,
                                 plasma_1D_data* pls_data) {
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

        real p1, p2;
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
