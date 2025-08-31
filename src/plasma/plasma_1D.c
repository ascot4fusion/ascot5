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
 * @param nrho number of rho grid points in the data
 * @param nion number of ion species
 * @param rho rho grid in which data is tabulated [1]
 * @param anum mass number of ions present in the plasma
 * @param znum charge number of ions present in the plasma
 * @param mass mass of ions present in the plasma [kg]
 * @param charge charge of ions present in the plasma [C]
 * @param Te electron temperature [J]
 * @param Ti ion temperature [J]
 * @param ne electron density [m^-3]
 * @param ni density of ion species [m^-3]
 * @param vtor toroidal rotation [rad/s]
 *
 * @return zero if initialization succeeded
 */
int plasma_1D_init(plasma_1D_data* data, int nrho, int nion, real* rho,
                   int* anum, int* znum, real* mass, real* charge,
                   real* Te, real* Ti, real* ne, real* ni, real* vtor) {

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
    data->vtor = (real*) malloc( nrho*sizeof(real) );
    data->temp = (real*) malloc( 2*nrho*sizeof(real) );
    data->dens = (real*) malloc( (nion+1)*nrho*sizeof(real) );
    for(int i = 0; i < data->n_rho; i++) {
        data->rho[i] = rho[i];
        data->vtor[i] = vtor[i];
        data->temp[i] = Te[i];
        data->temp[nrho + i] = Ti[i];
        data->dens[i] = ne[i];
        for(int j = 0; j < nion; j++) {
            data->dens[(j+1) * nrho + i] = ni[j*nrho + i];
        }
    }
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
        data->vtor[0:data->n_rho], \
        data->dens[0:data->n_rho*data->n_species]
    )
}

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param temp pointer to where evaluated temperature is stored [J]
 * @param rho radial coordinate [1]
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

        real p1 = pls_data->temp[(species>0)*pls_data->n_rho + i_rho];
        real p2 = pls_data->temp[(species>0)*pls_data->n_rho + i_rho+1];
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
 * @param dens pointer to where evaluated density is stored [m^-3]
 * @param rho radial coordinate [1]
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
 * @param dens pointer to where interpolated densities are stored [m^-3]
 * @param temp pointer to where interpolated temperatures are stored [J]
 * @param rho radial coordinate [1]
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

/**
 * @brief Evalate plasma flow along the field lines
 *
 * @param vflow pointer where the flow value is stored [m/s]
 * @param rho particle rho coordinate [1]
 * @param r particle R coordinate [m]
 * @param pls_data pointer to plasma data
 */
a5err plasma_1D_eval_flow(real* vflow, real rho, real r,
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
        *vflow = pls_data->vtor[i_rho];
    }
    *vflow *= r;
    return err;
}
