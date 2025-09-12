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
#include "plasma_1Dt.h"


/**
 * @brief Initialize 1Dt plasma data and check inputs
 *
 * @param data pointer to the data struct
 * @param nrho number of rho grid points in the data
 * @param nrho number of time grid points in the data
 * @param nion number of ion species
 * @param rho rho grid points in which data is tabulated [1]
 * @param time time grid points in which data is tabulated [s]
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
int plasma_1Dt_init(plasma_1Dt_data* data, int nrho, int ntime, int nion,
                    real* rho, real* time, int* anum, int* znum, real* mass,
                    real* charge, real* Te, real* Ti, real* ne, real* ni,
                    real* vtor) {

    data->n_rho = nrho;
    data->n_time = ntime;
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
    for(int i = 0; i < nrho; i++) {
        data->rho[i] = rho[i];
    }
    data->time = (real*) malloc( ntime*sizeof(real) );
    for(int i = 0; i < ntime; i++) {
        data->time[i] = time[i];
    }
    data->vtor = (real*) malloc( nrho*ntime*sizeof(real) );
    data->temp = (real*) malloc( 2*nrho*ntime*sizeof(real) );
    data->dens = (real*) malloc( (nion+1)*nrho*ntime*sizeof(real) );
    for(int i = 0; i < nrho; i++) {
        for(int j = 0; j < ntime; j++) {
            data->vtor[j*nrho + i] = vtor[j*nrho + i];
            data->temp[j*2*nrho + i] = Te[j*nrho + i];
            data->temp[(j*2+1)*nrho + i] = Ti[j*nrho + i];
            data->dens[j*nrho + i] = ne[j*nrho + i];
            for(int k = 0; k < nion; k++) {
                data->dens[(k+1)*nrho*ntime + j*nrho + i] =
                    ni[k*nrho*ntime + j*nrho + i];
            }
        }
    }
    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void plasma_1Dt_free(plasma_1Dt_data* data) {
    free(data->mass);
    free(data->charge);
    free(data->anum);
    free(data->znum);
    free(data->rho);
    free(data->time);
    free(data->temp);
    free(data->dens);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void plasma_1Dt_offload(plasma_1Dt_data* data) {
    GPU_MAP_TO_DEVICE(
        data->mass[0:data->n_species], data->charge[0:data->n_species], \
        data->anum[0:data->n_species-1], data->znum[0:data->n_species-1], \
        data->rho[0:data->n_rho], data->time[0:data->n_time], \
        data->vtor[0:data->n_time*data->n_rho], \
        data->temp[0:data->n_time*data->n_rho*2], \
        data->dens[0:data->n_rho*data->n_species*data->n_time]
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
 * @param t time instant [s]
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
 * @param dens pointer to where evaluated density is stored [m^-3]
 * @param rho radial coordinate [1]
 * @param t time instant [s]
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
 * @param dens pointer to where interpolated densities are stored [m^-3]
 * @param temp pointer to where interpolated temperatures are stored [J]
 * @param rho radial coordinate [1]
 * @param t time instant [s]
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
 * @brief Evalate plasma flow along the field lines
 *
 * @param vflow pointer where the flow value is stored [m/s]
 * @param rho particle rho coordinate [1]
 * @param t particle time coordinate [s]
 * @param r particle R coordinate [m]
 * @param pls_data pointer to plasma data
 */
a5err plasma_1Dt_eval_flow(real* vflow, real rho, real t, real r,
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

        real p11, p12, p21, p22, p1, p2;

        p11 = pls_data->vtor[i_time*pls_data->n_rho + i_rho];
        p12 = pls_data->vtor[(i_time+1)*pls_data->n_rho + i_rho];
        p21 = pls_data->vtor[i_time*pls_data->n_rho + i_rho + 1];
        p22 = pls_data->vtor[(i_time+1)*pls_data->n_rho + i_rho + 1];

        p1 = p11 + t_rho * (p12 - p11);
        p2 = p21 + t_rho * (p22 - p21);

        *vflow = p1 + t_time * (p2 - p1);
    }
    *vflow *= r;
    return err;
}
