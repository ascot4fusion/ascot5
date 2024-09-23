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
 * @param data pointer to the data struct
 *
 * @return zero if initialization succeeded
 */
int plasma_1DS_init(plasma_1DS_data* data, int nrho, real rhomin, real rhomax,
                    int nion, int* anum, int* znum, real* mass, real* charge,
                    real* Te, real* Ti, real* ne, real* ni) {
    int err = 0;
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

    real* Te_scaled = (real*) malloc( nrho*sizeof(real) );
    real* Ti_scaled = (real*) malloc( nrho*sizeof(real) );
    real* ne_scaled = (real*) malloc( nrho*sizeof(real) );
    real* ni_scaled = (real*) malloc( nion*nrho*sizeof(real) );
    for(int i=0; i < nrho; i++) {
#if PLASMA_1DS_NONEG == PLASMA_1DS_LOG
        Te_scaled[i] = log(Te[i]);
        Ti_scaled[i] = log(Ti[i]);
        ne_scaled[i] = log(ne[i]);
        for(int j = 0; j < nion; j++) {
            ni_scaled[j*nrho + i] = log(ni[j*nrho + i]);
        }
#elif PLASMA_1DS_NONEG == PLASMA_1DS_SQRT
        Te_scaled[i] = sqrt(Te[i]);
        Ti_scaled[i] = sqrt(Ti[i]);
        ne_scaled[i] = sqrt(ne[i]);
        for(int j = 0; j < nion; j++) {
            ni_scaled[j*nrho + i] = sqrt(ni[j*nrho + i]);
        }
#endif
    }
    err = interp1Dcomp_setup(&data->temp[0], Te_scaled, nrho, NATURALBC,
                             rhomin, rhomax);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        free(Te_scaled);
        free(Ti_scaled);
        free(ne_scaled);
        free(ni_scaled);
        return 1;
    }
    err = interp1Dcomp_setup(&data->temp[1], Ti_scaled, nrho, NATURALBC,
                             rhomin, rhomax);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        free(Te_scaled);
        free(Ti_scaled);
        free(ne_scaled);
        free(ni_scaled);
        return 1;
    }

    data->dens = (interp1D_data*) malloc(data->n_species*sizeof(interp1D_data));
    err = interp1Dcomp_setup(&data->dens[0], ne_scaled, nrho, NATURALBC,
                             rhomin, rhomax);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        free(Te_scaled);
        free(Ti_scaled);
        free(ne_scaled);
        free(ni_scaled);
        free(data->dens);
        return 1;
    }
    for(int i = 0; i < nion; i++) {
        err = interp1Dcomp_setup(&data->dens[i+1], &ni_scaled[i*nrho],
                                 nrho, NATURALBC, rhomin, rhomax);
        if(err) {
            print_err("Error: Failed to initialize splines.\n");
            free(Te_scaled);
            free(Ti_scaled);
            free(ne_scaled);
            free(ni_scaled);
            free(data->dens);
            return 1;
        }
    }
    free(Te_scaled);
    free(Ti_scaled);
    free(ne_scaled);
    free(ni_scaled);

    print_out(VERBOSE_IO, "\n1D plasma profiles (P_1DS)\n");
    print_out(VERBOSE_IO,
              "Min rho = %1.2le, Max rho = %1.2le,"
              " Number of rho grid points = %d,"
              " Number of ion species = %d\n",
              rhomin, rhomax, nrho, nion);
    print_out(VERBOSE_IO,
              "Species Z/A  charge [e]/mass [amu] Density [m^-3] at Min/Max rho"
              "    Temperature [eV] at Min/Max rho\n");
    real T0, T1, n0, n1;
    for(int i=0; i < nion; i++) {
        plasma_1DS_eval_temp(&T0, rhomin, i+1, data);
        plasma_1DS_eval_temp(&T1, rhomax, i+1, data);
        plasma_1DS_eval_dens(&n0, rhomin, i+1, data);
        plasma_1DS_eval_dens(&n1, rhomax, i+1, data);
        print_out(VERBOSE_IO,
                  " %3d  /%3d   %3d  /%7.3f             %1.2le/%1.2le     "
                  "           %1.2le/%1.2le       \n",
                  data->znum[i], data->anum[i],
                  (int)round(data->charge[i+1]/CONST_E),
                  data->mass[i+1]/CONST_U,
                  n0, n1, T0 / CONST_E, T1 / CONST_E);
    }

    plasma_1DS_eval_temp(&T0, rhomin, 0, data);
    plasma_1DS_eval_temp(&T1, rhomax, 0, data);
    plasma_1DS_eval_dens(&n0, rhomin, 0, data);
    plasma_1DS_eval_dens(&n1, rhomax, 0, data);
    print_out(VERBOSE_IO,
              "[electrons]  %3d  /%7.3f             %1.2le/%1.2le          "
              "      %1.2le/%1.2le       \n", -1, CONST_M_E/CONST_U,
              n0, n1, T0 / CONST_E, T1 / CONST_E);
    real quasineutrality = 0;
    for(int k = 0; k < nrho; k++) {
        real rho = rhomin + k * (rhomax - rhomin) / (nrho - 1);
        real ele_qdens;
        plasma_1DS_eval_dens(&ele_qdens, rho, 0, data);
        real ion_qdens = 0;
        for(int i=0; i < nion; i++) {
            plasma_1DS_eval_dens(&n0, rho, i+1, data);
            ion_qdens += n0;
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
void plasma_1DS_free(plasma_1DS_data* data) {
    free(data->mass);
    free(data->charge);
    free(data->anum);
    free(data->znum);
    for(int i = 0; i < data->n_species; i++) {
        free(data->dens[i].c);
    }
    free(data->dens);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void plasma_1DS_offload(plasma_1DS_data* data) {
    GPU_MAP_TO_DEVICE(
        data->mass[0:data->n_species], data->charge[0:data->n_species], \
        data->anum[0:data->n_species-1], data->znum[0:data->n_species-1], \
        data->temp[0:2], data->dens[0:data->n_species], \
        data->temp[0].c[0:data->temp[0].n_x*NSIZE_COMP1D], \
        data->temp[1].c[0:data->temp[1].n_x*NSIZE_COMP1D]
    )
    for (int i = 0; i < data->n_species; i++) {
        GPU_MAP_TO_DEVICE( data->dens[i].c[0:data->dens[i].n_x*NSIZE_COMP1D] )
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

#if PLASMA_1DS_NONEG == PLASMA_1DS_LOG
    *temp = exp(*temp);
#elif PLASMA_1DS_NONEG == PLASMA_1DS_SQRT
    *temp = (*temp) * (*temp);
#endif
    if(!err && *temp < 0){
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

#if PLASMA_1DS_NONEG == PLASMA_1DS_LOG
    *dens = exp(*dens);
#elif PLASMA_1DS_NONEG == PLASMA_1DS_SQRT
    *dens = (*dens) * (*dens);
#endif
    if(!err && *dens < 0){
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

#if PLASMA_1DS_NONEG == PLASMA_1DS_LOG
    for(int i=0; i<plasma_data->n_species; i++) {
        dens[i] = exp(dens[i]);
        temp[i] = exp(temp[i]);
    }
#elif PLASMA_1DS_NONEG == PLASMA_1DS_SQRT
    for(int i=0; i<plasma_data->n_species; i++) {
        dens[i] = dens[i]*dens[i];
        temp[i] = temp[i]*temp[i];
    }
#endif
    for(int i=0; i<plasma_data->n_species; i++) {
        if(!err && (dens[i] < 0 || temp[i] < 0) ) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__,
                                EF_PLASMA_1DS );
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
 * @param phi particle toroidal coordinate [rad]
 * @param z particle z coordinate [m]
 * @param t particle time coordinate [s]
 * @param pls_data pointer to plasma data
 */
a5err plasma_1DS_eval_flow(real* vflow, real rho, plasma_1DS_data* pls_data) {
    a5err err = 0;

    return err;
}
