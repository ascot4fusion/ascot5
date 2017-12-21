/**
 * @file plasma_1DS.c
 * @brief 1D spline plasma evaluation functions
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "plasma_1DS.h"
#include "consts.h"
#include "spline/interp1Dcomp.h"

/**
 * @brief Free offload array and reset parameters 
 *
 *This function deallocates the offload_array.
 
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void plasma_1DS_free(plasma_1DS_data* plasma_data) {
    interp1Dcomp_free(plasma_data->dens);
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
a5err plasma_1DS_init(plasma_1DS_data* plasma_data,
                    plasma_1DS_offload_data* offload_data,
                    real* offload_array) {

    a5err err = 0;
    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(offload_array[0]);
    real* temp_e = &(offload_array[offload_data->n_rho]);
    real* temp_i = &(offload_array[offload_data->n_rho*2]);
    real* dens = &(offload_array[offload_data->n_rho
                                  + offload_data->n_rho*offload_data->n_species]);

    plasma_data->n_rho = offload_data->n_rho;
    plasma_data->rho_min = rho[0];
    plasma_data->rho_max = rho[plasma_data->n_rho - 1];
    plasma_data->rho_grid = (plasma_data->rho_max - plasma_data->rho_min)/(plasma_data->n_rho - 1);
    plasma_data->n_species = offload_data->n_species;
    int i;
    for(i = 0; i < plasma_data->n_species; i++) {
        plasma_data->mass[i] = offload_data->mass[i];
        plasma_data->charge[i] = offload_data->charge[i];
    }
    err += interp1Dcomp_init(&plasma_data->temp[0], temp_e,
			     plasma_data->n_rho, plasma_data->rho_min, plasma_data->rho_max,
                             plasma_data->rho_grid);
    err += interp1Dcomp_init(&plasma_data->temp[1], temp_i,
			     plasma_data->n_rho, plasma_data->rho_min, plasma_data->rho_max,
                             plasma_data->rho_grid);
    plasma_data->dens = (interp1D_data*) malloc(plasma_data->n_species*sizeof(interp1D_data));
    for(i = 0; i < plasma_data->n_species; i++) {
        err += interp1Dcomp_init(&plasma_data->dens[i],
                                 dens + i*plasma_data->n_rho,
                                 plasma_data->n_rho, plasma_data->rho_min, plasma_data->rho_max,
                                 plasma_data->rho_grid);
    }
    return err;
}

/**
 * @brief Evaluate plasma temperature
 *
 * This function evaluates the temperature of a plasma species at the given
 * radial coordinate using linear interpolation.
 *
 * @param temp temperature value will be stored in temp[0]
 * @param rho radial coordinate
 * @param species index of plasma species
 * @param plasma_data pointer to plasma data struct
 */
a5err plasma_1DS_eval_temp(real temp[], real rho, int species, plasma_1DS_data* plasma_data) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    if (rho < plasma_data->rho_min){
        rho = plasma_data->rho_min;
    } else if (rho > plasma_data->rho_max) {
        rho = plasma_data->rho_max;        
    }

    if (species == 0) {
        interperr += interp1Dcomp_eval_B(&temp[0], &plasma_data->temp[0], rho);
    } else {
        interperr += interp1Dcomp_eval_B(&temp[0], &plasma_data->temp[1], rho);
    }

    if(interperr) {err = error_raise( ERR_OUTSIDE_PLASMA, __LINE__ );}
    
    return err;
}

/**
 * @brief Evaluate plasma density
 *
 * @see plasma_1DS_eval_temp
 */
a5err plasma_1DS_eval_dens(real dens[], real rho, int species, plasma_1DS_data* plasma_data) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    if (rho < plasma_data->rho_min){
        rho = plasma_data->rho_min;
    } else if (rho > plasma_data->rho_max) {
        rho = plasma_data->rho_max;        
    }

    interperr += interp1Dcomp_eval_B(&dens[0], &plasma_data->dens[species], rho);

    if(interperr) {err = error_raise( ERR_OUTSIDE_PLASMA, __LINE__ );}
    
    return err;
}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 */
a5err plasma_1DS_eval_densandtemp(real* dens, real* temp, real rho, plasma_1DS_data* plasma_data) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    if (rho < plasma_data->rho_min){
        rho = plasma_data->rho_min;
    } else if (rho > plasma_data->rho_max) {
        rho = plasma_data->rho_max;        
    }

    /* Electrons */
    interperr += interp1Dcomp_eval_B(&temp[0], &plasma_data->temp[0], rho);
    interperr += interp1Dcomp_eval_B(&dens[0], &plasma_data->dens[0], rho);
    /* Ions, all have the same temperature */
    real temp_temp_i;
    interperr += interp1Dcomp_eval_B(&temp_temp_i, &plasma_data->temp[1], rho);
    for(int i = 1; i < plasma_data->n_species; i++) {
        temp[i] = temp_temp_i;
        interperr += interp1Dcomp_eval_B(&dens[i], &plasma_data->dens[i], rho);
    }
    
    if(interperr) {err = error_raise( ERR_OUTSIDE_PLASMA, __LINE__ );}

    return err;
}
