/**
 * @file plasma.c
 * @brief Plasma interface
 */
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "plasma.h"
#include "plasma/plasma_1D.h"
#include "plasma/plasma_1DS.h"
#include "consts.h"

/**
 * @brief Load plasma data and prepare parameters
 *
 * The reading of the ASCOT4 plasma data has been moved into
 * ascot4_interface, so this function is now a dummy.
 *
 */
void plasma_init_offload(plasma_offload_data* offload_data,
                            real** offload_array) {
    // Dummy function
}


void plasma_free_offload(plasma_offload_data* offload_data,
                            real** offload_array) {
    switch(offload_data->type) {
        case plasma_type_1D:
	    plasma_1D_free_offload(&(offload_data->plasma_1D), offload_array);
	    break;
	    
        case plasma_type_1DS:
	    //plasma_1DS_free_offload(&(offload_data->plasma_1DS), offload_array);
	    break;
    }
}


int plasma_init(plasma_data* pls_data,
                    plasma_offload_data* offload_data,
                    real* offload_array) {
    int err = 0;
    switch(offload_data->type) {
        case plasma_type_1D:
	    plasma_1D_init(&(pls_data->plasma_1D), &(offload_data->plasma_1D), offload_array);
	    break;
	    
        case plasma_type_1DS:
	    //plasma_1DS_init(&(pls_data->plasma_1DS), &(offload_data->plasma_1DS), offload_array);
	    break;
    }
    pls_data->type = offload_data->type;
    
    return err;
}


real plasma_eval_temp(real rho, int species, plasma_data* pls_data) {
    real p = 0;
    
    switch(pls_data->type) {
        case plasma_type_1D:
	    p = plasma_1D_eval_temp(rho, species, &(pls_data->plasma_1D));
	    break;
	    
        case plasma_type_1DS:
	    //p = plasma_1DS_eval_temp(rho, species, &(pls_data->plasma_1DS));
	    break;
    }
    
    return p;
}

/**
 * @brief Evaluate plasma density
 *
 * @see plasma_1d_eval_temp
 */
real plasma_eval_dens(real rho, int species, plasma_data* pls_data) {
    real p = 0;
    switch(pls_data->type) {
        case plasma_type_1D:
	    p = plasma_1D_eval_dens(rho, species, &(pls_data->plasma_1D));
	    break;
	    
        case plasma_type_1DS:
	    //p = plasma_1DS_eval_dens(rho, species, &(pls_data->plasma_1DS));
	    break;
    }
    return p;
}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 */
a5err plasma_eval_densandtemp(real rho, plasma_data* pls_data, real* dens, real* temp) {
    a5err err = 0;
    
    switch(pls_data->type) {
        case plasma_type_1D:
	    plasma_1D_eval_densandtemp(rho, &(pls_data->plasma_1D), dens, temp);
	    break;
	    
        case plasma_type_1DS:
	    //plasma_1DS_eval_densandtemp(rho, &(pls_data->plasma_1DS), dens, temp);
	    break;
    }

    return err;
    
}

int plasma_get_n_species(plasma_data* pls_data) {
    int n = 0;
    switch(pls_data->type) {
        case plasma_type_1D:
	    n = plasma_1D_get_n_species(&(pls_data->plasma_1D));
	    break;
	    
        case plasma_type_1DS:
	    //n = plasma_1DS_get_n_species(&(pls_data->plasma_1DS));
	    break;
    }

    return n;
}


real* plasma_get_species_mass(plasma_data* pls_data) {
    real* mass = NULL;
    switch(pls_data->type) {
        case plasma_type_1D:
	    mass = plasma_1D_get_species_mass(&(pls_data->plasma_1D));
	    break;
	    
        case plasma_type_1DS:
	    //mass = plasma_1DS_get_species_mass(&(plasma_data->plasma_1DS));
	    break;
    }

    return mass;
}

real* plasma_get_species_charge(plasma_data* pls_data) {
    real* charge = NULL;
    switch(pls_data->type) {
        case plasma_type_1D:
	    charge = plasma_1D_get_species_charge(&(pls_data->plasma_1D));
	    break;
	    
        case plasma_type_1DS:
	    //charge = plasma_1DS_get_species_charge(&(plasma_data->plasma_1DS));
	    break;
    }

    return charge;
}
