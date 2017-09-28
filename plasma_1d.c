/**
 * @file plasma_1d.c
 * @brief 1D plasma evaluation functions
 *
 * Offloading is performed with a method similar to the magnetic field offload.
 * @see B_field.h
 */
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "plasma_1d.h"
#include "consts.h"

/**
 * @brief Load 1D plasma data and prepare parameters
 *
 * This function reads the 1D plasma data from input.plasma_1d, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @todo Move reading the file to ascot4_interface
 */
void plasma_1d_init_offload(plasma_1d_offload_data* offload_data,
                            real** offload_array) {
    int i, j;

    FILE* f = fopen("input.plasma_1d", "r");

    /* Skip comment lines */
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');

    /* Number of ion species */
    int n_rho, n_ions;
    fscanf(f, "%d %d", &n_rho, &n_ions);
    while(fgetc(f) != '\n');

    /* One extra species for the electrons; they will always be species 0 */
    offload_data->n_species = n_ions + 1;
    offload_data->n_rho = n_rho;

    /* Allocate space for rho + temperature for each ion species and electrons
     * + density for each ion species and electrons */
    offload_data->offload_array_length = n_rho + 2*(n_ions+1)*n_rho;
    *offload_array = (real*) malloc(sizeof(real)
                                    * offload_data->offload_array_length);
    
    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* temp_e = &(*offload_array)[n_rho];
    real* temp_i = &(*offload_array)[n_rho*2];
    real* dens_e = &(*offload_array)[n_rho*2 + n_rho*n_ions];
    real* dens_i = &(*offload_array)[n_rho*2 + n_rho*n_ions + n_rho];

    /* Read Znum and calculate charge */
    /* Electron charge -1 */
    offload_data->charge[0] = -1 * CONST_E;
    for(i = 0; i < n_ions; i++) {
        int Znum;
        fscanf(f, "%d", &Znum);
        offload_data->charge[i+1] = Znum * CONST_E;
    }
    while(fgetc(f) != '\n');

    /* Read Amass and calculate mass */
    offload_data->mass[0] = CONST_M_E;
    for(i = 0; i < n_ions; i++) {
        int Amass;
        fscanf(f, "%d", &Amass);
        offload_data->mass[i+1] = Amass * CONST_U;
    }
    while(fgetc(f) != '\n');

    /* Skip collision mode and header line */
    while(fgetc(f) != '\n');
    while(fgetc(f) != '\n');

    /* Read actual data into array */
    for(i = 0; i < n_rho; i++) {
        fscanf(f, "%lf", &rho[i]);
        real temp_temp_e;
        fscanf(f, "%lf", &temp_temp_e);
        temp_e[i] = temp_temp_e * CONST_E / CONST_KB;
        fscanf(f, "%lf", &dens_e[i]);
        /* Don't need Vtor_I */
        fscanf(f, "%*f");
        /* All ions have same temperature */
        real temp_temp_i;
        fscanf(f, "%lf", &temp_temp_i);
        for(j = 0; j < n_ions; j++) {
            temp_i[j*n_rho + i] = temp_temp_i *CONST_E / CONST_KB;
            fscanf(f, "%lf", &dens_i[j*n_rho + i]);
        }
    }
    
    fclose(f);
}

/**
 * @brief Free offload array and reset parameters 
 *
 *This function deallocates the offload_array.

 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
*/
void plasma_1d_free_offload(plasma_1d_offload_data* offload_data,
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
int plasma_1d_init(plasma_1d_data* plasma_data,
                    plasma_1d_offload_data* offload_data,
                    real* offload_array) {

    int err = 0;

    plasma_data->n_rho = offload_data->n_rho;
    plasma_data->n_species = offload_data->n_species;
    int i;
    for(i = 0; i < plasma_data->n_species; i++) {
        plasma_data->mass[i] = offload_data->mass[i];
        plasma_data->charge[i] = offload_data->charge[i];
    }
    plasma_data->rho = &offload_array[0];
    plasma_data->temp = &offload_array[plasma_data->n_rho];
    plasma_data->dens = &offload_array[plasma_data->n_rho
                                  + plasma_data->n_rho*plasma_data->n_species];

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
real plasma_1d_eval_temp(real rho, int species, plasma_1d_data* plasma_data) {
    /* As the plasma data may be provided at irregular intervals, we must
     * search for the correct grid index */
    /** @todo Implement a more efficient search algorithm */
    
    real p = 0;
    if(rho < plasma_data->rho[0]) {
	p = plasma_data->temp[species*plasma_data->n_rho];
    }
    else if(rho >= plasma_data->rho[plasma_data->n_rho-1]) {
	p = plasma_data->temp[species*plasma_data->n_rho + plasma_data->n_rho - 1];
    }
    else {
	int i_rho = 0;
	while(i_rho < plasma_data->n_rho - 1 && plasma_data->rho[i_rho] <= rho) {
	    i_rho++;
	}
	i_rho--;
	real t_rho = (rho - plasma_data->rho[i_rho])
                 / (plasma_data->rho[i_rho+1] - plasma_data->rho[i_rho]);

	real p1 = plasma_data->temp[species*plasma_data->n_rho + i_rho];
	real p2 = plasma_data->temp[species*plasma_data->n_rho + i_rho+1];
	p = p1 + t_rho * (p2 - p1);
    }
    
    return p;
}

/**
 * @brief Evaluate plasma density
 *
 * @see plasma_1d_eval_temp
 */
real plasma_1d_eval_dens(real rho, int species, plasma_1d_data* plasma_data) {
    real p = 0;
    if(rho < plasma_data->rho[0]) {
	p = plasma_data->dens[species*plasma_data->n_rho];
    }
    else if(rho >= plasma_data->rho[plasma_data->n_rho-1]) {
	p = plasma_data->dens[species*plasma_data->n_rho + plasma_data->n_rho - 1];
    }
    else {
	int i_rho = 0;
	while(i_rho < plasma_data->n_rho - 1 && plasma_data->rho[i_rho] <= rho) {
	    i_rho++;
	}
	i_rho--;
	real t_rho = (rho - plasma_data->rho[i_rho])
                 / (plasma_data->rho[i_rho+1] - plasma_data->rho[i_rho]);

	real p1 = plasma_data->dens[species*plasma_data->n_rho + i_rho];
	real p2 = plasma_data->dens[species*plasma_data->n_rho + i_rho+1];
	p = p1 + t_rho * (p2 - p1);
    }
    
    return p;
}

/**
 * @brief Evaluate plasma density and temperature for all species
 *
 */
a5err plasma_1d_eval_densandtemp(real rho, plasma_1d_data* plasma_data, real* dens, real* temp) {
    a5err err = 0;
    real p1, p2;
    if(rho < plasma_data->rho[0]) {
	for(int i = 0; i < plasma_data->n_species; i++) {
	    dens[i] = plasma_data->dens[i*plasma_data->n_rho];

	    if(i < 2) {
		/* Electron and ion temperature */
		 temp[i] = plasma_data->temp[i*plasma_data->n_rho];
	    }
	    else {
		/* Temperature is same for all ion species */
		temp[i] = temp[1];
	    }
	}
    }
    else if(rho >= plasma_data->rho[plasma_data->n_rho-1]) {
	for(int i = 0; i < plasma_data->n_species; i++) {
	    dens[i] = plasma_data->dens[i*plasma_data->n_rho + plasma_data->n_rho - 1];

	    if(i < 2) {
		/* Electron and ion temperature */
		 temp[i] = plasma_data->temp[i*plasma_data->n_rho + plasma_data->n_rho - 1];
	    }
	    else {
		/* Temperature is same for all ion species */
		temp[i] = temp[1];
	    }
	}
    }
    else {
	int i_rho = 0;
	while(i_rho < plasma_data->n_rho-1 && plasma_data->rho[i_rho] <= rho) {
	    i_rho++;
	}
	i_rho--;

	real t_rho = (rho - plasma_data->rho[i_rho])
                 / (plasma_data->rho[i_rho+1] - plasma_data->rho[i_rho]);

	for(int i = 0; i < plasma_data->n_species; i++) {
	    p1 = plasma_data->dens[i*plasma_data->n_rho + i_rho];
	    p2 = plasma_data->dens[i*plasma_data->n_rho + i_rho+1];
	    dens[i] = p1 + t_rho * (p2 - p1);

	    if(i < 2) {
		/* Electron and ion temperature */
		p1 = plasma_data->temp[i*plasma_data->n_rho + i_rho];
		p2 = plasma_data->temp[i*plasma_data->n_rho + i_rho+1];
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
