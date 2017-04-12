/**
 * @author Joona Kontula joona.kontula@aalto.fi
 * @file E_1D.c
 * @brief 1D electric field evaluation functions
 *
 * Offloading is performed with a method similar to the magnetic field offload.
 * @see B_field.h
 */
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "E_1D.h"

/**
 * @brief Load 1D electric field data and prepare parameters
 *
 * This function reads the 1D electric field data from input.erad, fills the
 * offload struct with parameters and allocates and fills the offload array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @todo read data from input.h5 file
 */
void E_1D_init_offload(E_1D_offload_data* offload_data, real** offload_array) {

    FILE* f = fopen("input.erad", "r");

    /* Skip comment line */
    while(fgetc(f) != '\n');

    /* Number of rho points */
    int n_rho;
    fscanf(f, "%d", &n_rho);

    /* Allocate n_rho space for both rho and erad */
    offload_data->offload_array_length = 2*n_rho;
    *offload_array = (real*) malloc(2*n_rho*sizeof(real));

    /* Pointers to beginning of different data series to make code more
     * readable */
    real* rho = &(*offload_array)[0];
    real* dV = &(*offload_array)[n_rho];

    /* For data in format dV/rho, we can ignore effective minor radius */
    real a = 1.0;
    
    /* Read actual data into array */
    for(int i = 0; i < n_rho; i++) {
        fscanf(f, "%lf %lf", &rho[i], &dV[i]);
        /* Scale derivatives by effective minor radius */
        dV[i] = a * dV[i];
    }
    
    fclose(f);
}

/**
 * @brief Free offload array and reset parameters
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void E_1D_free_offload(E_1D_offload_data* offload_data,
                       real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize 1D electric field data struct on target
 *
 * This function copies the 1D electric field parameters from the offload struct
 * to the struct on target and sets the 1D electric field data pointers to
 * correct offsets in the offload array.
 *
 * @param Edata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void E_1D_init(E_1D_data* Edata,
               E_1D_offload_data* offload_data,
               real* offload_array) {
    Edata->n_rho = offload_data->n_rho;
    Edata->rho = &offload_array[0];
    Edata->dV = &offload_array[Edata->n_rho];
}

/**
 * @brief Evaluate 1D radial electric field
 *
 * This function evaluates the 1D potential gradient of the plasma at the given
 * radial coordinate using linear interpolation, and then calculates the
 * radial electric field by multiplying that with the rho-gradient.
 *
 * @param E array where the electric field will be stored (E_r -> E[1],
 *        E_phi -> E[1], E_z -> E[2])
 * @param rho_drho array where rho and components of grad-rho a are stored
 *          (rho -> rho_drho[0], drho/dr -> rho_drho[1], drho/dphi -> rho_drho[2],
 *          drho/dz -> rho_drho[3])
 * @param Edata pointer to electric field data
 */
void E_1D_eval_E(real E[], real rho_drho[], E_1D_data* Edata) {
    /* As the erad data may be provided at irregular intervals, we must
     * search for the correct grid index */
    /** @todo Implement a more efficient search algorithm */
    int i_rho = 0;
    while(i_rho < Edata->n_rho && Edata->rho[i_rho] <= rho_drho[0]) {
        i_rho++;
    }
    i_rho--;
    
    real t_rho = (rho_drho[0] - Edata->rho[i_rho])
        / (Edata->rho[i_rho+1] - Edata->rho[i_rho]);

    real p1 = Edata->dV[i_rho];
    real p2 = Edata->dV[i_rho+1];
    real dV = p1 + t_rho * (p2 - p1);

    E[0] = dV * rho_drho[1];
    E[1] = dV * rho_drho[2];
    E[2] = dV * rho_drho[3];
}
