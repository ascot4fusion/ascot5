/**
 * @file B_field.c
 * @brief Magnetic field interface
 */
#include <stdio.h>
#include "B_field.h"
#include "B_GS.h"
#include "B_2D.h"
#include "B_3D.h"
#include "B_ST.h"

void B_field_init_offload(B_field_offload_data* offload_data,
                          real** offload_array) {
    FILE* f = fopen("input.magn_bkg", "r");

    if(f == NULL) {
        FILE* f = fopen("input.h5", "r");
        if(f == NULL) {
            /* no magnetic field input, use analytic field */
            offload_data->type = 1;
        } else {
            /* assuming input.h5 includes stellarator bfield */
            offload_data->type = 5;
        }
    } else {
        /* 2D if number of sectors 0 */
        int nsector;
        fscanf(f, "%*d %d", &nsector);
        
        if(nsector == 0) {
            offload_data->type = 3;
        } else {
            offload_data->type = 4;
        }
    }

    switch(offload_data->type) {
        case 1:
        B_GS_init_offload(&(offload_data->BGS), offload_array);
        offload_data->offload_array_length = offload_data->BGS.offload_array_length;
        break;

        case 3:
        B_2D_init_offload(&(offload_data->B2D), offload_array);
        offload_data->offload_array_length = offload_data->B2D.offload_array_length;
        break;

        case 4:
        B_3D_init_offload(&(offload_data->B3D), offload_array);
        offload_data->offload_array_length = offload_data->B3D.offload_array_length;
        break;

        case 5:
        B_ST_init_offload(&(offload_data->BST), offload_array);
        offload_data->offload_array_length = offload_data->BST.offload_array_length;
        break;
    }
}

void B_field_free_offload(B_field_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {
        case 1:
        B_GS_free_offload(&(offload_data->BGS), offload_array);
        break;

        case 3:
        B_2D_free_offload(&(offload_data->B2D), offload_array);
        break;

        case 4:
        B_3D_free_offload(&(offload_data->B3D), offload_array);
        break;

        case 5:
        B_ST_free_offload(&(offload_data->BST), offload_array);
        break;
    }
}

void B_field_init(B_field_data* Bdata, B_field_offload_data* offload_data,
                  real* offload_array) {
    switch(offload_data->type) {
        case 1:
        B_GS_init(&(Bdata->BGS), &(offload_data->BGS), offload_array);
        break;

        case 3:
        B_2D_init(&(Bdata->B2D), &(offload_data->B2D), offload_array);
        break;

        case 4:
        B_3D_init(&(Bdata->B3D), &(offload_data->B3D), offload_array);
        break;

        case 5:
        B_ST_init(&(Bdata->BST), &(offload_data->BST), offload_array);
        break;
    }
    Bdata->type = offload_data->type;
}

void B_field_eval_B(real B[], real r, real phi, real z, B_field_data* Bdata) {
    switch(Bdata->type) {
        case 1:
        B_GS_eval_B(B, r, phi, z, &(Bdata->BGS));
        break;

        case 3:
        B_2D_eval_B(B, r, phi, z, &(Bdata->B2D));
        break;

        case 4:
        B_3D_eval_B(B, r, phi, z, &(Bdata->B3D));
        break;

        case 5:
        B_ST_eval_B(B, r, phi, z, &(Bdata->BST));
        break;
    }
}

void B_field_eval_psi(real psi[], real r, real phi, real z,
                      B_field_data* Bdata) {
    switch(Bdata->type) {
        case 1:
        B_GS_eval_psi(psi, r, phi, z, &(Bdata->BGS));
        break;

        case 3:
        B_2D_eval_psi(psi, r, phi, z, &(Bdata->B2D));
        break;

        case 4:
        B_3D_eval_psi(psi, r, phi, z, &(Bdata->B3D));
        break;

        case 5:
        B_ST_eval_psi(psi, r, phi, z, &(Bdata->BST));
        break;
    }
}

void B_field_eval_rho(real rho[], real psi, B_field_data* Bdata) {
    switch(Bdata->type) {
        case 1:
        B_GS_eval_rho(rho, psi, &(Bdata->BGS));
        break;

        case 3:
        B_2D_eval_rho(rho, psi, &(Bdata->B2D));
        break;

        case 4:
        B_3D_eval_rho(rho, psi, &(Bdata->B3D));
        break;

        case 5:
        B_ST_eval_rho(rho, psi, &(Bdata->BST));
        break;
    }
}

void B_field_eval_B_dB(real B_dB[], real r, real phi, real z,
                       B_field_data* Bdata) {
    switch(Bdata->type) {
        case 1:
        B_GS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BGS));
        break;

        case 3:
        B_2D_eval_B_dB(B_dB, r, phi, z, &(Bdata->B2D));
        break;

        case 4:
        B_3D_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3D));
        break;

        case 5:
        B_ST_eval_B_dB(B_dB, r, phi, z, &(Bdata->BST));
        break;
    }
}

real B_field_get_axis_r(B_field_data* Bdata) {
    switch(Bdata->type) {
        case 1:
        B_GS_get_axis_r(&(Bdata->BGS));
        break;

        case 3:
        B_2D_get_axis_r(&(Bdata->B2D));
        break;

        case 4:
        B_3D_get_axis_r(&(Bdata->B3D));
        break;

        case 5:
        B_ST_get_axis_r(&(Bdata->BST));
        break;
    }
}

real B_field_get_axis_z(B_field_data* Bdata) {
    switch(Bdata->type) {
        case 1:
        B_GS_get_axis_z(&(Bdata->BGS));
        break;

        case 3:
        B_2D_get_axis_z(&(Bdata->B2D));
        break;

        case 4:
        B_3D_get_axis_z(&(Bdata->B3D));
        break;

        case 5:
        B_ST_get_axis_z(&(Bdata->BST));
        break;
    }
}
