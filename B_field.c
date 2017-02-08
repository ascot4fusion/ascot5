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
#include "B_TC.h"

void B_field_init_offload(B_field_offload_data* offload_data,
                          real** offload_array) {
    FILE* f = fopen("input.magn_bkg", "r");

    if(f == NULL) {
        FILE* f = fopen("input.h5", "r");
        if(f == NULL) {
            /* no magnetic field input, use analytic field */
            offload_data->type = B_field_type_GS;

	    /* trivial cartesian field is only for debugging */
	    /* offload_data->type = B_field_type_TC; */
        } else {
            /* assuming input.h5 includes stellarator bfield */
            offload_data->type = B_field_type_ST;
        }
    } else {
        /* 2D if number of sectors 0 */
        int nsector;
        fscanf(f, "%*d %d", &nsector);
        
        if(nsector == 0) {
            offload_data->type = B_field_type_2D;
        } else {
            offload_data->type = B_field_type_3D;
        }
    }

    switch(offload_data->type) {
        case B_field_type_GS:
        B_GS_init_offload(&(offload_data->BGS), offload_array);
        offload_data->offload_array_length = offload_data->BGS.offload_array_length;
        break;

        case B_field_type_2D:
        B_2D_init_offload(&(offload_data->B2D), offload_array);
        offload_data->offload_array_length = offload_data->B2D.offload_array_length;
        break;

        case B_field_type_3D:
        B_3D_init_offload(&(offload_data->B3D), offload_array);
        offload_data->offload_array_length = offload_data->B3D.offload_array_length;
        break;

        case B_field_type_ST:
        B_ST_init_offload(&(offload_data->BST), offload_array);
        offload_data->offload_array_length = offload_data->BST.offload_array_length;
        break;

	case B_field_type_TC:
        B_TC_init_offload(&(offload_data->BTC), offload_array);
        offload_data->offload_array_length = offload_data->BTC.offload_array_length;
        break;
    }
}

void B_field_free_offload(B_field_offload_data* offload_data,
                          real** offload_array) {
    switch(offload_data->type) {
        case B_field_type_GS:
        B_GS_free_offload(&(offload_data->BGS), offload_array);
        break;

        case B_field_type_2D:
        B_2D_free_offload(&(offload_data->B2D), offload_array);
        break;

        case B_field_type_3D:
        B_3D_free_offload(&(offload_data->B3D), offload_array);
        break;

        case B_field_type_ST:
        B_ST_free_offload(&(offload_data->BST), offload_array);
        break;

	case B_field_type_TC:
        B_TC_free_offload(&(offload_data->BTC), offload_array);
        break;
    }
}

void B_field_init(B_field_data* Bdata, B_field_offload_data* offload_data,
                  real* offload_array) {
    switch(offload_data->type) {
        case B_field_type_GS:
        B_GS_init(&(Bdata->BGS), &(offload_data->BGS), offload_array);
        break;

        case B_field_type_2D:
        B_2D_init(&(Bdata->B2D), &(offload_data->B2D), offload_array);
        break;

        case B_field_type_3D:
        B_3D_init(&(Bdata->B3D), &(offload_data->B3D), offload_array);
        break;

        case B_field_type_ST:
        B_ST_init(&(Bdata->BST), &(offload_data->BST), offload_array);
        break;

	case B_field_type_TC:
        B_TC_init(&(Bdata->BTC), &(offload_data->BTC), offload_array);
        break;
    }
    Bdata->type = offload_data->type;
}

void B_field_eval_B(real B[], real r, real phi, real z, B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_B(B, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_eval_B(B, r, phi, z, &(Bdata->B2D));
        break;

        case B_field_type_3D:
        B_3D_eval_B(B, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_ST:
        B_ST_eval_B(B, r, phi, z, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_B(B, r, phi, z, &(Bdata->BTC));
        break;
    }
}

void B_field_eval_psi(real psi[], real r, real phi, real z,
                      B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_psi(psi, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_eval_psi(psi, r, phi, z, &(Bdata->B2D));
        break;

        case B_field_type_3D:
        B_3D_eval_psi(psi, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_ST:
        B_ST_eval_psi(psi, r, phi, z, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_psi(psi, r, phi, z, &(Bdata->BTC));
        break;
    }
}

void B_field_eval_rho(real rho[], real psi, B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_rho(rho, psi, &(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_eval_rho(rho, psi, &(Bdata->B2D));
        break;

        case B_field_type_3D:
        B_3D_eval_rho(rho, psi, &(Bdata->B3D));
        break;

        case B_field_type_ST:
        B_ST_eval_rho(rho, psi, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_rho(rho, psi, &(Bdata->BTC));
        break;
    }
}


/* Br   -> B[0]    dBr/dr -> B[1]   dBr/dphi -> B[2]    dBr/dz -> B[3]
 * Bphi -> B[4]  dBphi/dr -> B[5] dBphi/dphi -> B[6]  dBphi/dz -> B[7]
 * Bz   -> B[8]    dBz/dr -> B[9]   dBz/dphi -> B[10]   dBz/dz -> B[11]
 */
void B_field_eval_B_dB(real B_dB[], real r, real phi, real z,
                       B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_eval_B_dB(B_dB, r, phi, z, &(Bdata->B2D));
        break;

        case B_field_type_3D:
        B_3D_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_ST:
        B_ST_eval_B_dB(B_dB, r, phi, z, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_B_dB(B_dB, r, phi, z, &(Bdata->BTC));
        break;
    }
}

real B_field_get_axis_r(B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_get_axis_r(&(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_get_axis_r(&(Bdata->B2D));
        break;

        case B_field_type_3D:
        B_3D_get_axis_r(&(Bdata->B3D));
        break;

        case B_field_type_ST:
        B_ST_get_axis_r(&(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_get_axis_r(&(Bdata->BTC));
        break;
    }
}

real B_field_get_axis_z(B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_get_axis_z(&(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_get_axis_z(&(Bdata->B2D));
        break;

        case B_field_type_3D:
        B_3D_get_axis_z(&(Bdata->B3D));
        break;

        case B_field_type_ST:
        B_ST_get_axis_z(&(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_get_axis_z(&(Bdata->BTC));
        break;
    }
}
