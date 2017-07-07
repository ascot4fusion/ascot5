/**
 * @file B_field.c
 * @brief Magnetic field interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "B_field.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2D.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3D.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_ST.h"
#include "Bfield/B_TC.h"

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
            offload_data->type = B_field_type_2DS;
        } else {
            offload_data->type = B_field_type_3DS;
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

        case B_field_type_2DS:
        B_2DS_init_offload(&(offload_data->B2DS), offload_array);
        offload_data->offload_array_length = offload_data->B2DS.offload_array_length;
        break;

        case B_field_type_3D:
        B_3D_init_offload(&(offload_data->B3D), offload_array);
        offload_data->offload_array_length = offload_data->B3D.offload_array_length;
        break;

        case B_field_type_3DS:
        B_3DS_init_offload(&(offload_data->B3DS), offload_array);
        offload_data->offload_array_length = offload_data->B3DS.offload_array_length;
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

        case B_field_type_2DS:
        B_2DS_free_offload(&(offload_data->B2DS), offload_array);
        break;

        case B_field_type_3D:
        B_3D_free_offload(&(offload_data->B3D), offload_array);
        break;

        case B_field_type_3DS:
        B_3DS_free_offload(&(offload_data->B3DS), offload_array);
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

        case B_field_type_2DS:
        B_2DS_init(&(Bdata->B2DS), &(offload_data->B2DS), offload_array);
        break;

        case B_field_type_3D:
        B_3D_init(&(Bdata->B3D), &(offload_data->B3D), offload_array);
        break;

        case B_field_type_3DS:
        B_3DS_init(&(Bdata->B3DS), &(offload_data->B3DS), offload_array);
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

        case B_field_type_2DS:
        B_2DS_eval_B(B, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3D:
        B_3D_eval_B(B, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_3DS:
        B_3DS_eval_B(B, r, phi, z, &(Bdata->B3DS));
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

        case B_field_type_2DS:
        B_2DS_eval_psi(psi, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3D:
        B_3D_eval_psi(psi, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_3DS:
        B_3DS_eval_psi(psi, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_ST:
        B_ST_eval_psi(psi, r, phi, z, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_psi(psi, r, phi, z, &(Bdata->BTC));
        break;
    }
}

/* psi   -> psi_dpsi[0]    dpsi/dr -> psi_dpsi[1]
 * dpsi/dphi -> psi_dpsi[2]    dpsi/dz -> psi_dpsi[3]
 */
void B_field_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z,
                       B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B2D));
        break;

        case B_field_type_2DS:
        B_2DS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3D:
        B_3D_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_3DS:
        B_3DS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_ST:
        B_ST_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BTC));
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

        case B_field_type_2DS:
        B_2DS_eval_rho(rho, psi, &(Bdata->B2DS));
        break;

        case B_field_type_3D:
        B_3D_eval_rho(rho, psi, &(Bdata->B3D));
        break;

        case B_field_type_3DS:
        B_3DS_eval_rho(rho, psi, &(Bdata->B3DS));
        break;

        case B_field_type_ST:
        B_ST_eval_rho(rho, psi, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_rho(rho, psi, &(Bdata->BTC));
        break;
    }
}

/* rho   -> rho_drho[0]    drho/dr -> rho_drho[1]
 * drho/dphi -> rho_drho[2]    drho/dz -> rho_drho[3]
 */
void B_field_eval_rho_drho(real rho_drho[], real r, real phi, real z,
                       B_field_data* Bdata) {
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2D:
        B_2D_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B2D));
        break;

        case B_field_type_2DS:
        B_2DS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3D:
        B_3D_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_3DS:
        B_3DS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_ST:
        B_ST_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BST));
        break;

	case B_field_type_TC:
        B_TC_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BTC));
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

        case B_field_type_2DS:
        B_2DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3D:
        B_3D_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3D));
        break;

        case B_field_type_3DS:
        B_3DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3DS));
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
    real axis_r = 0;
    switch(Bdata->type) {
        case B_field_type_GS:
        axis_r = B_GS_get_axis_r(&(Bdata->BGS));
        break;

        case B_field_type_2D:
        axis_r = B_2D_get_axis_r(&(Bdata->B2D));
        break;

        case B_field_type_2DS:
        axis_r = B_2DS_get_axis_r(&(Bdata->B2DS));
        break;

        case B_field_type_3D:
        axis_r = B_3D_get_axis_r(&(Bdata->B3D));
        break;

        case B_field_type_3DS:
        axis_r = B_3DS_get_axis_r(&(Bdata->B3DS));
        break;

        case B_field_type_ST:
        axis_r = B_ST_get_axis_r(&(Bdata->BST));
        break;

	case B_field_type_TC:
        axis_r = B_TC_get_axis_r(&(Bdata->BTC));
        break;
    }
    return axis_r;
}

real B_field_get_axis_z(B_field_data* Bdata) {
    real axis_z = 0;
    switch(Bdata->type) {
        case B_field_type_GS:
        axis_z = B_GS_get_axis_z(&(Bdata->BGS));
        break;

        case B_field_type_2D:
        axis_z = B_2D_get_axis_z(&(Bdata->B2D));
        break;

        case B_field_type_2DS:
        axis_z = B_2DS_get_axis_z(&(Bdata->B2DS));
        break;

        case B_field_type_3D:
        axis_z = B_3D_get_axis_z(&(Bdata->B3D));
        break;

        case B_field_type_3DS:
        axis_z = B_3DS_get_axis_z(&(Bdata->B3DS));
        break;

        case B_field_type_ST:
        axis_z = B_ST_get_axis_z(&(Bdata->BST));
        break;

	case B_field_type_TC:
        axis_z = B_TC_get_axis_z(&(Bdata->BTC));
        break;
    }
    return axis_z;
}
