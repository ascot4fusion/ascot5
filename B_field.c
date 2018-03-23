/**
 * @file B_field.c
 * @brief Magnetic field interface
 */
#include <stdio.h>
#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "Bfield/B_GS.h"
#include "Bfield/B_2DS.h"
#include "Bfield/B_3DS.h"
#include "Bfield/B_STS.h"
#include "Bfield/B_TC.h"

void B_field_init_offload(B_field_offload_data* offload_data,
                          real** offload_array) {
    
    switch(offload_data->type) {
        case B_field_type_GS:
        B_GS_init_offload(&(offload_data->BGS), offload_array);
        offload_data->offload_array_length = offload_data->BGS.offload_array_length;
        break;

        case B_field_type_2DS:
        B_2DS_init_offload(&(offload_data->B2DS), offload_array);
        offload_data->offload_array_length = offload_data->B2DS.offload_array_length;
        break;

        case B_field_type_3DS:
        B_3DS_init_offload(&(offload_data->B3DS), offload_array);
        offload_data->offload_array_length = offload_data->B3DS.offload_array_length;
        break;

        case B_field_type_STS:
        B_STS_init_offload(&(offload_data->BSTS), offload_array);
        offload_data->offload_array_length = offload_data->BSTS.offload_array_length;
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

        case B_field_type_2DS:
        B_2DS_free_offload(&(offload_data->B2DS), offload_array);
        break;

        case B_field_type_3DS:
        B_3DS_free_offload(&(offload_data->B3DS), offload_array);
        break;

        case B_field_type_STS:
        B_STS_free_offload(&(offload_data->BSTS), offload_array);
        break;

	case B_field_type_TC:
        B_TC_free_offload(&(offload_data->BTC), offload_array);
        break;
    }
}

int B_field_init(B_field_data* Bdata, B_field_offload_data* offload_data,
                  real* offload_array) {
    int err = 0;

    switch(offload_data->type) {
        case B_field_type_GS:
        B_GS_init(&(Bdata->BGS), &(offload_data->BGS), offload_array);
        break;

        case B_field_type_2DS:
	err = B_2DS_init(&(Bdata->B2DS), &(offload_data->B2DS), offload_array);
        break;

        case B_field_type_3DS:
        err = B_3DS_init(&(Bdata->B3DS), &(offload_data->B3DS), offload_array);
        break;

        case B_field_type_STS:
        B_STS_init(&(Bdata->BSTS), &(offload_data->BSTS), offload_array);
        break;

	case B_field_type_TC:
        B_TC_init(&(Bdata->BTC), &(offload_data->BTC), offload_array);
        break;
    }
    Bdata->type = offload_data->type;

    return err;
}

a5err B_field_eval_psi(real psi[], real r, real phi, real z,
                      B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_psi(psi, r, phi, z, &(Bdata->BGS));
        break;
	
        case B_field_type_2DS:
        err = B_2DS_eval_psi(psi, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        err = B_3DS_eval_psi(psi, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_psi(psi, r, phi, z, &(Bdata->BSTS));
        break;

	case B_field_type_TC:
        B_TC_eval_psi(psi, r, phi, z, &(Bdata->BTC));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	psi[0] = 1;
    }

    return err;
}


a5err B_field_eval_psi_SIMD(int i, real psi[NSIMD], real r, real phi, real z,
                      B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_psi_SIMD(i, psi, r, phi, z, &(Bdata->BGS));
        break;

	case B_field_type_2DS:
        B_2DS_eval_psi_SIMD(i, psi, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        B_3DS_eval_psi_SIMD(i, psi, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_psi_SIMD(i, psi, r, phi, z, &(Bdata->BSTS));
        break;

	default:
        // Do nothing
	break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	psi[i] = 1;
    }

    return err;
}

/* psi   -> psi_dpsi[0]    dpsi/dr -> psi_dpsi[1]
 * dpsi/dphi -> psi_dpsi[2]    dpsi/dz -> psi_dpsi[3]
 */
a5err B_field_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z,
                       B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2DS:
        err = B_2DS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        err = B_3DS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BSTS));
        break;

	case B_field_type_TC:
        B_TC_eval_psi_dpsi(psi_dpsi, r, phi, z, &(Bdata->BTC));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	psi_dpsi[0] = 1;
	for(int k=1; k<4; k++) {psi_dpsi[k] = 0;}
    }

    return err;
}

a5err B_field_eval_rho(real rho[], real psi, B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_rho(rho, psi, &(Bdata->BGS));
        break;

        case B_field_type_2DS:
        err = B_2DS_eval_rho(rho, psi, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        err = B_3DS_eval_rho(rho, psi, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_rho(rho, psi, &(Bdata->BSTS));
        break;

	case B_field_type_TC:
        B_TC_eval_rho(rho, psi, &(Bdata->BTC));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	rho[0] = 1;
    }

    return err;
}

a5err B_field_eval_rho_SIMD(int i, real rho[NSIMD], real psi, B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_rho_SIMD(i, rho, psi, &(Bdata->BGS));
        break;

        case B_field_type_2DS:
        B_2DS_eval_rho_SIMD(i, rho, psi, &(Bdata->B2DS));
        break;

	case B_field_type_3DS:
        B_3DS_eval_rho_SIMD(i, rho, psi, &(Bdata->B3DS));
        break;

    	case B_field_type_STS:
        B_STS_eval_rho_SIMD(i, rho, psi, &(Bdata->BSTS));
        break;

	default:
        // Do nothing
	break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	rho[i] = 1;
    }

    return err;
}

/* rho   -> rho_drho[0]    drho/dr -> rho_drho[1]
 * drho/dphi -> rho_drho[2]    drho/dz -> rho_drho[3]
 */
a5err B_field_eval_rho_drho(real rho_drho[], real r, real phi, real z,
                       B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2DS:
        err = B_2DS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        err = B_3DS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BSTS));
        break;

	case B_field_type_TC:
        B_TC_eval_rho_drho(rho_drho, r, phi, z, &(Bdata->BTC));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	rho_drho[0] = 1;
	for(int k=1; k<4; k++) {rho_drho[k] = 0;}
    }

    return err;
}

a5err B_field_eval_B(real B[], real r, real phi, real z, B_field_data* Bdata) {
    a5err err = 0;
    
    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_B(B, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2DS:
        err = B_2DS_eval_B(B, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        err = B_3DS_eval_B(B, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_B(B, r, phi, z, &(Bdata->BSTS));
        break;

	case B_field_type_TC:
        B_TC_eval_B(B, r, phi, z, &(Bdata->BTC));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	B[0] = 1;
	for(int k=1; k<3; k++) {B[k] = 0;}
    }

    return err;
}

/* Br   -> B[0]    dBr/dr -> B[1]   dBr/dphi -> B[2]    dBr/dz -> B[3]
 * Bphi -> B[4]  dBphi/dr -> B[5] dBphi/dphi -> B[6]  dBphi/dz -> B[7]
 * Bz   -> B[8]    dBz/dr -> B[9]   dBz/dphi -> B[10]   dBz/dz -> B[11]
 */
a5err B_field_eval_B_dB(real B_dB[], real r, real phi, real z,
			B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {
        case B_field_type_GS:
        B_GS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BGS));
        break;

        case B_field_type_2DS:
        err = B_2DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        err = B_3DS_eval_B_dB(B_dB, r, phi, z, &(Bdata->B3DS));
        break;

        case B_field_type_STS:
        B_STS_eval_B_dB(B_dB, r, phi, z, &(Bdata->BSTS));
        break;

	case B_field_type_TC:
        B_TC_eval_B_dB(B_dB, r, phi, z, &(Bdata->BTC));
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	B_dB[0] = 1;
	for(int k=1; k<12; k++) {B_dB[k] = 0;}
    }

    return err = 0;
}


/* Br   -> B[0]    dBr/dr -> B[1]   dBr/dphi -> B[2]    dBr/dz -> B[3]
 * Bphi -> B[4]  dBphi/dr -> B[5] dBphi/dphi -> B[6]  dBphi/dz -> B[7]
 * Bz   -> B[8]    dBz/dr -> B[9]   dBz/dphi -> B[10]   dBz/dz -> B[11]
 */
a5err B_field_eval_B_dB_SIMD(int i, real B_dB[12][NSIMD], real r, real phi, real z,
			B_field_data* Bdata) {
    a5err err = 0;

    switch(Bdata->type) {	
        case B_field_type_GS:
        B_GS_eval_B_dB_SIMD(i, B_dB, r, phi, z, &(Bdata->BGS));
        break;

	case B_field_type_2DS:
        B_2DS_eval_B_dB_SIMD(i, B_dB, r, phi, z, &(Bdata->B2DS));
        break;

	case B_field_type_3DS:
        B_3DS_eval_B_dB_SIMD(i, B_dB, r, phi, z, &(Bdata->B3DS));
        break;

	case B_field_type_STS:
        B_STS_eval_B_dB_SIMD(i, B_dB, r, phi, z, &(Bdata->BSTS));
        break;

        default:
        // Do nothing
	break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	B_dB[0][i] = 1;
	for(int k=1; k<12; k++) {B_dB[k][i] = 0;}
    }

    return err = 0;
}

real B_field_get_axis_r(B_field_data* Bdata) {
    real axis_r = 0;
    switch(Bdata->type) {
        case B_field_type_GS:
        axis_r = B_GS_get_axis_r(&(Bdata->BGS));
        break;
	
        case B_field_type_2DS:
        axis_r = B_2DS_get_axis_r(&(Bdata->B2DS));
        break;

        case B_field_type_3DS:
        axis_r = B_3DS_get_axis_r(&(Bdata->B3DS));
        break;

        case B_field_type_STS:
        axis_r = B_STS_get_axis_r(&(Bdata->BSTS));
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

        case B_field_type_2DS:
        axis_z = B_2DS_get_axis_z(&(Bdata->B2DS));
        break;
	
        case B_field_type_3DS:
        axis_z = B_3DS_get_axis_z(&(Bdata->B3DS));
        break;

        case B_field_type_STS:
        axis_z = B_STS_get_axis_z(&(Bdata->BSTS));
        break;

	case B_field_type_TC:
        axis_z = B_TC_get_axis_z(&(Bdata->BTC));
        break;
    }
    return axis_z;
}
