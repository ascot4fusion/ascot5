/**
 * @file B_STS.c
 * @brief Stellarator magnetic field with cubic spline interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "B_STS.h"
#include "../spline/interp2D.h" /* for 2D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */
#include "../spline/interp2Dcomp.h" 
#include "../spline/interp3Dcomp.h" 
#include "../spline/interp2Dexpl.h" 
#include "../spline/interp3Dexpl.h"

void B_STS_init_offload(B_STS_offload_data* offload_data, real** offload_array) {
    
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_STS_free_offload(B_STS_offload_data* offload_data, real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize magnetic field data struct on target 
 *
 * This function copies the magnetic field parameters from the offload struct
 * to the struct on target and sets the magnetic field data pointers to
 * correct offsets in the offload array.
 *
 * @param BData pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
int B_STS_init(B_STS_data* Bdata, B_STS_offload_data* offload_data,
               real* offload_array) {
    int err = 0;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;
    Bdata->periods = offload_data->periods;

    /* Spline initialization and storage. */
    
    #if INTERP_SPL_EXPL    
    err += interp3Dexpl_init(&Bdata->B_r, offload_array,
	offload_data->n_r, offload_data->n_phi, offload_data->n_z,
	offload_data->r_min, offload_data->r_max, offload_data->r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    
    err += interp3Dexpl_init(&Bdata->B_phi, offload_array
	+offload_data->n_phi*offload_data->n_z*offload_data->n_r,
	offload_data->n_r, offload_data->n_phi, offload_data->n_z,
	offload_data->r_min, offload_data->r_max, offload_data->r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    	
    err += interp3Dexpl_init(&Bdata->B_z, offload_array
	+2*offload_data->n_phi*offload_data->n_z*offload_data->n_r,
	offload_data->n_r, offload_data->n_phi, offload_data->n_z,
	offload_data->r_min, offload_data->r_max, offload_data->r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    
    err += interp3Dexpl_init(&Bdata->s, offload_array
        +3*offload_data->n_phi*offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_phi, offload_data->n_z,
	offload_data->r_min, offload_data->r_max, offload_data->r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    #else
    err += interp3Dcomp_init(&Bdata->B_r, offload_array,
			     offload_data->n_r, offload_data->n_phi, offload_data->n_z,
			     offload_data->r_min, offload_data->r_max, offload_data->r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    
    err += interp3Dcomp_init(&Bdata->B_phi, offload_array
			     +offload_data->n_phi*offload_data->n_z*offload_data->n_r,
			     offload_data->n_r, offload_data->n_phi, offload_data->n_z,
			     offload_data->r_min, offload_data->r_max, offload_data->r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    
    err += interp3Dcomp_init(&Bdata->B_z, offload_array
			     +2*offload_data->n_phi*offload_data->n_z*offload_data->n_r,
			     offload_data->n_r, offload_data->n_phi, offload_data->n_z,
			     offload_data->r_min, offload_data->r_max, offload_data->r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp3Dcomp_init(&Bdata->s, offload_array
			     +3*offload_data->n_phi*offload_data->n_z*offload_data->n_r,
			     offload_data->n_r, offload_data->n_phi, offload_data->n_z,
			     offload_data->r_min, offload_data->r_max, offload_data->r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->z_min, offload_data->z_max, offload_data->z_grid);
    #endif
    return err;
    
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * This function evaluates the poloidal flux psi at the given coordinates using
 * tricubic interpolation on the stellarator 3D s data. This is a SIMD
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param psi psi value will be stored in psi[0][i]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Change to a scalar elemental function and compare performance
 */
a5err B_STS_eval_psi(real psi[], real r, real phi, real z,
                   B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    #if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_B(&psi[0], &Bdata->s, r, phi, z);
    #else
    interperr += interp3Dcomp_eval_B(&psi[0], &Bdata->s, r, phi, z);
    #endif

    /* Test for psi interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

a5err B_STS_eval_psi_SIMD(int i, real psi[NSIMD], real r, real phi, real z,
                   B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    interperr += interp3Dcomp_eval_B_SIMD(i, &psi[0], &Bdata->s, r, phi, z);

    /* Test for psi interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and derivatives
 *
 * This function evaluates the poloidal flux psi and it's derivatives at the given 
 * coordinates using tricubic interpolation on the stellarator magnetic field data.
 * This is a SIMD function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param psi psi values (psi   -> psi_dpsi[0][i]    dpsi/dr -> psi_dpsi[1][i]
 *        dpsi/dphi -> psi_dpsi[2][i]    dpsi/dz -> psi_dpsi[3][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 */
a5err B_STS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    interperr += interp3Dcomp_eval_dB(psi_dpsi, &Bdata->s, r, phi, z);

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate radial coordinate rho
 *
 * This function evaluates the radial coordinate rho at the given psi value.
 * This is a SIMD function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param rho rho value will be stored in rho[0][i]
 * @param psi poloidal flux value 
 * @param Bdata pointer to magnetic field data struct
 *
 */
/* TODO error checking */
a5err B_STS_eval_rho(real rho[], real psi, B_STS_data* Bdata) {
    a5err err = 0;

    rho[0] = sqrt(psi);

    return err;
}

/* TODO error checking */
a5err B_STS_eval_rho_SIMD(int i, real rho[NSIMD], real psi, B_STS_data* Bdata) {
    a5err err = 0;

    rho[i] = sqrt(psi);

    return err;
}

/**
 * @brief Evaluate radial coordinate rho and its derivatives
 *
 * This function evaluates the radial coordinate rho and its derivatives
 * at the given coordinates using tricubic interpolation on the 
 * stellarator magnetic field data. This is a SIMD 
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param rho rho values (rho   -> rho_drho[0][i]    drho/dr -> rho_drho[1][i]
 *        drho/dphi -> rho_drho[2][i]    drho/dz -> rho_drho[3][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 */
/* TODO error checking */
a5err B_STS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    real rho;
    B_STS_eval_psi_dpsi(rho_drho, r, phi, z, Bdata);
    /* Convert: rho = sqrt(psi), drho = dpsi/(2 * sqrt(psi)) */
    rho = sqrt(rho_drho[0]);
    rho_drho[0] = rho;
    rho_drho[1] = rho_drho[1] / (2*rho);
    rho_drho[2] = rho_drho[2] / (2*rho);
    rho_drho[3] = rho_drho[3] / (2*rho);

    return err;
}

/**
 * @brief Evaluate magnetic field
 *
 * This function evaluates the magnetic field at the given coordinates using
 * tricubic interpolation on the 3D magnetic field data. This is a SIMD
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param B array where magnetic field values will be stored (Br -> B[0][i],
 *          Bphi -> B[1][i], Bz -> B[2][i]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 */
/* TODO check the psi section */
a5err B_STS_eval_B(real B[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    interperr += interp3Dcomp_eval_B(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[2], &Bdata->B_z, r, phi, z);

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;

}

a5err B_STS_eval_B_SIMD(int i, real B[3][NSIMD], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    interperr += interp3Dcomp_eval_B_SIMD(i, B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_B_SIMD(i, B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_B_SIMD(i, B[2], &Bdata->B_z, r, phi, z);

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0][i]*B[0][i] + B[1][i]*B[1][i] + B[2][i]*B[2][i]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate magnetic field and derivatives
 *
 * This function evaluates the magnetic field and it's derivatives at the given 
 * coordinates using bicubic interpolation on the stellarator magnetic field data.
 * This is a SIMD function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param B_dB array where magnetic field values will be stored (Br -> B[0][i],
 *          dBr/dr -> B[1][i], dBr/dphi -> B[2][i], dBr/dz -> B[3][i],
 *          Bphi -> B[4][i], dBphi/dr -> B[5][i], dBphi/dphi -> B[6][i],
 *          dBphi/dz -> B[7][i], Bz -> B[8][i], dBz/dr -> B[9][i],
 *          dBz/dphi -> B[10][i], dBz/dz -> B[11][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 */
a5err B_STS_eval_B_dB(real B_dB[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];
   
    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }

    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_r, r, phi, z);

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];
    
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_phi, r, phi, z);

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];
    
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_z, r, phi, z);

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;
}

a5err B_STS_eval_B_dB_SIMD(int i, real B_dB[12][NSIMD], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10][NSIMD];

    phi = fmod(phi, 2*math_pi/Bdata->periods);
    if(phi < 0) {
        phi += 2*math_pi/Bdata->periods;
    }
    
    interperr += interp3Dcomp_eval_dB_SIMD(i, B_dB_temp, &Bdata->B_r, r, phi, z);

    B_dB[0][i] = B_dB_temp[0][i];
    B_dB[1][i] = B_dB_temp[1][i];
    B_dB[2][i] = B_dB_temp[2][i];
    B_dB[3][i] = B_dB_temp[3][i];

    interperr += interp3Dcomp_eval_dB_SIMD(i, B_dB_temp, &Bdata->B_phi, r, phi, z);

    B_dB[4][i] = B_dB_temp[0][i];
    B_dB[5][i] = B_dB_temp[1][i];
    B_dB[6][i] = B_dB_temp[2][i];
    B_dB[7][i] = B_dB_temp[3][i];

    interperr += interp3Dcomp_eval_dB_SIMD(i, B_dB_temp, &Bdata->B_z, r, phi, z);

    B_dB[8][i] = B_dB_temp[0][i];
    B_dB[9][i] = B_dB_temp[1][i];
    B_dB[10][i] = B_dB_temp[2][i];
    B_dB[11][i] = B_dB_temp[3][i];

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0][i]*B_dB[0][i] + B_dB[4][i]*B_dB[4][i] + B_dB[8][i]*B_dB[8][i]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;
}

real B_STS_get_axis_r(B_STS_data* Bdata) {
    // 3D magnetic axis not implemented yet
    return 0;
}

real B_STS_get_axis_z(B_STS_data* Bdata) {
    // 3D magnetic axis not implemented yet
    return 0;
}
