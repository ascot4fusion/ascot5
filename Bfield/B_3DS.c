/**
 * @file B_3DS.c
 * @brief 3D magnetic field with cubic spline interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "B_3DS.h"
#include "../spline/interp2D.h" /* for 2D interpolation routines */
#include "../spline/interp3D.h" /* for 3D interpolation routines */
#include "../spline/interp2Dcomp.h" 
#include "../spline/interp3Dcomp.h" 
#include "../spline/interp2Dexpl.h" 
#include "../spline/interp3Dexpl.h"

/**
 * @brief Load magnetic field data and prepare parameters
 *
 * This function reads the magnetic field data from input.magn_bkg and 
 * input.magn_header files, fills the offload struct with parameters and 
 * allocates and fills the offload array.
 *
 * The reading of the ASCOT4 3D magnetic field has been moved into
 * ascot4_interface, so this function is now a dummy.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3DS_init_offload(B_3DS_offload_data* offload_data, real** offload_array) {
    // Dummy function
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3DS_free_offload(B_3DS_offload_data* offload_data, real** offload_array) {
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
int B_3DS_init(B_3DS_data* Bdata, B_3DS_offload_data* offload_data,
               real* offload_array) {
    int err = 0;
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;
    /* Spline initialization and storage. */

    int B_size = offload_data->n_phi*offload_data->Bgrid_n_z*offload_data->Bgrid_n_r;
    
    #if INTERP_SPL_EXPL
    err += interp2Dexpl_init(&Bdata->psi, offload_array + 3*B_size,
	offload_data->psigrid_n_r, offload_data->psigrid_n_z,
	offload_data->psigrid_r_min, offload_data->psigrid_r_max, offload_data->psigrid_r_grid,
	offload_data->psigrid_z_min, offload_data->psigrid_z_max, offload_data->psigrid_z_grid);
    
    err += interp3Dexpl_init(&Bdata->B_r, offload_array + 0*B_size,
	offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
	offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
    err += interp3Dexpl_init(&Bdata->B_phi, offload_array + 1*B_size,
	offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
	offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    	
    err += interp3Dexpl_init(&Bdata->B_z, offload_array + 2*B_size,
	offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
	offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
    #else
    err += interp2Dcomp_init(&Bdata->psi, offload_array + 3*B_size,
			     offload_data->psigrid_n_r, offload_data->psigrid_n_z,
			     offload_data->psigrid_r_min, offload_data->psigrid_r_max, offload_data->psigrid_r_grid,
			     offload_data->psigrid_z_min, offload_data->psigrid_z_max, offload_data->psigrid_z_grid);
    
    err += interp3Dcomp_init(&Bdata->B_r, offload_array + 0*B_size,
			     offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
			     offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
    err += interp3Dcomp_init(&Bdata->B_phi, offload_array + 1*B_size,
			     offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
			     offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
    err += interp3Dcomp_init(&Bdata->B_z, offload_array + 2*B_size,
			     offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
			     offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    #endif
    return err;
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * This function evaluates the poloidal flux psi at the given coordinates using
 * bicubic interpolation on the 3D magnetic field data. This is a SIMD
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
a5err B_3DS_eval_psi(real psi[], real r, real phi, real z,
                   B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    #if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_B(&psi[0], &Bdata->psi, r, z);
    #else
    interperr += interp2Dcomp_eval_B(&psi[0], &Bdata->psi, r, z);
    #endif

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * This function evaluates the poloidal flux psi at the given coordinates using
 * tricubic interpolation on the 3D magnetic field data. This is a SIMD
 * function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param psi psi values (psi   -> psi_dpsi[0][i]    dpsi/dr -> psi_dpsi[1][i]
 *        dpsi/dphi -> psi_dpsi[2][i]    dpsi/dz -> psi_dpsi[3][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Change to a scalar elemental function and compare performance
 */
a5err B_3DS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z,
                   B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp[6];
    #if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(psi_dpsi_temp, &Bdata->psi, r, z);
    #else
    interperr += interp2Dcomp_eval_dB(psi_dpsi_temp, &Bdata->psi, r, z);
    #endif
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate radial coordinate rho
 *
 * This function evaluates the radial coordinate rho at the given psi value
 * using linear interpolation. This is a SIMD function, so the values are 
 * placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param rho rho value will be stored in rho[0][i]
 * @param psi poloidal flux value 
 * @param Bdata pointer to magnetic field data struct
 */
a5err B_3DS_eval_rho(real rho[], real psi, B_3DS_data* Bdata) {
    a5err err = 0;

    /* Check that the values seem valid */
    //if( (psi - Bdata->psi0) < 0 ) {err = error_raise( ERR_UNPHYSICAL_PSI, __LINE__ );}
    //else {
    /* Normalize psi to get rho */
    rho[0] = sqrt(fabs( (psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0) ));
    //}

    return err;
}

/**
 * @brief Evaluate radial coordinate rho and its derivatives
 *
 * This function evaluates the radial coordinate rho and its derivatives
 * at the given coordinates using tricubic interpolation on the 
 * 3D magnetic field data. This is a SIMD function, so the values are
 * placed in an NSIMD length struct.
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
a5err B_3DS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi[6];
    #if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
    #else
    interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
    #endif

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    //if(!err && (psi_dpsi[0] - Bdata->psi0) < 0) {err = error_raise( ERR_UNPHYSICAL_PSI, __LINE__ );}
    if(!err) {
        /* Normalize psi to get rho */
        real delta = Bdata->psi1 - Bdata->psi0;
        rho_drho[0] = sqrt(fabs( (psi_dpsi[0] - Bdata->psi0) / delta ));
	
	if(rho_drho[0] == 0) {
            rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
	    rho_drho[2] = 0;
	    rho_drho[3] = psi_dpsi[2] / (2*delta*rho_drho[0]);
        }
	else {
            rho_drho[1] = 0;
            rho_drho[2] = 0;
	    rho_drho[3] = 0;
        }

    }

    return err;
}

/**
 * @brief Evaluate magnetic field
 *
 * This function evaluates the magnetic field at the given coordinates using
 * bicubic interpolation on the 3D magnetic field data. This is a SIMD
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
a5err B_3DS_eval_B(real B[], real r, real phi, real z, B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    #if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_B(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dexpl_eval_B(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dexpl_eval_B(&B[2], &Bdata->B_z, r, phi, z);
    #else
    interperr += interp3Dcomp_eval_B(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[2], &Bdata->B_z, r, phi, z);
    #endif

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

    #ifndef NOPSI
    if(!err) {
	real psi_dpsi[6];
	#if INTERP_SPL_EXPL
	interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
	#else
	interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
	#endif
	B[0] = B[0] - psi_dpsi[2]/r;
	B[2] = B[2] + psi_dpsi[1]/r;

	/* Test for psi interpolation error */
	if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}
    }
    #endif

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate magnetic field and derivatives
 *
 * This function evaluates the magnetic field and it's derivatives at the given 
 * coordinates using bicubic interpolation on the 3D magnetic field data. This 
 * is a SIMD function, so the values are placed in an NSIMD length struct.
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
a5err B_3DS_eval_B_dB(real B_dB[], real r, real phi, real z, B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];
    #if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->B_r, r, phi, z);
    #else
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_r, r, phi, z);
    #endif

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];


    #if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->B_phi, r, phi, z);
    #else
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_phi, r, phi, z);
    #endif

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];


    #if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->B_z, r, phi, z);
    #else
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_z, r, phi, z);
    #endif

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

    #ifndef NOPSI
    real psi_dpsi[6];
    if(!err) {
	#if INTERP_SPL_EXPL
	interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
	#else
	interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
	#endif

	B_dB[0] = B_dB[0] - psi_dpsi[2]/r;
	B_dB[1] = B_dB[1] + psi_dpsi[2]/(r*r)-psi_dpsi[5]/r;
	B_dB[3] = B_dB[3] - psi_dpsi[4]/r;
	B_dB[8] = B_dB[8] + psi_dpsi[1]/r;
	B_dB[9] = B_dB[9] - psi_dpsi[1]/(r*r) + psi_dpsi[3]/r;
	B_dB[11] = B_dB[11] + psi_dpsi[5]/r;

	/* Test for psi interpolation error */
	if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}
    }
    #endif

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;
}

real B_3DS_get_axis_r(B_3DS_data* Bdata) {
    return Bdata->axis_r;
}

real B_3DS_get_axis_z(B_3DS_data* Bdata) {
    return Bdata->axis_z;
}
