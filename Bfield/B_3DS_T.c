/**
 * @file B_3DS_T.c
 * @brief 3D magnetic field with cubic spline interpolation
 * time linear interpolation
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
#include "../linint/linint1D.h"

/**
 * @brief Load magnetic field data and prepare parameters
 *
 * This function reads the magnetic field data from input.magn_bkg and 
 * input.magn_header files, fills the offload struct with parameters and 
 * allocates and fills the offload array.
 *
 * @todo Error checking
 * @todo Move reading the file to ascot4_interface
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3DS_T_init_offload(B_3DS_T_offload_data* offload_data, real** offload_array) {
    
    // The contents of this function are located in hdf5_bfield.c
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3DS_T_free_offload(B_3DS_T_offload_data* offload_data, real** offload_array) {
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
int B_3DS_T_init(B_3DS_T_data* Bdata, B_3DS_T_offload_data* offload_data,
               real* offload_array) {
    int err = 0;
    int i;

    int B_size = offload_data->n_phi*offload_data->Bgrid_n_z*offload_data->Bgrid_n_r;
    int psi_size = offload_data->psigrid_n_z*offload_data->psigrid_n_r;
    
    Bdata->n_time = offload_data->n_time;

    for(i=0;i<offload_data->n_time;i++){

      Bdata->time[i] = offload_data->time[i];
      
      Bdata->Bslice[i].psi0 = offload_data->psi0[i];
      Bdata->Bslice[i].psi1 = offload_data->psi1[i];
      Bdata->Bslice[i].axis_r = offload_data->axis_r[i];
      Bdata->Bslice[i].axis_z = offload_data->axis_z[i];
      /* Spline initialization and storage. */

      #if INTERP_SPL_EXPL
      err += interp2Dexpl_init(&Bdata->Bslice[i].psi, offload_array + 3*B_size*offload_data->n_time + i*psi_size,
	offload_data->psigrid_n_r, offload_data->psigrid_n_z,
	offload_data->psigrid_r_min, offload_data->psigrid_r_max, offload_data->psigrid_r_grid,
	offload_data->psigrid_z_min, offload_data->psigrid_z_max, offload_data->psigrid_z_grid);
    
      err += interp3Dexpl_init(&Bdata->Bslice[i].B_r, offload_array + 0*B_size*offload_data->n_time + i*B_size,
	offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
	offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
      err += interp3Dexpl_init(&Bdata->Bslice[i].B_phi, offload_array + 1*B_size*offload_data->n_time + i*B_size,
	offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
	offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    	
      err += interp3Dexpl_init(&Bdata->Bslice[i].B_z, offload_array + 2*B_size*offload_data->n_time + i*B_size,
	offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
	offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
	offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
	offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
      #else
      err += interp2Dcomp_init(&Bdata->Bslice[i].psi, offload_array + 3*B_size*offload_data->n_time + i*psi_size,
			     offload_data->psigrid_n_r, offload_data->psigrid_n_z,
			     offload_data->psigrid_r_min, offload_data->psigrid_r_max, offload_data->psigrid_r_grid,
			     offload_data->psigrid_z_min, offload_data->psigrid_z_max, offload_data->psigrid_z_grid);
    
      err += interp3Dcomp_init(&Bdata->Bslice[i].B_r, offload_array + 0*B_size*offload_data->n_time + i*B_size,
			     offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
			     offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
      err += interp3Dcomp_init(&Bdata->Bslice[i].B_phi, offload_array + 1*B_size*offload_data->n_time + i*B_size,
			     offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
			     offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
    
      err += interp3Dcomp_init(&Bdata->Bslice[i].B_z, offload_array + 2*B_size*offload_data->n_time + i*B_size,
			     offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
			     offload_data->Bgrid_r_min, offload_data->Bgrid_r_max, offload_data->Bgrid_r_grid,
			     offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
			     offload_data->Bgrid_z_min, offload_data->Bgrid_z_max, offload_data->Bgrid_z_grid);
      #endif
    }

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
 * @param time time coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Change to a scalar elemental function and compare performance
 */
a5err B_3DS_T_eval_psi(real psi[], real r, real phi, real z, real time,
                   B_3DS_T_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    real psi_temp_i0;
    real psi_temp_i1;

    int i = 0;

    if(time <= Bdata->time[0]){      
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_B(&psi[0], &Bdata->Bslice[0].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_B(&psi[0], &Bdata->Bslice[0].psi, r, z);
    #endif
    }
    else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_B(&psi[0], &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_B(&psi[0], &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #endif
    }
    else{
      while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_B(&psi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
      interperr += interp2Dexpl_eval_B(&psi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_B(&psi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
      interperr += interp2Dcomp_eval_B(&psi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #endif

      interperr += linint1D_eval(&psi[0], Bdata->time[i-1], psi_temp_i0, Bdata->time[i], psi_temp_i1, time);
    }



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
a5err B_3DS_T_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z, real time,
                   B_3DS_T_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp_i0[6];
    real psi_dpsi_temp_i1[6];

    int i = 0;
 
    psi_dpsi[2] = 0;

    if(time <= Bdata->time[0]){
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #endif
      psi_dpsi[0] = psi_dpsi_temp_i0[0];
      psi_dpsi[1] = psi_dpsi_temp_i0[1];
      psi_dpsi[3] = psi_dpsi_temp_i0[2];
    }
    else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #endif
      psi_dpsi[0] = psi_dpsi_temp_i0[0];
      psi_dpsi[1] = psi_dpsi_temp_i0[1];
      psi_dpsi[3] = psi_dpsi_temp_i0[2];
    }
    else{
      while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #endif
      interperr += linint1D_eval(&psi_dpsi[0], Bdata->time[i-1], psi_dpsi_temp_i0[0], Bdata->time[i], psi_dpsi_temp_i1[0], time);
      interperr += linint1D_eval(&psi_dpsi[1], Bdata->time[i-1], psi_dpsi_temp_i0[1], Bdata->time[i], psi_dpsi_temp_i1[1], time);
      interperr += linint1D_eval(&psi_dpsi[3], Bdata->time[i-1], psi_dpsi_temp_i0[2], Bdata->time[i], psi_dpsi_temp_i1[2], time);
    }

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


a5err B_3DS_T_eval_rho(real rho[], real psi, real time, B_3DS_T_data* Bdata) {
    a5err err = 0;
    int interperr = 0;
    
    real psi0_temp;
    real psi1_temp;


    int i = 0;


    if(time <= Bdata->time[0]){
      psi0_temp = B_data->Bslice[0].psi0;
      psi1_temp = B_data->Bslice[0].psi1;
    }
    else if (time >= Bdata->time[Bdata->n_time-1]){
      psi0_temp = B_data->Bslice[Bdata->n_time-1].psi0;
      psi1_temp = B_data->Bslice[Bdata->n_time-1].psi1;
    }
    else{
      while(time > Bdata->time[i]) i++;
      interperr += linint1D_eval(&psi0_temp, Bdata->time[i-1], Bdata->Bslice[i-1].psi0, Bdata->time[i], Bdata->Bslice[i].psi0, time);
      interperr += linint1D_eval(&psi1_temp, Bdata->time[i-1], Bdata->Bslice[i-1].psi1, Bdata->time[i], Bdata->Bslice[i].psi1, time);
    }

    /* Check that the values seem valid */
    //if( (psi - Bdata->psi0) < 0 ) {err = error_raise( ERR_UNPHYSICAL_PSI, __LINE__ );}
    //else {
    /* Normalize psi to get rho */
    rho[0] = sqrt(fabs( (psi - psi0_temp) / (psi1_temp - psi0_temp) ));
    //}


    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

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
a5err B_3DS_T_eval_rho_drho(real rho_drho[], real r, real phi, real z, real time, B_3DS_T_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp_i0[6];
    real psi_dpsi_temp_i1[6];

    int i = 0;

    if(time <= Bdata->time[0]){
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #endif
    }
    else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #endif
    }
    else{
      while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
      interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #else
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
      interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #endif
      interperr += linint1D_eval(&psi_dpsi_temp_i0[0], Bdata->time[i-1], psi_dpsi_temp_i0[0], Bdata->time[i], psi_dpsi_temp_i1[0], time);
      interperr += linint1D_eval(&psi_dpsi_temp_i0[1], Bdata->time[i-1], psi_dpsi_temp_i0[1], Bdata->time[i], psi_dpsi_temp_i1[1], time);
      interperr += linint1D_eval(&psi_dpsi_temp_i0[2], Bdata->time[i-1], psi_dpsi_temp_i0[2], Bdata->time[i], psi_dpsi_temp_i1[2], time);
      interperr += linint1D_eval(&psi_dpsi_temp_i0[3], Bdata->time[i-1], psi_dpsi_temp_i0[3], Bdata->time[i], psi_dpsi_temp_i1[3], time);
      interperr += linint1D_eval(&psi_dpsi_temp_i0[4], Bdata->time[i-1], psi_dpsi_temp_i0[4], Bdata->time[i], psi_dpsi_temp_i1[4], time);
      interperr += linint1D_eval(&psi_dpsi_temp_i0[5], Bdata->time[i-1], psi_dpsi_temp_i0[5], Bdata->time[i], psi_dpsi_temp_i1[5], time);
    }

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}
    

    //if(!err && (psi_dpsi[0] - Bdata->psi0) < 0) {err = error_raise( ERR_UNPHYSICAL_PSI, __LINE__ );}
    if(!err) {
      
      real psi0_temp;
      real psi1_temp;

      i = 0;

      if(time <= Bdata->time[0]){
	psi0_temp = B_data->Bslice[0].psi0;
	psi1_temp = B_data->Bslice[0].psi1;
      }
      else if (time >= Bdata->time[Bdata->n_time-1]){
	psi0_temp = B_data->Bslice[Bdata->n_time-1].psi0;
	psi1_temp = B_data->Bslice[Bdata->n_time-1].psi1;
      }
      else{
	while(time > Bdata->time[i]) i++;
	interperr += linint1D_eval(&psi0_temp, Bdata->time[i-1], Bdata->Bslice[i-1].psi0, Bdata->time[i], Bdata->Bslice[i].psi0, time);
	interperr += linint1D_eval(&psi1_temp, Bdata->time[i-1], Bdata->Bslice[i-1].psi1, Bdata->time[i], Bdata->Bslice[i].psi1, time);
      }


        /* Normalize psi to get rho */
      real delta = psi1_temp - psi0_temp;
      rho_drho[0] = sqrt(fabs( (psi_dpsi_temp_i0[0] - psi0_temp) / delta ));
	
      if(rho_drho[0] == 0) {
	rho_drho[1] = psi_dpsi_temp_i0[1] / (2*delta*rho_drho[0]);
	rho_drho[2] = 0;
	rho_drho[3] = psi_dpsi_temp_i0[2] / (2*delta*rho_drho[0]);
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
a5err B_3DS_T_eval_B(real B[], real r, real phi, real z, real time, B_3DS_T_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_temp_i0[3];
    real B_temp_i1[3];


    int i = 0;

    if(time <= Bdata->time[0]){
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_B(&B[0], &Bdata->Bslice[0].B_r, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B[1], &Bdata->Bslice[0].B_phi, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B[2], &Bdata->Bslice[0].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_B(&B[0], &Bdata->Bslice[0].B_r, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B[1], &Bdata->Bslice[0].B_phi, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B[2], &Bdata->Bslice[0].B_z, r, phi, z);
    #endif
    }
    else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_B(&B[0], &Bdata->Bslice[Bdata->n_time-1].B_r, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B[1], &Bdata->Bslice[Bdata->n_time-1].B_phi, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B[2], &Bdata->Bslice[Bdata->n_time-1].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_B(&B[0], &Bdata->Bslice[Bdata->n_time-1].B_r, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B[1], &Bdata->Bslice[Bdata->n_time-1].B_phi, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B[2], &Bdata->Bslice[Bdata->n_time-1].B_z, r, phi, z);
    #endif
    }
    else{
      while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_B(&B_temp_i0[0], &Bdata->Bslice[i-1].B_r, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B_temp_i0[1], &Bdata->Bslice[i-1].B_phi, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B_temp_i0[2], &Bdata->Bslice[i-1].B_z, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B_temp_i1[0], &Bdata->Bslice[i].B_r, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B_temp_i1[1], &Bdata->Bslice[i].B_phi, r, phi, z);
      interperr += interp3Dexpl_eval_B(&B_temp_i1[2], &Bdata->Bslice[i].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_B(&B_temp_i0[0], &Bdata->Bslice[i-1].B_r, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B_temp_i0[1], &Bdata->Bslice[i-1].B_phi, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B_temp_i0[2], &Bdata->Bslice[i-1].B_z, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B_temp_i1[0], &Bdata->Bslice[i].B_r, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B_temp_i1[1], &Bdata->Bslice[i].B_phi, r, phi, z);
      interperr += interp3Dcomp_eval_B(&B_temp_i1[2], &Bdata->Bslice[i].B_z, r, phi, z);
    #endif

      interperr += linint1D_eval(&B[0], Bdata->time[i-1], B_temp_i0[0], Bdata->time[i], B_temp_i1[0], time);
      interperr += linint1D_eval(&B[1], Bdata->time[i-1], B_temp_i0[1], Bdata->time[i], B_temp_i1[1], time);
      interperr += linint1D_eval(&B[2], Bdata->time[i-1], B_temp_i0[2], Bdata->time[i], B_temp_i1[2], time);

    }
    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}


    #ifndef NOPSI
    if(!err) {
	real psi_dpsi_temp_i0[6];
	real psi_dpsi_temp_i1[6];

	i = 0;

	if(time <= Bdata->time[0]){
    #if INTERP_SPL_EXPL
	  interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #else
	  interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #endif
	}
	else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
	  interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #else
	  interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #endif
	}
	else{
	  while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
	  interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
	  interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #else
	  interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
	  interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #endif
	  interperr += linint1D_eval(&psi_dpsi_temp_i0[0], Bdata->time[i-1], psi_dpsi_temp_i0[0], Bdata->time[i], psi_dpsi_temp_i1[0], time);
	  interperr += linint1D_eval(&psi_dpsi_temp_i0[1], Bdata->time[i-1], psi_dpsi_temp_i0[1], Bdata->time[i], psi_dpsi_temp_i1[1], time);
	  interperr += linint1D_eval(&psi_dpsi_temp_i0[2], Bdata->time[i-1], psi_dpsi_temp_i0[2], Bdata->time[i], psi_dpsi_temp_i1[2], time);
	  interperr += linint1D_eval(&psi_dpsi_temp_i0[3], Bdata->time[i-1], psi_dpsi_temp_i0[3], Bdata->time[i], psi_dpsi_temp_i1[3], time);
	  interperr += linint1D_eval(&psi_dpsi_temp_i0[4], Bdata->time[i-1], psi_dpsi_temp_i0[4], Bdata->time[i], psi_dpsi_temp_i1[4], time);
	  interperr += linint1D_eval(&psi_dpsi_temp_i0[5], Bdata->time[i-1], psi_dpsi_temp_i0[5], Bdata->time[i], psi_dpsi_temp_i1[5], time);
	}

	if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

	B[0] = B[0] - psi_dpsi_temp_i0[2]/r;
	B[2] = B[2] + psi_dpsi_temp_i0[1]/r;
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
 * coordinates using bicubic interpolation on the 3D magnetic field data. Values
 * are linear interpolated between time slices. This is a SIMD function, so the 
 * values are placed in an NSIMD length struct.
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
 * @param time time coordinate
 * @param Bdata pointer to magnetic field data struct
 */
a5err B_3DS_T_eval_B_dB(real B_dB[], real r, real phi, real z, real time, B_3DS_T_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];
    real B_dB_temp_i0[12];
    real B_dB_temp_i1[12];

    int i = 0;
    int j = 0;

    if(time <= Bdata->time[0]){
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[0].B_r, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[0].B_r, r, phi, z);
    #endif
      B_dB[0] = B_dB_temp[0];
      B_dB[1] = B_dB_temp[1];
      B_dB[2] = B_dB_temp[2];
      B_dB[3] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[0].B_phi, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[0].B_phi, r, phi, z);
    #endif
      B_dB[4] = B_dB_temp[0];
      B_dB[5] = B_dB_temp[1];
      B_dB[6] = B_dB_temp[2];
      B_dB[7] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[0].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[0].B_z, r, phi, z);
    #endif
      B_dB[8] = B_dB_temp[0];
      B_dB[9] = B_dB_temp[1];
      B_dB[10] = B_dB_temp[2];
      B_dB[11] = B_dB_temp[3];
    }
    else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[Bdata->n_time-1].B_r, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[Bdata->n_time-1].B_r, r, phi, z);
    #endif
      B_dB[0] = B_dB_temp[0];
      B_dB[1] = B_dB_temp[1];
      B_dB[2] = B_dB_temp[2];
      B_dB[3] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[Bdata->n_time-1].B_phi, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[Bdata->n_time-1].B_phi, r, phi, z);
    #endif
      B_dB[4] = B_dB_temp[0];
      B_dB[5] = B_dB_temp[1];
      B_dB[6] = B_dB_temp[2];
      B_dB[7] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[Bdata->n_time-1].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[Bdata->n_time-1].B_z, r, phi, z);
    #endif
      B_dB[8] = B_dB_temp[0];
      B_dB[9] = B_dB_temp[1];
      B_dB[10] = B_dB_temp[2];
      B_dB[11] = B_dB_temp[3];
    }
    else{
      while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[i-1].B_r, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[i-1].B_r, r, phi, z);
    #endif
      B_dB_temp_i0[0] = B_dB_temp[0];
      B_dB_temp_i0[1] = B_dB_temp[1];
      B_dB_temp_i0[2] = B_dB_temp[2];
      B_dB_temp_i0[3] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[i-1].B_phi, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[i-1].B_phi, r, phi, z);
    #endif
      B_dB_temp_i0[4] = B_dB_temp[0];
      B_dB_temp_i0[5] = B_dB_temp[1];
      B_dB_temp_i0[6] = B_dB_temp[2];
      B_dB_temp_i0[7] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[i-1].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[i-1].B_z, r, phi, z);
    #endif
      B_dB_temp_i0[8] = B_dB_temp[0];
      B_dB_temp_i0[9] = B_dB_temp[1];
      B_dB_temp_i0[10] = B_dB_temp[2];
      B_dB_temp_i0[11] = B_dB_temp[3];
    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[i].B_r, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[i].B_r, r, phi, z);
    #endif
      B_dB_temp_i1[0] = B_dB_temp[0];
      B_dB_temp_i1[[1] = B_dB_temp[1];
      B_dB_temp_i1[[2] = B_dB_temp[2];
      B_dB_temp_i1[[3] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[i].B_phi, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[i].B_phi, r, phi, z);
    #endif
      B_dB_temp_i1[[4] = B_dB_temp[0];
      B_dB_temp_i1[[5] = B_dB_temp[1];
      B_dB_temp_i1[[6] = B_dB_temp[2];
      B_dB_temp_i1[[7] = B_dB_temp[3];

    #if INTERP_SPL_EXPL
      interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->Bslice[i].B_z, r, phi, z);
    #else
      interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->Bslice[i].B_z, r, phi, z);
    #endif
      B_dB_temp_i1[[8] = B_dB_temp[0];
      B_dB_temp_i1[[9] = B_dB_temp[1];
      B_dB_temp_i1[[10] = B_dB_temp[2];
      B_dB_temp_i1[[11] = B_dB_temp[3];
      
      for(j=0;j<=11;j++){
	interperr += linint1D_eval(&B_dB[j], Bdata->time[i-1], B_dB_temp_i0[j], Bdata->time[i], B_dB_temp_i1[j], time);
      }

    }
    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}


    #ifndef NOPSI

    if(!err) {

      real psi_dpsi_temp_i0[6];
      real psi_dpsi_temp_i1[6];

      i = 0;

      if(time <= Bdata->time[0]){
    #if INTERP_SPL_EXPL
	interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #else
	interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[0].psi, r, z);
    #endif
      }
      else if (time >= Bdata->time[Bdata->n_time-1]){
    #if INTERP_SPL_EXPL
	interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #else
	interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[Bdata->n_time-1].psi, r, z);
    #endif
      }
      else{
	while(time > Bdata->time[i]) i++;
    #if INTERP_SPL_EXPL
	interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
	interperr += interp2Dexpl_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #else
	interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i0, &Bdata->Bslice[i-1].psi, r, z);
	interperr += interp2Dcomp_eval_dB(psi_dpsi_temp_i1, &Bdata->Bslice[i].psi, r, z);
    #endif
	interperr += linint1D_eval(&psi_dpsi_temp_i0[0], Bdata->time[i-1], psi_dpsi_temp_i0[0], Bdata->time[i], psi_dpsi_temp_i1[0], time);
	interperr += linint1D_eval(&psi_dpsi_temp_i0[1], Bdata->time[i-1], psi_dpsi_temp_i0[1], Bdata->time[i], psi_dpsi_temp_i1[1], time);
	interperr += linint1D_eval(&psi_dpsi_temp_i0[2], Bdata->time[i-1], psi_dpsi_temp_i0[2], Bdata->time[i], psi_dpsi_temp_i1[2], time);
	interperr += linint1D_eval(&psi_dpsi_temp_i0[3], Bdata->time[i-1], psi_dpsi_temp_i0[3], Bdata->time[i], psi_dpsi_temp_i1[3], time);
	interperr += linint1D_eval(&psi_dpsi_temp_i0[4], Bdata->time[i-1], psi_dpsi_temp_i0[4], Bdata->time[i], psi_dpsi_temp_i1[4], time);
	interperr += linint1D_eval(&psi_dpsi_temp_i0[5], Bdata->time[i-1], psi_dpsi_temp_i0[5], Bdata->time[i], psi_dpsi_temp_i1[5], time);
      }

      B_dB[0] = B_dB[0] - psi_dpsi_temp_i0[2]/r;
      B_dB[1] = B_dB[1] + psi_dpsi_temp_i0[2]/(r*r)-psi_dpsi_temp_i0[5]/r;
      B_dB[3] = B_dB[3] - psi_dpsi_temp_i0[4]/r;
      B_dB[8] = B_dB[8] + psi_dpsi_temp_i0[1]/r;
      B_dB[9] = B_dB[9] - psi_dpsi_temp_i0[1]/(r*r) + psi_dpsi_temp_i0[3]/r;
      B_dB[11] = B_dB[11] + psi_dpsi_temp_i0[5]/r;

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


real B_3DS_get_axis_r(real time, B_3DS_data* Bdata) {
	a5err err;
	int interperr;
	int i = 0;
	real out;
	if(time <= Bdata->time[0]){ 
	  out = Bdata->Bslide[0].axis_r;
	}
	else if (time >= Bdata->time[Bdata->n_time-1]){
	  out = Bdata->Bslide[Bdata->n_time-1].axis_r;
	}
	else{
	  while(time > Bdata->time[i]) i++;
	  interperr += linint1D_eval(&out, Bdata->time[i-1], Bdata->Bslide[i-1].axis_r, Bdata->time[i], Bdata->Bslide[i].axis_r, time);
	}

	if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return out;
}

real B_3DS_get_axis_z(real time, B_3DS_data* Bdata) {
  a5err err;
  int interperr;
  int i = 0;
  real out;
  if(time <= Bdata->time[0]){
    out = Bdata->Bslide[0].axis_z;
  }
  else if (time >= Bdata->time[Bdata->n_time-1]){
    out = Bdata->Bslide[Bdata->n_time-1].axis_z;
  }
  else{
    while(time > Bdata->time[i]) i++;
    interperr += linint1D_eval(&out, Bdata->time[i-1], Bdata->Bslide[i-1].axis_z, Bdata->time[i], Bdata->Bslide[i].axis_z, time);
  }
  
  if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}
  
  return out;
}
