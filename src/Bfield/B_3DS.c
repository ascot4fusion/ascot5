/**
 * @file B_3DS.c
 * @brief 3D magnetic field with tricubic spline interpolation
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "B_3DS.h"
#include "../spline/interp.h"

/**
 * @brief Initialize magnetic field data
 *
 * The offload array must contain the following data:
 * - offload_array[                     j*Bn_r*Bn_phi + z*Bn_r + i]
 *   = B_R(R_i, phi_z, z_j)   [T]
 * - offload_array[  Bn_r*Bn_z*Bn_phi + j*Bn_r*Bn_phi + z*Bn_r + i]
 *   = B_phi(R_i, phi_z, z_j) [T]
 * - offload_array[2*Bn_r*Bn_z*Bn_phi + j*Bn_r*Bn_phi + z*Bn_r + i]
 *   = B_z(R_i, phi_z, z_j)   [T]
 * - offload_array[3*Bn_r*Bn_z*Bn_phi + j*n_r + i]
 *   = psi(R_i, z_j)   [V*s*m^-1]
 *
 * @param data pointer to the data struct
 * @param p_n_r number of r grid points in psi data
 * @param p_r_min minimum R coordinate in psi data grid [m]
 * @param p_r_max maximum R coordinate in psi data grid [m]
 * @param p_n_z number of z grid points in psi data
 * @param p_z_min minimum z coordinate in psi data grid [m]
 * @param p_z_max maximum z coordinate in psi data grid [m]
 * @param p_n_r number of r grid points in B data
 * @param b_r_min minimum R coordinate in B data grid [m]
 * @param b_r_max maximum R coordinate in B data grid [m]
 * @param b_n_phi number of phi grid points in B data
 * @param b_phi_min minimum phi coordinate in B data grid [rad]
 * @param b_phi_max maximum phi coordinate in B data grid [rad]
 * @param b_n_z number of z grid points in B data
 * @param b_z_min minimum z coordinate in B data grid [m]
 * @param b_z_max maximum z coordinate in B data grid [m]
 * @param axis_r R coordinate of magnetic axis [m]
 * @param axis_z z coordinate of magnetic axis [m]
 * @param psi0 poloidal flux at magnetic axis [Vs/m]
 * @param psi1 poloidal flux at separatrix [Vs/m]
 * @param psi poloidal flux psi(R_i,z_j) = arr[j*n_r + i] [Vs/m]
 * @param B_r Magnetic field R component
 *        B_r(R_i,phi_j,z_k) = arr[k*b_n_r*b_n_phi + j*b_n_r + i] [T]
 * @param B_phi Magnetic field phi component
 *        B_phi(R_i,phi_j,z_k) = arr[k*b_n_r*b_n_phi + j*b_n_r + i] [T]
 * @param B_z Magnetic field z component
 *        B_z(R_i,phi_j,z_k) = arr[k*b_n_r*b_n_phi + j*b_n_r + i] [T]
 *
 * @return zero if initialization succeeded
 */
int B_3DS_init(B_3DS_data* data,
               int p_n_r, real p_r_min, real p_r_max,
               int p_n_z, real p_z_min, real p_z_max,
               int b_n_r, real b_r_min, real b_r_max,
               int b_n_phi, real b_phi_min, real b_phi_max,
               int b_n_z, real b_z_min, real b_z_max,
               real axis_r, real axis_z, real psi0, real psi1,
               real* psi, real* B_r, real* B_phi, real* B_z) {

    int err = 0;
    data->psi0   = psi0;
    data->psi1   = psi1;
    data->axis_r = axis_r;
    data->axis_z = axis_z;

    /* Set up the splines */
    err += interp2Dcomp_setup(&data->psi, psi, p_n_r, p_n_z,
                              NATURALBC, NATURALBC,
                              p_r_min, p_r_max, p_z_min, p_z_max);
    err += interp3Dcomp_setup(&data->B_r, B_r, b_n_r, b_n_phi, b_n_z,
                              NATURALBC, PERIODICBC, NATURALBC,
                              b_r_min, b_r_max, b_phi_min, b_phi_max,
                              b_z_min, b_z_max);
    err += interp3Dcomp_setup(&data->B_phi, B_phi, b_n_r, b_n_phi, b_n_z,
                              NATURALBC, PERIODICBC, NATURALBC,
                              b_r_min, b_r_max, b_phi_min, b_phi_max,
                              b_z_min, b_z_max);
    err += interp3Dcomp_setup(&data->B_z, B_z, b_n_r, b_n_phi, b_n_z,
                              NATURALBC, PERIODICBC, NATURALBC,
                              b_r_min, b_r_max, b_phi_min, b_phi_max,
                              b_z_min, b_z_max);

    /* Evaluate psi and magnetic field on axis for checks */
    real psival[1], Bval[3];
    err += B_3DS_eval_psi(psival, data->axis_r, 0, data->axis_z, data);
    err += B_3DS_eval_B(Bval, data->axis_r, 0, data->axis_z, data);
    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void B_3DS_free(B_3DS_data* data) {
    free(data->psi.c);
    free(data->B_r.c);
    free(data->B_phi.c);
    free(data->B_z.c);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void B_3DS_offload(B_3DS_data* data) {
    GPU_MAP_TO_DEVICE(
        data->psi, data->B_r, data->B_phi, data->B_z, \
        data->psi.c[0:data->psi.n_x*data->psi.n_y*NSIZE_COMP2D], \
        data->B_r.c[0:data->B_r.n_x*data->B_r.n_y*data->B_r.n_z*NSIZE_COMP3D], \
        data->B_phi.c[0:data->B_phi.n_x*data->B_phi.n_y*data->B_phi.n_z*NSIZE_COMP3D], \
        data->B_z.c[0:data->B_z.n_x*data->B_z.n_y*data->B_z.n_z*NSIZE_COMP3D]
    )
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * @param psi pointer where psi [V*s*m^-1] value will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_psi(real* psi, real r, real phi, real z,
                   B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += interp2Dcomp_eval_f(&psi[0], &Bdata->psi, r, z);

    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
    }

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * @param psi_dpsi pointer for storing psi [V*s*m^-1] and its derivatives
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                   B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0;
    real psi_dpsi_temp[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi_temp, &Bdata->psi, r, z);

    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
    }

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho and its derivatives
 *
 * @param rho_drho pointer where rho and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_3DS_data* Bdata) {
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (psi_dpsi[0] - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DS );
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt(fabs((psi_dpsi[0] - Bdata->psi0) / delta));

    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[2] / (2*delta*rho_drho[0]);

    return 0;
}

/**
 * @brief Evaluate magnetic field
 *
 * @param B pointer to array where magnetic field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_B(real B[3], real r, real phi, real z, B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0;

    interperr += interp3Dcomp_eval_f(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_f(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_f(&B[2], &Bdata->B_z, r, phi, z);

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
    }

    if(!err) {
        real psi_dpsi[6];
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

        B[0] = B[0] - psi_dpsi[2]/r;
        B[2] = B[2] + psi_dpsi[1]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DS );
    }

    return err;
}

/**
 * @brief Evaluate magnetic field and its derivatives
 *
 * @param B_dB pointer to array where the field and its derivatives are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_r, r, phi, z);
    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_phi, r, phi, z);
    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_z, r, phi, z);
    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
    }

    if(!err) {
        real psi_dpsi[6];
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

        B_dB[0] = B_dB[0] - psi_dpsi[2]/r;
        B_dB[1] = B_dB[1] + psi_dpsi[2]/(r*r)-psi_dpsi[5]/r;
        B_dB[3] = B_dB[3] - psi_dpsi[4]/r;
        B_dB[8] = B_dB[8] + psi_dpsi[1]/r;
        B_dB[9] = B_dB[9] - psi_dpsi[1]/(r*r) + psi_dpsi[3]/r;
        B_dB[11] = B_dB[11] + psi_dpsi[5]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DS );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DS );
    }

    return err;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param rz pointer where axis R and z [m] values will be stored
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Zero a5err value as this function can't fail.
 */
a5err B_3DS_get_axis_rz(real rz[2], B_3DS_data* Bdata) {
    a5err err = 0;
    rz[0] = Bdata->axis_r;
    rz[1] = Bdata->axis_z;
    return err;
}
