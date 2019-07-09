/**
 * @file B_3DST.c
 * @brief Time dependent 3D magnetic field with 4D spline interpolation
 *
 * This module represents a magnetic field where data is given in \f$R\phi z t\f$-
 * grid from which it is interpolated with 4-parameters cubic splines.
 *
 * The magnetic field is evaluated from magnetic field strength \f$\mathbf{B}\f$
 * which may not be divergence free. However, \f$B_R\f$ and \f$B_z\f$ components
 * are also evaluated from poloidal magnetic flux \f$\psi(R,z)\f$ as
 *
 * \f{align*}{
 * B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
 * B_z &= \frac{1}{R}\frac{\partial\psi}{\partial R}
 * \f}
 *
 * The total field is then a sum of components interpolated directly from
 * \f$\mathbf{B}\f$ and components calculated via interpolated \f$\psi\f$.
 * Note that \f$\psi\f$ is assumed to be axisymmetric and is interpolated with
 * bicubic splines. \f$\psi(R,z)\f$ is considered constant in the time domain.
 * \f$\psi\f$ and \f$\mathbf{B}\f$ are given in separate grids.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown. For \f$\phi\f$-grid, periodic boundary
 * conditions are used but it is user's responsibility to provide input
 * whose \f$\phi\f$-grid makes sense (in that it actually represents a periodic
 * field), i.e., \f$\phi_\mathrm{max}-\phi_\mathrm{min} = 2\pi/(N+1)\f$.
 * However, do note that in this module \f$\phi_\mathrm{max}\f$ is not the
 * "last" grid point but the second last, e.g. if \f$\phi_\mathrm{min}=0\f$
 * and \f$n_\phi = 360\f$, then \f$\phi_\mathrm{max}=359\f$ if periodicity is
 * \f$N=0\f$. This module does not extrapolate in the time domain so an error
 * will be thrown or the marker will be stopped.
 *
 * @see B_field.c
 */
#include <stdlib.h>
#include <math.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "B_3DST.h"
#include "../spline/interp.h"

/**
 * @brief Initialize magnetic field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct and to construct
 * cubic splines whose coefficients are stored in re-allocated offload array.
 *
 * The offload data struct must have the following fields initialized:
 * - B_3DST_offload_data.psigrid_n_r
 * - B_3DST_offload_data.psigrid_n_z
 * - B_3DST_offload_data.psigrid_r_min
 * - B_3DST_offload_data.psigrid_r_max
 * - B_3DST_offload_data.psigrid_z_min
 * - B_3DST_offload_data.psigrid_z_max
 *
 * - B_3DST_offload_data.Bgrid_n_r
 * - B_3DST_offload_data.Bgrid_n_z
 * - B_3DST_offload_data.Bgrid_r_min
 * - B_3DST_offload_data.Bgrid_r_max
 * - B_3DST_offload_data.Bgrid_z_min
 * - B_3DST_offload_data.Bgrid_z_max
 * - B_3DST_offload_data.Bgrid_n_phi
 * - B_3DST_offload_data.Bgrid_phi_min
 * - B_3DST_offload_data.Bgrid_phi_max
 * - B_3DST_offload_data.Bgrid_n_t
 * - B_3DST_offload_data.Bgrid_t_min
 * - B_3DST_offload_data.Bgrid_t_max
 *
 * - B_3DST_offload_data.psi0
 * - B_3DST_offload_data.psi1
 * - B_3DST_offload_data.axis_r
 * - B_3DST_offload_data.axis_z
 *
 * B_3DST_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[                           m*Bn_r*Bn_z*Bn_phi + k*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_R(R_i, phi_k, z_j, t_m)   [T]
 * - offload_array[Bn_r*Bn_z*Bn_phi*Bn_t + m*Bn_r*Bn_z*Bn_phi + k*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_phi(R_i, phi_k, z_j, t_m)   [T]
 * - offload_array[2*Bn_r*Bn_z*Bn_phi*Bn_t + m*Bn_r*Bn_z*Bn_phi + k*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_z(R_i, phi_k, z_j, t_m)   [T]
 * - offload_array[3*Bn_r*Bn_z*Bn_phi*Bn_t + j*n_r + i]
 *   = psi(R_i, z_j)   [V*s*m^-1]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int B_3DST_init_offload(B_3DST_offload_data* offload_data, real** offload_array) {

    /* Spline initialization. */
    int err = 0;
    int psi_size = offload_data->psigrid_n_r * offload_data->psigrid_n_z;
    int B_size   = offload_data->Bgrid_n_r   * offload_data->Bgrid_n_z
                   * offload_data->Bgrid_n_phi * offload_data->Bgrid_n_t;

    /* Allocate enough space to store three 3D arrays and one 2D array */
    real* coeff_array = (real*) malloc( (3*NSIZE_COMP4D*B_size
                                         + NSIZE_COMP2D*psi_size)*sizeof(real));
    real* B_r   = &(coeff_array[0*B_size*NSIZE_COMP4D]);
    real* B_phi = &(coeff_array[1*B_size*NSIZE_COMP4D]);
    real* B_z   = &(coeff_array[2*B_size*NSIZE_COMP4D]);
    real* psi   = &(coeff_array[3*B_size*NSIZE_COMP4D]);

    err += interp2Dcomp_init_coeff(
        psi, *offload_array + 3*B_size,
        offload_data->psigrid_n_r, offload_data->psigrid_n_z,
        NATURALBC, NATURALBC,
        offload_data->psigrid_r_min, offload_data->psigrid_r_max,
        offload_data->psigrid_z_min, offload_data->psigrid_z_max);

    err += interp4Dcomp_init_coeff(
        B_r, *offload_array + 0*B_size,
        offload_data->Bgrid_n_r, offload_data->Bgrid_n_phi,
        offload_data->Bgrid_n_z, offload_data->Bgrid_n_t,
        NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
        offload_data->Bgrid_r_min,   offload_data->Bgrid_r_max,
        offload_data->Bgrid_phi_min, offload_data->Bgrid_phi_max,
        offload_data->Bgrid_z_min,   offload_data->Bgrid_z_max,
        offload_data->Bgrid_t_min,   offload_data->Bgrid_t_max);

    err += interp4Dcomp_init_coeff(
        B_phi, *offload_array + 1*B_size,
        offload_data->Bgrid_n_r, offload_data->Bgrid_n_phi,
        offload_data->Bgrid_n_z, offload_data->Bgrid_n_t,
        NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
        offload_data->Bgrid_r_min,   offload_data->Bgrid_r_max,
        offload_data->Bgrid_phi_min, offload_data->Bgrid_phi_max,
        offload_data->Bgrid_z_min,   offload_data->Bgrid_z_max,
        offload_data->Bgrid_t_min,   offload_data->Bgrid_t_max);

    err += interp4Dcomp_init_coeff(
        B_z, *offload_array + 2*B_size,
        offload_data->Bgrid_n_r, offload_data->Bgrid_n_phi,
        offload_data->Bgrid_n_z, offload_data->Bgrid_n_t,
        NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
        offload_data->Bgrid_r_min,   offload_data->Bgrid_r_max,
        offload_data->Bgrid_phi_min, offload_data->Bgrid_phi_max,
        offload_data->Bgrid_z_min,   offload_data->Bgrid_z_max,
        offload_data->Bgrid_t_min,   offload_data->Bgrid_t_max);

    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = NSIZE_COMP2D*psi_size
                                       + NSIZE_COMP4D*B_size*3;

    /* Evaluate psi and magnetic field on axis for checks */

    B_3DST_data Bdata;
    B_3DST_init(&Bdata, offload_data, *offload_array);
    real psival[1], Bval[3];
    err = B_3DST_eval_psi(psival, offload_data->axis_r, 0, offload_data->axis_z,
                          &Bdata);
    err = B_3DST_eval_B(Bval, offload_data->axis_r, 0, offload_data->axis_z,
                        offload_data->Bgrid_t_min, &Bdata);
    if(err) {
        print_err("Error: Initialization failed.\n");
        return err;
    }

    /* Print some sanity check on data */
    print_out(VERBOSE_IO,
              "\nTime-dependent 3D magnetic field (B_3DST)\n"
              "Psi-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n"
              "          nz = %4.d zmin = %3.3f m zmax = %3.3f m\n"
	      "B-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n"
	      "        nz = %4.d zmin = %3.3f m zmax = %3.3f m\n"
	      "        nphi = %4.d phimin = %3.3f deg phimax = %3.3f deg\n"
	      "        ntime = %4.d tmin = %3.3f s tmax = %3.3f s\n"
	      "Psi at magnetic axis (%1.3f m, %1.3f m)\n"
              "%3.3f (evaluated)\n%3.3f (given)\n"
              "Magnetic field on axis:\n"
              "B_R = %3.3f B_phi = %3.3f B_z = %3.3f\n",
              offload_data->psigrid_n_r,
              offload_data->psigrid_r_min, offload_data->psigrid_r_max,
              offload_data->psigrid_n_z,
              offload_data->psigrid_z_min, offload_data->psigrid_z_max,
              offload_data->Bgrid_n_r,
              offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
              offload_data->Bgrid_n_z,
              offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
              offload_data->Bgrid_n_phi,
              math_rad2deg(offload_data->Bgrid_phi_min),
	      math_rad2deg(offload_data->Bgrid_phi_max),
              offload_data->Bgrid_n_t,
              offload_data->Bgrid_t_min, offload_data->Bgrid_t_max,
              offload_data->axis_r, offload_data->axis_z,
              psival[0], offload_data->psi0,
              Bval[0], Bval[1], Bval[2]);

    return err;
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3DST_free_offload(B_3DST_offload_data* offload_data,
                        real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize magnetic field data struct on target
 *
 * @param Bdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array offload array
 */
void B_3DST_init(B_3DST_data* Bdata, B_3DST_offload_data* offload_data,
                real* offload_array) {

    int B_size = NSIZE_COMP4D * offload_data->Bgrid_n_r
        * offload_data->Bgrid_n_z * offload_data->Bgrid_n_phi * offload_data->Bgrid_n_t;

    /* Initialize target data struct */
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;

    /* Initialize spline structs from the coefficients */
    interp4Dcomp_init_spline(&Bdata->B_r, &(offload_array[0*B_size]),
                             offload_data->Bgrid_n_r,
                             offload_data->Bgrid_n_phi,
                             offload_data->Bgrid_n_z,
                             offload_data->Bgrid_n_t,
                             NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
                             offload_data->Bgrid_r_min,
                             offload_data->Bgrid_r_max,
                             offload_data->Bgrid_phi_min,
                             offload_data->Bgrid_phi_max,
                             offload_data->Bgrid_z_min,
                             offload_data->Bgrid_z_max,
                             offload_data->Bgrid_t_min,
                             offload_data->Bgrid_t_max);

    interp4Dcomp_init_spline(&Bdata->B_phi, &(offload_array[1*B_size]),
                             offload_data->Bgrid_n_r,
                             offload_data->Bgrid_n_phi,
                             offload_data->Bgrid_n_z,
                             offload_data->Bgrid_n_t,
                             NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
                             offload_data->Bgrid_r_min,
                             offload_data->Bgrid_r_max,
                             offload_data->Bgrid_phi_min,
                             offload_data->Bgrid_phi_max,
                             offload_data->Bgrid_z_min,
                             offload_data->Bgrid_z_max,
                             offload_data->Bgrid_t_min,
                             offload_data->Bgrid_t_max);

    interp4Dcomp_init_spline(&Bdata->B_z, &(offload_array[2*B_size]),
                             offload_data->Bgrid_n_r,
                             offload_data->Bgrid_n_phi,
                             offload_data->Bgrid_n_z,
                             offload_data->Bgrid_n_t,
                             NATURALBC, PERIODICBC, NATURALBC, NATURALBC,
                             offload_data->Bgrid_r_min,
                             offload_data->Bgrid_r_max,
                             offload_data->Bgrid_phi_min,
                             offload_data->Bgrid_phi_max,
                             offload_data->Bgrid_z_min,
                             offload_data->Bgrid_z_max,
                             offload_data->Bgrid_t_min,
                             offload_data->Bgrid_t_max);

    interp2Dcomp_init_spline(&Bdata->psi, &(offload_array[3*B_size]),
                             offload_data->psigrid_n_r,
                             offload_data->psigrid_n_z,
                             NATURALBC, NATURALBC,
                             offload_data->psigrid_r_min,
                             offload_data->psigrid_r_max,
                             offload_data->psigrid_z_min,
                             offload_data->psigrid_z_max);

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
a5err B_3DST_eval_psi(real psi[1], real r, real phi, real z,
                   B_3DST_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    interperr += interp2Dcomp_eval_f(&psi[0], &Bdata->psi, r, z);

    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
    }

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * @param psi pointer where psi [V*s*m^-1] and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DST_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                   B_3DST_data* Bdata) {
    a5err err = 0;
    int interperr = 0;
    real psi_dpsi_temp[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi_temp, &Bdata->psi, r, z);

    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
    }

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho
 *
 * @param rho pointer where rho value will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DST_eval_rho(real rho[1], real psi, B_3DST_data* Bdata) {
    a5err err = 0;

    /* Check that the values seem valid */
    real delta = (Bdata->psi1 - Bdata->psi0);
    if( (psi - Bdata->psi0) / delta < 0 ) {
         err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DST );
    }

    if(!err) {
        rho[0] = sqrt( (psi - Bdata->psi0) / delta );
    }
    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho and its derivatives
 *
 * @param rho pointer where rho and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DST_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_3DST_data* Bdata) {
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (psi_dpsi[0] - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DST );
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
a5err B_3DST_eval_B(real B[3], real r, real phi, real z, real t, B_3DST_data* Bdata) {
    a5err err = 0;
    int interperr = 0;

    interperr += interp4Dcomp_eval_f(&B[0], &Bdata->B_r, r, phi, z, t);
    interperr += interp4Dcomp_eval_f(&B[1], &Bdata->B_phi, r, phi, z, t);
    interperr += interp4Dcomp_eval_f(&B[2], &Bdata->B_z, r, phi, z, t);

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
    }

    if(!err) {
        real psi_dpsi[6];
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

        B[0] = B[0] - psi_dpsi[2]/r;
        B[2] = B[2] + psi_dpsi[1]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DST );
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
a5err B_3DST_eval_B_dB(real B_dB[], real r, real phi, real z, real t,
                      B_3DST_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];

    interperr += interp4Dcomp_eval_df(B_dB_temp, &Bdata->B_r, r, phi, z, t);
    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];

    interperr += interp4Dcomp_eval_df(B_dB_temp, &Bdata->B_phi, r, phi, z, t);
    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];

    interperr += interp4Dcomp_eval_df(B_dB_temp, &Bdata->B_z, r, phi, z, t);
    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
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
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_3DST );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_3DST );
    }

    return err;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis R-coordinate [m]
 */
real B_3DST_get_axis_r(B_3DST_data* Bdata) {
    return Bdata->axis_r;
}

/**
 * @brief Return magnetic axis z-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis z-coordinate [m]
 */
real B_3DST_get_axis_z(B_3DST_data* Bdata) {
    return Bdata->axis_z;
}
