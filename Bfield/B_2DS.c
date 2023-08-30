/**
 * @file B_2DS.c
 * @brief 2D magnetic field with bicubic spline interpolation
 *
 * This module represents a magnetic field where data is given in \f$Rz\f$-grid
 * from which it is interpolated with bicubic splines. The field is
 * axisymmetric.
 *
 * The magnetic field is evaluated from magnetic field strength \f$\mathbf{B}\f$
 * which may not be divergence free. However, \f$B_R\f$ and \f$B_z\f$ components
 * are also evaluated from poloidal magnetic flux \f$\psi\f$ as
 *
 * \f{align*}{
 * B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
 * B_z &= \frac{1}{R}\frac{\partial\psi}{\partial R}
 * \f}
 *
 * The total field is then a sum of components interpolated directly from
 * \f$\mathbf{B}\f$ and components calculated via interpolated \f$\psi\f$.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown.
 *
 * @see B_field.c
 */
#include <stdlib.h>
#include <math.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "B_2DS.h"
#include "../spline/interp.h"

/**
 * @brief Initialize magnetic field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct and to construct
 * bicubic splines whose coefficients are stored in re-allocated offload array.
 *
 * The offload data struct must have the following fields initialized:
 * - B_2DS_offload_data.n_r
 * - B_2DS_offload_data.n_z
 * - B_2DS_offload_data.r_min
 * - B_2DS_offload_data.r_max
 * - B_2DS_offload_data.z_min
 * - B_2DS_offload_data.z_max
 * - B_2DS_offload_data.psi0
 * - B_2DS_offload_data.psi1
 * - B_2DS_offload_data.axis_r
 * - B_2DS_offload_data.axis_z
 *
 * B_2DS_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[            j*n_r + i] = psi(R_i, z_j)   [V*s*m^-1]
 * - offload_array[  n_r*n_z + j*n_r + i] = B_R(R_i, z_j)   [T]
 * - offload_array[2*n_r*n_z + j*n_r + i] = B_phi(R_i, z_j) [T]
 * - offload_array[3*n_r*n_z + j*n_r + i] = B_z(R_i, z_j)   [T]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int B_2DS_init_offload(B_2DS_offload_data* offload_data,
                       real** offload_array) {

    /* Spline initialization */
    int err = 0;
    int datasize = offload_data->n_r*offload_data->n_z;

    /* Allocate enough space to store four 2D arrays */
    real* coeff_array = (real*) malloc(4*NSIZE_COMP2D*datasize*sizeof(real));
    real* psi   = &(coeff_array[0*datasize*NSIZE_COMP2D]);
    real* B_r   = &(coeff_array[1*datasize*NSIZE_COMP2D]);
    real* B_phi = &(coeff_array[2*datasize*NSIZE_COMP2D]);
    real* B_z   = &(coeff_array[3*datasize*NSIZE_COMP2D]);

    /* Evaluate spline coefficients */
    err += interp2Dcomp_init_coeff(
        psi, *offload_array + 0*datasize,
        offload_data->n_r, offload_data->n_z,
        NATURALBC, NATURALBC,
        offload_data->r_min, offload_data->r_max,
        offload_data->z_min, offload_data->z_max);

    err += interp2Dcomp_init_coeff(
        B_r, *offload_array + 1*datasize,
        offload_data->n_r, offload_data->n_z,
        NATURALBC, NATURALBC,
        offload_data->r_min, offload_data->r_max,
        offload_data->z_min, offload_data->z_max);

    err += interp2Dcomp_init_coeff(
        B_phi, *offload_array + 2*datasize,
        offload_data->n_r, offload_data->n_z,
        NATURALBC, NATURALBC,
        offload_data->r_min, offload_data->r_max,
        offload_data->z_min, offload_data->z_max);

    err += interp2Dcomp_init_coeff(
        B_z, *offload_array + 3*datasize,
        offload_data->n_r, offload_data->n_z,
        NATURALBC, NATURALBC,
        offload_data->r_min, offload_data->r_max,
        offload_data->z_min, offload_data->z_max);

    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Free offload array and and replace it with the coefficient array */
    free(*offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = 4*NSIZE_COMP2D*datasize;

    /* Initialization complete. Check that the data seem valid. */

    /* Evaluate psi and magnetic field on axis for checks */
    B_2DS_data Bdata;
    B_2DS_init(&Bdata, offload_data, *offload_array);
    real psival[1], Bval[3];
    err = B_2DS_eval_psi(psival, offload_data->axis_r, 0, offload_data->axis_z,
                         &Bdata);
    err = B_2DS_eval_B(Bval, offload_data->axis_r, 0, offload_data->axis_z,
                       &Bdata);
    if(err) {
        print_err("Error: Initialization failed.\n");
        return err;
    }

    /* Print some sanity check on data */
    print_out(VERBOSE_IO,
              "\n2D magnetic field (B_2DS)\n"
              "Grid: nR = %4.d Rmin = %3.3f Rmax = %3.3f\n"
              "      nz = %4.d zmin = %3.3f zmax = %3.3f\n"
              "Psi at magnetic axis (%1.3f m, %1.3f m)\n"
              "%3.3f (evaluated)\n%3.3f (given)\n"
              "Magnetic field on axis:\n"
              "B_R = %3.3f B_phi = %3.3f B_z = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max,
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max,
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
void B_2DS_free_offload(B_2DS_offload_data* offload_data,
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
void B_2DS_init(B_2DS_data* Bdata, B_2DS_offload_data* offload_data,
                real* offload_array) {

    int splinesize = NSIZE_COMP2D * offload_data->n_r * offload_data->n_z;

    /* Initialize target data struct */
    Bdata->psi0   = offload_data->psi0;
    Bdata->psi1   = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;

    /* Copy parameters and assign pointers to offload array to initialize the
       spline structs */
    interp2Dcomp_init_spline(&Bdata->psi, &(offload_array[0*splinesize]),
                             offload_data->n_r,
                             offload_data->n_z,
                             NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->z_min,
                             offload_data->z_max);

    interp2Dcomp_init_spline(&Bdata->B_r, &(offload_array[1*splinesize]),
                             offload_data->n_r,
                             offload_data->n_z,
                             NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->z_min,
                             offload_data->z_max);

    interp2Dcomp_init_spline(&Bdata->B_phi, &(offload_array[2*splinesize]),
                             offload_data->n_r,
                             offload_data->n_z,
                             NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->z_min,
                             offload_data->z_max);

    interp2Dcomp_init_spline(&Bdata->B_z, &(offload_array[3*splinesize]),
                             offload_data->n_r,
                             offload_data->n_z,
                             NATURALBC, NATURALBC,
                             offload_data->r_min,
                             offload_data->r_max,
                             offload_data->z_min,
                             offload_data->z_max);
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
a5err B_2DS_eval_psi(real* psi, real r, real phi, real z, B_2DS_data* Bdata) {

    int interperr = 0;
    interperr += interp2Dcomp_eval_f(&psi[0], &Bdata->psi, r, z);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
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
a5err B_2DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_2DS_data* Bdata) {

    int interperr = 0;
    real psi_dpsi_temp[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi_temp, &Bdata->psi, r, z);
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

     a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
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
a5err B_2DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_2DS_data* Bdata) {
    int interperr = 0;
    real psi_dpsi[6];

    interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

    a5err err = 0;
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( !err && (psi_dpsi[0] - Bdata->psi0) / delta < 0 ) {
         err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt(fabs((psi_dpsi[0] - Bdata->psi0) / delta));

    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[2] / (2*delta*rho_drho[0]);

    return err;
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
a5err B_2DS_eval_B(real B[3], real r, real phi, real z, B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0;

    interperr += interp2Dcomp_eval_f(&B[0], &Bdata->B_r, r, z);
    interperr += interp2Dcomp_eval_f(&B[1], &Bdata->B_phi, r, z);
    interperr += interp2Dcomp_eval_f(&B[2], &Bdata->B_z, r, z);

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }

    if(!err) {
        real psi_dpsi[6];
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);
        B[0] = B[0] - psi_dpsi[2]/r;
        B[2] = B[2] + psi_dpsi[1]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
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
a5err B_2DS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0;
    real B_dB_temp[6];

    interperr += interp2Dcomp_eval_df(B_dB_temp, &Bdata->B_r, r, z);

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = 0;
    B_dB[3] = B_dB_temp[2];

    interperr += interp2Dcomp_eval_df(B_dB_temp, &Bdata->B_phi, r, z);

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = 0;
    B_dB[7] = B_dB_temp[2];

    interperr += interp2Dcomp_eval_df(B_dB_temp, &Bdata->B_z, r, z);

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = 0;
    B_dB[11] = B_dB_temp[2];

    /* Test for B field interpolation error */
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }


    real psi_dpsi[6];

    if(!err) {
        interperr += interp2Dcomp_eval_df(psi_dpsi, &Bdata->psi, r, z);

        B_dB[0] = B_dB[0] - psi_dpsi[2]/r;
        B_dB[1] = B_dB[1] + psi_dpsi[2]/(r*r)-psi_dpsi[5]/r;
        B_dB[3] = B_dB[3] - psi_dpsi[4]/r;
        B_dB[8] = B_dB[8] + psi_dpsi[1]/r;
        B_dB[9] = B_dB[9] - psi_dpsi[1]/(r*r) + psi_dpsi[3]/r;
        B_dB[11] = B_dB[11] + psi_dpsi[5]/r;

        /* Test for psi interpolation error */
        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
        }
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {
        err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
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
a5err B_2DS_get_axis_rz(real rz[2], B_2DS_data* Bdata) {
    a5err err = 0;
    rz[0] = Bdata->axis_r;
    rz[1] = Bdata->axis_z;
    return err;
}
