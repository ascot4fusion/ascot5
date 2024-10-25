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
 * @brief Initialize magnetic field data
 *
 * @param n_r number of r grid points
 * @param r_min minimum R coordinate in the grid [m]
 * @param r_max maximum R coordinate in the grid [m]
 * @param n_z number of z grid points
 * @param z_min minimum z coordinate in the grid [m]
 * @param z_max maximum z coordinate in the grid [m]
 * @param axis_r R coordinate of magnetic axis [m]
 * @param axis_z z coordinate of magnetic axis [m]
 * @param psi0 poloidal flux at magnetic axis [Vs/m]
 * @param psi1 poloidal flux at separatrix [Vs/m]
 * @param psi poloidal flux psi(R_i,z_j) = arr[j*n_r + i] [Vs/m]
 * @param B_r Magnetic field R component B_r(R_i,z_j) = arr[j*n_r + i] [T]
 * @param B_phi Magnetic field phi component B_phi(R_i,z_j) = arr[j*n_r + i] [T]
 * @param B_z Magnetic field z component B_z(R_i,z_j) = arr[j*n_r + i] [T]
 *
 * @return zero if initialization succeeded
 */
int B_2DS_init(B_2DS_data* data,
               int n_r, real r_min, real r_max,
               int n_z, real z_min, real z_max,
               real axis_r, real axis_z, real psi0, real psi1,
               real* psi, real* B_r, real* B_phi, real* B_z) {

    int err = 0;
    data->psi0   = psi0;
    data->psi1   = psi1;
    data->axis_r = axis_r;
    data->axis_z = axis_z;

    /* Set up the splines */
    err = interp2Dcomp_setup(&data->psi, psi, n_r, n_z, NATURALBC, NATURALBC,
                             r_min, r_max, z_min, z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }
    err = interp2Dcomp_setup(&data->B_r, B_r, n_r, n_z, NATURALBC, NATURALBC,
                             r_min, r_max, z_min, z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }
    err = interp2Dcomp_setup(&data->B_phi, B_phi, n_r, n_z,
                             NATURALBC, NATURALBC, r_min, r_max, z_min, z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }
    err = interp2Dcomp_setup(&data->B_z, B_z, n_r, n_z, NATURALBC, NATURALBC,
                             r_min, r_max, z_min, z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }

    /* Evaluate psi and magnetic field on axis for checks */
    real psival[1], Bval[3];
    err = B_2DS_eval_psi(psival, data->axis_r, 0, data->axis_z, &data);
    err = B_2DS_eval_B(Bval, data->axis_r, 0, data->axis_z, &data);
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
              n_r, r_min, r_max,
              n_z, z_min, z_max,
              data->axis_r, data->axis_z,
              data[0], data->psi0,
              Bval[0], Bval[1], Bval[2]);

    return err;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to data struct
 */
void B_2DS_free(B_2DS_data* data) {
    free(data->psi.c);
    free(data->B_r.c);
    free(data->B_phi.c);
    free(data->B_z.c);
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
    rho_drho[0] = sqrt((psi_dpsi[0] - Bdata->psi0) / delta);
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
