/**
 * @file B_TC.c
 * @brief Trivial Cartesian magnetic field
 *
 * Magnetic field whose Cartesian components \f$(B_x, B_y, B_z)\f$ have
 * user-defined values at the origo \f$(R,z) = (0,0)\f$, and constant Jacobian
 *
 * \f[
 * \begin{bmatrix}
 * \frac{\partial B_x}{\partial x}, &\frac{\partial B_x}{\partial y},
 * &\frac{\partial B_x}{\partial z} \\
 * \frac{\partial B_y}{\partial x}, &\frac{\partial B_y}{\partial y},
 * &\frac{\partial B_y}{\partial z} \\
 * \frac{\partial B_z}{\partial x}, &\frac{\partial B_z}{\partial y},
 * &\frac{\partial B_z}{\partial z}
 * \end{bmatrix}
 * \f]
 *
 * This magnetic field is intended to serve debugging and testing purposes.
 * Other values \f$\psi\f$, \f$\rho\f$ and magnetic axis coordinates have
 * fixed user-defined values.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"
#include "../math.h"
#include "../error.h"
#include "../print.h"
#include "B_TC.h"

/**
 * @brief Initialize magnetic field data
 *
 * @param data pointer to the data struct
 * @param axisr value returned when quering magnetic axis R coordinate [m]
 * @param axisz value returned when quering magnetic axis z coordinate [m]
 * @param psival value returned when quering magnetic flux [Vs/m]
 * @param rhoval value returned when quering normalized poloidal flux [1]
 * @param B magnetic field at origo [B_x, B_y, B_z] [T]
 * @param dB magnetic field Jacobian [dB_x/dx, dB_x/dy, dB_x/dz, dB_y/dx,
 *        dB_y/dy, dB_y/dz, dB_z/dx, dB_z/dy, dB_z/dz] [T/m]
 *
 * @return zero to indicate success
 */
int B_TC_init(B_TC_data* data, real axisr, real axisz, real psival, real rhoval,
              real B[3], real dB[9]) {

    data->axisr = axisr;
    data->axisz = axisz;
    data->psival = psival;
    data->rhoval = rhoval;
    data->B[0] = B[0];
    data->B[1] = B[1];
    data->B[2] = B[2];
    for(int i = 0; i < 9; i++) {
        data->dB[i] = dB[i];
    }

    print_out(VERBOSE_IO,
              "\nTrivial cartesian magnetic field (B_TC)\n"
              "Magnetic axis at (R,z) = (%.1f,%.1f)\n"
              "psi = %.1f, rho = %.1f\n"
              "Magnetic field at origo\n"
              "B_x = %.3f B_y = %.3f B_z = %.3f\n"
              "Magnetic field gradient\n"
              "dB_x/dx = %.3f dB_x/dy = %.3f B_x/dz = %.3f\n"
              "dB_y/dx = %.3f dB_y/dy = %.3f B_y/dz = %.3f\n"
              "dB_z/dx = %.3f dB_z/dy = %.3f B_z/dz = %.3f\n",
              data->axisr, data->axisz,
              data->psival, data->rhoval,
              data->B[0],  data->B[1],  data->B[2],
              data->dB[0], data->dB[1], data->dB[2],
              data->dB[3], data->dB[4], data->dB[5],
              data->dB[6], data->dB[7], data->dB[8]);

    return 0;
}

/**
 * @brief Free allocated resources
 */
void B_TC_free(B_TC_data* data) {
    // No resources were dynamically allocated
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void B_TC_offload(B_TC_data* data) {
    GPU_MAP_TO_DEVICE( data->B[0:3], data->dB[0:9] )
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
 * @return zero to indicate success
 */
a5err B_TC_eval_psi(real* psi, real r, real phi, real z,
                   B_TC_data* Bdata) {
    psi[0] = Bdata->psival;

    return 0;
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
 * @return zero to indicate success
 */
a5err B_TC_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                         B_TC_data* Bdata) {
    psi_dpsi[0] = Bdata->psival;
    psi_dpsi[1] = 0;
    psi_dpsi[2] = 0;
    psi_dpsi[3] = 0;

    return 0;
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
 * @return zero to indicate success
 */
a5err B_TC_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                         B_TC_data* Bdata) {
    rho_drho[0] = Bdata->rhoval;

    rho_drho[1] = 0;
    rho_drho[2] = 0;
    rho_drho[3] = 0;

    return 0;
}

/**
 * @brief Evaluate magnetic field
 *
 * @param B pointer to array where magnetic field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return zero to indicate success
 */
a5err B_TC_eval_B(real B[3], real r, real phi,
                 real z, B_TC_data* Bdata) {
    /* Find the Cartesian position and evaluate the field there */
    real xyz[3];
    real rpz[3] = {r, phi, z};
    math_rpz2xyz(rpz, xyz);

    real Bxyz[3];
    Bxyz[0] = Bdata->B[0] + Bdata->dB[0]*xyz[0] + Bdata->dB[1]*xyz[1]
        + Bdata->dB[2]*xyz[2];
    Bxyz[1] = Bdata->B[1] + Bdata->dB[3]*xyz[0] + Bdata->dB[4]*xyz[1]
        + Bdata->dB[5]*xyz[2];
    Bxyz[2] = Bdata->B[2] + Bdata->dB[6]*xyz[0] + Bdata->dB[7]*xyz[1]
        + Bdata->dB[8]*xyz[2];

    /* Transform the Cartesian field vector to cylindrical coordinates */
    math_vec_xyz2rpz(Bxyz, B, phi);

    return 0;
}

/**
 * @brief Evaluate magnetic field and its derivatives
 *
 * @param B_dB pointer to array where the field and its derivatives are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return zero to indicate success
 */
a5err B_TC_eval_B_dB(real B_dB[12], real r, real phi, real z,
                     B_TC_data* Bdata) {
    /* Find the Cartesian position and evaluate the field there */
    real xyz[3];
    real rpz[3] = {r, phi, z};
    math_rpz2xyz(rpz, xyz);

    real Bxyz[3];
    Bxyz[0] = Bdata->B[0] + Bdata->dB[0]*xyz[0] + Bdata->dB[1]*xyz[1]
        + Bdata->dB[2]*xyz[2];
    Bxyz[1] = Bdata->B[1] + Bdata->dB[3]*xyz[0] + Bdata->dB[4]*xyz[1]
        + Bdata->dB[5]*xyz[2];
    Bxyz[2] = Bdata->B[2] + Bdata->dB[6]*xyz[0] + Bdata->dB[7]*xyz[1]
        + Bdata->dB[8]*xyz[2];

    /* Transform the Cartesian field vector and Jacobian
       to cylindrical coordinates */
    real Brpz[3];
    math_vec_xyz2rpz(Bxyz, Brpz, phi);
    B_dB[0] = Brpz[0];
    B_dB[4] = Brpz[1];
    B_dB[8] = Brpz[2];

    real B_dBxyz[12] = {Bdata->B[0], Bdata->dB[0], Bdata->dB[1], Bdata->dB[2],
                        Bdata->B[1], Bdata->dB[3], Bdata->dB[4], Bdata->dB[5],
                        Bdata->B[2], Bdata->dB[6], Bdata->dB[7], Bdata->dB[8]};
    math_jac_xyz2rpz(B_dBxyz, B_dB, r, phi);

    return 0;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param rz pointer where axis R and z [m] values will be stored
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Zero a5err value as this function can't fail.
 */
a5err B_TC_get_axis_rz(real rz[2], B_TC_data* Bdata) {
    a5err err = 0;
    rz[0] = Bdata->axisr;
    rz[1] = Bdata->axisz;
    return err;
}
