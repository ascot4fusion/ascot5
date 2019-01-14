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
 * @brief Initialize magnetic field offload data
 *
 * The offload data struct is very simple and contains only the necessary fields
 * which are all initialized when the data is read (except
 * B_TC.offload_array_length). The offload array is not required at all.
 * Therefore, this function only sets the offload_array_length to zero, assigns
 * a NULL pointer to offload_array, and prints sanity checks so that user may
 * verify that data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @return zero to indicate success
 */
int B_TC_init_offload(B_TC_offload_data* offload_data, real** offload_array) {
    // No initialization done
    offload_data->offload_array_length = 0;
    *offload_array = NULL;

    print_out(VERBOSE_IO, "\nTrivial cartesian magnetic field (B_TC)\n");
    print_out(VERBOSE_IO, "Magnetic axis at (R,z) = (%.1f,%.1f)\n",
              offload_data->axisr, offload_data->axisz);
    print_out(VERBOSE_IO, "psi = %.1f, rho = %.1f\n",
              offload_data->psival, offload_data->rhoval);
    print_out(VERBOSE_IO, "Magnetic field at origo\n"
              "B_x = %.3f B_y = %.3f B_z = %.3f\n",
              offload_data->B[0], offload_data->B[1], offload_data->B[2]);
    print_out(VERBOSE_IO, "Magnetic field gradient\n"
              "dB_x/dx = %.3f dB_x/dy = %.3f B_x/dz = %.3f\n",
              offload_data->dB[0], offload_data->dB[1], offload_data->dB[2]);
    print_out(VERBOSE_IO,
              "dB_y/dx = %.3f dB_y/dy = %.3f B_y/dz = %.3f\n",
              offload_data->dB[3], offload_data->dB[4], offload_data->dB[5]);
    print_out(VERBOSE_IO,
              "dB_z/dx = %.3f dB_z/dy = %.3f B_z/dz = %.3f\n",
              offload_data->dB[6], offload_data->dB[7], offload_data->dB[8]);

    return 0;
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_TC_free_offload(B_TC_offload_data* offload_data,
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
void B_TC_init(B_TC_data* Bdata, B_TC_offload_data* offload_data,
               real* offload_array) {
    Bdata->rhoval = offload_data->rhoval;
    Bdata->psival = offload_data->psival;
    Bdata->axisr = offload_data->axisr;
    Bdata->axisz = offload_data->axisz;

    Bdata->B  = offload_data->B;
    Bdata->dB = offload_data->dB;

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
 * @param psi pointer where psi [V*s*m^-1] and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return zero to indicate success
 */
a5err B_TC_eval_psi_dpsi(real* psi, real r, real phi, real z,
                   B_TC_data* Bdata) {
    psi[0] = Bdata->psival;
    psi[1] = 0;
    psi[2] = 0;
    psi[3] = 0;

    return 0;
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
 * @return zero to indicate success
 */
a5err B_TC_eval_rho(real* rho, real psi, B_TC_data* Bdata) {
    rho[0] = Bdata->rhoval;

    return 0;
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
 * @return zero to indicate success
 */
a5err B_TC_eval_rho_drho(real* rho, real r, real phi, real z,
                         B_TC_data* Bdata) {
    rho[0] = Bdata->rhoval;

    rho[1] = 0;
    rho[2] = 0;
    rho[3] = 0;

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
 * @return zero to indicate success
 */
a5err B_TC_eval_B(real* B, real r, real phi,
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
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return zero to indicate success
 */
a5err B_TC_eval_B_dB(real* B_dB, real r, real phi, real z,
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
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis R-coordinate [m]
 */
real B_TC_get_axis_r(B_TC_data* Bdata) {
    return Bdata->axisr;
}

/**
 * @brief Return magnetic axis z-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis z-coordinate [m]
 */
real B_TC_get_axis_z(B_TC_data* Bdata) {
    return Bdata->axisz;
}
