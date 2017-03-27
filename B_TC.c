/**
 * @author Konsta Sarkimaki konsta.sarkimaki@aalto.fi
 * @file B_TC.c
 * @brief Trivial Cartesian magnetic field
 *
 * !!THIS MODULE IS INTENDED PURELY FOR TESTING AND DEBUGGING INTEGRATORS!!
 * 
 * Representation of a magnetic field in Cartesian system.
 * The magnetic field is defined everywhere with user-defined B0 being
 * in the origo. The gradient of the magnetic field is constant everywhere.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "math.h"
#include "B_TC.h"

/* Define the magnetic field here. 
 * Note that these are in Cartesian basis. 
 * The default is B_x(0,0,0) = 1 T dB_x/dy = 1 T/m.
 */

#pragma omp declare target
static const int B_TC_B[3] = {1.0, 0.0, 0.0};
static const int B_TC_dB[9] = {0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
#pragma omp end declare target

/**
 * @brief Dummy function
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_TC_init_offload(B_TC_offload_data* offload_data, real** offload_array) {
    
}

/**
 * @brief Dummy function
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_TC_free_offload(B_TC_offload_data* offload_data,
                       real** offload_array) {
    
}

/**
 * @brief Initialize magnetic field data struct
 *
 * @param Bdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void B_TC_init(B_TC_data* Bdata, B_TC_offload_data* offload_data,
               real* offload_array) {
    Bdata->B[0] = B_TC_B[0];
    Bdata->B[1] = B_TC_B[1];
    Bdata->B[2] = B_TC_B[2];

    Bdata->dB[0] = B_TC_dB[0];
    Bdata->dB[1] = B_TC_dB[1];
    Bdata->dB[2] = B_TC_dB[2];
    Bdata->dB[3] = B_TC_dB[3];
    Bdata->dB[4] = B_TC_dB[4];
    Bdata->dB[5] = B_TC_dB[5];
    Bdata->dB[6] = B_TC_dB[6];
    Bdata->dB[7] = B_TC_dB[7];
    Bdata->dB[8] = B_TC_dB[8];
}

/**
 * @brief Evaluate magnetic field 
 *
 * This function evaluates the analytic magnetic field at the given coordinates.
 * This is a SIMD function, so the values are placed in an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated 
 * @param B array where magnetic field values will be stored (Br -> B[0][i],
 *          Bphi -> B[1][i], Bz -> B[2][i]
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Replace the divisions by powers of R0 by setting r/=R0, z/=R0 and 
 *       psi_mult/=R0 to reduce the number of computations.
 */
void B_TC_eval_B(real* B, real r, real phi,
                 real z, B_TC_data* Bdata) {
    /* Find the Cartesian position and evaluate the field there */
    real xyz[3];
    real rpz[3] = {r, phi, z};
    math_rpz2xyz(rpz, xyz);

    real Bxyz[3];
    Bxyz[0] = Bdata->B[0] + Bdata->dB[0]*xyz[0] + Bdata->dB[1]*xyz[1] + Bdata->dB[2]*xyz[2];
    Bxyz[1] = Bdata->B[1] + Bdata->dB[3]*xyz[0] + Bdata->dB[4]*xyz[1] + Bdata->dB[5]*xyz[2];
    Bxyz[2] = Bdata->B[2] + Bdata->dB[6]*xyz[0] + Bdata->dB[7]*xyz[1] + Bdata->dB[8]*xyz[2];

    /* Transform the Cartesian field vector to cylindrical coordinates */
    math_vec_xyz2rpz(Bxyz, B, phi);

}

/**
 * @brief Evaluate poloidal flux psi
 * 
 * This function evaluates the poloidal flux psi at the given coordinates from
 * the analytic solution. This is a SIMD function, so the values are placed in 
 * an NSIMD length struct.
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
void B_TC_eval_psi(real* psi, real r, real phi, real z,
                   B_TC_data* Bdata) {
    printf("There is no psi in this field.");
    exit(EXIT_FAILURE);
}

/**
 * @brief Evaluate poloidal flux psi and derivatives
 *
 */
void B_TC_eval_psi_dpsi(real* psi, real r, real phi, real z,
                   B_TC_data* Bdata) {
    printf("There is no psi in this field.");
    exit(EXIT_FAILURE);
}

/**
 * @brief Evaluate radial coordinate rho
 *
 * This function is identical to B_2D_eval_rho.
 *
 * @see B_2D_eval_rho
 */
void B_TC_eval_rho(real* rho, real psi, B_TC_data* Bdata) {
    printf("There is no rho in this field.");
    exit(EXIT_FAILURE);
}

/**
 * @brief Evaluate radial coordinate rho and derivatives
 *
 */
void B_TC_eval_rho_drho(real* rho, real r, real phi, real z, B_TC_data* Bdata) {
    printf("There is no rho in this field.");
    exit(EXIT_FAILURE);
}


/**
 * @brief Evaluate magnetic field 
 *
 * This function evaluates the analytic magnetic field and it's derivatives at 
 * the given coordinates. This is a SIMD function, so the values are placed in 
 * an NSIMD length struct.
 *
 * @param i index in the NSIMD struct that will be populated
 * @param B_dB array where magnetic field values will be stored (Br -> B[0][i],
 *        dBr/dr -> B[1][i], dBr/dphi -> B[2][i], dBr/dz -> B[3][i],
 *        Bphi -> B[4][i], dBphi/dr -> B[5][i], dBphi/dphi -> B[6][i],
 *        dBphi/dz -> B[7][i], Bz -> B[8][i], dBz/dr -> B[9][i],
 *        dBz/dphi -> B[10][i], dBz/dz -> B[11][i])
 * @param r r coordinate
 * @param phi phi coordinate
 * @param z z coordinate
 * @param Bdata pointer to magnetic field data struct
 *
 * @todo Replace the divisions by powers of R0 by setting r/=R0, z/=R0 and 
 *       psi_mult/=R0 to reduce the number of computations.
 */
void B_TC_eval_B_dB(real* B_dB, real r, real phi, real z,
                    B_TC_data* Bdata) {
    /* Find the Cartesian position and evaluate the field there */
    real xyz[3];
    real rpz[3] = {r, phi, z};
    math_rpz2xyz(rpz, xyz);

    real Bxyz[3];
    Bxyz[0] = Bdata->B[0] + Bdata->dB[0]*xyz[0] + Bdata->dB[1]*xyz[1] + Bdata->dB[2]*xyz[2];
    Bxyz[1] = Bdata->B[1] + Bdata->dB[3]*xyz[0] + Bdata->dB[4]*xyz[1] + Bdata->dB[5]*xyz[2];
    Bxyz[2] = Bdata->B[2] + Bdata->dB[6]*xyz[0] + Bdata->dB[7]*xyz[1] + Bdata->dB[8]*xyz[2];

    /* Transform the Cartesian field vector and Jacobian to cylindrical coordinates */
    real Brpz[3];
    math_vec_xyz2rpz(Bxyz, Brpz, phi);
    B_dB[0] = Brpz[0];
    B_dB[4] = Brpz[1];
    B_dB[8] = Brpz[2];

    real c = cos(phi);
    real s = sin(phi);

    B_dB[1] = c*Bdata->dB[0]+s*Bdata->dB[1];
    B_dB[2] = -r*s*Bdata->dB[0]+r*c*Bdata->dB[1];
    B_dB[3] = Bdata->dB[2];

    B_dB[5] = c*Bdata->dB[3]+s*Bdata->dB[4];
    B_dB[6] = -r*s*Bdata->dB[3]+r*c*Bdata->dB[4];
    B_dB[7] = Bdata->dB[5];

    B_dB[9]  = c*Bdata->dB[6]+s*Bdata->dB[7];
    B_dB[11] = -r*s*Bdata->dB[6]+r*c*Bdata->dB[7];
    B_dB[10] = Bdata->dB[8];
}

real B_TC_get_axis_r(B_TC_data* Bdata) {
    printf("There is no axis in this field.");
    exit(EXIT_FAILURE);
    return 0.0;
}

real B_TC_get_axis_z(B_TC_data* Bdata) {
    printf("There is no axis in this field.");
    exit(EXIT_FAILURE);
    return 0.0;
}
