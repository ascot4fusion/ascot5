/**
 * @file B_GS.c
 * @brief Analytic magnetic field
 *
 * This module implements a toroidal magnetic field based on an analytical
 * solution to the Grad-Shafranov equation [1]. The field is scaled in
 * dimensions and field strength to approximate the Iter field at 9 MA.
 *
 * [1] A.J. Cerfon, J.P. Freidberg. "One size fits all" analytic solutions to
 *     the Grad-Shafranov equation. Physics of Plasmas 17 (3) (2010) 032502. 
 *     http://scitation.aip.org/content/aip/journal/pop/17/3/10.1063/1.3328818
 */
#include <stdlib.h>
#include <math.h>
#include "B_GS.h"

#define R0 6.2 /**< @brief Major radius of the torus */
#define B_phi0 5 /**< @brief Toroidal field at the magnetic axis */
#define psi_mult 138 /**< @brief Psi multiplier yielding 9 MA plasma current */

/**
 * @brief Load magnetic field data and prepare parameters
 *
 * This function allocates the offload array and fills it with hardcoded
 * magnetic field parameters.
 * 
 * @todo Read parameters from a file
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 * @param filename file in input.magn_bkg format
 */
void B_GS_init_offload(B_GS_offload_data* offload_data, real** offload_array) {
    offload_data->psi0 = -0.0365;
    offload_data->psi1 = 0.0;
    offload_data->offload_array_length = 13;

    *offload_array = (real*) malloc(13 * sizeof(real));
    (*offload_array)[0] = 0.0862949108578025;   /* C1 ... */
    (*offload_array)[1] = 0.327930658772432;
    (*offload_array)[2] = 0.526867770124117;
    (*offload_array)[3] = -0.236620894691229;
    (*offload_array)[4] = 0.382582676559332;
    (*offload_array)[5] = -0.357315314775458;
    (*offload_array)[6] = -0.0148416683303735;
    (*offload_array)[7] = 0.150604594328656;
    (*offload_array)[8] = 0.742822645941527;
    (*offload_array)[9] = -0.444715310510484;
    (*offload_array)[10] = -0.108464039573686;
    (*offload_array)[11] = 0.0128159923595110;  /* ... C12 */
    (*offload_array)[12] = -0.1550;             /* and A */
}

/**
 * @brief Free offload array and reset parameters 
 *
 * This function deallocates the offload_array.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_GS_free_offload(B_GS_offload_data* offload_data,
                       real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize magnetic field data struct on target 
 *
 * This function copies the magnetic field parameters from the offload struct
 * to the struct on target and sets the magnetic field data pointer.
 *
 * @param BData pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void B_GS_init(B_GS_data* Bdata, B_GS_offload_data* offload_data,
               real* offload_array) {
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->psi_coeff = offload_array;
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
void B_GS_eval_B(real B[], real r, real phi,
                 real z, B_GS_data* Bdata) {
    /* Normalize the coordinates */
    r /= R0;
    z /= R0;

    real* C = Bdata->psi_coeff;

    /* Precalculated terms used in the components */
    real r2 = r*r;
    real r3 = r2*r;
    real r4 = r3*r;
    real r5 = r4*r;
    real z2 = z*z;
    real z3 = z2*z;
    real z4 = z3*z;
    real z5 = z4*z;
    real logr = log(r);

    /* r field component */
    B[0] =   C[2]  * (-2*z) 
              + C[3]  * (-8*r2*z) 
              + C[4]  * (-18*r2*z - 24*r2*logr*z + 8*z3) 
              + C[5]  * (-24*r4*z + 32*r2*z3) 
              + C[6]  * (48*z5 - 560*r2*z3 - 480*r2*logr*z3 +360*r4*logr*z
                         + 150*r4*z) 
              + C[7]  * (1) 
              + C[8]  * (r2) 
              + C[9]  * (3*z2 - 3*r2*logr) 
              + C[10] * (3*r4 - 12*z2*r2) 
              + C[11] * (40*z4 - 45*r4 - 240*z2*r2*logr + 60*r4*logr);
    B[0] *= -psi_mult / (r * R0 * R0);

    /* phi field component */
    B[1] = B_phi0 / r;

    /* z field component */
    B[2] = (1-C[12]) * (r3/2)
              + C[12] * (r/2 + r*logr)
              + C[1]  * (2*r)
              + C[2]  * (2*r*logr + r)
              + C[3]  * (4*r3 - 8*r*z2)
              + C[4]  * (12*r3*logr + 3*r3 - 30*r*z2 - 24*r*logr*z2)
              + C[5]  * (6*r5 - 48*r3*z2 + 16*r*z4)
              + C[6]  * (-400*r*z4 -240*r*logr*z4 + 720*r3*logr*z2 + 480*r3*z2
                         -90*r5*logr - 15*r5)
              + C[8]  * (2*z*r)
              + C[9]  * (-6*z*r*logr - 3*z*r)
              + C[10] * (12*z*r3 - 8*z3*r)
              + C[11] * (-120*z*r3-160*z3*r*logr-80*z3*r+240*z*r3*logr);
    B[2] *= psi_mult / (r * R0 * R0);
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
void B_GS_eval_psi(real psi[], real r, real phi, real z,
                   B_GS_data* Bdata) {
    /* Normalize the coordinates */
    r /= R0;
    z /= R0;

    real* C = Bdata->psi_coeff;

    /* Precalculated terms used in the components */
    double r2 = r*r;
    double r3 = r2*r;
    double r4 = r3*r;
    double r5 = r4*r;
    double r6 = r5*r;
    double z2 = z*z;
    double z3 = z2*z;
    double z4 = z3*z;
    double z5 = z4*z;
    double z6 = z5*z;
    double logr = log(r);

    psi[0] = (1-C[12]) * (r4/8)
                + C[12] * (r2*logr/2)
                + C[0]  * (1)
                + C[1]  * (r2)
                + C[2]  * (r2*logr - z2)
                + C[3]  * (r4 - 4*r2*z2)
                + C[4]  * (3*r4*logr - 9*r2*z2 - 12*r2*logr*z2 + 2*z4)
                + C[5]  * (r6 - 12*r4*z2 + 8*r2*z4)
                + C[6]  * (8*z6 - 140*r2*z4 - 120*r2*logr*z4 + 180*r4*logr*z2
                           + 75*r4*z2 - 15*r6*logr)
                + C[7]  * (z)
                + C[8]  * (z*r2)
                + C[9]  * (z3 - 3*z*r2*logr)
                + C[10] * (3*z*r4 - 4*z3*r2)
                + C[11] * (8*z5 - 45*z*r4 - 80*z3*r2*logr + 60*z*r4*logr);
}

/**
 * @brief Evaluate radial coordinate rho
 *
 * This function is identical to B_2D_eval_rho.
 *
 * @see B_2D_eval_rho
 */
void B_GS_eval_rho(real rho[], real psi, B_GS_data* Bdata) {
    if(psi - Bdata->psi0 < 0) {
        rho[0] = 0;
    }
    else {
        rho[0] = sqrt((psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0));
    }
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
void B_GS_eval_B_dB(real B_dB[], real r, real phi, real z,
                    B_GS_data* Bdata) {
    real* C = Bdata->psi_coeff;

    /* Precalculated terms used in the components */
    real r2 = r*r;
    real r3 = r2*r;
    real r4 = r3*r;
    real r5 = r4*r;
    real z2 = z*z;
    real z3 = z2*z;
    real z4 = z3*z;
    real z5 = z4*z;
    real logr = log(r/R0);
    real R02 = R0*R0;
    real R03 = R02*R0;
    real R04 = R03*R0;
    real R05 = R04*R0;
    real R06 = R05*R0;

    /* r field component */
    B_dB[0] =   C[2]  * (-2*z) / R02
              + C[3]  * (-8*r2*z) / R04
              + C[4]  * (-18*r2*z - 24*r2*logr*z + 8*z3) / R04
              + C[5]  * (-24*r4*z + 32*r2*z3) / R06
              + C[6]  * (48*z5 - 560*r2*z3 - 480*r2*logr*z3 +360*r4*logr*z
                         + 150*r4*z) / R06
              + C[7]  * (1) / R0
              + C[8]  * (r2) / R03
              + C[9]  * (3*z2 - 3*r2*logr) / R03
              + C[10] * (3*r4 - 12*z2*r2) / R05
              + C[11] * (40*z4 - 45*r4 - 240*z2*r2*logr + 60*r4*logr) / R05;
    B_dB[0] *= -psi_mult / r;
    B_dB[1] =   C[3]  * (-16*r*z) / R04
                 + C[4]  * (-60*r*z - 48*r*logr*z) / R04
                 + C[5]  * (-96*r3*z  + 64*r*z3) / R06
                 + C[6]  * (-1600*r*z3 - 960*r*logr*z3 +1440*r3*logr*z
                            + 960*r3*z) / R06
                 + C[8]  * (2*r) / R03
                 + C[9]  * (-6*r*logr - 3*r) / R03
                 + C[10] * (12*r3 - 24*z2*r) / R05
                 + C[11] * (-120*r3 - 480*z2*r*logr -240*z2*r +240*r3*logr)/R05;
    B_dB[1] = -B_dB[0] / r - B_dB[1] * psi_mult / r;
    B_dB[2] = 0;
    B_dB[3] =   C[2]  * (-2) / R02
                 + C[3]  * (-8*r2) / R04
                 + C[4]  * (-18*r2 - 24*r2*logr + 24*z2) / R04
                 + C[5]  * (-24*r4 + 96*r2*z2) / R06
                 + C[6]  * (240*z4 - 1680*r2*z2 - 1440*r2*logr*z2 + 360*r4*logr
                            + 150*r4) / R06
                 + C[9]  * (6*z) / R03
                 + C[10] * (-24*z*r2) / R05
                 + C[11] * (160*z3 - 480*z*r2*logr) / R05;
    B_dB[3] *= -psi_mult / r;

    /* phi field component */
    B_dB[4] = B_phi0 * R0 / r;
    B_dB[5] = -B_phi0 * R0 / r2;
    B_dB[6] = 0;
    B_dB[7] = 0;

    /* z field component */
    B_dB[8] = (1-C[12]) * (r3/2) / R04
              + C[12] * (r/2 + r*logr) / R02
              + C[1]  * (2*r) / R02
              + C[2]  * (2*r*logr + r) / R02
              + C[3]  * (4*r3 - 8*r*z2) / R04
              + C[4]  * (12*r3*logr + 3*r3 - 30*r*z2 - 24*r*logr*z2) / R04
              + C[5]  * (6*r5 - 48*r3*z2 + 16*r*z4) / R06
              + C[6]  * (-400*r*z4 -240*r*logr*z4 + 720*r3*logr*z2 + 480*r3*z2
                         -90*r5*logr - 15*r5) / R06
              + C[8]  * (2*z*r) / R03
              + C[9]  * (-6*z*r*logr - 3*z*r) / R03
              + C[10] * (12*z*r3 - 8*z3*r) / R05
              + C[11] * (-120*z*r3-160*z3*r*logr-80*z3*r+240*z*r3*logr) / R05;
    B_dB[8] *= psi_mult / r;
    B_dB[9] = (1-C[12]) * (3*r2/2) / R04
                 + C[12] * (1.5 + logr) / R02
                 + C[1]  * (2) / R02
                 + C[2]  * (2*logr + 3) / R02
                 + C[3]  * (12*r2 - 8*z2) / R04
                 + C[4]  * (36*r2*logr + 21*r2 - 54*z2 - 24*logr*z2) / R04
                 + C[5]  * (30*r4 - 144*r2*z2 + 16*z4) / R06
                 + C[6]  * (-640*z4 - 240*logr*z4 + 2160*r2*logr*z2 + 2160*r2*z2
                            -450*r4*logr -165*r4) / R06
                 + C[8]  * (2*z) / R03
                 + C[9]  * (-6*z*logr - 9*z) / R03
                 + C[10] * (36*z*r2 - 8*z3) / R05
                 + C[11] * (-120*z*r2 - 160*z3*logr -240*z3 +720*z*r2*logr)/R05;
    B_dB[9] = B_dB[9] * psi_mult / r - B_dB[8] / r;
    B_dB[10] = 0;
    B_dB[11] = -B_dB[1] - B_dB[0] / r;
}

real B_GS_get_axis_r(B_GS_data* Bdata) {
    return 6.52306293;
}

real B_GS_get_axis_z(B_GS_data* Bdata) {
    return 0.17144154;
}
