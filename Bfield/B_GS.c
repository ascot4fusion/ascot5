/**
 * @file B_GS.c
 * @brief Analytic magnetic field
 *
 * This module implements a toroidal magnetic field based on an analytical
 * solution to the Grad-Shafranov equation [1].
 *
 * In this model, the poloidal flux psi is calculated as
 *
 * \f{align*}{
   \psi(R,Z) =
   &\psi_c [(1-A) (r^4/8)                                         \\
   &+ A (r^2\log(r)/2)                                            \\
   &+ c_0                                                         \\
   &+ c_1 (r^2)                                                   \\
   &+ c_2 (r^2\log(r) - z^2)                                      \\
   &+ c_3 (r^4 - 4r^2z^2)                                         \\
   &+ c_4 (3r^4\log(r) - 9r^2z^2 - 12r^2\log(r)z^2 + 2z^4)        \\
   &+ c_5 (r^6 - 12r^4z^2 + 8r^2z^4)                              \\
   &+ c_6 (8z^6 - 140r^2z^4 - 120r^2\log(r)z^4 + 180r^4\log(r)z^2
                         + 75r^4z^2 - 15r^6\log(r))               \\
   &+ c_7  (z)                                                    \\
   &+ c_8  (zr^2)                                                 \\
   &+ c_9  (z^3 - 3zr^2\log(r))                                   \\
   &+ c_{10} (3zr^4 - 4z^3r^2)                                    \\
   &+ c_{12} (8z^5 - 45zr^4 - 80z^3r^2\log(r) + 60zr^4\log(r)) ]
   \f}
 *
 * where \f$c_i\f$ and \f$A\f$ are pre-calculated coefficients which can be
 * chosen so that realistic equilibria resembling different machines are
 * produced. The equilibrium can be non-symmetric with respect to magnetic
 * plane, and can have zero, one, or two X-points. \f$\psi_c\f$ is a scaling
 * constant, and \f$r = R/R_0\f$ \f$z = Z/R_0\f$. From \f$\psi\f$ the poloidal
 * magnetic field components can be evaluated from Grad-Shafranov relations
 *
 * \f{align*}{
   B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
   B_z &= \frac{1}{R}\frac{\partial\psi}{\partial R}
   \f}
 *
 * and toroidal field is evaluated as
 *
 * \f{equation*}{
   B_\phi = \frac{B_0R_0}{R}
   \f}
 *
 * This module also includes the possibility to have an analytical model for
 * toroidal field ripple, which is used if ripple period \f$N>0\f$. The rippled
 * toroidal field is
 *
 * \f{equation*}{
   \tilde{B}_\phi = B_\phi( 1 + \delta\cos(N\phi) )
   \f}
 *
 * where
 * \f{equation*}{
   \delta = \delta_0 \frac{r'}{a_0}^{\alpha_0}  e^{-\theta^2}
   \f}
 *
 * and \f$\theta = \arctan(Z - z_0, R- R_0)\f$ and
 * \f$r' = \sqrt{(R-R_0)^2 + (Z-Z_0)^2}\f$.
 *
 * [1] A.J. Cerfon, J.P. Freidberg. "One size fits all" analytic solutions to
 *     the Grad-Shafranov equation. Physics of Plasmas 17 (3) (2010) 032502.
 *     http://scitation.aip.org/content/aip/journal/pop/17/3/10.1063/1.3328818
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../ascot5.h"
#include "../consts.h"
#include "../print.h"
#include "../error.h"
#include "B_GS.h"

/**
 * @brief Initialize magnetic field offload data
 *
 * The offload data struct is very simple and contains only the necessary fields
 * which are all initialized when the data is read (except
 * B_GS.offload_array_length). The offload array is not required at all.
 * Therefore, this function only sets the offload_array_length to zero, assigns
 * a NULL pointer to offload_array, and prints sanity checks so that user may
 * verify that data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @return zero to indicate success
 */
int B_GS_init_offload(B_GS_offload_data* offload_data, real** offload_array) {

    offload_data->offload_array_length = 0;

    /* Evaluate psi and magnetic field on axis for checks */
    a5err err = 0;
    B_GS_data Bdata;
    B_GS_init(&Bdata, offload_data, *offload_array);
    real psival[1], Bval[3];
    err = B_GS_eval_psi(psival, offload_data->R0, 0, offload_data->z0,
                        &Bdata);
    err = B_GS_eval_B(Bval, offload_data->R0, 0, offload_data->z0,
                      &Bdata);
    if(err) {
        print_err("Error: Initialization failed.\n");
        return err;
    }

    /* Print some sanity check on data */
    print_out(VERBOSE_IO,
              "\nAnalytical tokamak magnetic field (B_GS)\n"
              "Psi at magnetic axis (%1.3f m, %1.3f m)\n"
              "%3.3f (evaluated)\n%3.3f (given)\n"
              "Magnetic field on axis:\n"
              "B_R = %3.3f B_phi = %3.3f B_z = %3.3f\n"
              "Number of toroidal field coils %d\n",
              offload_data->R0, offload_data->z0,
              psival[0], offload_data->psi0,
              Bval[0], Bval[1], Bval[2],
              offload_data->Nripple);

    return 0;
}

/**
 * @brief Free offload array
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
 * @param Bdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array offload array
 */
void B_GS_init(B_GS_data* Bdata, B_GS_offload_data* offload_data,
               real* offload_array) {
    Bdata->R0        = offload_data->R0;
    Bdata->z0        = offload_data->z0;
    Bdata->B_phi0    = offload_data->B_phi0;
    Bdata->psi0      = offload_data->psi0;
    Bdata->psi1      = offload_data->psi1;
    Bdata->psi_mult  = offload_data->psi_mult;

    Bdata->Nripple   = offload_data->Nripple;
    Bdata->a0        = offload_data->a0;
    Bdata->delta0    = offload_data->delta0;
    Bdata->alpha0    = offload_data->alpha0;

    Bdata->psi_coeff = offload_data->psi_coeff;
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
a5err B_GS_eval_psi(real* psi, real r, real phi, real z,
                    B_GS_data* Bdata) {
    /* Normalize the coordinates */
    z -= Bdata->z0;
    r /= Bdata->R0;
    z /= Bdata->R0;

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

    psi[0] = Bdata->psi_mult * (
             (1-C[12]) * (r4/8)
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
              + C[11] * (8*z5 - 45*z*r4 - 80*z3*r2*logr + 60*z*r4*logr) );

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
a5err B_GS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                         B_GS_data* Bdata) {

    // The psi value itself we can get from here
    B_GS_eval_psi(&(psi_dpsi[0]), r, phi, z, Bdata);

    /* Normalize the coordinates */
    z -= Bdata->z0;
    r /= Bdata->R0;
    z /= Bdata->R0;

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

    psi_dpsi[1] = ( Bdata->psi_mult / (r * Bdata->R0) ) *
        (1-C[12]) * (r3/2)
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

    psi_dpsi[2] = 0;

    psi_dpsi[3] =
        ( Bdata->psi_mult / (r * Bdata->R0) ) *
        ( C[2]  * (-2*z)
          + C[3]  * (-8*r2*z)
          + C[4]  * (-18*r2*z - 24*r2*logr*z + 8*z3)
          + C[5]  * (-24*r4*z + 32*r2*z3)
          + C[6]  * (48*z5 - 560*r2*z3 - 480*r2*logr*z3 +360*r4*logr*z
                     + 150*r4*z)
          + C[7]  * (1)
          + C[8]  * (r2)
          + C[9]  * (3*z2 - 3*r2*logr)
          + C[10] * (3*r4 - 12*z2*r2)
          + C[11] * (40*z4 - 45*r4 - 240*z2*r2*logr + 60*r4*logr) );

    return 0;
}

/**
 * @brief Evaluate normalized poloidal flux rho
 *
 * @param rho pointer where rho value will be stored
 * @param psi poloidal flux from which rho is evaluated
 * @param Bdata pointer to magnetic field data struct
 *
 * @return zero to indicate success
 */
a5err B_GS_eval_rho(real* rho, real psi, B_GS_data* Bdata) {

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (psi - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_GS );
    }

    rho[0] = sqrt( (psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0) );

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
a5err B_GS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                         B_GS_data* Bdata) {
    real psi_dpsi[4];

    B_GS_eval_psi_dpsi(psi_dpsi, r, phi, z, Bdata);

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (psi_dpsi[0] - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_GS );
    }

    rho_drho[0] = sqrt( (psi_dpsi[0] - Bdata->psi0) / delta );
    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[3] / (2*delta*rho_drho[0]);

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
a5err B_GS_eval_B(real B[3], real r, real phi,
                 real z, B_GS_data* Bdata) {
    /* Normalize the coordinates */
    z -= Bdata->z0;
    r /= Bdata->R0;
    z /= Bdata->R0;

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
    B[0] *= -Bdata->psi_mult / (r * Bdata->R0 * Bdata->R0);

    /* phi field component */
    B[1] = Bdata->B_phi0 / r;

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
    B[2] *= Bdata->psi_mult / (r * Bdata->R0 * Bdata->R0);

    /* Ripple */
    if(Bdata->Nripple > 0) {
        r *= Bdata->R0;
        z *= Bdata->R0;
        real radius = sqrt( ( r - Bdata->R0 ) * ( r - Bdata->R0 )
                      + ( z - Bdata->z0 ) * ( z - Bdata->z0 ));
        real theta = atan2( z - Bdata->z0, r - Bdata->R0 );
        real delta = Bdata->delta0 * exp(-0.5*theta*theta)
                     * pow( radius / Bdata->a0, Bdata->alpha0 );
        B[1] = B[1] * ( 1 + delta * cos(Bdata->Nripple * phi) );
    }

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
a5err B_GS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                     B_GS_data* Bdata) {

    real* C = Bdata->psi_coeff;

    z -= Bdata->z0;
    real R0 = Bdata->R0;
    real z0 = Bdata->z0;
    real B_phi0 = Bdata->B_phi0;
    real psi_mult = Bdata->psi_mult;

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
    real B0 =   C[2]  * (-2*z) / R02
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
    B0 *= -psi_mult / r;
    real B1 =   C[3]  * (-16*r*z) / R04
              + C[4]  * (-60*r*z - 48*r*logr*z) / R04
              + C[5]  * (-96*r3*z  + 64*r*z3) / R06
              + C[6]  * (-1600*r*z3 - 960*r*logr*z3 +1440*r3*logr*z
                            + 960*r3*z) / R06
              + C[8]  * (2*r) / R03
              + C[9]  * (-6*r*logr - 3*r) / R03
              + C[10] * (12*r3 - 24*z2*r) / R05
              + C[11] * (-120*r3 - 480*z2*r*logr -240*z2*r +240*r3*logr)/R05;
    B1 = -B0 / r - B1 * psi_mult / r;
    real B2 = 0;
    real B3 =   C[2]  * (-2) / R02
              + C[3]  * (-8*r2) / R04
              + C[4]  * (-18*r2 - 24*r2*logr + 24*z2) / R04
              + C[5]  * (-24*r4 + 96*r2*z2) / R06
              + C[6]  * (240*z4 - 1680*r2*z2 - 1440*r2*logr*z2 + 360*r4*logr
                            + 150*r4) / R06
              + C[9]  * (6*z) / R03
              + C[10] * (-24*z*r2) / R05
              + C[11] * (160*z3 - 480*z*r2*logr) / R05;
    B3 *= -psi_mult / r;

    /* phi field component */
    real B4 = B_phi0 * R0 / r;
    real B5 = -B_phi0 * R0 / r2;
    real B6 = 0;
    real B7 = 0;

    /* z field component */
    real B8 = (1-C[12]) * (r3/2) / R04
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
    B8 *= psi_mult / r;
    real B9 = (1-C[12]) * (3*r2/2) / R04
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
    B9 = B9 * psi_mult / r - B8 / r;
    real B10 = 0;
    real B11 = -B1 - B0 / r;

    /* Ripple */
    if(Bdata->Nripple > 0) {
        real radius = sqrt( ( r - R0 ) * ( r - R0 )
                      + ( z - z0 ) * ( z - z0 ));
        real theta = atan2( z - z0, r - R0 );
        real delta = Bdata->delta0 * exp(-0.5*theta*theta)
                     * pow( radius / Bdata->a0, Bdata->alpha0 );

        real Bphi = B4;
        real Bpert = Bphi * delta * cos(Bdata->Nripple * phi);
        B4 += Bpert;
        B6 += - Bphi * delta * Bdata->Nripple * sin(Bdata->Nripple * phi);

        real dBpertdR = Bpert * (
                ( (r - R0) /radius) * ( Bdata->alpha0 / radius )
                + ( (z - z0) /(radius*radius) ) * theta
            );

        real dBpertdz = Bpert * (
                ( (z - z0) /radius) * ( Bdata->alpha0 / radius )
                - ( (r - R0) /(radius*radius) ) * theta
            );

        B5 += B5 * Bpert / Bphi + dBpertdR;
        B7 += dBpertdz;
    }

    B_dB[0] = B0;
    B_dB[1] = B1;
    B_dB[2] = B2;
    B_dB[3] = B3;
    B_dB[4] = B4;
    B_dB[5] = B5;
    B_dB[6] = B6;
    B_dB[7] = B7;
    B_dB[8] = B8;
    B_dB[9] = B9;
    B_dB[10] = B10;
    B_dB[11] = B11;

    return 0;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis R-coordinate [m]
 */
real B_GS_get_axis_r(B_GS_data* Bdata) {
    return Bdata->R0;
}

/**
 * @brief Return magnetic axis z-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis z-coordinate [m]
 */
real B_GS_get_axis_z(B_GS_data* Bdata) {
    return Bdata->z0;
}
