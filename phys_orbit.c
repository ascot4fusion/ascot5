/**
 * @file gctransform.c
 * @brief Module for performing guiding center transformations
 *
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "math.h"
#include "consts.h"
#include "phys_orbit.h"

/**
 * @brief Transform particle to guiding center  phase space
 *
 * The transformation is done from coordinates [r, phi, z, vr , vphi, vz] to
 * [R, Phi, Z, vpar, mu]. The transformation is done to first order
 *
 * @param mass   mass [kg]
 * @param charge charge [C]
 * @param B_dB   magnetic field and Jacobian as given by B_field_eval_BdB()
 * @param r      particle R coordinate [m]
 * @param phi    particle phi coordinate [rad]
 * @param z      particle z coordinate [m]
 * @param vr     particle velocity R component [m/s]
 * @param vphi   particle velocity phi component [m/s]
 * @param vz     particle velocity z component [m/s]
 * @param R      pointer to guiding center R coordinate [m]
 * @param Phi    pointer to guiding center Phi coordinate [rad]
 * @param Z      pointer to guiding center Z coordinate [m]
 * @param vpar   pointer to guiding center parallel velocity [m/s]
 * @param mu     pointer to guiding center magnetic moment [J/T]
 * @param theta  pointer to guiding center gyroangle [rad]
 */
void gctransform_prt2gc(real mass, real charge, real* B_dB,
                        real r, real phi, real z, real vr, real vphi, real vz,
                        real* R, real* Phi, real* Z, real* vpar, real* mu, real* theta);
    /* Temporary re-used variables */
    real TEMP_S1;
    real TEMP_V1[3];

    /* Helper variables */
    real p_vec [3] = {p_r, p_phi, p_z};
    real p_unit[3];
    math_unit(p_vec, p_unit);
    real momentum = math_norm(p_vec);
    real gamma = phys_gammaprtp(mass, momentum);

    real B_vec[3]  = {B_dB[0], B_dB[4], B_dB[8]};
    real B_norm = math_norm(B_vec);
    real B_unit[3];
    math_unit(B_vec, B_unit);

    real pitch = math_dot(p_unit,B_unit);

    /* Magnetic field jacobian, gradient and curl in cylindrical coordinates */
    real jacB[9] = {B_dB[1], B_dB[2]/r, B_dB[3],
                    B_dB[5], B_dB[6]/r, B_dB[7],
                    B_dB[9], B_dB[10]/r, B_dB[11]};

    real gradB[3];
    gradB[0] = B_unit[0]*jacB[0] + B_unit[0]*jacB[3] + B_unit[0]*jacB[6];
    gradB[1] = B_unit[1]*jacB[1] + B_unit[1]*jacB[4] + B_unit[1]*jacB[7];
    gradB[2] = B_unit[2]*jacB[2] + B_unit[2]*jacB[5] + B_unit[2]*jacB[8];
    math_matmul(jacB,B_unit,3,3,1,gradB);

    real curlB[3] = {jacB[7]-jacB[5], jacB[2]-jacB[6], jacB[3]+B_vec[1]/r-jacB[1]};

    real tau_B = math_dot(B_unit, curlB)/B_norm;
    real nablabhat[9];
    nablabhat[0]=(jacB[0] - gradB[0]*B_unit[0])/B_norm;
    nablabhat[1]=(jacB[1] - gradB[0]*B_unit[1])/B_norm;
    nablabhat[2]=(jacB[2] - gradB[0]*B_unit[2])/B_norm;
    nablabhat[3]=(jacB[3] - gradB[1]*B_unit[0])/B_norm;
    nablabhat[4]=(jacB[4] - gradB[1]*B_unit[1])/B_norm;
    nablabhat[5]=(jacB[5] - gradB[1]*B_unit[2])/B_norm;
    nablabhat[6]=(jacB[6] - gradB[2]*B_unit[0])/B_norm;
    nablabhat[7]=(jacB[7] - gradB[2]*B_unit[1])/B_norm;
    nablabhat[8]=(jacB[8] - gradB[2]*B_unit[2])/B_norm;

    real kappa[3];
    math_matmul(nablabhat,B_unit,3,3,1,kappa);

    /* Zeroth order momentum terms */
    real p_para0 = pitch*momentum;
    real mu_0 = ( 1 - pow(pitch,2) )*pow(momentum,2)/(2*mass*B_norm);

    /* Make the spatial transformation */
    real rho[3];
    math_cross(B_unit,p_unit,rho);
    math_prod(rho,momentum/(charge*B_norm));

    real rho_unit[3];
    math_unit(rho,rho_unit);
    real rho_norm = math_norm(rho);

    gcpos[0] = r - rho[0];
    gcpos[1] = phi - rho[1];
    gcpos[2] = z - rho[2];

    /* First order momentum terms */
    real perphat[3];
    math_cross(rho_unit, B_unit, perphat); /* Does this hold for ions and elecs? */

    real a1ddotgradb = -0.5*(2*(rho_unit[0]*perphat[0]*nablabhat[0]+
                                rho_unit[1]*perphat[1]*nablabhat[4]+
                                rho_unit[2]*perphat[2]*nablabhat[8])
                             +(rho_unit[0]*perphat[1]+rho_unit[1]*perphat[0])*
                              (nablabhat[1]+nablabhat[3])
                             +(rho_unit[0]*perphat[2]+rho_unit[2]*perphat[0])*
                              (nablabhat[2]+nablabhat[6])
                             +(rho_unit[1]*perphat[2]+rho_unit[2]*perphat[1])*
                              (nablabhat[5]+nablabhat[7]));

    real p_para1 = -p_para0*math_dot(rho,kappa)+((mass*mu_0)/charge)*(tau_B+a1ddotgradb);
    TEMP_S1 = pow(p_para0,2)/(mass*B_norm);
    TEMP_V1[0] = TEMP_S1 * kappa[0] + mu_0*gradB[0]/B_norm;
    TEMP_V1[1] = TEMP_S1 * kappa[1] + mu_0*gradB[1]/B_norm;
    TEMP_V1[2] = TEMP_S1 * kappa[2] + mu_0*gradB[2]/B_norm;
    real mu_1 = math_dot(rho,TEMP_V1)-((mu_0*p_para0)/(charge*B_norm))*(tau_B+a1ddotgradb);

    /* Make the momentum transformation */
    gcpos[3] = p_para0 + p_para1;
    gcpos[4] = mu_0 + mu_1;

    /* Calculate gyroangle */
    real a1[3];
    real z_unit[3];
    z_unit[0] = 0.0;
    z_unit[1] = 0.0;
    z_unit[2] = 1.0;
    math_cross(B_unit,z_unit,a1);
    math_unit(a1,a1);

    real a2[3];
    math_cross(a1,B_unit,a2);
    math_unit(a2,a2);

    gcpos[5] = atan2(math_dot(rho_unit,a2),math_dot(rho_unit,a1));
    gcpos[5] = fmod(gcpos[5] + CONST_2PI, CONST_2PI); /* theta is in interval 0 to 2PI */
}

/**
 * @brief Transform guiding center to particle phase space
 *
 * The transformation is done from coordinates [R, Phi, Z, vpar, mu] to
 * [r, phi, z, vr , vphi, vz]. The transformation is done to first order.
 *
 * @param mass   mass [kg]
 * @param charge charge [C]
 * @param B_dB   magnetic field and Jacobian as given by B_field_eval_BdB()
 * @param R      guiding center R coordinate [m]
 * @param Phi    guiding center Phi coordinate [rad]
 * @param Z      guiding center Z coordinate [m]
 * @param vpar   guiding center parallel velocity [m/s]
 * @param mu     guiding center magnetic moment [J/T]
 * @param theta  guiding center gyroangle [rad]
 * @param r      pointer to particle R coordinate [m]
 * @param phi    pointer to particle phi coordinate [rad]
 * @param z      pointer to particle z coordinate [m]
 * @param vr     pointer to particle velocity R component [m/s]
 * @param vphi   pointer to particle velocity phi component [m/s]
 * @param vz     pointer to particle velocity z component [m/s]
 */
void gctransform_gc2prt(real mass, real charge, real* B_dB,
                        real R, real Phi, real Z, real vpar, real mu, real theta,
                        real* r, real* phi, real* z, real* vr, real* vphi, real* vz);
    /* Temporary re-used variables */
    real TEMP_S1;
    real TEMP_V1[3];

    /* Helper variables */
    real B_vec[3]  = {B_dB[0], B_dB[4], B_dB[8]};
    real B_norm = math_norm(B_vec);
    real B_unit[3];
    math_unit(B_vec, B_unit);

    real gamma = phys_gammagcv(mass,v_para,mu);

    /* Magnetic field jacobian, gradient and curl in cylindrical coordinates */
    real jacB[9] = {B_dB[1], B_dB[2]/R, B_dB[3],
                    B_dB[5], B_dB[6]/R, B_dB[7],
                    B_dB[9], B_dB[10]/R, B_dB[11]};

    real gradB[3];
    gradB[0] = B_unit[0]*jacB[0] + B_unit[0]*jacB[3] + B_unit[0]*jacB[6];
    gradB[1] = B_unit[1]*jacB[1] + B_unit[1]*jacB[4] + B_unit[1]*jacB[7];
    gradB[2] = B_unit[2]*jacB[2] + B_unit[2]*jacB[5] + B_unit[2]*jacB[8];
    math_matmul(jacB,B_unit,3,3,1,gradB);

    real curlB[3] = {jacB[7]-jacB[5], jacB[2]-jacB[6], jacB[3]+B_vec[1]/R-jacB[1]};

    real tau_B = math_dot(B_unit, curlB)/B_norm;
    real nablabhat[9];
    nablabhat[0]=(jacB[0] - gradB[0]*B_unit[0])/B_norm;
    nablabhat[1]=(jacB[1] - gradB[0]*B_unit[1])/B_norm;
    nablabhat[2]=(jacB[2] - gradB[0]*B_unit[2])/B_norm;
    nablabhat[3]=(jacB[3] - gradB[1]*B_unit[0])/B_norm;
    nablabhat[4]=(jacB[4] - gradB[1]*B_unit[1])/B_norm;
    nablabhat[5]=(jacB[5] - gradB[1]*B_unit[2])/B_norm;
    nablabhat[6]=(jacB[6] - gradB[2]*B_unit[0])/B_norm;
    nablabhat[7]=(jacB[7] - gradB[2]*B_unit[1])/B_norm;
    nablabhat[8]=(jacB[8] - gradB[2]*B_unit[2])/B_norm;

    real kappa[3];
    math_matmul(nablabhat,B_unit,3,3,1,kappa);

    /* Zeroth order momentum terms */
    real p_para0 = gamma*mass*v_para;
    real mu_0 = mu;

    /* Calculate gyroradius and its unit vector. Choose basis vectors a1,
       which is perpendicular to B and z (since we always have Bphi),
       and a2 which is perpendicular to a1 and B */
    real a1[3];
    real z_unit[3];
    z_unit[0] = 0.0;
    z_unit[1] = 0.0;
    z_unit[2] = 1.0;
    math_cross(B_unit,z_unit,a1);
    math_unit(a1,a1);

    real a2[3];
    math_cross(a1,B_unit,a2);
    math_unit(a2,a2);

    real rho_unit[3];
    rho_unit[0] = cos(theta)*a1[0]+sin(theta)*a2[0]; /* theta + for ions, - for elecs */
    rho_unit[1] = cos(theta)*a1[1]+sin(theta)*a2[1];
    rho_unit[2] = cos(theta)*a1[2]+sin(theta)*a2[2];

    real rho[3];
    rho[0] = rho_unit[0];
    rho[1] = rho_unit[1];
    rho[2] = rho_unit[2];
    math_prod(rho,sqrt(2*pow(gamma,2)*mass*mu_0/(pow(charge,2)*B_norm))); /* Is the Lorentz
                                                               factor used correctly? */

    /* Make the spatial transformation */
    prtpos[0] = R+rho[0];
    prtpos[1] = Phi+rho[1];
    prtpos[2] = Z+rho[2];

    /* First order momentum terms */
    real perphat[3];
    math_cross(rho_unit,B_unit,perphat); /* Does this hold for ions and elecs? */

    real a1ddotgradb = -0.5*(2*(rho_unit[0]*perphat[0]*nablabhat[0]+
                                rho_unit[1]*perphat[1]*nablabhat[4]+
                                rho_unit[2]*perphat[2]*nablabhat[8])
                             +(rho_unit[0]*perphat[1]+rho_unit[1]*perphat[0])*
                              (nablabhat[1]+nablabhat[3])
                             +(rho_unit[0]*perphat[2]+rho_unit[2]*perphat[0])*
                              (nablabhat[2]+nablabhat[6])
                             +(rho_unit[1]*perphat[2]+rho_unit[2]*perphat[1])*
                              (nablabhat[5]+nablabhat[7]));

    real p_para1 = -p_para0*math_dot(rho,kappa)+((mass*mu_0)/charge)*(tau_B+a1ddotgradb);
    TEMP_S1 = pow(p_para0,2)/(mass*B_norm);
    TEMP_V1[0] = TEMP_S1 * kappa[0] + mu_0*gradB[0]/B_norm;
    TEMP_V1[1] = TEMP_S1 * kappa[1] + mu_0*gradB[1]/B_norm;
    TEMP_V1[2] = TEMP_S1 * kappa[2] + mu_0*gradB[2]/B_norm;
    real mu_1 = math_dot(rho,TEMP_V1)-((mu_0*p_para0)/(charge*B_norm))*(tau_B+a1ddotgradb);

    /* Calculate the parallel momentum vector */
    real p_para[3];
    p_para[0] = B_unit[0];
    p_para[1] = B_unit[1];
    p_para[2] = B_unit[2];
    math_prod(p_para,p_para0-p_para1);

    /* Calculate perpendicular momentum vector. */
    real p_perp[3];
    p_perp[0] = perphat[0];
    p_perp[1] = perphat[1];
    p_perp[2] = perphat[2];
    math_prod(p_perp,sqrt(2*mass*(mu_0-mu_1)*B_norm));

    /* Calculate the momentum vector */
    prtpos[3] = p_para[0]+p_perp[0];
    prtpos[4] = p_para[1]+p_perp[1];
    prtpos[5] = p_para[2]+p_perp[2];
}