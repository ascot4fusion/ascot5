/**
 * @file phys_orbit.h
 * @brief Header file for phys_orbit.c
 */
#ifndef PHYS_ORBIT_H
#define PHYS_ORBIT_H

#include "ascot5.h"
#include "consts.h"

#pragma omp declare target
/**
 * @brief Lorentz factor in particle phase space using momentum
 * gamma(mass,p) = sqrt( 1 + (p/mc)^2 )
 */
#define phys_gammaprtp(m,p) ( sqrt(1 + pow(p/(m*CONST_C), 2) ) );

/**
 * @brief Lorentz factor in particle phase space using velocity
 * gamma(v) = 1/sqrt( 1 - v^2/c^2)
 */
#define phys_gammaprtv(v) ( 1/sqrt(1 - v*v/CONST_C2 ) );

/**
 * @brief Lorentz factor in guiding center phase space using momentum
 * gamma(mass,ppara,mu) = sqrt( 1 + 2 * mu / (mass * c^2) +(ppara/m * c)^2)
 */
#define phys_gammagcp(m,ppara,mu) ( sqrt(1 + 2 * mu / ( m * CONST_C2 ) + pow(ppara/( m*CONST_C ), 2) ) );

/**
 * @brief Lorentz factor in guiding center phase space using velocity
 * gamma(mass,vpara,mu) = sqrt( ( 1 + 2 * mu / m * c^2 )/( 1-vpara^2/c^2 ))
 */
#define phys_gammagcv(m,vpara,mu) ( sqrt( ( 1 + 2 * mu / ( m * CONST_C2 ) ) / ( 1 - vpara * vpara / CONST_C2 ) ) );

/**
 * @brief Gyro length (m) in particle phase space (mass, charge, momentum[3], magnetic field [3])
 *        r_g = p_perp / |q| * |B|
 */
#define phys_gyrolengthprt(m, q, p, B) ( math_dot(p,B) / (fabs(q) * ( B[0]*B[0] + B[1]*B[1] + B[2]*B[2] ) ) );

/**
 * @brief Gyro length (m) in guiding center phase space
 */
#define phys_gyrolengthgc() ();

/**
 * @brief Gyro frequency (rad/s) in particle phase space (mass, charge, momentum, magnetic field)
 *        w_g = |q| * B / gamma * m
 */
#define phys_gyrofreqprt(m, q, p, B) ( fabs(q) * math_norm(B) / ( m * phys_gammaprtp(m,p) ) );

/**
 * @brief Gyro frequency (rad/s) in guiding center phase space
 */
#define phys_gyrofreqgc() ();

/**
 * @brief Kinetic energy to magnitude of velocity
 * v(m,E) = c*sqrt( 1 - 1 /(1 + E/(mc^2))^2 )
 */
#define phys_Ekintovelocity(m,E) ( CONST_C * sqrt( 1 - 1 / pow( 1 + E/(m * CONST_C2), 2) ) );

/**
 * @brief Magnetic moment from v_para and v_perp
 * mu(v_para,v_perp,m,E) = (gamma*v_perp*m)^2 / (2 * m * B)
 */
#define phys_mu(vpa,vpe,m,B) ( pow(phys_gammaprtv( sqrt(vpa*vpa+vpe*vpe) )*vpe*m, 2) / (2*m*B) );

#pragma omp declare simd
void phys_prttogc(real mass, real charge, real r, real phi, real z, 
		  real p_r, real p_phi, real p_z, real* B_dB, real* gcpos);
#pragma omp declare simd
void phys_gctoprt(real R, real Phi, real Z, real v_para, real mu, 
		  real* B_dB, real* prtpos);
#pragma omp declare simd
void phys_eomprt();
#pragma omp declare simd
void phys_eomgc(real* ydot, real* y, real mass, real charge, real* B_dB, real* E);
#pragma omp end declare target

#endif
