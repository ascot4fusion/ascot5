/**
 * @file physlib.h
 * @brief Methods to evaluate elementary physical quantities
 */
#ifndef PHYSLIB_H
#define PHYSLIB_H

#include <math.h>
#include "ascot5.h"
#include "const.h"

/**
* @brief Evaluate Lorentz factor from velocity norm
*
* \$f\gamma = \sqrt{\frac{1}{1-v^2/c^2}}\$f
*
* where
*
* - \$fv\$f is velocity norm [m/s]
*/
#define physlib_gamma_vnorm(v) ( \
        sqrt( 1.0 / ( (1.0 - v / CONST_C) * (1.0 + v / CONST_C) ) ) )

/**
* @brief Evaluate Lorentz factor from momentum norm
*
* \$f\gamma = \sqrt{1 + \left(\frac{p}{mc}\right)^2}\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$fp\$f is momentum norm [kg m/s]
*/
#define physlib_gamma_pnorm(mass, p) ( \
        sqrt(1.0 + ( p * p ) / ( m * m * CONST_C2 ) ) )

/**
* @brief Evaluate Lorentz factor from parallel velocity
*
* \$f\gamma = \sqrt{\frac{1 + (2\mu B/mc^2)}{1 - v_\parallel^2/c^2}}\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$f\mu\$f is magnetic moment [J/T]
* - \$fv_\parallel\$f is parallel velocity [m/s]
* - \$fB\$f is magnetic field norm [T]
*/
#define physlib_gamma_vpar(m, mu, vpar, B) ( \
        sqrt( ( 1.0 + 2 * mu * B / ( m * CONST_C2 ) ) / \
        ( (1.0 - vpar / CONST_C) * (1 + vpar / CONST_C) ) ) )

/**
* @brief Evaluate Lorentz factor from parallel momentum
*
* \$f\gamma = \sqrt{1 + 2\mu B/mc^2 + (p_\parallel/mc)^2}\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$f\mu\$f is magnetic moment [J/T]
* - \$fp_\parallel\$f is parallel momentum [kg m/s]
* - \$fB\$f is magnetic field norm [T]
*/
#define physlib_gamma_ppar(m, mu, ppar, B) ( \
        sqrt( 1.0 + 2\mu B / ( m * CONST_C2 ) +      \
        ppar * ppar / ( m * m * CONST_C2 ) ) )

/**
* @brief Evaluate kinetic energy [J] from velocity norm
*
* \$fE_\mathrm{kin}=(\gamma(v) - 1) * m c^2\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$fv\$f is velocity norm [kg m/s]
*/
#define physlib_Ekin_vnorm(m, v) ( \
        ( physlib_gamma_vpar(v) - 1.0 ) * m * CONST_C2 )

/**
* @brief Evaluate kinetic energy [J] from parallel velocity
*
* \$fE_\mathrm{kin}=(\gamma(m, \mu, v_\parallel, B) - 1) * m c^2\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$fv\$f is velocity norm [kg m/s]
*/
#define physlib_Ekin_vpar(m, mu, vpar, B) ( \
        ( physlib_relfactor_vpar(m, mu, vpar, B) - 1.0 ) * m * CONST_C2 )

/**
* @brief Evaluate guiding center parallel velocity [m/s] from velocity norm and pitch
*
* \$fv_\parallel = \xi v_\parallel\$f
*
* where
*
* - \$fv\$f is velocity norm [m/s]
* - \$f\xi\$f is pitch
*/
#define physlib_gc_vpar(v, xi) ( v*\xi ) 

/**
* @brief Evaluate guiding center magnetic moment [J/T] from velocity norm and pitch
*
* \$fv_\parallel = \xi v_\parallel\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$fv\$f is velocity norm [m/s]
* - \$f\xi\$f is pitch
* - \$fB\$f is magnetic field [T]
*/
#define physlib_gc_mu(m, v, xi, B) ( m * pow(physlib_gamma_vnorm(v) * v, 2) * \
        ( 1.0 - xi * xi ) / ( 2 * B ) )

/**
* @brief Evaluate guiding center velocity norm [m/s] from parallel velocity and magnetic moment
*
* \$fv = \sqrt{v_\parallel^2 + 2\mu B / m\gamma(m,\mu,v_\parallel,B)}\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$f\mu\$f is magnetic moment [J/T]
* - \$fv_\parallel\$f is parallel velocity
* - \$fB\$f is magnetic field [T]
*/
#define physlib_gc_v(m, mu, vpar, B) ( sqrt( vpar * vpar + \
        2 * mu * B / ( m * physlib_gamma_vpar(m, mu, vpar, B) ) ) )

/**
* @brief Evaluate guiding center pitch from parallel velocity and magnetic moment
*
* \$f\xi = v_\parallel\sqrt{v_\parallel^2 + 2\mu B 
*          / m\gamma(m,\mu,v_\parallel,B)}\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$f\mu\$f is magnetic moment [J/T]
* - \$fv_\parallel\$f is parallel velocity
* - \$fB\$f is magnetic field [T]
*/
#define physlib_gc_xi(m, mu, vpar, B) ( vpar / sqrt( vpar * vpar + \
        2 * mu * B / ( m * physlib_gamma_vpar(m, mu, vpar, B) ) ) )

/**
 * @brief Gyro length (m) in particle phase space (mass, charge, momentum[3], magnetic field [3])
 *        r_g = p_perp / |q| * |B|
 */
#define physlib_gyrolength_vnorm(m, q, p, B) ( math_dot(p,B) / (fabs(q) * ( B[0]*B[0] + B[1]*B[1] + B[2]*B[2] ) ) );

/**
 * @brief Gyro length (m) in guiding center phase space
 */
#define phys_gyrolength_() ();

/**
 * @brief Gyro frequency (rad/s) in particle phase space (mass, charge, momentum, magnetic field)
 *        w_g = |q| * B / gamma * m
 */
#define phys_gyrofreqprt(m, q, p, B) ( fabs(q) * math_norm(B) / ( m * phys_gammaprtp(m,p) ) );

/**
 * @brief Gyro frequency (rad/s) in guiding center phase space
 */
#define phys_gyrofreqgc() ();

#endif
