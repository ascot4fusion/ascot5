/**
 * @file physlib.h
 * @brief Methods to evaluate elementary physical quantities
 */
#ifndef PHYSLIB_H
#define PHYSLIB_H

#include <math.h>
#include "math.h"
#include "consts.h"

/**
* @brief Evaluate Lorentz factor from velocity norm
*
* \$f\gamma = \sqrt{\frac{1}{1-v^2/c^2}}\$f
*
* where
*
* - \$fv\$f is velocity norm [m/s]
*/
#define physlib_gamma_vnorm(v) (                                        \
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
#define physlib_gamma_pnorm(m, p) (                  \
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
#define physlib_gamma_vpar(m, mu, vpar, B) (                            \
        sqrt( ( 1.0 + (2.0 * mu * B) / ( m * CONST_C2 ) ) /             \
              ( (1.0 - vpar / CONST_C) * (1.0 + vpar / CONST_C) ) ) )

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
#define physlib_gamma_ppar(m, mu, ppar, B) (             \
        sqrt( 1.0 + 2 * mu * B / ( m * CONST_C2 ) +      \
              ppar * ppar / ( m * m * CONST_C2 ) ) )

/**
* @brief Evaluate kinetic energy [J] from velocity norm
*
* \$fE_\mathrm{kin}=(\gamma(v) - 1) * m c^2\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$fv\$f is velocity norm [m/s]
*/
#define physlib_Ekin_vnorm(m, v) (                              \
        ( physlib_gamma_vpar(v) - 1.0 ) * m * CONST_C2 )

/**
* @brief Evaluate kinetic energy [J] from parallel velocity
*
* \$fE_\mathrm{kin}=(\gamma(m, \mu, v_\parallel, B) - 1) * m c^2\$f
*
* where
*
* - \$fm\$f is mass [kg]
* - \$fv\$f is velocity norm [m/s]
*/
#define physlib_Ekin_vpar(m, mu, vpar, B) (                             \
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
#define physlib_gc_vpar(v, xi) ( v * xi )

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
* - \$fB\$f is magnetic field norm [T]
*/
#define physlib_gc_mu(m, v, xi, B) (                    \
        m * pow(physlib_gamma_vnorm(v) * v, 2) *        \
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
* - \$fv_\parallel\$f is parallel velocity [m/s]
* - \$fB\$f is magnetic field norm [T]
*/
#define physlib_gc_v(m, mu, vpar, B) ( sqrt( vpar * vpar +              \
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
* - \$fv_\parallel\$f is parallel velocity [m/s]
* - \$fB\$f is magnetic field [T]
*/
#define physlib_gc_xi(m, mu, vpar, B) (                                 \
        vpar / sqrt( vpar * vpar +                                      \
                     2 * mu * B / ( m * physlib_gamma_vpar(m, mu, vpar, B) ) ) )

/**
 * @brief Evaluate gyroradius [m] from velocity norm
 *
 * \$f\rho_g = \frac{m\gamma(v)\mathbf{v}\cdot\mathbf{B}}{|q|B^2}\$f
 *
 * where
 *
 * - \$fm\$f is mass [kg]
 * - \$fq\$f is charge [C]
 * - \$f\$f is velocity vector [m/s]
 * - \$f\mathbf{B}\$f is magnetic field vector [T]
 */
#define physlib_gyrolength_v(m, q, v, B) (                              \
        m * physlib_gamma_vnorm( math_norm(v) ) * math_dot(v, B) /       \
        ( fabs(q) * math_dot(B, B) ) )

/**
 * @brief Evaluate gyroradius [m] from parallel velocity and magnetic moment
 *
 * \$f\rho_g = \frac{1}{|q|}\sqrt{\frac{2\gamma(m, \mu, v_\parallel, B) m \mu}
 *             {B}}\$f
 *
 * where
 *
 * - \$fm\$f is mass [kg]
 * - \$fq\$f is charge [C]
 * - \$f\mu\$f is magnetic moment [J/T]
 * - \$fv_\parallel\$f is parallel velocity [m/s]
 * - \$fB\$f is magnetic field norm [T]
 */
#define phys_gyrolength_vpar(m, q, mu, vpar, B) (                   \
        sqrt( 2 * m * mu *                                          \
              physlib_gamma_vpar(m, mu, vpar, B) / B ) / fabs(q) )

/**
 * @brief Evaluate gyrofrequency [rad/s] from velocity norm
 *
 * \$f\omega_g = \frac{q B}{\gamma(v) m}$f
 *
 * where
 *
 * - \$fm\$f is mass [kg]
 * - \$fq\$f is charge [C]
 * - \$fv\$f is velocity norm [m/s]
 * - \$fB\$f is magnetic field norm [T]
 */
#define phys_gyrofreq_vnorm(m, q, v, B) (               \
        fabs(q) * B / ( m * physlib_gamma_vnorm(v) ) )

/**
 * @brief Evaluate gyrofrequency [rad/s] from parallel velocity and magnetic moment
 *
 * \$f\omega_g = \frac{q B}{\gamma(m, \mu, v_\parallel, B) m}$f
 *
 * where
 *
 * - \$fm\$f is mass [kg]
 * - \$fq\$f is charge [C]
 * - \$f\mu\$f is magnetic moment [J/T]
 * - \$fv_\parallel\$f is parallel velocity [m/s]
 * - \$fB\$f is magnetic field norm [T]
 */
#define phys_gyrofreq_vpar(m, q, mu, vpar, B) (                 \
        fabs(q) * B / ( m * physlib_gamma_vpar(m, mu, vpar, B) ) )

#endif
