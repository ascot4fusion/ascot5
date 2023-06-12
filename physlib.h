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
* \f$ \gamma = \sqrt{\frac{1}{1-v^2/c^2}}\f$
*
* where
*
* - \f$v\f$ is velocity norm [m/s]
*/
#define physlib_gamma_vnorm(v) (                                        \
        sqrt( 1.0 / ( (1.0 - v / CONST_C) * (1.0 + v / CONST_C) ) ) )

/**
* @brief Evaluate Lorentz factor from momentum norm
*
* \f$\gamma = \sqrt{1 + \left(\frac{p}{mc}\right)^2}\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$p\f$ is momentum norm [kg m/s]
*/
#define physlib_gamma_pnorm(m, p) (                  \
        sqrt(1.0 + ( p * p ) / ( m * m * CONST_C2 ) ) )

/**
* @brief Evaluate Lorentz factor from parallel velocity
*
* \f$\gamma = \sqrt{\frac{1 + (2\mu B/mc^2)}{1 - v_\parallel^2/c^2}}\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$\mu\f$ is magnetic moment [J/T]
* - \f$v_\parallel\f$ is parallel velocity [m/s]
* - \f$B\f$ is magnetic field norm [T]
*/
#define physlib_gamma_vpar(m, mu, vpar, B) (                            \
        sqrt( ( 1.0 + (2.0 * mu * B) / ( m * CONST_C2 ) ) /             \
              ( (1.0 - vpar / CONST_C) * (1.0 + vpar / CONST_C) ) ) )

/**
* @brief Evaluate Lorentz factor from parallel momentum
*
* \f$\gamma = \sqrt{1 + 2\mu B/mc^2 + (p_\parallel/mc)^2}\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$\mu\f$ is magnetic moment [J/T]
* - \f$p_\parallel\f$ is parallel momentum [kg m/s]
* - \f$B\f$ is magnetic field norm [T]
*/
#define physlib_gamma_ppar(m, mu, ppar, B) (             \
        sqrt( 1.0 + 2 * mu * B / ( m * CONST_C2 ) +      \
              ppar * ppar / ( m * m * CONST_C2 ) ) )

/**
* @brief Evaluate kinetic energy [J] from momentum norm
*
* \f$E_\mathrm{kin}=(\gamma(p) - 1) * m c^2\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$p\f$ is momentum norm [kg m/s]
*/
#define physlib_Ekin_pnorm(m, p) (                              \
        ( physlib_gamma_pnorm(m, p) - 1.0 ) * m * CONST_C2 )

/**
* @brief Evaluate kinetic energy [J] from parallel momentum
*
* \f$E_\mathrm{kin}=(\gamma(m, \mu, p_\parallel, B) - 1) * m c^2\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$p_\parallel\f$ is parallel momentum [kg m/s]
*/
#define physlib_Ekin_ppar(m, mu, ppar, B) (                             \
        ( physlib_gamma_ppar(m, mu, ppar, B) - 1.0 ) * m * CONST_C2 )

/**
* @brief Evaluate velocity norm [m/s] from momentum norm
*
* \f$v = p/\gamma(p)m\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$p\f$ is momentum norm [kg m/s]
*/
#define physlib_vnorm_pnorm(m, p) (                  \
        p / sqrt(m * m + ( p * p ) / CONST_C2 ) )

/**
* @brief Evaluate momentum norm [kg m/s] from velocity norm
*
* \f$p = \gamma(v)mv\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$v\f$ is velocity norm [m/s]
*/
#define physlib_pnorm_vnorm(m, v) (                  \
        m * v * physlib_gamma_vnorm(v) )

/**
* @brief Evaluate guiding center parallel momentum [kg m/s] from momentum norm and pitch
*
* \f$p_\parallel = \xi p\f$
*
* where
*
* - \f$p\f$ is momentum norm [kg m/s]
* - \f$\xi\f$ is pitch
*/
#define physlib_gc_ppar(p, xi) ( p * xi )

/**
* @brief Evaluate guiding center magnetic moment [J/T] from momentum norm and pitch
*
* \f$\mu = (1-\xi^2)p/(2mB)\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$p\f$ is momentum norm [kg m/s]
* - \f$\xi\f$ is pitch
* - \f$B\f$ is magnetic field norm [T]
*/
#define physlib_gc_mu(m, p, xi, B) (                    \
        p * p * ( 1.0 - xi * xi ) / ( 2 * B * m ) )

/**
* @brief Evaluate guiding center momentum norm [kg m/s] from parallel momentum and magnetic moment
*
* \f$p = \sqrt{\gamma(m,\mu,p_\parallel,B)^2 - 1} m c\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$\mu\f$ is magnetic moment [J/T]
* - \f$p_\parallel\f$ is parallel momentum [kg m/s]
* - \f$B\f$ is magnetic field norm [T]
*/
#define physlib_gc_p(m, mu, ppar, B) ( \
        m * CONST_C * sqrt( pow( physlib_gamma_ppar(m, mu, ppar, B), 2) - 1 ) )

/**
* @brief Evaluate guiding center pitch from parallel momentum and magnetic moment
*
* \f$\xi = p_\parallel / p(m,\mu,p_\parallel,B) \f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$\mu\f$ is magnetic moment [J/T]
* - \f$p_\parallel\f$ is parallel momentum [kg m/s]
* - \f$B\f$ is magnetic field [T]
*/
#define physlib_gc_xi(m, mu, ppar, B) (                                 \
        ppar / physlib_gc_p(m, mu, ppar, B) )

/**
 * @brief Evaluate gyroradius [m] from momentum vector
 *
 * \f$\rho_g = \frac{\mathbf{p}\cdot\mathbf{B}}{|q|B^2}\f$
 *
 * where
 *
 * - \f$q\f$ is charge [C]
 * - \f$\mathbf{v}\f$ is momentum vector [kg m/s]
 * - \f$\mathbf{B}\f$ is magnetic field vector [T]
 */
#define physlib_gyrolength_p(q, p, B) (                              \
        math_dot(p, B) / ( fabs(q) * math_dot(B, B) ) )

/**
 * @brief Evaluate gyroradius [m] from parallel momentum and magnetic moment
 *
 * \f$\rho_g = \frac{1}{|q|}\sqrt{\frac{2\gamma(m, \mu, p_\parallel, B) m \mu}
 *             {B}}\f$
 *
 * where
 *
 * - \f$m\f$ is mass [kg]
 * - \f$q\f$ is charge [C]
 * - \f$\mu\f$ is magnetic moment [J/T]
 * - \f$p_\parallel\f$ is parallel momentum [kg m/s]
 * - \f$B\f$ is magnetic field norm [T]
 */
#define phys_gyrolength_ppar(m, q, mu, ppar, B) (                   \
        sqrt( 2 * m * mu *                                          \
              physlib_gamma_ppar(m, mu, ppar, B) / B ) / fabs(q) )

/**
 * @brief Evaluate gyrofrequency [rad/s] from momentum norm
 *
 * \f$\omega_g = \frac{q B}{\gamma(p) m}\f$
 *
 * where
 *
 * - \f$m\f$ is mass [kg]
 * - \f$q\f$ is charge [C]
 * - \f$p\f$ is momentum norm [kg m/s]
 * - \f$B\f$ is magnetic field norm [T]
 */
#define phys_gyrofreq_pnorm(m, q, p, B) (               \
        fabs(q) * B / ( m * physlib_gamma_pnorm(m, p) ) )

/**
 * @brief Evaluate gyrofrequency [rad/s] from parallel momentum and magnetic moment
 *
 * \f$\omega_g = \frac{q B}{\gamma(m, \mu, p_\parallel, B) m}\f$
 *
 * where
 *
 * - \f$m\f$ is mass [kg]
 * - \f$q\f$ is charge [C]
 * - \f$\mu\f$ is magnetic moment [J/T]
 * - \f$p_\parallel\f$ is parallel momentum [kg m/s]
 * - \f$B\f$ is magnetic field norm [T]
 */
#define phys_gyrofreq_ppar(m, q, mu, ppar, B) (                 \
        fabs(q) * B / ( m * physlib_gamma_ppar(m, mu, ppar, B) ) )


/**
 * P_Tor for guiding center 

 */
#define phys_ptoroid_gc()



/**
 * P_Tor for particle 

 */
#define phys_ptoroid_prt()


#endif
