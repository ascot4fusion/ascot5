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
 * @brief Evaluate velocity norm from Lorentz factor
 *
 * \f$ v = \sqrt{1 - \frac{1}{\gamma^2}}c\f$
 *
 * where
 *
 * - \f$v\f$ is velocity norm [m/s]
*/
#define physlib_vnorm_gamma(gamma) (                                   \
        sqrt( 1.0 - 1.0 / ( gamma * gamma ) ) * CONST_C )

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
* @brief Evaluate kinetic energy [J] from Lorentz factor
*
* \f$E_\mathrm{kin}=(\gamma - 1) * m c^2\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$\gamma\f$ is the Lorentz factor
*/
#define physlib_Ekin_gamma(m, gamma) ( ( gamma - 1.0 ) * m * CONST_C2 )

/**
* @brief Evaluate Lorentz factor from kinetic energy [J]
*
* \f$\gamma = \frac{E_\mathrm{kin}}{m c^2} + 1\f$
*
* where
*
* - \f$m\f$ is mass [kg]
* - \f$\gamma\f$ is the Lorentz factor
*/
#define physlib_gamma_Ekin(m, ekin) ( ekin / ( m * CONST_C2 ) + 1.0 )

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
#define physlib_pnorm_vnorm(m, v) ( m * v * physlib_gamma_vnorm(v) )

/**
* @brief Evaluate guiding center parallel momentum [kg m/s] from momentum norm
* and pitch
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
* @brief Evaluate guiding center magnetic moment [J/T] from momentum norm and
* pitch
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
* @brief Evaluate guiding center momentum norm [kg m/s] from parallel momentum
* and magnetic moment
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
* @brief Evaluate guiding center pitch from parallel momentum and magnetic
* moment
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
 * @brief Evaluate gyrofrequency [rad/s] from parallel momentum and magnetic
 * moment
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
 * @brief Evaluate toroidal canonical momentum for particle
 *
 * \f$\Pi_{\phi,fo} = Rp_{\phi} + q\psi\f$
 *
 * where
 *
 * \f$\Pi_{\phi,fo}\f$ is the canonical momentum conjugate to (the toroidal
 *                     angle) \f$\phi\f$ [J s/rad]
 * \f$R\f$ is the major radius [m]
 * \f$p_{\phi}\f$ is the toroidal momentum [kg m/s]
 * \fq\f$ is the charge of the particle [C]
 * \f$\psi\f$ is the poloidal magnetic flux [Tm^2]
 */
#define phys_ptoroid_fo(q, R, pphi, psi) ( \
        R * pphi  +  q * psi )

/**
 * @brief Evaluate toroidal canonical momentum for guiding center
 *
 * \f$\Pi_{\phi,gc} = p_{\parallel}R\frac{B_{\phi}}{B} + q\psi\f$
 *
 * where
 *
 * \f$\Pi_{\phi,gc}\f$ is the canonical momentum conjugate to (the toroidal
 *                     angle) \f$\phi\f$ [J s/rad]
 * \f$p_{\parallel}\f$ is the parallel momentum [kg m/s]
 * \f$R\f$ is the major radius [m]
 * \f$B_{\phi}\f$ is the toroidal component of the magnetic flux density [T]
 * \f$B\f$ is the magnitude of the magnetic flux density [T]
 * \fq\f$ is the charge of the particle [C]
 * \f$\psi\f$ is the poloidal magnetic flux [Tm^2]
 */
#define phys_ptoroid_gc(q, R, ppar, psi, B, Bphi) (             \
        ppar * R * (Bphi / B)  +  q * psi )

/**
 * @brief Evaluate perpendicular speed from velocity norm and parallel
 * component of velocity for guiding centre.
 *
 * \f$v_{\perp} = \sqrt{v^2 - v_{\parallel}^2}\f$
 *
 * where
 *
 * - \f$v_{\perp}\f$ is perpendicular speed [m/s]
 * - \f$v\f$ is speed [m/s]
 * - \f$v_{\parallel}\f$ is parallel velocity component [m/s]
 */
#define phys_vperp_gc(v,vpar) sqrt(pow(v,2) - pow(vpar,2))

/**
 * @brief Evaluate perpendicular momentum from momentum norm and parallel
 * component of momentum for guiding centre.
 *
 * \f$p_{\perp} = \sqrt{p^2 - p_{\parallel}^2}\f$
 *
 * where
 *
 * - \f$p_{\perp}\f$ is perpendicular momentum [kg m/s]
 * - \f$p\f$ is momentum [kg m/s]
 * - \f$p_{\parallel}\f$ is parallel momentum component [kg m/s]
 */
#define phys_pperp_gc(p,ppar) sqrt(pow(p,2) - pow(ppar,2))

/**
 * @brief Evaluate magnitude of parallel momentum from kinetic energy
 *
 * \f$abs(p_{\parallel}) = \sqrt{(\gamma^2 - 1)m^2c^2 - 2m\mu B}\f$
 *
 * where
 *
 * \f$abs(p_{\parallel})\f$ is the magnitude of the parallel component of
 *                           momentum
 * \f$\gammaf$ is the Lorentz factor [-]
 * \f$mf$ is the mass of the particle [kg]
 * \f$c\f$ is the speed of light in vacuum [m/s]
 * \f$\mu\f$ is the magnetic moment of the particle [J/T]
 * \f$B\f$ is the magnitude of magnetic flux density [T]
 */
#define phys_ppar_Ekin(m, ekin, mu, B) (sqrt((physlib_gamma_Ekin(m, ekin)*\
        physlib_gamma_Ekin(m, ekin) - 1.0)*m*m*CONST_C2 - 2.0*m*mu*B))

/**
 * @brief Evaluates magnitude of perpendicular velocity based on the magnetic
 * moment
 *
 * \f$v_{\perp} = \sqrt{\frac{2\mu B}{m}}\f$
 *
 * where
 *
 * \f$v_{\perp}\f$ is the perpendicular speed [m/S]
 * \f$\mu\f$ is the magnetic moment [J/T]
 * \fB\f$ is the magnitude of the magnetic flux density [T]
 * \f$\m\f$ is the mass of the particle [kg]
 */
#define phys_vperp_mu(m, mu, B) sqrt(2*mu*B/m)

/**
* @brief  Evaluates magnitude of parallel momentum from magnetic moment and
* gyrofrequency.
*
* \f$abs(p_{\parallel}) = mc\sqrt{\left(\frac{qB}{\Omega m}\right)^2 - \frac{2
*                         \mu B}{mc^2} - 1}\f$
*
* where
*
* \f$abs(p_{\parallel})\f$ is the magnitude of the parallel component of
*                           momentum
* \fm\f$ is the mass of the particle [kg]
* \f$c\f$ is the speed of light in vacuum [m/s]
* \f$q\f$ is the charge of the particle [C]
* \f$B\f$ is the magnitude of the magnetic flux density [T]
* \f$\Omega\f$ is the gyrofrequency [rad/s]
* \f$\mu\f$ is the magnetic moment [J/T]
*/
#define phys_ppar_mu(m, mu, B, q, gyrof) (m*CONST_C*sqrt((q*B/(gyrof*m)) * \
        (q*B/(gyrof*m)) - 2*mu*B/(m*CONST_C2) - 1))


// ISSUE: SEEMS TO CONTAIN A BUG!
/**
 * @brief Evaluate parallel momentum from canonical momentum conjugate to phi.
 *
 * \f$p_{\parallel} = \frac{B}{R B_{\phi}}\left(\Pi_{\phi} - q\psi\right)\f$
 *
 * where
 *
 * \f$p_{\parallel}\f$ is the parallel component of the momentum [kg m/s]
 * \fB\f$ is the magnitude of the magnetic flux density [T]
 * \fR\f$ is the major radius [m]
 * \fB_{\phi}\f$ is the toroidal component of the magnetic flux density [T]
 * \f\Pi_{\phi}\f$ is the canonical momentum conjugate to (the toroidal angle)
 *  \f$\phi\f$ [J s/rad]
 * \fq\f$ is the charge of the particle [C]
 * \f$\psi\f$ is the poloidal magnetic flux [Tm^2]
 */
#define phys_ppar_pphi(B, R, B_phi, p_phi, q, psi) (B/(R*B_phi)*(p_phi - q*psi))
#endif
