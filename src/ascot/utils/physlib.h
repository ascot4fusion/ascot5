/**
 * @file physlib.h
 * Methods to evaluate elementary physical quantities
 */
#ifndef PHYSLIB_H
#define PHYSLIB_H

#include "defines.h"
#include "parallel.h"
#include "consts.h"
#include "mathlib.h"
#include <math.h>

/**
 * @brief Evaluate Lorentz factor from velocity norm
 *
 * \f$ \gamma = \sqrt{\frac{1}{1-v^2/c^2}}\f$
 *
 * where
 *
 * - \f$v\f$ is velocity norm [m/s]
 */
#define physlib_gamma_vnorm(v)                                                 \
    (sqrt(1.0 / ((1.0 - v / CONST_C) * (1.0 + v / CONST_C))))

/**
 * @brief Evaluate velocity norm from Lorentz factor
 *
 * \f$ v = \sqrt{1 - \frac{1}{\gamma^2}}c\f$
 *
 * where
 *
 * - \f$v\f$ is velocity norm [m/s]
 */
#define physlib_vnorm_gamma(gamma) (sqrt(1.0 - 1.0 / (gamma * gamma)) * CONST_C)

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
#define physlib_gamma_pnorm(m, p) (sqrt(1.0 + (p * p) / (m * m * CONST_C2)))

/**
 * @brief Evaluate momentum norm from Lorentz factor
 *
 */
#define physlib_pnorm_gamma(m, gamma)                                          \
    (sqrt((gamma * gamma - 1.0)) * m * CONST_C)

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
#define physlib_gamma_vpar(m, mu, vpar, B)                                     \
    (sqrt(                                                                     \
        (1.0 + (2.0 * mu * B) / (m * CONST_C2)) /                              \
        ((1.0 - vpar / CONST_C) * (1.0 + vpar / CONST_C))))

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
#define physlib_gamma_ppar(m, mu, ppar, B)                                     \
    (sqrt(1.0 + 2 * mu * B / (m * CONST_C2) + ppar * ppar / (m * m * CONST_C2)))

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
#define physlib_Ekin_gamma(m, gamma) ((gamma - 1.0) * m * CONST_C2)

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
#define physlib_gamma_Ekin(m, ekin) (ekin / (m * CONST_C2) + 1.0)

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
#define physlib_Ekin_pnorm(m, p)                                               \
    ((physlib_gamma_pnorm(m, p) - 1.0) * m * CONST_C2)

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
#define physlib_Ekin_ppar(m, mu, ppar, B)                                      \
    ((physlib_gamma_ppar(m, mu, ppar, B) - 1.0) * m * CONST_C2)

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
#define physlib_vnorm_pnorm(m, p) (p / sqrt(m * m + (p * p) / CONST_C2))

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
#define physlib_pnorm_vnorm(m, v) (m * v * physlib_gamma_vnorm(v))

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
#define physlib_gc_ppar(p, xi) (p * xi)

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
#define physlib_gc_mu(m, p, xi, B) (p * p * (1.0 - xi * xi) / (2 * B * m))

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
#define physlib_gc_p(m, mu, ppar, B)                                           \
    (m * CONST_C * sqrt(pow(physlib_gamma_ppar(m, mu, ppar, B), 2) - 1))

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
#define physlib_gc_xi(m, mu, ppar, B) (ppar / physlib_gc_p(m, mu, ppar, B))

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
#define physlib_gyrolength_p(q, p, B)                                          \
    (math_dot(p, B) / (fabs(q) * math_dot(B, B)))

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
#define phys_gyrolength_ppar(m, q, mu, ppar, B)                                \
    (sqrt(2 * m * mu * physlib_gamma_ppar(m, mu, ppar, B) / B) / fabs(q))

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
#define phys_gyrofreq_pnorm(m, q, p, B)                                        \
    (fabs(q) * B / (m * physlib_gamma_pnorm(m, p)))

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
#define phys_gyrofreq_ppar(m, q, mu, ppar, B)                                  \
    (fabs(q) * B / (m * physlib_gamma_ppar(m, mu, ppar, B)))

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
#define phys_ptoroid_fo(q, R, pphi, psi) (R * pphi + q * psi)

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
#define phys_ptoroid_gc(q, R, ppar, psi, B, Bphi)                              \
    (ppar * R * (Bphi / B) + q * psi)

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
#define phys_vperp_gc(v, vpar) sqrt(pow(v, 2) - pow(vpar, 2))

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
#define phys_pperp_gc(p, ppar) sqrt(pow(p, 2) - pow(ppar, 2))

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
#define phys_ppar_Ekin(m, ekin, mu, B)                                         \
    (sqrt(                                                                     \
        (physlib_gamma_Ekin(m, ekin) * physlib_gamma_Ekin(m, ekin) - 1.0) *    \
            m * m * CONST_C2 -                                                 \
        2.0 * m * mu * B))

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
#define phys_vperp_mu(m, mu, B) sqrt(2 * mu * B / m)

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
#define phys_ppar_mu(m, mu, B, q, gyrof)                                       \
    (m * CONST_C *                                                             \
     sqrt(                                                                     \
         (q * B / (gyrof * m)) * (q * B / (gyrof * m)) -                       \
         2 * mu * B / (m * CONST_C2) - 1))

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
#define phys_ppar_pphi(B, R, B_phi, p_phi, q, psi)                             \
    (B / (R * B_phi) * (p_phi - q * psi))

/**
 * @brief Characteristic time for the radiation reaction force
 *
 * From E. Hirvijoki et al. 2015 Guiding-center transformation of
 * the Abrahams-Lorentz-Dirac radiation reaction force
 * http://de.arxiv.org/abs/1412.1966
 */
#define phys_ald_force_chartime(q, m, B, gamma)                                \
    ((q * q * q * q) * (B * B) /                                               \
     (6.0 * CONST_PI * CONST_E0 * gamma * (m * m * m) * CONST_C3))

/**
 * @brief Evaluate collision parameter [kg^2 m^3 / s^4]
 *
 *\f$c_{ab} = \frac{n_b q_a^2q_b^2 \ln\Lambda_{ab}}{4\pi\epsilon_0^2}\f$
 *
 * where
 *
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_cab(qa, qb, nb, clogab) (                            \
        nb * qa*qa * qb*qb * clogab / ( 4 * CONST_PI * CONST_E0*CONST_E0 ) )

/**
 * @brief Evaluate non-relativistic drag coefficient [m/s^2]
 *
 *\f$Q =-c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0(v_a/v_b) / (m_a m_b v_b^2)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$m_b\f$ is plasma species mass [kg]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_Q(ma, qa, mb, qb, nb, vb, clogab, mu0) (           \
        -mccc_coefs_cab(qa, qb, nb, clogab) * mu0 / ( ma * mb * vb*vb ) )

/**
 * @brief Evaluate derivative of non-relativistic drag coefficient [m/s^2]
 *
 *\f$Q'=-c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0'(v_a/v_b) / (m_a m_b v_b^2)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$m_b\f$ is plasma species mass [kg]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_dQ(ma, qa, mb, qb, nb, vb, clogab, dmu0) (           \
        -mccc_coefs_cab(qa, qb, nb, clogab) * dmu0 / ( ma * mb * vb*vb*vb ) )

/**
 * @brief Evaluate non-relativistic friction coefficient [m/s^2]
 *
 *\f$F = -c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\left(m_a^{-1} + m_b^{-1}\right)
  \mu0(v_a/v_b) / (m_a v_b^2)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$m_b\f$ is plasma species mass [kg]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_F(ma, qa, mb, qb, nb, vb, clogab, mu0) (             \
        -( 1/mb + 1/ma ) * mccc_coefs_cab(qa, qb, nb, clogab) * mu0     \
        / ( ma * vb*vb ) )

/**
 * @brief Evaluate non-relativistic parallel diffusion coefficient [m^2/s^3]
 *
 *\f$D_\parallel=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0(v_a/v_b)
  /(2m_a^2v_a)\f$
 *
 * or \f$D_\parallel = 4c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})/(6\sqrt{\pi}m_a^2vb)
      \;\mathrm{when}\; v_a=0\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_Dpara(ma, qa, va, qb, nb, vb, clogab, mu0) (         \
        ( va > 0 ) ?                                                    \
        mccc_coefs_cab(qa, qb, nb, clogab) * mu0 / ( 2 * ma*ma * va ) : \
        mccc_coefs_cab(qa, qb, nb, clogab) * 4                          \
        / ( 6 * CONST_SQRTPI * ma*ma * vb ) )

/**
 * @brief Evaluate derivative of non-relativistic parallel diffusion
 * coefficient [m/s^2]
 *
 * \f$D_\parallel'=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})(\mu_0'(v_a/v_b)/v_b
  - \mu(v_a/v_b)/v_a)/(2m_a^2v_a)\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_dDpara(ma, qa, va, qb, nb, vb, clogab, mu0, dmu0) (  \
        mccc_coefs_cab(qa, qb, nb, clogab) * ( dmu0/vb - mu0/va )       \
        / ( 2 * ma*ma * va ) )

/**
 * @brief Evaluate non-relativistic perpendicular diffusion
 * coefficient [m^2/s^3]
 *
 *\f$D_\perp=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_1(v_a/v_b)
  /(2m_a^2v_a)\f$
 *
 * or \f$D_\perp = 4c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})/(6\sqrt{\pi}m_a^2vb)
      \;\mathrm{when}\; v_a=0\f$
 *
 * where
 *
 * - \f$m_a\f$ is test particle mass [kg]
 * - \f$q_a\f$ is test particle charge [C]
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$q_b\f$ is plasma species charge [C]
 * - \f$n_b\f$ is plasma species density [m^-3]
 * - \f$v_b\f$ is plasma species thermal velocity [m/s]
 * - \f$\ln\Lambda_{ab}\f$ is Coulomb logarithm.
 */
#define mccc_coefs_Dperp(ma, qa, va, qb, nb, vb, clogab, mu1) (         \
        ( va > 0 ) ?                                                    \
        mccc_coefs_cab(qa, qb, nb, clogab) * mu1 / ( 2 * ma*ma * va ) : \
        mccc_coefs_cab(qa, qb, nb, clogab) * 4                          \
        / ( 6 * CONST_SQRTPI * ma*ma * vb ) )

/**
 * @brief Evaluate guiding center drag coefficient [m/s^2]
 *
 *\f$K = Q + D_\parallel' + 2D_\parallel/va\f$
 *
 * where
 *
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$D_parallel\f$ is  [1/s]
 * - \f$Q\f$ is  [m/s^2]
 */
#define mccc_coefs_K(va, Dpara, dDpara, Q) (    \
        Q + dDpara + 2*Dpara / va )

/**
 * @brief Evaluate pitch collision frequency [1/s]
 *
 *\f$\nu = 2D_\perp/v_a^2\f$
 *
 * where
 *
 * - \f$v_a\f$ is test particle velocity [m/s]
 * - \f$D_perp\f$ is  []
 */
 #define mccc_coefs_nu(va, Dperp) (             \
        2 * Dperp / ( va * va ) )

/**
 * @brief Evaluate spatial diffusion coefficient [m^2/s]
 *
 *\f$D_X = (\frac{1}{2}(D_\parallel-D_\perp)(1-\xi^2) + D_\perp)/\omega_g^2\f$
 *
 * where
 *
 * - \f$xi\f$ is test particle pitch
 * - \f$D_\parallel\f$ is parallel diffusion coefficient [m^2/s^3]
 * - \f$D_perp\f$ is perpendicular diffusion coefficient [m^2/s^3]
 * - \f$\omega_g\f$ is gyrofrequency [1/s]
 */
#define mccc_coefs_DX(xi, Dpara, Dperp, gyrofreq) (              \
        ( 0.5 * ( Dpara - Dperp ) * ( 1 - xi*xi ) + Dperp )      \
        / (gyrofreq*gyrofreq) )

GPU_DECLARE_TARGET_SIMD_UNIFORM(mdata)
static void mccc_coefs_mufun(real mufun[3], real x, mccc_data* mdata);
DECLARE_TARGET_END

GPU_DECLARE_TARGET_SIMD_UNIFORM(nspec, mb, qb, nb, Tb)
inline static void mccc_coefs_clog(
        real* clogab, real ma, real qa, real va, int nspec,
        const real* mb, const real* qb, const real* nb, const real* Tb);
DECLARE_TARGET_END

/**
 * @brief Evaluate Coulomb logarithm.
 *
 * Coulomb logarithm is evaluated separately with respect to each plasma
 * species. It is calculated as a logarithm of the ratio of maximum and
 * minimum impact parameters. Maximum impact parameter is the Debye length
 * and minimum impact parameter is either classical particle radius or
 * inverse of De Broglie wavelength.
 *
 * @param clogab array where evaluated values for Coulomb logarithm are stored.
 * @param ma test particle mass [kg]
 * @param qa test particle charge [C]
 * @param va test particle velocity [m/s]
 * @param nspec number of plasma species
 * @param mb plasma species masses [kg]
 * @param qb plasma species charges [C]
 * @param nb plasma species densities [m^-3]
 * @param Tb plasma species temperatures [J]
 */
inline static void mccc_coefs_clog(
        real* clogab, real ma, real qa, real va, int nspec, const real* mb,
        const real* qb, const real* nb, const real* Tb) {

    /* Evaluate Debye length */
    real sum = 0;
    GPU_SEQUENTIAL_LOOP
    for(int i = 0; i < nspec; i++){
        sum += nb[i] * qb[i] * qb[i] / Tb[i];
    }
    real debyeLength = sqrt(CONST_E0/sum);

    /* Evaluate classical and quantum mechanical impact parameter. The one *
     * that is larger is used to evaluate Coulomb logarithm.               */
    GPU_SEQUENTIAL_LOOP
    for(int i=0; i < nspec; i++){
        real vbar = va * va + 2 * Tb[i] / mb[i];
        real mr   = ma * mb[i] / ( ma + mb[i] );
        real bcl  = fabs( qa * qb[i] / ( 4*CONST_PI*CONST_E0 * mr * vbar ) );
        real bqm  = fabs( CONST_HBAR / ( 2 * mr * sqrt( vbar ) ) );

        if(bcl > bqm){
            clogab[i] = log( debyeLength / bcl );
        }
        else{
            clogab[i] = log( debyeLength / bqm );
        }
    }
}

/**
 * @brief Evaluate special functions needed by collision coefficients
 *
 * This function either evaluates the special functions directly or interpolates
 * them from look-up table which should be initialized with mccc_init() before
 * calling this function.
 *
 * Special functions are
 *
 * - mufun[0] = \f$\mu_0(x) = (\mathrm{erf}(x) -2 x \pi^{-1/2} e^{-x^2})/x^2\f$
 * - mufun[1] = \f$\mu_1(x) = \mathrm{erf}(x) - \frac{1}{2}\mu_0(x)\f$
 * - mufun[2] = \f$\mu_0'(x)\f$
 *
 * @param mufun pointer to array where values are stored
 * @param x argument for the special functions
 * @param mdata pointer to mccc data
 */
inline static void mccc_coefs_mufun(real mufun[3], real x, mccc_data* mdata) {

    if(!mdata->usetabulated && x!= 0) {
        real expm2x = exp(-x*x);
        real erfx   = erf(x);

        mufun[0] = ( erfx - 2 * x * expm2x / CONST_SQRTPI ) / (x*x);
        mufun[1] = erfx - 0.5 * mufun[0];
        mufun[2] = 4 * expm2x / CONST_SQRTPI - 2 * mufun[0] / x;
    }
    else if(mdata->usetabulated && x != 0) {
        // TODO implement me
    }
    else {
        mufun[0] = 0;
        mufun[1] = 0;
        mufun[2] = 4 / ( 3 * CONST_SQRTPI );
    }

}

/**
 * @brief Calculate guiding center equations of motion for a single particle
 *
 * @param ydot output right hand side of the equations of motion in a
 *             6-length array (rdot, phidot, zdot, ppardot, mudot, zetadot)
 * @param y input coordinates in a 6-length array (r, phi, z, vpar, mu, zeta)
 * @param mass mass [kg]
 * @param charge charge [C]
 * @param B_dB magnetic field and derivatives at the guiding center location
 * @param E electric field at the guiding center location
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
DECLARE_TARGET_SIMD
inline static void step_gceom(real* ydot, real* y, real mass, real charge,
                              real* B_dB, real* E, int aldforce) {

    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[4];
    B[2] = B_dB[8];

    real normB = sqrt(math_dot(B, B));
    real gamma = physlib_gamma_ppar(mass, y[4], y[3], normB);

    real gradB[3];
    gradB[0] = (B[0]*B_dB[1] + B[1]*B_dB[5] + B[2]*B_dB[9]) / normB;
    gradB[1] = (B[0]*B_dB[2] + B[1]*B_dB[6] + B[2]*B_dB[10])
               / (normB * y[0]);
    gradB[2] = (B[0]*B_dB[3] + B[1]*B_dB[7] + B[2]*B_dB[11]) / normB;

    real gradBcrossB[3];
    math_cross(gradB, B, gradBcrossB);

    real curlB[3];
    curlB[0] = B_dB[10] / y[0] - B_dB[7];
    curlB[1] = B_dB[3] - B_dB[9];
    curlB[2] = (B[1] - B_dB[2]) / y[0] + B_dB[5];

    real Bstar[3];
    Bstar[0] = B[0] + (y[3] / charge)
                      * (curlB[0] / normB - gradBcrossB[0] / (normB*normB));
    Bstar[1] = B[1] + (y[3] / charge)
                      * (curlB[1] / normB - gradBcrossB[1] / (normB*normB));
    Bstar[2] = B[2] + (y[3] / charge)
                      * (curlB[2] / normB - gradBcrossB[2] / (normB*normB));

    real Estar[3];
    Estar[0] = E[0] - y[4] * gradB[0] / ( charge * gamma );
    Estar[1] = E[1] - y[4] * gradB[1] / ( charge * gamma );
    Estar[2] = E[2] - y[4] * gradB[2] / ( charge * gamma );

    real Bhat[3];
    Bhat[0] = B[0] / normB;
    Bhat[1] = B[1] / normB;
    Bhat[2] = B[2] / normB;

    real BhatDotBstar = math_dot(Bhat, Bstar);

    real EstarcrossBhat[3];
    math_cross(Estar, Bhat, EstarcrossBhat);

    ydot[0] = ( y[3] * Bstar[0] / ( gamma * mass )
                + EstarcrossBhat[0] ) / BhatDotBstar;
    ydot[1] = ( y[3] * Bstar[1] / ( gamma * mass )
                + EstarcrossBhat[1] ) / ( y[0]*BhatDotBstar );
    ydot[2] = ( y[3] * Bstar[2] / ( gamma * mass )
                + EstarcrossBhat[2] ) / BhatDotBstar;
    ydot[3] = charge * math_dot(Bstar,Estar) / BhatDotBstar;
    ydot[4] = 0;
    ydot[5] = charge * normB / ( gamma * mass );

    real t_ald = phys_ald_force_chartime(charge, mass, normB, gamma) * aldforce;
    real C = 2 * y[4] * normB / (mass * CONST_C2);
    ydot[3] += -t_ald * y[3] * C;
    ydot[4] += -2 * t_ald * y[4] * (1 + C);
}


/**
 * @brief Calculate guiding center equations of motion for a single particle
 *
 * @param ydot output right hand side of the equations of motion in a
 *             6-length array (rdot, phidot, zdot, ppardot, mudot, chidot)
 * @param y input coordinates in a 5-length array (r, phi, z, rhopar, mu)
 * @param mass mass
 * @param charge charge
 * @param B_dB magnetic field and derivatives at the guiding center location
 * @param E electric field at the guiding center location
 * @param mhd_dmhd mhd perturbation information evaluated by mhd.c
 * @param aldforce indicates whether Abraham-Lorentz-Dirac force is enabled
 */
DECLARE_TARGET_SIMD
inline static void step_gceom_mhd(
    real* ydot, real* y, real mass, real charge, real* B_dB, real* E,
    real alpha[5], real Phi[5], int aldforce) {

    real B[3];
    B[0] = B_dB[0];
    B[1] = B_dB[4];
    B[2] = B_dB[8];

    real normB = sqrt(math_dot(B, B));
    real gamma = physlib_gamma_ppar(mass, y[4], y[3], normB);

    real gradB[3];
    gradB[0] = (B[0]*B_dB[1] + B[1]*B_dB[5] + B[2]*B_dB[9]) / normB;
    gradB[1] = (B[0]*B_dB[2] + B[1]*B_dB[6] + B[2]*B_dB[10])
               / (normB * y[0]);
    gradB[2] = (B[0]*B_dB[3] + B[1]*B_dB[7] + B[2]*B_dB[11]) / normB;

    real gradBcrossB[3];
    math_cross(gradB, B, gradBcrossB);

    real curlB[3];
    curlB[0] = B_dB[10] / y[0] - B_dB[7];
    curlB[1] = B_dB[3] - B_dB[9];
    curlB[2] = (B[1] - B_dB[2]) / y[0] + B_dB[5];

    real* gradalpha = &alpha[2];

    real gradalphacrossB[3];
    math_cross(gradalpha, B, gradalphacrossB);

    real Bstar[3];
    Bstar[0] = B[0] + gradalphacrossB[0] + alpha[0] * curlB[0]
        + ( y[3] / charge ) * ( curlB[0] / normB
                                - gradBcrossB[0] / (normB*normB));
    Bstar[1] = B[1] + gradalphacrossB[1] + alpha[0] * curlB[1]
        + ( y[3] / charge ) * ( curlB[1] / normB
                                - gradBcrossB[0] / (normB*normB));
    Bstar[2] = B[2] + gradalphacrossB[2] + alpha[0] * curlB[2]
        + ( y[3] / charge ) * ( curlB[2] / normB
                                - gradBcrossB[2] / (normB*normB));

    real Estar[3];
    Estar[0] = E[0] - Phi[2] - y[4] * gradB[0] / ( charge * gamma )
        - B[0] * alpha[1];
    Estar[1] = E[1] - Phi[3] - y[4] * gradB[1] / ( charge * gamma )
        - B[1] * alpha[1];
    Estar[2] = E[2] - Phi[4] - y[4] * gradB[2] / ( charge * gamma )
        - B[2] * alpha[1];

    real Bhat[3];
    Bhat[0] = B[0] / normB;
    Bhat[1] = B[1] / normB;
    Bhat[2] = B[2] / normB;

    real BhatDotBstar = math_dot(Bhat, Bstar);

    real EstarcrossBhat[3];
    math_cross(Estar, Bhat, EstarcrossBhat);

    ydot[0] = (y[3] * Bstar[0] / ( gamma * mass )
               + EstarcrossBhat[0]) / BhatDotBstar;
    ydot[1] = (y[3] * Bstar[1] / ( gamma * mass )
               + EstarcrossBhat[1]) / (y[0]*BhatDotBstar);
    ydot[2] = (y[3] * Bstar[2] / ( gamma * mass )
               + EstarcrossBhat[2]) / BhatDotBstar;
    ydot[3] = charge * math_dot(Bstar,Estar) / BhatDotBstar;
    ydot[4] = 0;
    ydot[5] = charge * normB / (gamma*mass);

    real t_ald = phys_ald_force_chartime(charge, mass, normB, gamma) * aldforce;
    real C = 2 * y[4] * normB / (mass * CONST_C2);
    ydot[3] += -t_ald * y[3] * C;
    ydot[4] += -2 * t_ald * y[4] * (1 + C);
}

#endif
