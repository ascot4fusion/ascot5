/**
 * @file mccc_coefficients.h
 * @brief Routines to evaluate coefficients needed to evaluate collisions.
 */
#ifndef MCCC_COEFFICIENTS_H
#define MCCC_COEFFICIENTS_H

#include <math.h>
#include "../../ascot5.h"

/**
 * @brief Evaluate collision parameter
 *
 *\$fc_{ab} = \frac{n_b q_a^2q_b^2 \ln\Lambda_{ab}}{4\pi\epsilon_0^2}\$f
 *
 * where
 *
 * - \$fq_a\$f is test particle charge [C]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_cab(qa, qb, nb, clogab) (                            \
        nb * qa*qa * qb*qb * clogab / ( 4 * CONST_PI * CONST_E0*CONST_E0 ) )

/**
 * @brief Evaluate non-relativistic drag coefficient []
 *
 *\$fQ =-c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0(v_a/v_b) / (m_a m_b v_b^2)\$f
 *
 * where
 *
 * - \$fm_a\$f is test particle mass [kg]
 * - \$fq_a\$f is test particle charge [C]
 * - \$fm_b\$f is plasma species mass [kg]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$fv_b\$f is plasma species thermal velocity [m/s]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_Q(ma, qa, mb, qb, nb, vb, clogab, mu0) (           \
        -mccc_coefs_cab(qa, qb, nb, clogab) * mu0 / ( ma * mb * vb*vb ) )

/**
 * @brief Evaluate derivative of non-relativistic drag coefficient []
 *
 *\$fQ'=-c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0'(v_a/v_b) / (m_a m_b v_b^2)\$f
 *
 * where
 *
 * - \$fm_a\$f is test particle mass [kg]
 * - \$fq_a\$f is test particle charge [C]
 * - \$fm_b\$f is plasma species mass [kg]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$fv_b\$f is plasma species thermal velocity [m/s]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_dQ(ma, qa, mb, qb, nb, vb, clogab, dmu0) (           \
        -mccc_coefs_cab(qa, qb, nb, clogab) * dmu0 / ( ma * mb * vth*vth ) )

/**
 * @brief Evaluate non-relativistic friction coefficient
 *
 *\$fF = -c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\left(m_a^{-1} + m_b^{-1}\right)
  \mu0(v_a/v_b) / (m_a v_b^2)\$f
 *
 * where
 *
 * - \$fm_a\$f is test particle mass [kg]
 * - \$fq_a\$f is test particle charge [C]
 * - \$fm_b\$f is plasma species mass [kg]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$fv_b\$f is plasma species thermal velocity [m/s]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_F(ma, qa, mb, qb, nb, vb, clogab, mu0) (             \
        -( 1/mb + 1/ma ) * mccc_coefs_cab(qa, qb, nb, clogab) * mu0     \
        / ( ma * vb*vb ) )

/**
 * @brief Evaluate non-relativistic parallel diffusion coefficient []
 *
 *\$fD_\parallel=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_0(v_a/v_b)
  /(2m_a^2v_a)\$f
 *
 * or $fD_\parallel = 4c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})/(6\sqrt{\pi}m_a^2vb)
      \;\mathrm{when}\; v_a=0\$f
 *
 * where
 *
 * - \$fm_a\$f is test particle mass [kg]
 * - \$fq_a\$f is test particle charge [C]
 * - \$fv_a\$f is test particle velocity [m/s]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$fv_b\$f is plasma species thermal velocity [m/s]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_Dpara(ma, qa, va, qb, nb, vb, clogab, mu0) (         \
        ( va > 0 ) ?                                                    \
        mccc_coefs_cab(qa, qb, nb, clogab) * mu0 / ( 2 * ma*ma * va ) : \
        mccc_coefs_cab(qa, qb, nb, clogab) * 4                          \
        / ( 6 * CONST_SQRTPI * ma*ma * vb ) )

/**
 * @brief Evaluate non-relativistic parallel diffusion coefficient []
 *
 * \$fD_\parallel'=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})(\mu_0'(v_a/v_b)/v_b
  - \mu(v_a/v_b)/v_a)/(2m_a^2v_a)\$f
 *
 * where
 *
 * - \$fm_a\$f is test particle mass [kg]
 * - \$fq_a\$f is test particle charge [C]
 * - \$fv_a\$f is test particle velocity [m/s]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$fv_b\$f is plasma species thermal velocity [m/s]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_dDpara(ma, qa, va, qb, nb, vb, clogab, mu0, dmu0) (  \
        mccc_coefs_cab(qa, qb, nb, clogab) * ( dmu0/vb - mu0/va )       \
        / ( 2 * ma*ma * va ) )

/**
 * @brief Evaluate non-relativistic perpendicular diffusion coefficient []
 *
 *\$fD_\perp=c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})\mu_1(v_a/v_b)
  /(2m_a^2v_a)\$f
 *
 * or $fD_\perp = 4c_{ab}(q_a,q_b,n_b,\ln\Lambda_{ab})/(6\sqrt{\pi}m_a^2vb)
      \;\mathrm{when}\; v_a=0\$f
 *
 * where
 *
 * - \$fm_a\$f is test particle mass [kg]
 * - \$fq_a\$f is test particle charge [C]
 * - \$fv_a\$f is test particle velocity [m/s]
 * - \$fq_b\$f is plasma species charge [C]
 * - \$fn_b\$f is plasma species density [m^-3]
 * - \$fv_b\$f is plasma species thermal velocity [m/s]
 * - \$f\ln\Lambda_{ab}\$f is Coulomb logarithm.
 */
#define mccc_coefs_Dperp(ma, qa, va, qb, nb, vb, clogab, mu1) (         \
        ( va > 0 ) ?                                                    \
        mccc_coefs_cab(qa, qb, nb, clogab) * mu1 / ( 2 * ma*ma * va ) : \
        mccc_coefs_cab(qa, qb, nb, clogab) * 4                          \
        / ( 6 * CONST_SQRTPI * ma*ma * vb ) )

/**
 * @brief Evaluate guiding center drag coefficient []
 *
 *\$fK = Q + D_\parallel' + 2D_\parallel/va\$f
 *
 * where
 *
 * - \$fv_a\$f is test particle velocity [m/s]
 * - \$fQ\$f is  [C]
 * - \$D_parallel\$f is  [m^-3]
 */
#define mccc_coefs_K(va, Dpara, dDpara, Q) (    \
        Q + dDpara + 2*Dpara / va )

/**
 * @brief Evaluate pitch collision frequency
 */
 #define mccc_coefs_nu(va, Dperp) (             \
        2 * Dperp / ( va * va ) )

/**
 * @brief Evaluate spatial diffusion coefficient
 */
#define mccc_coefs_Dx(xi, Dpara, Dperp, gyrofreq) (              \
        ( 0.5 * ( Dpara - Dperp ) * ( 1 - xi*xi ) + Dperp )      \
        / (gyrofreq*gyrofreq) )

#pragma omp declare target

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
#pragma omp declare simd
static void mccc_coefs_clog(real* clogab, real ma, real qa, real va, int nspec,
                            const real* mb, const real* qb, const real* nb,
                            const real* Tb) {

    /* Evaluate Debye length */
    real sum = 0;
    for(int i = 0; i < nspec; i++){
        sum += (nb[i] * qb[i] * qb[i])/(Tb[i]);
    }
    real debyeLength = sqrt(CONST_E0/sum);

    /* Evaluate classical and quantum mechanical impact parameter. The one *
     * that is larger is used to evaluate Coulomb logarithm.               */
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

#pragma omp declare simd
static void mccc_coefs_mufun(real mufun[3], real x, real* coldata) {

    if(coldata == NULL && x!= 0) {
        real expm2x = exp(-x*x);
        real erfx   = erf(x);

        mufun[0] = ( erfx - 2 * x * expm2x / CONST_SQRTPI ) / (x*x);
        mufun[1] = erfx - 0.5 * mufun[0];
        mufun[2] = 4 * expm2x / CONST_SQRTPI - 2 * mufun[0] / x;
    }
    else if(coldata != NULL && x != 0) {
        // TODO implement me
    }
    else {
        mufun[0] = 0;
        mufun[1] = 0;
        mufun[2] = 4 / ( 3 * CONST_SQRTPI );
    }

}

#pragma omp end declare target

#endif
