/**
 * @file interact.c
 * @brief Functions for Coulomb collision evaluation
 */
#include <math.h>
#include <stdlib.h>
//#include <mkl_vsl.h>
#include "ascot5.h"
#include "interact.h"
#include "plasma_1d.h"
#include "B_field.h"
#include "math.h"
#include "particle.h"

/**
 * @brief Integrate a full orbit step of Coulomb interactions
 *
 * This function integrates the Coulomb drag and diffusion terms in the
 * differential equation (31) in the ASCOT paper. The deterministic term is
 * integrated with Euler method and the stochastic terms with Euler-Maruyama.
 *
 * @param p particle struct that will be updated
 * @param t time
 * @param h length of time step
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to 1D plasma data
 * 
 * @todo Move deterministic term into the orbit integration function
 */
void interact_step_fo_euler(particle_simd_fo* p, real t, real h,
                            B_field_data* Bdata, plasma_1d_data* pdata) {

    int i;
    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real v[3];
            real vprev[3];
            real psi[1];
            real rho[1];

            vprev[0] = p->rdot[i] * cos(p->phi[i])
	    	- (p->phidot[i]*p->r[i]) * sin(p->phi[i]);
            vprev[1] = p->rdot[i] * sin(p->phi[i])
	    	+ (p->phidot[i]*p->r[i]) * cos(p->phi[i]);
            vprev[2] = p->zdot[i];

            real absv = sqrt(vprev[0]*vprev[0]+vprev[1]*vprev[1]
                        +vprev[2]*vprev[2]);
            real absv2 = absv*absv;

            B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);
            B_field_eval_rho(rho, psi[0], Bdata);

            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];

            int j;
            for(j = 0; j < pdata->n_species; j++) {
                temp[j] = plasma_1d_eval_temp(rho[0], j, pdata);
                dens[j] = plasma_1d_eval_dens(rho[0], j, pdata);
            }

            real nu = interact_nu_fo(absv, p->mass[i], p->charge[i], temp, dens,
                                  pdata->charge, pdata->mass, pdata->n_species);

            real D_par = interact_D_par(absv, p->mass[i], p->charge[i], temp,
                                        dens, pdata->charge, pdata->mass,
                                        pdata->n_species);
            real D_perp = interact_D_par(absv, p->mass[i], p->charge[i], temp,
                                         dens, pdata->charge, pdata->mass,
                                         pdata->n_species);

            v[0] = vprev[0];
            v[1] = vprev[1];
            v[2] = vprev[2];

            /* Deterministic step */
            v[0] -= nu * vprev[0] * h;
            v[1] -= nu * vprev[1] * h;
            v[2] -= nu * vprev[2] * h;

            /* Stochastic step */
            real A = sqrt(D_par);
            real B = sqrt(D_perp);
            real w1 = drand48();
            if(w1 > 0.5)
                w1 = sqrt(h);
            else
                w1 = -sqrt(h);
            real w2 = drand48();
            if(w2 > 0.5)
                w2 = sqrt(h);
            else
                w2 = -sqrt(h);
            real w3 = drand48();
            if(w3 > 0.5)
                w3 = sqrt(h);
            else
                w3 = -sqrt(h);
            v[0] +=   w1 * (B+(A-B)*vprev[0]*vprev[0]/absv2)
                       + w2 * ((A-B)*vprev[0]*vprev[1]/absv2)
                       + w3 * ((A-B)*vprev[0]*vprev[2]/absv2);
            v[1] +=   w1 * ((A-B)*vprev[1]*vprev[0]/absv2)
                       + w2 * (B+(A-B)*vprev[1]*vprev[1]/absv2)
                       + w3 * ((A-B)*vprev[1]*vprev[2]/absv2);
            v[2] +=   w1 * ((A-B)*vprev[2]*vprev[0]/absv2)
                       + w2 * ((A-B)*vprev[2]*vprev[1]/absv2)
                       + w3 * (B+(A-B)*vprev[2]*vprev[2]/absv2);

            p->rdot[i] = v[0] * cos(p->phi[i]) + v[1] * sin(p->phi[i]);
            p->phidot[i] = (-v[0] * sin(p->phi[i]) + v[1] * cos(p->phi[i]) ) / p->r[i];
            p->zdot[i] = v[2];
        }
    }
}

void interact_step_gc_euler(particle_simd_gc* p, real t, real h,
                            B_field_data* Bdata, plasma_1d_data* pdata) {
    real w1[1][NSIMD];
    real w2[1][NSIMD];
/*    VSLStreamStatePtr r;
    vslNewStream(&r, VSL_BRNG_SFMT19937, 123);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, NSIMD, &w1[0][0], 0.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, NSIMD, &w2[0][0], 0.0, 1.0);*/

    int i;
    for(i = 0; i < NSIMD; i++) {
        w1[0][i] = drand48();
        w2[0][i] = drand48();
    }

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real v[3];
            real vprev[3];
            real B[3];
            real psi[1];
            real rho[1];

            real absB = sqrt(p->B_r[i]*p->B_r[i] + p->B_phi[i]*p->B_phi[i]
                             + p->B_z[i]*p->B_z[i]);

            real vperp2 = 2 * absB * p->mu[i] / p->mass[i];
            real absv2 = p->vpar[i]*p->vpar[i] + vperp2;
            real absv = sqrt(absv2);

            real pitch = p->vpar[i] / sqrt(absv2);

            B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);
            B_field_eval_rho(rho, psi[0], Bdata);

            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];

            int j;
            for(j = 0; j < pdata->n_species; j++) {
                temp[j] = plasma_1d_eval_temp(rho[0], j, pdata);
                dens[j] = plasma_1d_eval_dens(rho[0], j, pdata);
            }

            real nu = interact_nu_gc(absv, p->mass[i], p->charge[i], temp, dens,
                                  pdata->charge, pdata->mass, pdata->n_species);

            real D_par, dD_par;
            interact_D_par_dD_par(&D_par, &dD_par, absv,
                           p->mass[i], p->charge[i], 
                           temp, dens, pdata->charge,
                           pdata->mass, pdata->n_species);
            real D_perp = interact_D_perp(absv, p->mass[i], p->charge[i], temp,
                                         dens, pdata->charge, pdata->mass,
                                         pdata->n_species);

            real dv_det = (-nu + dD_par + 2 * D_par / (p->mass[i] * absv))/p->mass[i];
            real dv_sto = sqrt(2*D_par) / p->mass[i];

            if(w1[0][i] > 0.5)
                w1[0][i] = sqrt(h);
            else
                w1[0][i] = -sqrt(h);

            absv = absv + dv_det * h + dv_sto * w1[0][i];
            
            real dp_det = -pitch * 2 * D_perp / (p->mass[i]*p->mass[i]*absv2);
            real dp_sto = sqrt((1-pitch*pitch)*2*D_perp / (p->mass[i]*p->mass[i]*absv2));

            if(w2[0][i] > 0.5)
                w2[0][i] = sqrt(h);
            else
                w2[0][i] = -sqrt(h);

            pitch = pitch + dp_det * h + dp_sto * w2[0][i];

            if(pitch > 1) {
                pitch = 2 - pitch;
            } else if(pitch < -1) {
                pitch = -2 - pitch;
            }

            absv2 = absv * absv; 
            p->mu[i] = (1-pitch*pitch) * absv2 * p->mass[i] / (2 * absB);
            p->vpar[i] = pitch * sqrt(absv2);
        }
    }
}

void interact_step_gc_euler_ascot4(particle_simd_gc* p, real t, real h,
                            B_field_data* Bdata, plasma_1d_data* pdata) {
    real w1[1][NSIMD];
    real w2[1][NSIMD];
/*    VSLStreamStatePtr r;
    vslNewStream(&r, VSL_BRNG_SFMT19937, 123);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, NSIMD, &w1[0][0], 0.0, 1.0);
    vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, r, NSIMD, &w2[0][0], 0.0, 1.0);*/

    int i;
    for(i = 0; i < NSIMD; i++) {
        w1[0][i] = drand48();
        w2[0][i] = drand48();
    }

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real v[3];
            real vprev[3];
            real B[3];
            real psi[1];
            real rho[1];

            real prevmu = p->mu[i];
            real prevvpar = p->vpar[i];
            B_field_eval_B(B, p->r[i], p->phi[i], p->z[i], Bdata);
            real absB = sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2]);

            real vperp2 = 2 * absB * p->mu[i] / p->mass[i];
            real absv2 = p->vpar[i]*p->vpar[i] + vperp2;
            real absv = sqrt(absv2);
            real Ekin = 0.5 * p->mass[i] * absv2;

            real pitch = p->vpar[i] / sqrt(absv2);

            B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);
            B_field_eval_rho(rho, psi[0], Bdata);

            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];

            int j;
            for(j = 0; j < pdata->n_species; j++) {
                temp[j] = plasma_1d_eval_temp(rho[0], j, pdata);
                dens[j] = plasma_1d_eval_dens(rho[0], j, pdata);
            }

            real debye_length2 = interact_debye_length2(temp, dens,
                pdata->charge, pdata->n_species);
            real a1 = 2 / (p->mass[i] * absv);
            real a2 = 2 * Ekin / (p->mass[i] * absv);
            real a3 = 1 / p->mass[i];

            real de1 = 0;
            real de2 = 0;
            real freq = 0;

            for(j = 0; j < pdata->n_species; j++) {
                real coulomb_log = interact_coulomb_log(absv, p->mass[i],
                    p->charge[i], temp[j], dens[j], pdata->charge[j],
                    pdata->mass[j], debye_length2);
                real dcoulomb_log = interact_dcoulomb_log(absv, p->mass[i],
                    p->charge[i], temp[j], dens[j], pdata->charge[j],
                    pdata->mass[j], debye_length2);
                real gs = pdata->charge[j] * pdata->charge[j] * dens[j]
                          * coulomb_log;
                real as = sqrt(pdata->mass[j] / (2 * CONST_KB * temp[j]));
                real x = absv * as;
                de1 += gs*(interact_psi(x)*(a1 - a2/(CONST_KB * temp[j])
                                            + a3 * dcoulomb_log/coulomb_log)
                           + a3*as*interact_dpsi(x));
                de2 += gs*interact_psi(x);
                freq += gs*as*interact_nu_ascot4(x);
            }

            de1 *= p->charge[i] * p->charge[i] / (4 * math_pi * CONST_E0
                                                  * CONST_E0);
            de2 *= Ekin * p->charge[i] * p->charge[i]
                   / (p->mass[i] * absv * math_pi * CONST_E0 * CONST_E0);
            freq *= p->charge[i] * p->charge[i]
                    / (4 * math_pi * CONST_E0 * CONST_E0 * absv * absv
                       * p->mass[i] * p->mass[i]);

            if(w1[0][i] > 0.5)
                w1[0][i] = sqrt(h);
            else
                w1[0][i] = -sqrt(h);

            Ekin = Ekin + de1 * h + sqrt(de2) * w1[0][i];

            if(w2[0][i] > 0.5)
                w2[0][i] = sqrt(h);
            else
                w2[0][i] = -sqrt(h);

            pitch = pitch - freq * pitch * h
                + w2[0][i] * sqrt((1-pitch*pitch)*freq);

            if(pitch > 1) {
                pitch = 2 - pitch;
            } else if(pitch < -1) {
                pitch = -2 - pitch;
            }

            absv2 = 2 * Ekin / p->mass[i]; 
            p->mu[i] = (1-pitch*pitch) * absv2 * p->mass[i] / (2 * absB);
            p->vpar[i] = pitch * sqrt(absv2);
        }
    }
}

/**
 * @brief Compute the parallel diffusion coefficient
 *
 * This function computes the parallel diffusion coefficient according to
 * Hinton (90).
 *
 * @param v velocity of the test particle (m/s)
 * @param m mass of the test particle (kg)
 * @param q charge of the test particle (C)
 * @param temp plasma species temperatures (K)
 * @param dens plasma species densities (1/m3)
 * @param charge plasma species charges (C)
 * @param mass plasma species masses (kg)
 * @param n_species number of plasma species
 */
real interact_D_par(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species) {
    int i;
    real debye_length2 = interact_debye_length2(temp, dens, charge, n_species);
    real D_par = 0;
    for(i = 0; i < n_species; i++) {
        real x = v / sqrt(2 * CONST_KB * temp[i] / mass[i]);
        real coulomb_log = interact_coulomb_log(v, m, q, temp[i], dens[i],
                                             charge[i], mass[i], debye_length2);
        D_par += dens[i] * charge[i] * charge[i] * coulomb_log
                 * interact_psi(x);
    }
    D_par *= q * q / (4 * math_pi * CONST_E0 * CONST_E0 * v);
    return D_par;
}

/**
 * @brief Compute the parallel diffusion coefficient and derivative
 *
 * This function computes the parallel diffusion coefficient according to
 * Hinton (90).
 *
 * @param v velocity of the test particle (m/s)
 * @param m mass of the test particle (kg)
 * @param q charge of the test particle (C)
 * @param temp plasma species temperatures (K)
 * @param dens plasma species densities (1/m3)
 * @param charge plasma species charges (C)
 * @param mass plasma species masses (kg)
 * @param n_species number of plasma species
 */
void interact_D_par_dD_par(real* D_par, real* dD_par,
                          real v, real m, real q, real temp[MAX_SPECIES],
                          real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                          real mass[MAX_SPECIES], int n_species) {
    int i;
    real debye_length2 = interact_debye_length2(temp, dens, charge, n_species);
    real sum_D_par = 0;
    real sum_dD_par = 0;
    for(i = 0; i < n_species; i++) {
        real vb = sqrt(2 * CONST_KB * temp[i] / mass[i]);
        real x = v / vb;
        real coulomb_log = interact_coulomb_log(v, m, q, temp[i], dens[i],
                                             charge[i], mass[i], debye_length2);
        real temp_sum = dens[i] * charge[i] * charge[i] * coulomb_log
                        * interact_psi(x);
        sum_D_par += temp_sum;
        sum_dD_par += -temp_sum/v + dens[i] * charge[i] * charge[i] 
                      * coulomb_log * interact_dpsi(x) / vb;
    }
    real coeff = q * q / (4 * math_pi * CONST_E0 * CONST_E0 * v);
    *D_par = coeff * sum_D_par;
    *dD_par = coeff * sum_dD_par;
}

/**
 * @brief Compute the parallel diffusion coefficient
 *
 * This function computes the perpendicular diffusion coefficient according to
 * Hinton (89).
 *
 * @see interact_D_par
 */
real interact_D_perp(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species) {
    int i;
    real debye_length2 = interact_debye_length2(temp, dens, charge, n_species);
    real D_perp = 0;
    for(i = 0; i < n_species; i++) {
        real x = v / sqrt(2 * CONST_KB * temp[i] / mass[i]);
        real coulomb_log = interact_coulomb_log(v, m, q, temp[i], dens[i],
                                             charge[i], mass[i], debye_length2);
        D_perp += dens[i] * charge[i] * charge[i] * coulomb_log
                  * (erf(x) - interact_psi(x));
    }
    D_perp *= q * q / (8 * math_pi * CONST_E0 * CONST_E0 * v);
    return D_perp;
}

/**
 * @brief Compute the slowing-down rate for full orbit
 *
 * This function computes the slowing-down rate according to Hinton (88).
 *
 * @see interact_D_par
 */
real interact_nu_fo(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species) {
    int i;
    real debye_length2 = interact_debye_length2(temp, dens, charge, n_species);
    real nu = 0;
    for(i = 0; i < n_species; i++) {
        real x = v / sqrt(2 * CONST_KB * temp[i] / mass[i]);
        real coulomb_log = interact_coulomb_log(v, m, q, temp[i], dens[i],
                                             charge[i], mass[i], debye_length2);
        nu += dens[i] * charge[i] * charge[i] * coulomb_log 
              / temp[i]
              * (1 + m / mass[i]) * interact_psi(x);
    }
    nu *= q * q / (4 * math_pi * CONST_E0 * CONST_E0 * CONST_KB);
    return nu;
}

/**
 * @brief Compute the slowing-down rate for guiding center
 *
 * This function computes the slowing-down rate according to ASCOT4 paper (28).
 *
 * @see interact_D_par
 */
real interact_nu_gc(real v, real m, real q, real temp[MAX_SPECIES],
                    real dens[MAX_SPECIES], real charge[MAX_SPECIES],
                    real mass[MAX_SPECIES], int n_species) {
    int i;
    real debye_length2 = interact_debye_length2(temp, dens, charge, n_species);
    real nu = 0;
    for(i = 0; i < n_species; i++) {
        real x = v / sqrt(2 * CONST_KB * temp[i] / mass[i]);
        real coulomb_log = interact_coulomb_log(v, m, q, temp[i], dens[i],
                                             charge[i], mass[i], debye_length2);
        nu += dens[i] * charge[i] * charge[i] * coulomb_log 
              / temp[i]
              * interact_psi(x);
    }
    nu *= q * q / (4 * math_pi * CONST_E0 * CONST_E0 * CONST_KB);
    return nu;
}

/**
 * @brief Compute the Coulomb logarithm 
 *
 * This function computes the Coulomb logarithm (formulas from ASCOT4).
 *
 * @param v velocity of the test particle (m/s)
 * @param m mass of the test particle (kg)
 * @param q charge of the test particle (C)
 * @param temp_s single plasma species temperature (K)
 * @param dens_s single plasma species density (1/m3)
 * @param charge_s single plasma species charge (C)
 * @param mass_s single plasma species mass (kg)
 * @param debye_length2 Debye length squared
 */
real interact_coulomb_log(real v, real m, real q, real temp_s, real dens_s,
                          real charge_s, real mass_s, real debye_length2) {
    real m_r = m * mass_s / (m + mass_s);
    real u2 = v * v + 2 * CONST_KB * temp_s / mass_s;
    real lambda_db2 = pow(CONST_HBAR / m_r, 2) / u2;
    real b0_2 = pow(q * charge_s / (4*math_pi*CONST_E0*m_r*u2), 2);
    if(b0_2 > lambda_db2) {
        return 0.5 * log(debye_length2 / b0_2);
    }
    else {
        return 0.5 * log(debye_length2 / lambda_db2);
    }
}

real interact_dcoulomb_log(real v, real m, real q, real temp_s, real dens_s,
                          real charge_s, real mass_s, real debye_length2) {
    real m_r = m * mass_s / (m + mass_s);
    real u2 = v * v + 2 * CONST_KB * temp_s / mass_s;
    real lambda_db2 = pow(CONST_HBAR / m_r, 2) / u2;
    real b0_2 = pow(q * charge_s / (4*math_pi*CONST_E0*m_r*u2), 2);
    if(b0_2 > lambda_db2) {
        return 2 * v / u2;
    }
    else {
        return v / u2;
    }
}

/**
 * @brief Compute the Debye length
 *
 * This function computes the Debye length for the plasma
 * (formulas from ASCOT4).
 *
 * @param temp plasma species temperatures (K)
 * @param dens plasma species densities (1/m3)
 * @param charge plasma species charges (C)
 * @param n_species number of plasma species
 */
real interact_debye_length2(real temp[MAX_SPECIES], real dens[MAX_SPECIES],
                            real charge[MAX_SPECIES], int n_species) {
    int i;
    real sum = 0;
    for(i = 0; i < n_species; i++) {
        sum += dens[i] * charge[i] * charge[i] / (CONST_KB * temp[i]);
    }
    return CONST_E0 / sum;
}

/**
 * @brief Evaluate the Chandrasekhar function
 *
 * This function evaluates the Chandrasekhar function in Hinton (91).
 */
real interact_psi(real x) {
    real x2 = x*x;
    return (erf(x) - x * 2 * exp(-x2) / sqrt(math_pi)) / (2 * x2);
}

/**
 * @brief Evaluate the derivative of the Chandrasekhar function
 *
 * This function evaluates the derivative of the Chandrasekhar function.
 */
real interact_dpsi(real x) {
    real x2 = x*x;
    real x3 = x*x2;

    return 1/x3*(2/sqrt(math_pi)*exp(-x2)*(x3+x)-erf(x));
}

real interact_nu_ascot4(real x) {
    real x2 = x*x;
    real x3 = x*x2;

    return ((2*x2-1)*erf(x) + 2/sqrt(math_pi)*x*exp(-x2))/(2*x3);
}
