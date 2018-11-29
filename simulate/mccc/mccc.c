/**
 * @file mccc.c
 * @brief Interface for using mccc package within ascot5
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../error.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../physlib.h"
#include "mccc.h"
#include "mccc_wiener.h"
#include "mccc_push.h"
#include "mccc_coefs.h"

/**
 * @brief Evaluates collision coefficients in fo picture
 *
 * Finds the rho coordinate first and uses it to evaluate plasma parameters
 * that are then used to evaluate Coulomb logarithm and collision coefficients.
 *
 * The coefficients are returned in arrays whose format is D_is = D[i*s] where
 * i is the particle SIMD position and s is species (maximum is MAX_SPECIES
 * defined in ascot5.h).
 *
 * Operations are performed simultaneously for NSIMD markers.
 *
 * @param p pointer to SIMD_fo struct containing the markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param clogab pointer to array storing the evaluated Coulomb logarithms
 * @param F pointer to array storing the evaluated Fokker-Planck coefficient [m/s^2]
 * @param Dpara pointer to array storing the evaluated parallel diffusion coefficient [m^2/s^3]
 * @param Dperp pointer to array storing the evaluated perpendicular diffusion coefficient [m²/s^3]
 * @param K pointer to array storing the evaluated friction coefficient [m/s^2]
 * @param nu pointer to array storing the evaluated pitch collision frequency [1/s]
 */
void mccc_update_fo(particle_simd_fo* p, B_field_data* Bdata, plasma_data* pdata, real* coldata, 
                    real* clogab, real* F, real* Dpara, real* Dperp, real* K, real* nu){
    int i;
    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            /* Update background data */
            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];

            int n_species = plasma_get_n_species(pdata);
            real* q_species = plasma_get_species_charge(pdata);
            real* m_species = plasma_get_species_mass(pdata);

            // Electron and ion temperature
            temp[0] = plasma_eval_temp(p->rho[i], 0, pdata)*CONST_KB;
            temp[1] = plasma_eval_temp(p->rho[i], 1, pdata)*CONST_KB;

            // Electron density
            dens[0] = plasma_eval_dens(p->rho[i], 0, pdata);

            // Ion densities (and temperatures)
            int j;
            for(j = 1; j < n_species; j++) {
                dens[j] = plasma_eval_dens(p->rho[i], j, pdata);
                temp[j] = temp[1];
            }

            /* Evaluate coefficients */
            real va = sqrt(p->rdot[i]*p->rdot[i] + (p->r[i]*p->phidot[i])*(p->r[i]*p->phidot[i]) + p->zdot[i]*p->zdot[i]);
            mccc_coefs_clog(p->mass[i], p->charge[i], va, m_species, q_species, dens, temp, &clogab[i*MAX_SPECIES], n_species);
            mccc_coefs_fo(p->mass[i], p->charge[i], va, m_species, q_species, dens, temp, &clogab[i*MAX_SPECIES], n_species, coldata,
                          &F[i*MAX_SPECIES], &Dpara[i*MAX_SPECIES], &Dperp[i*MAX_SPECIES], &K[i*MAX_SPECIES], &nu[i*MAX_SPECIES]);
        }
    }
}

/**
 * @brief Evaluates collision frequency in gc picture
 *
 * Finds the rho coordinate first and uses it to evaluate plasma parameters
 * that are then used to evaluate Coulomb logarithm and collision coefficients.
 *
 * @param p pointer to SIMD_gc struct containing the markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param nu pointer to pitch collision frequency
 * @param i index of the marker in simd array
 */
void mccc_collfreq_gc(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, real* coldata, 
                      real* nu, int i){

    /* Update background data */
    real B[3];
    real xi;

    B[0] = p->B_r[i];
    B[1] = p->B_phi[i];
    B[2] = p->B_z[i];

    real temp[MAX_SPECIES];
    real dens[MAX_SPECIES];

    int n_species = plasma_get_n_species(pdata);
    real* q_species = plasma_get_species_charge(pdata);
    real* m_species = plasma_get_species_mass(pdata);

    // Electron and ion temperature
    temp[0] = plasma_eval_temp(p->rho[i], 0, pdata)*CONST_KB;
    temp[1] = plasma_eval_temp(p->rho[i], 1, pdata)*CONST_KB;

    // Electron density
    dens[0] = plasma_eval_dens(p->rho[i], 0, pdata);

    // Ion densities (and temperatures)
    int j;
    for(j = 1; j < n_species; j++) {
        dens[j] = plasma_eval_dens(p->rho[i], j, pdata);
        temp[j] = temp[1];
    }

    /* Evaluate coefficients */
    real Bnorm = math_norm(B);
    real t = 2*p->mu[i]*Bnorm*p->mass[i];
    real va = sqrt(p->vpar[i]*p->vpar[i] + t*t);
    xi = p->vpar[i]/va;

    real clogab[MAX_SPECIES];
    real Dparab[MAX_SPECIES];
    real Kb[MAX_SPECIES];
    real nub[MAX_SPECIES];
    real DXb[MAX_SPECIES];
    mccc_coefs_clog(p->mass[i], p->charge[i], va, m_species, q_species, dens, temp, clogab, n_species);
    mccc_coefs_gcfixed(p->mass[i], p->charge[i], va, xi, m_species, q_species, dens, temp, Bnorm, clogab, n_species, coldata,
                       Dparab,DXb,Kb,nub);

    *nu = 0;
    for(j = 0; j < n_species; j++) {
        *nu += nub[j];
    }


}

/**
 * @brief Evaluates collision coefficients in gc picture
 *
 * Finds the rho coordinate first and uses it to evaluate plasma parameters
 * that are then used to evaluate Coulomb logarithm and collision coefficients.
 *
 * The coefficients are returned in arrays whose format is D_is = D[i*s] where
 * i is the particle SIMD position and s is species (maximum is MAX_SPECIES
 * defined in ascot5.h).
 *
 * Operations are performed simultaneously for NSIMD markers.
 *
 * @param p pointer to SIMD_gc struct containing the markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param clogab pointer to array storing the evaluated Coulomb logarithms
 * @param Dpara pointer to array the evaluated storing parallel diffusion coefficient [m^2/s^3]
 * @param DX pointer to array storing the classical diffusion coefficient [m^2/s]
 * @param K pointer to array storing the evaluated friction coefficient [m/s^2]
 * @param nu pointer to array storing the evaluated pitch collision frequency [1/s]
 * @param dQ pointer to array storing the evaluated derivative dQ/dv [1/s]
 * @param dDpara pointer to array storing the evaluated derivative dDpara/dv [m/s^2]
 *
 */
void mccc_update_gc(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, real* coldata, 
                    real* clogab, real* Dpara, real* DX, real* K, real* nu, real* dQ, real* dDpara){
    int i;
    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {

            /* Update background data */
            real B[3];
            real xi;

            B[0] = p->B_r[i];
            B[1] = p->B_phi[i];
            B[2] = p->B_z[i];

            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];
            int n_species = plasma_get_n_species(pdata);
            real* q_species = plasma_get_species_charge(pdata);
            real* m_species = plasma_get_species_mass(pdata);

            // Electron and ion temperature
            temp[0] = plasma_eval_temp(p->rho[i], 0, pdata)*CONST_KB;
            temp[1] = plasma_eval_temp(p->rho[i], 1, pdata)*CONST_KB;

            // Electron density
            dens[0] = plasma_eval_dens(p->rho[i], 0, pdata);

            // Ion densities (and temperatures)
            int j;
            for(j = 1; j < n_species; j++) {
                dens[j] = plasma_eval_dens(p->rho[i], j, pdata);
                temp[j] = temp[1];
            }

            /* Evaluate coefficients */
            real Bnorm = math_norm(B);
            real t = 2*p->mu[i]*Bnorm*p->mass[i];
            real va = sqrt(p->vpar[i]*p->vpar[i] + t*t);
            xi = p->vpar[i]/va;

            mccc_coefs_clog(p->mass[i], p->charge[i], va, m_species, q_species, dens, temp, &clogab[i*MAX_SPECIES], n_species);
            mccc_coefs_gcadaptive(p->mass[i], p->charge[i], va, xi, m_species, q_species, dens, temp, Bnorm, &clogab[i*MAX_SPECIES], n_species, coldata,
            &Dpara[i*MAX_SPECIES], &DX[i*MAX_SPECIES], &K[i*MAX_SPECIES], &nu[i*MAX_SPECIES], &dQ[i*MAX_SPECIES], &dDpara[i*MAX_SPECIES]);
        }
    }
}

/**
 * @brief Evaluates collisions in fo picture with fixed time-step
 *
 * This function first evaluates collision coefficients (see mccc_update_fo)
 * and then evaluates collisions using Euler-Maruyama method and updates
 * marker state.
 *
 * Operations are performed simultaneously for NSIMD markers.
 *
 * @param p pointer to SIMD_fo struct containing the markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param h pointer to time step values [s]
 */
void mccc_step_fo_fixed(particle_simd_fo* p, B_field_data* Bdata, plasma_data* pdata, random_data* rdata, real* coldata, real* h){
    int i;
    real rnd[3*NSIMD];
    random_normal_simd(rdata, 3*NSIMD, rnd);

    int n_species = plasma_get_n_species(pdata);
    real* q_species = plasma_get_species_charge(pdata);
    real* m_species = plasma_get_species_mass(pdata);

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            /* Gather relevant particle data */
            real vin[3], vout[3], va;
            vin[0] = p->rdot[i] * cos(p->phi[i]) - (p->phidot[i]*p->r[i]) * sin(p->phi[i]);
            vin[1] = p->rdot[i] * sin(p->phi[i]) + (p->phidot[i]*p->r[i]) * cos(p->phi[i]);
            vin[2] = p->zdot[i];
            va = math_norm(vin);

            /* Evaluate density and temperature */
            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];
            if(!errflag) {errflag = plasma_eval_densandtemp(p->rho[i], pdata, dens, temp);}
            for(int j = 0; j < n_species; j++) {
                temp[j] = temp[j]*CONST_KB;
            }

            /* Evaluate coefficients */
            real clogab[MAX_SPECIES];
            real Kb[MAX_SPECIES];
            real nub[MAX_SPECIES];
            real Fb[MAX_SPECIES],     F     = 0;
            real Dparab[MAX_SPECIES], Dpara = 0;
            real Dperpb[MAX_SPECIES], Dperp = 0;
            if(!errflag) {
                errflag = mccc_coefs_clog(p->mass[i], p->charge[i], va, 
                                          m_species, q_species, dens, temp, clogab, n_species);
            }
            if(!errflag) {
                errflag = mccc_coefs_fo(p->mass[i], p->charge[i], va,
                                        m_species, q_species, dens, temp, clogab, n_species, coldata,
                Fb, Dparab, Dperpb, Kb, nub);
            }
            for(int j = 0; j < n_species; j=j+1){
                F = F + Fb[j];
                Dpara = Dpara + Dparab[j];
                Dperp = Dperp + Dperpb[j];
            }

            /* Evaluate collisions */
            if(!errflag) {errflag = mccc_push_foEM(F, Dpara, Dperp, h[i], &rnd[i*3], vin, vout);}

            /* Update particle */
            real vnorm = va/(math_norm(vout));
#if A5_CCOL_NOENERGY
            vout[0] *= vnorm;
            vout[1] *= vnorm;
            vout[2] *= vnorm;
#endif
#if A5_CCOL_NOPITCH
            vout[0] = vin[0] * vnorm;
            vout[1] = vin[1] * vnorm;
            vout[2] = vin[2] * vnorm;
#endif

            if(!errflag) {
                p->rdot[i] = vout[0] * cos(p->phi[i]) + vout[1] * sin(p->phi[i]);
                p->phidot[i] = (-vout[0] * sin(p->phi[i]) + vout[1] * cos(p->phi[i]) ) / p->r[i];
                p->zdot[i] = vout[2];
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}

/**
 * @brief Evaluates collisions in gc picture with fixed time-step
 *
 * This function first evaluates collision coefficients (see mccc_update_gc)
 * and then evaluates collisions using Euler-Maruyama method and updates
 * marker state.
 *
 * Operations are performed simultaneously for NSIMD markers.
 *
 * @param p pointer to SIMD_gc struct containing the markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param h pointer to time step values [s]
 */
void mccc_step_gc_fixed(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, random_data* rdata, real* coldata, real* h){
    int i;
    real rnd[5*NSIMD];
    random_normal_simd(rdata, 5*NSIMD, rnd);

    int n_species = plasma_get_n_species(pdata);
    real* q_species = plasma_get_species_charge(pdata);
    real* m_species = plasma_get_species_mass(pdata);

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            /* Gather relevant gc information */
            real Brpz[3] = {p->B_r[i], p->B_phi[i], p->B_z[i]};
            real B[3];
            math_vec_rpz2xyz(Brpz, B, p->phi[i]);
            real vin, xiin, Xin[3];
            real vout, xiout, Xout[3];
            real Bnorm = math_norm(B);
            vin = physlib_gc_v(p->mass[i], p->mu[i], p->vpar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], p->vpar[i], Bnorm);

            real R0   = p->r[i];
            real z0   = p->z[i];
            Xin[0] = p->r[i]*cos(p->phi[i]);
            Xin[1] = p->r[i]*sin(p->phi[i]);
            Xin[2] = p->z[i];

            /* Evaluate density and temperature */
            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];
            if(!errflag) {errflag = plasma_eval_densandtemp(p->rho[i], pdata, dens, temp);}

            for(int j = 0; j < n_species; j++) {
                temp[j] = temp[j]*CONST_KB;
            }
            real cutoff = 0.1*sqrt(temp[0]/p->mass[i]);

            /* Evaluate coefficients */
            real clogab[MAX_SPECIES];
            real Dparab[MAX_SPECIES], Dpara = 0;
            real Kb[MAX_SPECIES],     K     = 0;
            real nub[MAX_SPECIES],    nu    = 0;
            real DXb[MAX_SPECIES],    DX    = 0;
            if(!errflag) {
                errflag = mccc_coefs_clog(p->mass[i], p->charge[i], vin,
                        m_species, q_species, dens, temp, clogab, n_species);
            }
            if(!errflag) {
                errflag = mccc_coefs_gcfixed(p->mass[i], p->charge[i], vin, xiin,
                            m_species, q_species, dens, temp, Bnorm, clogab, n_species, coldata,
                            Dparab,DXb,Kb,nub);
            }

            for(int j = 0; j < n_species; j=j+1){
                Dpara = Dpara + Dparab[j];
                K = K + Kb[j];
                nu = nu + nub[j];
                DX = DX + DXb[j];
            }

            /* Evaluate collisions */
            if(!errflag) {
                errflag = mccc_push_gcEM(K,nu,Dpara,DX,B,h[i],&rnd[i*5],
                        vin,&vout,xiin,&xiout,Xin,Xout,cutoff);
            }

            /* Update particle */
#if A5_CCOL_NOENERGY
            vout = vin;
#endif
#if A5_CCOL_NOPITCH
            xiout = xiin;
#endif
#if A5_CCOL_NOGCDIFF
            Xout[0] = Xin[0];
            Xout[1] = Xin[1];
            Xout[2] = Xin[2];
#endif

            if(!errflag) {
                p->r[i] = sqrt(Xout[0]*Xout[0] + Xout[1]*Xout[1]);
                p->z[i] = Xout[2];
                /* Evaluate phi and pol angles so that they are cumulative */
                real axis_r = B_field_get_axis_r(Bdata, p->phi[i]);
                real axis_z = B_field_get_axis_z(Bdata, p->phi[i]);
                p->pol[i] += atan2( (R0-axis_r) * (p->z[i]-axis_z) - (z0-axis_z) * (p->r[i]-axis_r), 
                        (R0-axis_r) * (p->r[i]-axis_r) + (z0-axis_z) * (p->z[i]-axis_z) );
                p->phi[i] += atan2( Xin[0] * Xout[1] - Xin[1] * Xout[0],
                    Xin[0] * Xout[0] + Xin[1] * Xout[1] );
            }

            /* Evaluate magnetic field (and gradient) and rho at new position */
            real B_dB[12], psi[1], rho[1];
            if(!errflag) {errflag = B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i], Bdata);}
            if(!errflag) {errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);}
            if(!errflag) {errflag = B_field_eval_rho(rho, psi[0], Bdata);}

            if(!errflag) {
                p->B_r[i]        = B_dB[0];
                p->B_r_dr[i]     = B_dB[1];
                p->B_r_dphi[i]   = B_dB[2];
                p->B_r_dz[i]     = B_dB[3];

                p->B_phi[i]      = B_dB[4];
                p->B_phi_dr[i]   = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i]   = B_dB[7];

                p->B_z[i]        = B_dB[8];
                p->B_z_dr[i]     = B_dB[9];
                p->B_z_dphi[i]   = B_dB[10];
                p->B_z_dz[i]     = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
                p->vpar[i] = physlib_gc_vpar(vout, xiout);
                p->mu[i]   = physlib_gc_mu(p->mass[i], vout, xiout, Bnorm);
            }

            /* Error handling */
            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
            }
        }
    }
}

/**
 * @brief Evaluates collisions in gc picture with adaptive time-step
 *
 * This function first evaluates collision coefficients (see mccc_update_gc)
 * and then evaluates collisions using Milstein method and updates
 * marker state irrespective whether time-step was accepted. Returns suggestion
 * for the next time-step which has minus sign if the time-step was rejected.
 *
 * Operations are performed simultaneously for NSIMD markers.
 *
 * @param p pointer to SIMD_gc struct containing the markers
 * @param Bdata pointer to magnetic field data
 * @param pdata pointer to plasma data
 * @param hin pointer to values for current time-step [s]
 * @param hout pointer to values for next time-step [s]
 * @param w NSIMD array of pointers to Wiener structs
 * @param tol integration error tolerance
 *
 * @todo There is no need to check for error when collision frequency is low and other processes determine adaptive time-step
 */
void mccc_step_gc_adaptive(particle_simd_gc* p, B_field_data* Bdata, plasma_data* pdata, random_data* rdata, real* coldata, real* hin, real* hout, mccc_wienarr** w, real tol){
    int i;
    real rand5[5*NSIMD];
    random_normal_simd(rdata, 5*NSIMD, rand5);

    int n_species = plasma_get_n_species(pdata);
    real* q_species = plasma_get_species_charge(pdata);
    real* m_species = plasma_get_species_mass(pdata);

    /* Error estimates */
    real kappa_k[NSIMD], kappa_d0[NSIMD], kappa_d1[NSIMD];
    //real dWopt0[NSIMD], dWopt1[NSIMD], alpha[NSIMD]; Needed only if the accurate error-checking is used
    int tindex[NSIMD];

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            /* Gather relevant gc information */
            real Brpz[3] = {p->B_r[i], p->B_phi[i], p->B_z[i]};
            real B[3];
            math_vec_rpz2xyz(Brpz, B, p->phi[i]);
            real vin, xiin, Xin[3];
            real vout, xiout, Xout[3];
            real Bnorm = math_norm(B);
            vin = physlib_gc_v(p->mass[i], p->mu[i], p->vpar[i], Bnorm);
            xiin = physlib_gc_xi(p->mass[i], p->mu[i], p->vpar[i], Bnorm);

            real R0   = p->r[i];
            real z0   = p->z[i];
            Xin[0] = p->r[i]*cos(p->phi[i]);
            Xin[1] = p->r[i]*sin(p->phi[i]);
            Xin[2] = p->z[i];

            real dW[5];
            if(!errflag) {errflag = mccc_wiener_generate(w[i], w[i]->time[0]+hin[i], &tindex[i], &rand5[i*MCCC_NDIM]);}
            if(!errflag) {
                dW[0] = w[i]->wiener[tindex[i]*MCCC_NDIM + 0] - w[i]->wiener[0];
                dW[1] = w[i]->wiener[tindex[i]*MCCC_NDIM + 1] - w[i]->wiener[1];
                dW[2] = w[i]->wiener[tindex[i]*MCCC_NDIM + 2] - w[i]->wiener[2];
                dW[3] = w[i]->wiener[tindex[i]*MCCC_NDIM + 3] - w[i]->wiener[3];
                dW[4] = w[i]->wiener[tindex[i]*MCCC_NDIM + 4] - w[i]->wiener[4];
            }

            /* Evaluate density and temperature */
            real temp[MAX_SPECIES];
            real dens[MAX_SPECIES];
            if(!errflag) {errflag = plasma_eval_densandtemp(p->rho[i], pdata, dens, temp);}
            for(int j = 0; j < n_species; j++) {
                temp[j] = temp[j]*CONST_KB;
            }
            real cutoff = 0.1*sqrt(temp[0]/p->mass[i]);

            /* Evaluate coefficients */
            real clogab[MAX_SPECIES];
            real dQb[MAX_SPECIES],     dQ     = 0;
            real dDparab[MAX_SPECIES], dDpara = 0;
            real Dparab[MAX_SPECIES],  Dpara  = 0;
            real Kb[MAX_SPECIES],      K      = 0;
            real nub[MAX_SPECIES],     nu     = 0;
            real DXb[MAX_SPECIES],     DX     = 0;
            if(!errflag) {
                errflag = mccc_coefs_clog(p->mass[i], p->charge[i], vin, m_species, q_species, 
                        dens, temp, clogab, n_species);
            }
            if(!errflag) {
                errflag = mccc_coefs_gcadaptive(p->mass[i], p->charge[i], vin, xiin, m_species, q_species,
                            dens, temp, Bnorm, clogab, n_species, coldata,
                            Dparab, DXb, Kb, nub, dQb, dDparab);
            }
            for(int j = 0; j < n_species; j=j+1){
                Dpara = Dpara + Dparab[j];
                K = K + Kb[j];
                nu = nu + nub[j];
                DX = DX + DXb[j];
                dQ = dQ + dQb[j];
                dDpara = dDpara + dDparab[j];
            }

            /* Evaluate collisions */
            if(!errflag) {
                errflag = mccc_push_gcMI(K,nu,Dpara,DX,B,hin[i],dW,dQ,dDpara,
                        vin,&vout,xiin,&xiout,Xin,Xout,cutoff,tol,
                        &kappa_k[i], &kappa_d0[i], &kappa_d1[i]);
            }

            /* Needed for finding the next time step (for the old method, see end of function)*/
            /*dWopt0[i] = 0.9*fabs(dW[3])*pow(kappa_d0[i],-1.0/3);
            dWopt1[i] = 0.9*fabs(dW[4])*pow(kappa_d1[i],-1.0/3);
            alpha[i]  = fabs(dW[3]);
            if(alpha[i] < fabs(dW[4])) {
            alpha[i] = fabs(dW[4]);
            }
            alpha[i] = alpha[i]/sqrt(hin[i]);*/

            /* Update particle */
#if A5_CCOL_NOENERGY
            vout = vin;
#endif
#if A5_CCOL_NOPITCH
            xiout = xiin;
#endif
#if A5_CCOL_NOGCDIFF
            Xout[0] = Xin[0];
            Xout[1] = Xin[1];
            Xout[2] = Xin[2];
#endif

            if(!errflag) {
                p->r[i] = sqrt(Xout[0]*Xout[0] + Xout[1]*Xout[1]);
                p->z[i] = Xout[2];

                /* Evaluate phi and pol angles so that they are cumulative */
                real axis_r = B_field_get_axis_r(Bdata, p->phi[i]);
                real axis_z = B_field_get_axis_z(Bdata, p->phi[i]);
                p->pol[i] += atan2( (R0-axis_r) * (p->z[i]-axis_z) - (z0-axis_z) * (p->r[i]-axis_r), 
                        (R0-axis_r) * (p->r[i]-axis_r) + (z0-axis_z) * (p->z[i]-axis_z) );
                p->phi[i] += atan2( Xin[0] * Xout[1] - Xin[1] * Xout[0],
                        Xin[0] * Xout[0] + Xin[1] * Xout[1] );
            }

            /* Evaluate magnetic field (and gradient) at new position */
            real B_dB[12], psi[1], rho[1];
            if(!errflag) {errflag = B_field_eval_B_dB(B_dB, p->r[i], p->phi[i], p->z[i], Bdata);}
            if(!errflag) {errflag = B_field_eval_psi(psi, p->r[i], p->phi[i], p->z[i], Bdata);}
            if(!errflag) {errflag = B_field_eval_rho(rho, psi[0], Bdata);}

            if(!errflag) {
                p->B_r[i]        = B_dB[0];
                p->B_r_dr[i]     = B_dB[1];
                p->B_r_dphi[i]   = B_dB[2];
                p->B_r_dz[i]     = B_dB[3];

                p->B_phi[i]      = B_dB[4];
                p->B_phi_dr[i]   = B_dB[5];
                p->B_phi_dphi[i] = B_dB[6];
                p->B_phi_dz[i]   = B_dB[7];

                p->B_z[i]        = B_dB[8];
                p->B_z_dr[i]     = B_dB[9];
                p->B_z_dphi[i]   = B_dB[10];
                p->B_z_dz[i]     = B_dB[11];

                p->rho[i] = rho[0];

                Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
                p->vpar[i] = physlib_gc_vpar(vout, xiout);
                p->mu[i]   = physlib_gc_mu(p->mass[i], vout, xiout, Bnorm);

                /* Test whether timestep was rejected and suggest next time step */
                int rejected = 0;
                if(kappa_k[i] > 1 || kappa_d0[i] > 1 || kappa_d1[i] > 1) {
                    rejected = 1;
                    tindex[i]=0;
                }

                if(kappa_k[i] >= kappa_d0[i] && kappa_k[i] >= kappa_d1[i]) {
                    hout[i] = 0.8*hin[i]/sqrt(kappa_k[i]);
                }
                else if(kappa_d0[i] >= kappa_k[i] && kappa_d0[i] >= kappa_d1[i]) {
                    hout[i] = fabs(hin[i])*pow(0.9*pow(kappa_d0[i],-1.0/3.0),2.0);
                }
                else {
                    hout[i] = fabs(hin[i])*pow(0.9*pow(kappa_d1[i],-1.0/3.0),2.0);
                }

                /* Negative value indicates time step was rejected*/
                if(rejected){
                    hout[i] = -hout[i];
                }
                else if(hout[i] > 1.5*hin[i]) {
                    /* Make sure we don't increase time step too much */
                    hout[i] = 1.5*hin[i];
                }
            }

            /* Error handling */
            if(!errflag && (fabs(hout[i]) < A5_EXTREMELY_SMALL_TIMESTEP))      {errflag = error_raise(ERR_INVALID_TIMESTEP, __LINE__, EF_MCCC);}

            if(errflag) {
                p->err[i]     = errflag;
                p->running[i] = 0;
                hout[i]       = hin[i];
            }
        }
    }

    // This is not currently defined anywhere...
#ifdef MCCC_MORE_ACCURATE_DT_GUESS
    // More accurate but probably less efficient method below
    /* Choose next time step (This loop can be vectorized if there is a
       suitable tool for drawing random numbers) */
    for(i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            real t = w[i]->time[0];

            /* Check whether time step was accepted and find value for the next time-step */
            int rejected = 0;
            if(kappa_k[i] > 1 || kappa_d0[i] > 1 || kappa_d1[i] > 1) {
                rejected = 1;
                tindex[i]=0;
            }

            /* Different time step estimates are used depending which error estimate dominates
             * This scheme automatically takes care of time step reduction (increase) when 
             * time step is rejected (accepted) */
            int ki, kmax;
            int windex;
            real dW[2];

            if(kappa_k[i] > kappa_d0[i] || kappa_k[i] > kappa_d1[i]) {
                real dti = 0.8*hin[i]/sqrt(kappa_k[i]);
                if(1.5*hin[i] < dti){dti = 1.5*hin[i];}
                kmax = 4;
                for(ki=1; ki < kmax; ki=ki+1){
                    random_normal_simd(rdata, MCCC_NDIM, &rand5[i*MCCC_NDIM]);
                    mccc_wiener_generate(w[i], t+ki*dti/3, &windex, &rand5[i*MCCC_NDIM]);
                    dW[0] = fabs(w[i]->wiener[3 + windex*MCCC_NDIM]
                                 - w[i]->wiener[3 + tindex[i]*MCCC_NDIM]);
                    if(dW[0] > dWopt0[i]) {
                        kmax = 0; // Exit loop
                    }
                    else {
                        dW[1] = fabs(w[i]->wiener[4 + windex*MCCC_NDIM]
                                     - w[i]->wiener[4 + tindex[i]*MCCC_NDIM]);
                        if(dW[1] > dWopt1[i]) {
                            kmax = 0; // Exit loop
                        }
                    }
                }
                if(ki == 1){
                    hout[i] = (dti/3);
                }
                else{
                    hout[i] = (ki-1)*(dti/3);
                }
            }
            else{
                kmax = 6;
                if (rejected) {
                    kmax = 2;
                }
                else if (alpha[i] > 2) {
                    kmax = 4;
                }

                for(ki=1; ki < kmax; ki=ki+1){
                    random_normal_simd(rdata, MCCC_NDIM, &rand5[i*MCCC_NDIM]);
                    mccc_wiener_generate(w[i], t+ki*hin[i]/3, &windex, &rand5[i*MCCC_NDIM]);
                    dW[0] = abs(w[i]->wiener[3 + windex*MCCC_NDIM] - w[i]->wiener[3 + tindex[i]*MCCC_NDIM]);
                    if(dW[0] > dWopt0[i]) {
                        kmax = 0; // Exit loop
                    }
                    else{
                        dW[1] = abs(w[i]->wiener[4 + windex*MCCC_NDIM] - w[i]->wiener[4 + tindex[i]*MCCC_NDIM]);
                        if(dW[1] > dWopt1[i]) {
                            kmax = 0; // Exit loop
                        }
                    }
                }
                if(ki == 1){
                    hout[i] = (hin[i]/3);
                }
                else{
                    hout[i] = (ki-1)*(hin[i]/3);
                }
            }

            /* Negative value indicates time step was rejected*/
            if(rejected){
                hout[i] = -hout[i];
            }
        }
    }
#endif
}
