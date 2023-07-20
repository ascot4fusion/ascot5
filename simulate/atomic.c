/**
 * @file atomic.c
 * @brief Atomic reactions model
 *
 * The atomic module models atomic reactions for the simulated test particles.
 * Using atomic reaction data supplied by the asigma (atomicsigma) helper
 * module, reaction probabilities are calculated and tested against a random
 * number. If a reaction occurs, the particle charge state is changed.
 *
 * The terms ionization and recombination are used loosely in variable names.
 * Ionization (ion) stands for all charge-increasing reactions, and
 * recombination (rec) for all charge-decreasing reactions.
 */
#include <stdlib.h>
#include <stdio.h>
#include "../ascot5.h"
#include "../math.h"
#include "../physlib.h"
#include "../consts.h"
#include "../error.h"
#include "../print.h"
#include "../particle.h"
#include "../plasma.h"
#include "../random.h"
#include "../asigma.h"
#include "atomic.h"

#pragma omp declare target
#pragma omp declare simd uniform(asgm_data)
a5err atomic_rates(
    real* rate_eff_ion, real* rate_eff_rec, int z_1, int a_1, real m_1,
    const int* z_2, const int* a_2, const real* m_2, asigma_data* asgm_data,
    int q, real E, int N_pls_spec, int N_ntl_spec, real* T, real* T_0,
    real* n, real* n_0, int* enable_atomic);
#pragma omp declare simd
a5err atomic_react(
    int* q, real dt, real rate_eff_ion, real rate_eff_rec, int z_1, real rnd);
#pragma omp end declare target

/**
 * @brief Determine if atomic reactions occur during time-step and change charge
 *
 * The terms ionization and recombination are used loosely in variable names.
 * Ionization (ion) stands for all charge-increasing reactions, and
 * recombination (rec) for all charge-decreasing reactions.
 *
 * @param p fo struct
 * @param h time-steps from NSIMD markers
 * @param p_data pointer to plasma data
 * @param n_data pointer to neutral data
 * @param r_data pointer to random-generator data
 * @param asgm_data pointer to atomic reaction data
 * @param enable_atomic pointer to atomic enable and functionality flag
 */
void atomic_fo(particle_simd_fo* p, real* h,
               plasma_data* p_data, neutral_data* n_data,
               random_data* r_data, asigma_data* asigmadata,
               int* enable_atomic) {

    /* Generate random numbers and get plasma information before going to the *
     * SIMD loop                                                              */
    real rnd[NSIMD];
    random_uniform_simd(r_data, NSIMD, rnd);

    int N_pls_spec  = plasma_get_n_species(p_data);
    int N_ntl_spec  = neutral_get_n_species(n_data);
    const real* m_2 = plasma_get_species_mass(p_data);
    const int* z_2  = plasma_get_species_znum(p_data);
    const int* a_2  = plasma_get_species_anum(p_data);

    #pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p->running[i]) {
            a5err errflag = 0;

            /* Calculate kinetic energy of test particle */
            real p_comps[3];
            p_comps[0] = p->p_r[i];
            p_comps[1] = p->p_phi[i];
            p_comps[2] = p->p_z[i];
            real p_norm   = math_norm(p_comps);
            real E = physlib_Ekin_pnorm(p->mass[i], p_norm);

            /* Evaluate plasma density and temperature */
            real n_2[MAX_SPECIES], T_2[MAX_SPECIES];
            if(!errflag) {
                errflag = plasma_eval_densandtemp(n_2, T_2, p->rho[i],
                                                  p->r[i], p->phi[i], p->z[i],
                                                  p->time[i], p_data);
            }

            /* Evaluate neutral density and temperature */
            real n_0[MAX_SPECIES], T_0[MAX_SPECIES];
            if(!errflag) {
                errflag = neutral_eval_n0(n_0, p->rho[i],
                                          p->r[i], p->phi[i], p->z[i],
                                          p->time[i], n_data);
            }
            if(!errflag) {
                errflag = neutral_eval_t0(T_0, p->rho[i],
                                          p->r[i], p->phi[i], p->z[i],
                                          p->time[i], n_data);
            }

            /* Evaluate the reaction rates for ionizing (charge-increasing) *
               and recombining (charge-decreasing) reactions                */
            int q = (int)round(p->charge[i]/CONST_E);
            real rate_eff_ion, rate_eff_rec;
            if(!errflag) {
                errflag = atomic_rates(&rate_eff_ion, &rate_eff_rec,
                                       p->znum[i], p->anum[i], p->mass[i],
                                       z_2, a_2, m_2,
                                       asigmadata,
                                       q, E,
                                       N_pls_spec, N_ntl_spec,
                                       T_2, T_0,
                                       n_2, n_0,
                                       enable_atomic);
            }

            /* Determine if an atomic reaction occurs */
            if(!errflag) {
                int q_prev = q;
                errflag = atomic_react(&q, h[i],
                                       rate_eff_ion, rate_eff_rec,
                                       p->znum[i], rnd[i]);
                if(q != q_prev) {
                    /* A reaction has occured, change particle charge */
                    p->charge[i] = q*CONST_E;
                }
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
 * @brief Determines atomic reaction rates
 *
 * This function checks the particle charge state, and determines the
 * reaction rates for the different possible atomic reactions.
 *
 * The terms ionization and recombination are used loosely in variable names.
 * Ionization (ion) stands for all charge-increasing reactions, and
 * recombination (rec) for all charge-decreasing reactions.
 *
 * @param rate_eff_ion pointer to evaluated reaction rate for ionization
 * @param rate_eff_rec pointer to evaluated reaction rate for recombination
 * @param z_1 atomic number of fast particle
 * @param a_1 atomic mass number of fast particle
 * @param m_1 mass of fast particle
 * @param z_2 atomic numbers of bulk particle species
 * @param a_2 atomic mass numbers of bulk particle species
 * @param m_2 masses of bulk particle species
 * @param asigmadata pointer to atomic data struct
 * @param q charge of fast particle
 * @param E energy of fast particle
 * @param N_pls_spec number of species in bulk plasma
 * @param N_ntl_spec number of species in bulk neutrals
 * @param T temperature of bulk plasma species
 * @param T_0 temperature of bulk neutral species
 * @param n density of bulk plasma species
 * @param n_0 density of bulk neutral species
 * @param enable_atomic pointer to atomic enable and functionality flag
 *
 * @return zero if evaluation succeeded
 *
 * @todo Implement a more general algorithm for the determination of reaction
 *   rates based on charge state!
 */
a5err atomic_rates(
    real* rate_eff_ion, real* rate_eff_rec, int z_1, int a_1, real m_1,
    const int* z_2, const int* a_2, const real* m_2, asigma_data* asigmadata,
    int q, real E, int N_pls_spec, int N_ntl_spec, real* T, real* T_0,
    real* n, real* n_0, int* enable_atomic) {
    a5err err = 0;

    /* Define a helper variable for storing rate coefficients, and
       initialize the reaction rate variables */
    real sigmav;
    *rate_eff_rec = 0;
    *rate_eff_ion = 0;

    /* Based on particle charge state, reaction rates are calculated for
       possible reactions using the asigma helper module. Typically, asigma
       does not return a reaction rate. Instead, it might for example return
       a rate coefficient. A rate coefficient must be multiplied by the
       density of the reaction counterpart to obtain the reaction rate.
       NOTE: Currently, only two charge states and a limited range of
       reaction types are implemented. */
    if(q == 1) {
        /* Test particle has charge state +1. Evaluate the rate coefficient
           and multiply it by the reaction counterpart density to find the
           reaction rate for a recombining (charge-decreasing) reaction. */
        /* Charge exchange (CX)
           NOTE: Ion density is not used in asigma_eval_sigmav() for
           reac_type_sigmav_CX, so just pass density of main bulk ion species
           as the dummy parameter. */
        for(int i_spec = 0; i_spec < N_ntl_spec; i_spec++) {
            err = asigma_eval_sigmav(&sigmav,
                                     z_1, a_1, m_1,
                                     z_2[i_spec], a_2[i_spec],
                                     reac_type_sigmav_CX,
                                     asigmadata,
                                     E,
                                     T[0], T_0[i_spec], n[1],
                                     enable_atomic);
            *rate_eff_rec += sigmav*n_0[i_spec];
        }
    } else if(q == 0) {
        /* Test particle is neutral. Evaluate the beam-stopping (BMS)
           coefficient. Because of how BMS coefficients are defined, to
           determine the reaction rate, the BMS coefficient returned by
           asigma_eval_sigmav() is multiplied by the electron density
           corresponding to the ionic reaction counterpart density. The
           species index maximum excludes electrons.
           NOTE: Fully ionized ions are assumed when multiplying the effective
           rate coefficient by the electron density corresponding to the ionic
           reaction counterpart density.
           NOTE: Neutral temperature is not used in asigma_eval_sigmav() for
           reac_type_BMS_sigmav, so just pass temperature of main bulk neutral
           species as the dummy parameter. */
        for(int i_spec = 0; i_spec < (N_pls_spec-1); i_spec++) {
            err = asigma_eval_sigmav(&sigmav,
                                     z_1, a_1, m_1,
                                     z_2[i_spec], a_2[i_spec],
                                     reac_type_BMS_sigmav,
                                     asigmadata,
                                     E,
                                     T[0], T_0[0], n[i_spec+1],
                                     enable_atomic);
            *rate_eff_ion += sigmav*z_2[i_spec]*n[i_spec+1];
        }
    } else {
        print_err("Warning: Test particle has a charge state for which "
                  "determination of atomic reaction rates has not "
                  "been implemented.\n");
        err = error_raise( ERR_ATOMIC_EVALUATION, __LINE__, EF_ATOMIC );
    }

    return err;
}

/**
 * @brief Determines if an atomic reaction occurs during one time step
 *
 * The terms ionization and recombination are used loosely in variable names.
 * Ionization (ion) stands for all charge-increasing reactions, and
 * recombination (rec) for all charge-decreasing reactions.
 *
 * @param q charge of fast particle
 * @param dt time-step
 * @param rate_eff_ion reaction rate for ionization
 * @param rate_eff_rec reaction rate for recombination
 * @param z_1 atomic number of fast particle
 * @param rnd random number
 *
 * @return zero if evaluation succeeded
 */
a5err atomic_react(
    int* q, real dt, real rate_eff_ion, real rate_eff_rec, int z_1, real rnd) {
    a5err err = 0;

    /* Calculate the reaction probabilities for ionizing (charge-increasing)
       and recombining (charge-decreasing) atomic reactions. But first,
       a fail-safe against zero rate values. */
    real prob_eff_ion;
    real prob_eff_rec;
    if(rate_eff_ion == 0.0) {
        prob_eff_ion = 0.0;
    } else {
        prob_eff_ion = rate_eff_ion/(rate_eff_ion+rate_eff_rec)
                       *(1.0-exp(-(rate_eff_ion+rate_eff_rec)*dt));
    }
    if(rate_eff_rec == 0.0) {
        prob_eff_rec = 0.0;
    } else {
        prob_eff_rec = rate_eff_rec/(rate_eff_ion+rate_eff_rec)
                       *(1.0-exp(-(rate_eff_ion+rate_eff_rec)*dt));
    }

    /* Check if reaction occurs, and update charge state accordingly */
    if(*q < z_1 && rnd < prob_eff_ion) {
        *q += 1;
    } else if(*q > 0 && rnd < (prob_eff_ion+prob_eff_rec)) {
        *q -= 1;
    }

    return err;
}
