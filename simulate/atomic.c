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
    int q, real E, int N_pls_spec, real* T, real T_0, real* n, real n_0);
#pragma omp declare simd
a5err atomic_react(
    int* q, real dt, real rate_eff_ion, real rate_eff_rec, int z_1, real rnd);
#pragma omp end declare target

/**
 * @brief Determine if atomic reactions occur during one time-step
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
 * @param asigma_data pointer to atomic reaction data
 */
void atomic_fo(particle_simd_fo* p, real* h,
               plasma_data* p_data, neutral_data* n_data,
               random_data* r_data, asigma_data* asigmadata) {

    /* Generate random numbers and get plasma information before going to the *
     * SIMD loop                                                              */
    real rnd[NSIMD];
    random_uniform_simd(r_data, NSIMD, rnd);

    int N_pls_spec  = plasma_get_n_species(p_data);
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
            real E = p_norm*p_norm/(2*p->mass[i]);

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
                                       N_pls_spec,
                                       T_2, T_0[0],
                                       n_2, n_0[0]);
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
 * @param N_pls_spec of species in bulk plasma
 * @param T temperature of bulk plasma species
 * @param T_0 temperature of bulk neutrals
 * @param n density of bulk plasma species
 * @param n_0 density of bulk neutrals
 *
 * @return zero if evaluation succeeded
 *
 * @todo
 * - Implement a more general algorithm for the determination of reaction
 * rates based on charge state!
 * - Implement multiple neutral species (in ASCOT5 as a whole), and
 * account for this in the call to asigma_eval_sigmav()!
 * - Any error raised by the asigma_eval_sigmav() call to determine the rate
 * coefficient for recombinaton to a free electron is overwritten by the
 * asigma_eval_sigmav() call to determine the rate coefficient for CX.
 */
a5err atomic_rates(
    real* rate_eff_ion, real* rate_eff_rec, int z_1, int a_1, real m_1,
    const int* z_2, const int* a_2, const real* m_2, asigma_data* asigmadata,
    int q, real E, int N_pls_spec, real* T, real T_0, real* n, real n_0) {
    a5err err = 0;

    /* Define a helper variable for storing rate coefficients, and
       initialize the reaction rate variables */
    real sigmav;
    *rate_eff_rec = 0;
    *rate_eff_ion = 0;

    /* Based on particle charge state, reaction rates are calculated for
       possible reactions using the asigma helper module. Asigma
       typically does not return a reaction rate. Instead, it might for
       example return a rate coefficient. A rate coefficient must be
       multiplied by the density of the reaction counterpart to obtain
       the reaction rate.
       NOTE: Currently, only two charge states and a limited range of
       reaction types are implemented. */
    if(q == 1) {
        /* Test particle has charge state +1. Several possible recombining
           (charge-decreasing) reactions exist. */
        /* Recombination to free electron */
        err = asigma_eval_sigmav(&sigmav,
                                 z_1, a_1,
                                 z_2[0], a_2[0],
                                 0,//reac_type_sigmav_rec,
                                 asigmadata,
                                 E,
                                 T[0], &T[1], T_0,
                                 n[0], &n[1], 0);
        *rate_eff_rec += sigmav*n[0];
        /* Charge exchange (CX).
           NOTE: The loop over neutral species is currently only a placeholder.
           Note that the upper loop limit is a hardcoded 1! The reason is that
           the current implementation allows for only one neutral species.
           Note that T_0 and n_0 are scalar variables! */
        for(int i_spec = 0; i_spec < 1; i_spec++) {
            err = asigma_eval_sigmav(&sigmav,
                                     z_1, a_1,
                                     z_2[i_spec], a_2[i_spec],
                                     0,//reac_type_sigmav_CX,
                                     asigmadata,
                                     E,
                                     T[0], &T[i_spec+1], T_0,
                                     n[0], &n[i_spec+1], i_spec);
            *rate_eff_rec += sigmav*n_0;
        }
    } else if(q == 0) {
        /* Test particle is neutral. We use beam-stopping (BMS) coefficients
           to calculate the ionization probability. The BMS coefficient
           returned by asigma_eval_sigmav() is multiplied by the electron
           density corresponding to the reaction counterpart density because
           of how BMS coefficients are defined. The species index maximum
           excludes electrons.
           NOTE: The below iteration over ion species and cumulation of
           the reaction rate does not work with the Suzuki model, which
           currently might be called inside asigma_loc_eval_sigmav(...), which
           is called inside the below asigma_eval_sigmav(...). However, the
           best fix is probably to move Suzuki away from
           asigma_loc_eval_sigmav(...), where it never belonged in the first
           place since the "loc" submodule is for when local-file atomic data
           is used.
           NOTE: Fully ionized ions are assumed when multiplying the effective
           rate coefficient by the electron density corresponding to the
           reaction counterpart density. */
        for(int i_spec = 0; i_spec < (N_pls_spec-1); i_spec++) {
            err = asigma_eval_sigmav(&sigmav,
                                     z_1, a_1,
                                     z_2[i_spec], a_2[i_spec],
                                     0,//reac_type_BMS_sigmav,
                                     asigmadata,
                                     E,
                                     T[0], &T[i_spec+1], T_0,
                                     n[0], &n[i_spec+1], i_spec);
            *rate_eff_ion += sigmav*z_2[i_spec]*n[i_spec+1];
        }
    } else {
        print_err("Warning: Test particle has a charge state for which "
                  "determination of atomic reaction rates has not "
                  "been implemented\n");
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
