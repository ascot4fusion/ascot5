/**
 * Implements "atomic.h".
 */
#include "atomic_reactions.h"
#include "consts.h"
#include "data/atomic.h"
#include "data/marker.h"
#include "data/neutral.h"
#include "data/plasma.h"
#include "defines.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "utils/random.h"
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

void atomic_reactions_eval_go(
    MarkerGyroOrbit *mrk, real *h, Plasma *plasma, Neutral *neutral,
    Atomic *atomic)
{
    int N_pls_spec = Plasma_get_n_species(plasma);
    int N_ntl_spec = Neutral_get_n_species(neutral);
    const int *z_2 = Plasma_get_species_znum(plasma);
    const int *a_2 = Plasma_get_species_anum(plasma);

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < mrk->size; i++)
    {
        if (mrk->running[i])
        {
            err_t errflag = 0;

            /* Calculate kinetic energy of test particle */
            real p_comps[3];
            p_comps[0] = mrk->p_r[i];
            p_comps[1] = mrk->p_phi[i];
            p_comps[2] = mrk->p_z[i];
            real p_norm = math_norm(p_comps);
            real E = physlib_Ekin_pnorm(mrk->mass, p_norm);

            real n_2[MAX_SPECIES], T_2[MAX_SPECIES];
            if (!errflag)
            {
                errflag = Plasma_eval_nT(
                    n_2, T_2, mrk->rho[i], mrk->r[i], mrk->phi[i], mrk->z[i],
                    mrk->time[i], plasma);
            }
            real n_0[MAX_SPECIES], T_0[MAX_SPECIES];
            if (!errflag)
            {
                errflag = Neutral_eval_nT(
                    n_0, T_0, mrk->rho[i], mrk->r[i], mrk->phi[i], mrk->z[i],
                    mrk->time[i], neutral);
            }

            /* Evaluate the reaction rates for ionizing (charge-increasing) *
               and recombining (charge-decreasing) reactions                */
            real rate_eff_ion, rate_eff_rec;
            if (!errflag)
            {
                errflag = atomic_reactions_rates(
                    &rate_eff_ion, &rate_eff_rec, mrk->znum, mrk->anum,
                    mrk->mass, z_2, a_2, atomic, mrk->charge[i], E, N_pls_spec,
                    N_ntl_spec, T_2, T_0, n_2, n_0);
            }

            /* Determine if an atomic reaction occurs */
            real rnd[2];
            random_uniform(rnd, 2, mrk->id, mrk->rng_count, i);
            if (!errflag)
            {
                atomic_reactions_apply(
                    &mrk->charge[i], h[i], rate_eff_ion, rate_eff_rec,
                    mrk->znum, rnd);
            }

            /* Error handling */
            if (errflag)
            {
                mrk->err[i] = errflag;
                mrk->running[i] = 0;
            }
        }
    }
}

err_t atomic_reactions_rates(
    real *rate_eff_ion, real *rate_eff_rec, int z_1, int a_1, real m_1,
    const int *z_2, const int *a_2, Atomic *atomic, int q, real E,
    int N_pls_spec, int N_ntl_spec, real *T, real *T_0, real *n, real *n_0)
{
    err_t err = 0;
    real rate;
    *rate_eff_rec = 0;
    *rate_eff_ion = 0;

    // ERROR_CHECK(err, condition, type, file)
    Atomic_eval_bms(
        &rate, E, m_1, (N_pls_spec - 1), z_2, a_2, T[0], n[0], &n[1], 0,
        atomic);
    *rate_eff_ion += rate * n[0];
    err = Atomic_eval_cx(
        &rate, q, E, N_ntl_spec, z_2, a_2, T_0[0], n_0, 0, atomic);
    *rate_eff_rec += rate;
    err = Atomic_eval_ei(&rate, z_1, a_1, q, T[0], n[0], 0, atomic);
    *rate_eff_ion += rate;

    return err;
}

void atomic_reactions_apply(
    int *charge, int znum, real dt, real rate_ionization,
    real rate_recombination, real rnd[2])
{
    /* Calculate the reaction probabilities for ionizing (charge-increasing)
       and recombining (charge-decreasing) atomic reactions. But first,
       a fail-safe against zero rate values. */
    int no_reaction = rate_ionization == 0.0 && rate_recombination == 0.0;
    int ionization = (rate_ionization == 0.0 || *charge == znum) && *charge > 0;
    int recombination = (rate_recombination == 0.0 || *charge == 0);
    int eitherway = !no_reaction && !ionization && !recombination;

    int success = 1.0 - ionization * exp(-rate_recombination * dt) -
                  recombination * exp(-rate_ionization * dt) -
                  eitherway * exp(-(rate_ionization + rate_recombination) * dt);

    *charge = *charge +
              success * (ionization - recombination) * (1 - eitherway) +
              success * eitherway;
    *charge =
        *charge -
        2 * success * eitherway *
            (rnd[1] < rate_recombination /
                          (rate_ionization + rate_recombination + DBL_MIN));
}
