/**
 * Implements endcond.h.
 */
#include "endcond.h"
#include "consts.h"
#include "math.h"
#include "options.h"
#include "particle.h"
#include "physlib.h"
#include "plasma.h"
#include "simulate.h"
#include <math.h>

const unsigned int endcond_tlim = (1u << 0);   /**< Simulation time limit     */
const unsigned int endcond_emin = (1u << 1);   /**< Minimum energy            */
const unsigned int endcond_therm = (1u << 2);  /**< Thermalized               */
const unsigned int endcond_wall = (1u << 3);   /**< Wall collision            */
const unsigned int endcond_rhomin = (1u << 4); /**< Minimum rho               */
const unsigned int endcond_rhomax = (1u << 5); /**< Maximum rho               */
const unsigned int endcond_polmax = (1u << 6); /**< Poloidal limit            */
const unsigned int endcond_tormax = (1u << 7); /**< Toroidal limit            */
const unsigned int endcond_cpumax = (1u << 8); /**< Wall time exceeded        */
const unsigned int endcond_hybrid = (1u << 9); /**< Hybrid mode condition     */
const unsigned int endcond_neutr = (1u << 10); /**< Neutralized               */
const unsigned int endcond_ioniz = (1u << 11); /**< Ionized                   */

void endcond_check_fo(
    particle_simd_fo *p_f, particle_simd_fo *p_i, sim_data *sim)
{
    sim_parameters *params = sim->params;

    /* Note which end conditions are set as active.
       Only these ones are checked */
    int active_tlim = params->endcond_active & endcond_tlim;
    int active_wall = params->endcond_active & endcond_wall;
    int active_emin = params->endcond_active & endcond_emin;
    int active_therm = params->endcond_active & endcond_therm;
    int active_rhomax = params->endcond_active & endcond_rhomax;
    int active_rhomin = params->endcond_active & endcond_rhomin;
    int active_polmax = params->endcond_active & endcond_polmax;
    int active_tormax = params->endcond_active & endcond_tormax;
    int active_cpumax = params->endcond_active & endcond_cpumax;
    int active_neutr = params->endcond_active & endcond_neutr;
    int active_ioniz = params->endcond_active & endcond_ioniz;

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (int i = 0; i < p_f->n_mrk; i++)
    {
        if (p_f->running[i])
        {

            /* Update bounces if pitch changed sign */
            if ((p_i->p_r[i] * p_i->B_r[i] + p_i->p_phi[i] * p_i->B_phi[i] +
                 p_i->p_z[i] * p_i->B_z[i]) *
                    (p_f->p_r[i] * p_f->B_r[i] + p_f->p_phi[i] * p_f->B_phi[i] +
                     p_f->p_z[i] * p_f->B_z[i]) <
                0)
            {
                if (p_f->bounces[i] > 0)
                {
                    /* Half bounce */
                    p_f->bounces[i] *= -1;
                }
                else if (p_f->bounces[i] < 0)
                {
                    /* Bounce complete */
                    p_f->bounces[i] *= -1;
                    p_f->bounces[i] += 1;
                }
                else
                {
                    /* Initial bounce */
                    p_f->bounces[i] += 1;
                }
            }

            /* Check if the marker time exceeds simulation time */
            if (active_tlim)
            {
                if (!params->reverse_time &&
                    p_f->time[i] > params->lab_time_limit)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if (params->reverse_time &&
                    p_f->time[i] < params->lab_time_limit)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if (p_f->mileage[i] > params->max_mileage)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if (active_wall)
            {
                real w_coll = 0;
                int tile = wall_hit_wall(
                    p_i->r[i], p_i->phi[i], p_i->z[i], p_f->r[i], p_f->phi[i],
                    p_f->z[i], &sim->wall_data, &w_coll);
                if (tile > 0)
                {
                    real w = w_coll;
                    p_f->time[i] =
                        p_i->time[i] + w * (p_f->time[i] - p_i->time[i]);
                    p_f->r[i] = p_i->r[i] + w * (p_f->r[i] - p_i->r[i]);
                    p_f->phi[i] = p_i->phi[i] + w * (p_f->phi[i] - p_i->phi[i]);
                    p_f->z[i] = p_i->z[i] + w * (p_f->z[i] - p_i->z[i]);

                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Evaluate marker energy, and check if it is below the minimum
             * energy limit or local thermal energy limit */
            if (active_emin || active_therm)
            {
                real pnorm =
                    math_normc(p_f->p_r[i], p_f->p_phi[i], p_f->p_z[i]);
                real ekin = physlib_Ekin_pnorm(p_f->mass[i], pnorm);

                real Ti;
                a5err errflag = plasma_eval_temp(
                    &Ti, p_f->rho[i], p_f->r[i], p_f->phi[i], p_f->z[i],
                    p_f->time[i], 1, &sim->plasma_data);

                /* Error handling */
                if (errflag)
                {
                    p_f->err[i] = errflag;
                    p_f->running[i] = 0;
                    Ti = 0;
                }

                if (active_emin && (ekin < params->min_energy))
                {
                    p_f->endcond[i] |= endcond_emin;
                    p_f->running[i] = 0;
                }
                if (active_therm &&
                    (ekin < (params->min_local_thermal_energy * Ti)))
                {
                    p_f->endcond[i] |= endcond_therm;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if (active_rhomax)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
            }
            if (active_rhomin)
            {
                if (p_f->rho[i] < params->rho_coordinate_limits[0])
                {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if (active_tormax)
            {
                if (fabs(p_f->phi[i]) > params->max_number_of_toroidal_orbits)
                {
                    maxorb |= endcond_tormax;
                }
            }
            if (active_polmax)
            {
                if (fabs(p_f->theta[i]) > params->max_number_of_poloidal_orbits)
                {
                    maxorb |= endcond_polmax;
                }
                else if (
                    p_f->bounces[i] - 1 >=
                    (int)(params->max_number_of_poloidal_orbits / CONST_2PI))
                {
                    maxorb |= endcond_polmax;
                }
            }
            if (params->require_both_tor_and_pol && maxorb & endcond_tormax &&
                maxorb & endcond_polmax)
            {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            else if (!params->require_both_tor_and_pol && maxorb)
            {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if (active_cpumax)
            {
                if (p_f->cputime[i] > params->max_real_time)
                {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }
            /* Check if the particle has been neutralized */
            if (active_neutr)
            {
                if (p_i->charge[i] != 0.0 && p_f->charge[i] == 0.0)
                {
                    p_f->endcond[i] |= endcond_neutr;
                    p_f->running[i] = 0;
                }
            }

            /* Check if the particle has been ionized */
            if (active_ioniz)
            {
                if (p_i->charge[i] == 0.0 && p_f->charge[i] != 0.0)
                {
                    p_f->endcond[i] |= endcond_ioniz;
                    p_f->running[i] = 0;
                }
            }

            /* Zero end condition if error happened in this function */
            if (p_f->err[i])
            {
                p_f->endcond[i] = 0;
            }
        }
    }
}

void endcond_check_gc(
    particle_simd_gc *p_f, particle_simd_gc *p_i, sim_data *sim)
{
    sim_parameters *params = sim->params;

    int active_tlim = params->endcond_active & endcond_tlim;
    int active_wall = params->endcond_active & endcond_wall;
    int active_emin = params->endcond_active & endcond_emin;
    int active_therm = params->endcond_active & endcond_therm;
    int active_rhomax = params->endcond_active & endcond_rhomax;
    int active_rhomin = params->endcond_active & endcond_rhomin;
    int active_polmax = params->endcond_active & endcond_polmax;
    int active_tormax = params->endcond_active & endcond_tormax;
    int active_cpumax = params->endcond_active & endcond_cpumax;

#pragma omp simd
    for (int i = 0; i < NSIMD; i++)
    {
        if (p_f->running[i])
        {
            /* Update bounces if pitch changed sign */
            if (p_i->ppar[i] * p_f->ppar[i] < 0)
            {
                if (p_f->bounces[i] > 0)
                {
                    /* Half bounce */
                    p_f->bounces[i] *= -1;
                }
                else if (p_f->bounces[i] < 0)
                {
                    /* Bounce complete */
                    p_f->bounces[i] *= -1;
                    p_f->bounces[i] += 1;
                }
                else
                {
                    /* Initial bounce */
                    p_f->bounces[i] += 1;
                }
            }

            /* Check if the marker time exceeds simulation time */
            if (active_tlim)
            {
                if (!params->reverse_time &&
                    p_f->time[i] > params->lab_time_limit)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if (params->reverse_time &&
                    p_f->time[i] < params->lab_time_limit)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if (p_f->mileage[i] > params->max_mileage)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if (active_wall)
            {
                real w_coll = 0;
                int tile = wall_hit_wall(
                    p_i->r[i], p_i->phi[i], p_i->z[i], p_f->r[i], p_f->phi[i],
                    p_f->z[i], &sim->wall_data, &w_coll);
                if (tile > 0)
                {
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Evaluate marker energy, and check if it is below the minimum
             * energy limit or local thermal energy limit */
            if (active_emin || active_therm)
            {
                real Bnorm =
                    math_normc(p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]);
                real ekin = physlib_Ekin_ppar(
                    p_f->mass[i], p_f->mu[i], p_f->ppar[i], Bnorm);

                real Ti;
                a5err errflag = plasma_eval_temp(
                    &Ti, p_f->rho[i], p_f->r[i], p_f->phi[i], p_f->z[i],
                    p_f->time[i], 1, &sim->plasma_data);

                /* Error handling */
                if (errflag)
                {
                    p_f->err[i] = errflag;
                    p_f->running[i] = 0;
                    Ti = 0;
                }

                if (active_emin && (ekin < params->min_energy))
                {
                    p_f->endcond[i] |= endcond_emin;
                    p_f->running[i] = 0;
                }
                if (active_therm &&
                    (ekin < (params->min_local_thermal_energy * Ti)))
                {
                    p_f->endcond[i] |= endcond_therm;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if (active_rhomax)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
            }
            if (active_rhomin)
            {
                if (p_f->rho[i] < params->rho_coordinate_limits[0])
                {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if (active_tormax)
            {
                if (fabs(p_f->phi[i]) > params->max_number_of_toroidal_orbits)
                {
                    maxorb |= endcond_tormax;
                }
            }
            if (active_polmax)
            {
                if (fabs(p_f->theta[i]) > params->max_number_of_poloidal_orbits)
                {
                    maxorb |= endcond_polmax;
                }
                else if (
                    p_f->bounces[i] - 1 >=
                    (int)(params->max_number_of_poloidal_orbits / CONST_2PI))
                {
                    maxorb |= endcond_polmax;
                }
            }
            if (params->require_both_tor_and_pol && maxorb & endcond_tormax &&
                maxorb & endcond_polmax)
            {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            else if (!params->require_both_tor_and_pol && maxorb)
            {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }

            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if (active_cpumax)
            {
                if (p_f->cputime[i] > params->max_real_time)
                {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }

            /* If hybrid mode is used, check whether this marker meets the
             * hybrid condition. */
            if (params->simulation_mode == simulate_mode_hybrid)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= endcond_hybrid;
                    p_f->running[i] = 0;
                }
            }

            /* Zero end condition if error happened in this function */
            if (p_f->err[i])
            {
                p_f->endcond[i] = 0;
            }
        }
    }
}

void endcond_check_ml(
    particle_simd_ml *p_f, particle_simd_ml *p_i, sim_data *sim)
{
    sim_parameters *params = sim->params;
    int active_tlim = params->endcond_active & endcond_tlim;
    int active_wall = params->endcond_active & endcond_wall;
    int active_rhomax = params->endcond_active & endcond_rhomax;
    int active_rhomin = params->endcond_active & endcond_rhomin;
    int active_polmax = params->endcond_active & endcond_polmax;
    int active_tormax = params->endcond_active & endcond_tormax;
    int active_cpumax = params->endcond_active & endcond_cpumax;

#pragma omp simd
    for (int i = 0; i < NSIMD; i++)
    {
        if (p_f->running[i])
        {
            /* Check if the marker time exceeds simulation time */
            if (active_tlim)
            {
                if (!params->reverse_time &&
                    p_f->time[i] > params->lab_time_limit)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if (params->reverse_time &&
                    p_f->time[i] < params->lab_time_limit)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
                if (p_f->mileage[i] > params->max_mileage)
                {
                    p_f->endcond[i] |= endcond_tlim;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if (active_wall)
            {
                real w_coll = 0;
                int tile = wall_hit_wall(
                    p_i->r[i], p_i->phi[i], p_i->z[i], p_f->r[i], p_f->phi[i],
                    p_f->z[i], &sim->wall_data, &w_coll);
                if (tile > 0)
                {
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if (active_rhomax)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
            }
            if (active_rhomin)
            {
                if (p_f->rho[i] < params->rho_coordinate_limits[0])
                {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if (active_tormax)
            {
                if (fabs(p_f->phi[i]) > params->max_number_of_toroidal_orbits)
                {
                    maxorb |= endcond_tormax;
                }
            }
            if (active_polmax)
            {
                if (fabs(p_f->theta[i]) > params->max_number_of_poloidal_orbits)
                {
                    maxorb |= endcond_polmax;
                }
            }
            if (params->require_both_tor_and_pol && maxorb & endcond_tormax &&
                maxorb & endcond_polmax)
            {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }
            else if (!params->require_both_tor_and_pol && maxorb)
            {
                p_f->endcond[i] |= maxorb;
                p_f->running[i] = 0;
            }

            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if (active_cpumax)
            {
                if (p_f->cputime[i] > params->max_real_time)
                {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }
        }
    }
}
