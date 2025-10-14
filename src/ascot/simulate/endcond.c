/**
 * Implements endcond.h.
 */
#include "endcond.h"
#include "consts.h"
#include "data/marker.h"
#include "data/plasma.h"
#include "data/wall.h"
#include "datatypes.h"
#include "options.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include <math.h>
#include <stdint.h>

void endcond_check_go(
    MarkerGyroOrbit *p_f, MarkerGyroOrbit *p_i, Simulation *sim)
{
    Options *params = sim->options;

    /* Note which end conditions are set as active.
       Only these ones are checked */
    endcond_t active_tlim = params->endcond_active & ENDCOND_TLIM;
    endcond_t active_wall = params->endcond_active & ENDCOND_WALL;
    endcond_t active_emin = params->endcond_active & ENDCOND_EMIN;
    endcond_t active_therm = params->endcond_active & ENDCOND_THERM;
    endcond_t active_rhomax = params->endcond_active & ENDCOND_RHOMAX;
    endcond_t active_rhomin = params->endcond_active & ENDCOND_RHOMIN;
    endcond_t active_polmax = params->endcond_active & ENDCOND_POLMAX;
    endcond_t active_tormax = params->endcond_active & ENDCOND_TORMAX;
    endcond_t active_cpumax = params->endcond_active & ENDCOND_CPUMAX;
    endcond_t active_neutr = params->endcond_active & ENDCOND_NEUTR;
    endcond_t active_ioniz = params->endcond_active & ENDCOND_IONIZ;

    GPU_PARALLEL_LOOP_ALL_LEVELS
    for (size_t i = 0; i < p_f->n_mrk; i++)
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
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
                if (params->reverse_time &&
                    p_f->time[i] < params->lab_time_limit)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
                if (p_f->mileage[i] > params->max_mileage)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if (active_wall)
            {
                real w_coll = 0;
                int tile = Wall_eval_intersection(
                    &w_coll, p_i->r[i], p_i->phi[i], p_i->z[i], p_f->r[i],
                    p_f->phi[i], p_f->z[i], &sim->wall);
                if (tile > 0)
                {
                    real w = w_coll;
                    p_f->time[i] =
                        p_i->time[i] + w * (p_f->time[i] - p_i->time[i]);
                    p_f->r[i] = p_i->r[i] + w * (p_f->r[i] - p_i->r[i]);
                    p_f->phi[i] = p_i->phi[i] + w * (p_f->phi[i] - p_i->phi[i]);
                    p_f->z[i] = p_i->z[i] + w * (p_f->z[i] - p_i->z[i]);

                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= ENDCOND_WALL;
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
                err_t errflag = Plasma_eval_temperature(
                    &Ti, p_f->rho[i], p_f->r[i], p_f->phi[i], p_f->z[i],
                    p_f->time[i], 1, &sim->plasma);

                /* Error handling */
                if (errflag)
                {
                    p_f->err[i] = errflag;
                    p_f->running[i] = 0;
                    Ti = 0;
                }

                if (active_emin && (ekin < params->min_energy))
                {
                    p_f->endcond[i] |= ENDCOND_EMIN;
                    p_f->running[i] = 0;
                }
                if (active_therm &&
                    (ekin < (params->min_local_thermal_energy * Ti)))
                {
                    p_f->endcond[i] |= ENDCOND_THERM;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if (active_rhomax)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= ENDCOND_RHOMAX;
                    p_f->running[i] = 0;
                }
            }
            if (active_rhomin)
            {
                if (p_f->rho[i] < params->rho_coordinate_limits[0])
                {
                    p_f->endcond[i] |= ENDCOND_RHOMIN;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if (active_tormax)
            {
                if (fabs(p_f->phi[i]) > params->max_number_of_toroidal_orbits)
                {
                    maxorb |= ENDCOND_TORMAX;
                }
            }
            if (active_polmax)
            {
                if (fabs(p_f->theta[i]) > params->max_number_of_poloidal_orbits)
                {
                    maxorb |= ENDCOND_POLMAX;
                }
                else if (
                    p_f->bounces[i] - 1 >=
                    (int)(params->max_number_of_poloidal_orbits / CONST_2PI))
                {
                    maxorb |= ENDCOND_POLMAX;
                }
            }
            if (params->require_both_tor_and_pol && maxorb & ENDCOND_TORMAX &&
                maxorb & ENDCOND_POLMAX)
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
                    p_f->endcond[i] |= ENDCOND_CPUMAX;
                    p_f->running[i] = 0;
                }
            }
            /* Check if the particle has been neutralized */
            if (active_neutr)
            {
                if (p_i->charge[i] != 0.0 && p_f->charge[i] == 0.0)
                {
                    p_f->endcond[i] |= ENDCOND_NEUTR;
                    p_f->running[i] = 0;
                }
            }

            /* Check if the particle has been ionized */
            if (active_ioniz)
            {
                if (p_i->charge[i] == 0.0 && p_f->charge[i] != 0.0)
                {
                    p_f->endcond[i] |= ENDCOND_IONIZ;
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
    MarkerGuidingCenter *p_f, MarkerGuidingCenter *p_i, Simulation *sim)
{
    Options *params = sim->options;

    endcond_t active_tlim = params->endcond_active & ENDCOND_TLIM;
    endcond_t active_wall = params->endcond_active & ENDCOND_WALL;
    endcond_t active_emin = params->endcond_active & ENDCOND_EMIN;
    endcond_t active_therm = params->endcond_active & ENDCOND_THERM;
    endcond_t active_rhomax = params->endcond_active & ENDCOND_RHOMAX;
    endcond_t active_rhomin = params->endcond_active & ENDCOND_RHOMIN;
    endcond_t active_polmax = params->endcond_active & ENDCOND_POLMAX;
    endcond_t active_tormax = params->endcond_active & ENDCOND_TORMAX;
    endcond_t active_cpumax = params->endcond_active & ENDCOND_CPUMAX;

#pragma omp simd
    for (size_t i = 0; i < NSIMD; i++)
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
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
                if (params->reverse_time &&
                    p_f->time[i] < params->lab_time_limit)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
                if (p_f->mileage[i] > params->max_mileage)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if (active_wall)
            {
                real w_coll = 0;
                int tile = Wall_eval_intersection(
                    &w_coll, p_i->r[i], p_i->phi[i], p_i->z[i], p_f->r[i],
                    p_f->phi[i], p_f->z[i], &sim->wall);
                if (tile > 0)
                {
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= ENDCOND_WALL;
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
                err_t errflag = Plasma_eval_temperature(
                    &Ti, p_f->rho[i], p_f->r[i], p_f->phi[i], p_f->z[i],
                    p_f->time[i], 1, &sim->plasma);

                /* Error handling */
                if (errflag)
                {
                    p_f->err[i] = errflag;
                    p_f->running[i] = 0;
                    Ti = 0;
                }

                if (active_emin && (ekin < params->min_energy))
                {
                    p_f->endcond[i] |= ENDCOND_EMIN;
                    p_f->running[i] = 0;
                }
                if (active_therm &&
                    (ekin < (params->min_local_thermal_energy * Ti)))
                {
                    p_f->endcond[i] |= ENDCOND_THERM;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if (active_rhomax)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= ENDCOND_RHOMAX;
                    p_f->running[i] = 0;
                }
            }
            if (active_rhomin)
            {
                if (p_f->rho[i] < params->rho_coordinate_limits[0])
                {
                    p_f->endcond[i] |= ENDCOND_RHOMIN;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if (active_tormax)
            {
                if (fabs(p_f->phi[i]) > params->max_number_of_toroidal_orbits)
                {
                    maxorb |= ENDCOND_TORMAX;
                }
            }
            if (active_polmax)
            {
                if (fabs(p_f->theta[i]) > params->max_number_of_poloidal_orbits)
                {
                    maxorb |= ENDCOND_POLMAX;
                }
                else if (
                    p_f->bounces[i] - 1 >=
                    (int)(params->max_number_of_poloidal_orbits / CONST_2PI))
                {
                    maxorb |= ENDCOND_POLMAX;
                }
            }
            if (params->require_both_tor_and_pol && maxorb & ENDCOND_TORMAX &&
                maxorb & ENDCOND_POLMAX)
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
                    p_f->endcond[i] |= ENDCOND_CPUMAX;
                    p_f->running[i] = 0;
                }
            }

            /* If hybrid mode is used, check whether this marker meets the
             * hybrid condition. */
            if (params->simulation_mode == simulate_mode_hybrid)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= ENDCOND_HYBRID;
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

void endcond_check_fl(
    MarkerFieldLine *p_f, MarkerFieldLine *p_i, Simulation *sim)
{
    Options *params = sim->options;
    endcond_t active_tlim = params->endcond_active & ENDCOND_TLIM;
    endcond_t active_wall = params->endcond_active & ENDCOND_WALL;
    endcond_t active_rhomax = params->endcond_active & ENDCOND_RHOMAX;
    endcond_t active_rhomin = params->endcond_active & ENDCOND_RHOMIN;
    endcond_t active_polmax = params->endcond_active & ENDCOND_POLMAX;
    endcond_t active_tormax = params->endcond_active & ENDCOND_TORMAX;
    endcond_t active_cpumax = params->endcond_active & ENDCOND_CPUMAX;

#pragma omp simd
    for (size_t i = 0; i < NSIMD; i++)
    {
        if (p_f->running[i])
        {
            /* Check if the marker time exceeds simulation time */
            if (active_tlim)
            {
                if (!params->reverse_time &&
                    p_f->time[i] > params->lab_time_limit)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
                if (params->reverse_time &&
                    p_f->time[i] < params->lab_time_limit)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
                if (p_f->mileage[i] > params->max_mileage)
                {
                    p_f->endcond[i] |= ENDCOND_TLIM;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if (active_wall)
            {
                real w_coll = 0;
                int tile = Wall_eval_intersection(
                    &w_coll, p_i->r[i], p_i->phi[i], p_i->z[i], p_f->r[i],
                    p_f->phi[i], p_f->z[i], &sim->wall);
                if (tile > 0)
                {
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= ENDCOND_WALL;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if (active_rhomax)
            {
                if (p_f->rho[i] > params->rho_coordinate_limits[1])
                {
                    p_f->endcond[i] |= ENDCOND_RHOMAX;
                    p_f->running[i] = 0;
                }
            }
            if (active_rhomin)
            {
                if (p_f->rho[i] < params->rho_coordinate_limits[0])
                {
                    p_f->endcond[i] |= ENDCOND_RHOMIN;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            int maxorb = 0;
            if (active_tormax)
            {
                if (fabs(p_f->phi[i]) > params->max_number_of_toroidal_orbits)
                {
                    maxorb |= ENDCOND_TORMAX;
                }
            }
            if (active_polmax)
            {
                if (fabs(p_f->theta[i]) > params->max_number_of_poloidal_orbits)
                {
                    maxorb |= ENDCOND_POLMAX;
                }
            }
            if (params->require_both_tor_and_pol && maxorb & ENDCOND_TORMAX &&
                maxorb & ENDCOND_POLMAX)
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
                    p_f->endcond[i] |= ENDCOND_CPUMAX;
                    p_f->running[i] = 0;
                }
            }
        }
    }
}
