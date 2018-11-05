/**
 * @file endcond.c
 * @brief Marker simulation end conditions
 *
 * In the absence of errors, marker simulation is ended when marker meets even
 * one of the active end conditions. User can choose which end conditions are
 * active.
 *
 * The end conditions are:
 * - tmax: Marker time exceeds maximum simulation time
 *
 * - emin: Marker energy is below minimum value
 *
 * - therm: energy is below value derived from local thermal electron energy
 *
 * - wall: Marker has intersected wall
 *
 * - rhomin: Marker has reached minimum rho (normalized poloidal flux) value
 *
 * - rhomax: has reached maximum rho value
 *
 * - polmax: The total cumulative distance marker has travelled poloidally
 *   exceeds maximum value
 *
 * - tormax: The total cumulative distance marker has travelled toroidally
 *   exceeds maximum value
 *
 * - cpumax: Marker simulation has exceeded maximum wall time
 *
 * - hybrid: Not an end condition per se but used to notate that the guiding
 *   center simulation will be resumed as a gyro-orbit simulation
 *
 * Note: emin and therm, rhomin and rhomax, and polmax and tormax are pairs
 * where both members are either active at the same time or neither is i.e.,
 * they cannot be toggled separately.
 *
 * As magnetic field lines have no energy, emin and therm are never checked for
 * them. Guiding centers are the only markers for which hybrid is checked.
 *
 * In the code, the end conditions are represented as bit arrays with each bit
 * corresponding to a specific end condition. Each marker has a field "endcond",
 * and when marker meets an end condition, the corresponding bit is flagged.
 * This way if marker simultaneously meets several end conditions, all can be
 * flagged.
 *
 * Additionally, when marker meets an end condition, its running state is set to
 * False which notates its simulation should be discontinued. If the end
 * condition is wall collision, the ID of the wall element the marker collided
 * with is stored in the marker fields.
 */
#include <math.h>
#include "endcond.h"
#include "particle.h"
#include "simulate.h"
#include "physlib.h"
#include "consts.h"
#include "math.h"
#include "plasma.h"

/**
 * @brief Check end conditions for FO markers
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void endcond_check_fo(particle_simd_fo* p_f, particle_simd_fo* p_i,
                      sim_data* sim) {

    /* Note which end conditions are set as active.
       Only these ones are checked */
    int active_tmax      = sim->endcond_active & endcond_tmax;
    int active_wall      = sim->endcond_active & endcond_wall;
    int active_emin      = sim->endcond_active & endcond_emin;
    int active_rholim    = sim->endcond_active & endcond_rhomax;
    int active_orbitlim  = sim->endcond_active & endcond_polmax;
    int active_cpumax    = sim->endcond_active & endcond_cpumax;

#pragma omp simd
    for(int i = 0; i < NSIMD; i++) {
        if(p_f->running[i]) {

            /* Check if the marker time exceeds simulation time */
            if(active_tmax) {
                if(p_f->time[i] > sim->endcond_maxSimTime) {
                    p_f->endcond[i] |= endcond_tmax;
                    p_f->running[i] = 0;
                }
            }

            /* Check, using the wall collision module, whether marker hit wall
             * during this time-step. Store the wall element ID if it did. */
            if(active_wall) {
                int tile = wall_hit_wall(
                    p_i->r[i], p_i->phi[i], p_i->z[i],
                    p_f->r[i], p_f->phi[i], p_f->z[i], &sim->wall_data);
                if(tile > 0) {
                    p_f->walltile[i] = tile;
                    p_f->endcond[i] |= endcond_wall;
                    p_f->running[i] = 0;
                }
            }

            /* Evaluate marker energy, and check if it is below the minimum
             * energy limit or local thermal energy limit */
            if(active_emin) {
                real vnorm = math_normc(
                    p_f->rdot[i], p_f->phidot[i] * p_f->r[i], p_f->zdot[i]);
                real gamma = physlib_relfactorv_fo(vnorm);
                real ekin = CONST_C2 * p_f->mass[i] * (gamma - 1);
                real Te = plasma_eval_temp(p_f->rho[i], 0, &sim->plasma_data)
                    * CONST_KB;

                if(ekin < sim->endcond_minEkin) {
                    p_f->endcond[i] |= endcond_emin;
                    p_f->running[i] = 0;
                }
                if(ekin < (sim->endcond_minEkinPerTe * Te)) {
                    p_f->endcond[i] |= endcond_therm;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker is not within the rho limits */
            if(active_rholim) {
                if(p_f->rho[i] > sim->endcond_maxRho) {
                    p_f->endcond[i] |= endcond_rhomax;
                    p_f->running[i] = 0;
                }
                if(p_f->rho[i] < sim->endcond_minRho) {
                    p_f->endcond[i] |= endcond_rhomin;
                    p_f->running[i] = 0;
                }
            }

            /* Check if marker exceeds toroidal or poloidal limits */
            if(active_orbitlim) {
                if(fabs(p_f->phi[i]) > sim->endcond_maxTorOrb) {
                    p_f->endcond[i] |= endcond_tormax;
                    p_f->running[i] = 0;
                }
                if(fabs(p_f->pol[i]) > sim->endcond_maxPolOrb) {
                    p_f->endcond[i] |= endcond_polmax;
                    p_f->running[i] = 0;
                }
            }

            /* Check if the time spent simulating this marker exceeds the
             * given limit*/
            if(active_cpumax) {
                if(p_f->cputime[i] > sim->endcond_maxCpuTime) {
                    p_f->endcond[i] |= endcond_cpumax;
                    p_f->running[i] = 0;
                }
            }
        }
    }
}

/**
 * @brief Check end conditions for GC markers
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 *
 * @todo Hybrid condition checks whether marker is over maximum rho limit. More
 *       smarter check is required.
 */
void endcond_check_gc(particle_simd_gc* p_f, particle_simd_gc* p_i,
                      sim_data* sim) {
    int i;

    int active_tmax      = sim->endcond_active & endcond_tmax;
    int active_wall      = sim->endcond_active & endcond_wall;
    int active_emin      = sim->endcond_active & endcond_emin;
    int active_rholim    = sim->endcond_active & endcond_rhomax;
    int active_orbitlim  = sim->endcond_active & endcond_polmax;
    int active_cpumax    = sim->endcond_active & endcond_cpumax;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {

        /* Check if the marker time exceeds simulation time */
        if(active_tmax) {
            if(p_f->time[i] > sim->endcond_maxSimTime) {
                p_f->endcond[i] |= endcond_tmax;
                p_f->running[i] = 0;
            }
        }

        /* Check, using the wall collision module, whether marker hit wall
         * during this time-step. Store the wall element ID if it did. */
        if(active_wall) {
            int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                                     p_f->r[i], p_f->phi[i], p_f->z[i],
                                     &sim->wall_data);
            if(tile > 0) {
                p_f->walltile[i] = tile;
                p_f->endcond[i] |= endcond_wall;
                p_f->running[i] = 0;
            }
        }

        /* Evaluate marker energy, and check if it is below the minimum
         * energy limit or local thermal energy limit */
        if(active_emin) {
            real Bnorm = math_normc(p_f->B_r[i], p_f->B_phi[i], p_f->B_z[i]);
            real gamma = physlib_relfactorv_gc(
                p_f->mass[i], p_f->mu[i], p_f->vpar[i], Bnorm);
            real ekin = CONST_C2 * p_f->mass[i] * (gamma - 1);
            real Te = plasma_eval_temp(p_f->rho[i], 0, &sim->plasma_data)
                * CONST_KB;

            if(ekin < sim->endcond_minEkin) {
                p_f->endcond[i] |= endcond_emin;
                p_f->running[i] = 0;
            }
            if(ekin < (sim->endcond_minEkinPerTe * Te)) {
                p_f->endcond[i] |= endcond_therm;
                p_f->running[i] = 0;
            }
        }

        /* Check if marker is not within the rho limits */
        if(active_rholim) {
            if(p_f->rho[i] > sim->endcond_maxRho) {
                p_f->endcond[i] |= endcond_rhomax;
                p_f->running[i] = 0;
            }
            if(p_f->rho[i] < sim->endcond_minRho) {
                p_f->endcond[i] |= endcond_rhomin;
                p_f->running[i] = 0;
            }
        }

        /* Check if marker exceeds toroidal or poloidal limits */
        if(active_orbitlim) {
            if(fabs(p_f->phi[i]) > sim->endcond_maxTorOrb) {
                p_f->endcond[i] |= endcond_tormax;
                p_f->running[i] = 0;
            }
            if(fabs(p_f->pol[i]) > sim->endcond_maxPolOrb) {
                p_f->endcond[i] |= endcond_polmax;
                p_f->running[i] = 0;
            }
        }

        /* Check if the time spent simulating this marker exceeds the
         * given limit*/
        if(active_cpumax) {
            if(p_f->cputime[i] > sim->endcond_maxCpuTime) {
                p_f->endcond[i] |= endcond_cpumax;
                p_f->running[i] = 0;
            }
        }

        /* If hybrid mode is used, check whether this marker meets the hybrid
         * condition. */
        if(sim->sim_mode == 3) {
            if(p_f->rho[i] > sim->endcond_maxRho) {
                p_f->endcond[i] |= endcond_hybrid;
                p_f->running[i] = 0;
            }
        }
    }
}

/**
 * @brief Check end conditions for ML markers
 *
 * The end conditions are checked for all markers within the SIMD marker struct.
 *
 * @param p_f pointer to SIMD struct storing marker states at the end of current
 *        time-step
 * @param p_i pointer to SIMD struct storing marker states at the beginning of
 *        current time-step
 */
void endcond_check_ml(particle_simd_ml* p_f, particle_simd_ml* p_i,
                      sim_data* sim) {
    int i;

    int active_tmax      = sim->endcond_active & endcond_tmax;
    int active_wall      = sim->endcond_active & endcond_wall;
    int active_rholim    = sim->endcond_active & endcond_rhomax;
    int active_orbitlim  = sim->endcond_active & endcond_polmax;
    int active_cpumax    = sim->endcond_active & endcond_cpumax;

    #pragma omp simd
    for(i = 0; i < NSIMD; i++) {

        /* Check if the marker time exceeds simulation time */
        if(active_tmax) {
            if(p_f->time[i] > sim->endcond_maxSimTime) {
                p_f->endcond[i] |= endcond_tmax;
                p_f->running[i] = 0;
            }
        }

        /* Check, using the wall collision module, whether marker hit wall
         * during this time-step. Store the wall element ID if it did. */
        if(active_wall) {
            int tile = wall_hit_wall(p_i->r[i], p_i->phi[i], p_i->z[i],
                                     p_f->r[i], p_f->phi[i], p_f->z[i],
                                     &sim->wall_data);
            if(tile > 0) {
                p_f->walltile[i] = tile;
                p_f->endcond[i] |= endcond_wall;
                p_f->running[i] = 0;
            }
        }

        /* Check if marker is not within the rho limits */
        if(active_rholim) {
            if(p_f->rho[i] > sim->endcond_maxRho) {
                p_f->endcond[i] |= endcond_rhomax;
                p_f->running[i] = 0;
            }
            if(p_f->rho[i] < sim->endcond_minRho) {
                p_f->endcond[i] |= endcond_rhomin;
                p_f->running[i] = 0;
            }
        }

        /* Check if marker exceeds toroidal or poloidal limits */
        if(active_orbitlim) {
            if(fabs(p_f->phi[i]) > sim->endcond_maxTorOrb) {
                p_f->endcond[i] |= endcond_tormax;
                p_f->running[i] = 0;
            }
            if(fabs(p_f->pol[i]) > sim->endcond_maxPolOrb) {
                p_f->endcond[i] |= endcond_polmax;
                p_f->running[i] = 0;
            }
        }

        /* Check if the time spent simulating this marker exceeds the
         * given limit*/
        if(active_cpumax) {
            if(p_f->cputime[i] > sim->endcond_maxCpuTime) {
                p_f->endcond[i] |= endcond_cpumax;
                p_f->running[i] = 0;
            }
        }
    }
}
