/**
 * @file bmc_push.c
 * @brief Methods to calculate push-matrix for the backward Monte Carlo
 *
 * The push-matrix is used to update the probability matrix in the BMC scheme.
 * Here it is constructed as follows:
 *
 *   1. A marker is initialized on each node in the mesh.
 *   2. Each marker is traced for a given amount of time (without collisions)
 *   3. The marker final position is recorded and also if it hit the wall
 *      (or FILD).
 *   4. If there was no wall hit, each marker is split into HERMITE_KNOT number
 *      of markers that experience collisions for the same time as
 *      the orbit-following lasted.
 *   5. The final positions are updated and returned (this is the push-matrix).
 */
#include <stdlib.h>
#include "../ascot5.h"
#include "../simulate.h"
#include "../particle.h"
#include "../physlib.h"
#include "../consts.h"
#include "../bmc_mesh.h"
#include "step/step_gc_rk4.h"
#include "mccc/mccc.h"
#include "simulate_bmc.h"

/**
 * @brief Generate push-matrix by initializing guiding center markers on mesh
 * nodes.
 *
 * A subsection of the mesh can be pushed by specifying indices start and stop,
 * in which case the push matrix is generated for mesh elements [start,stop).
 * This makes it possible to divide this operation among MPI processes.
 *
 * The generated matrix is given as arrays with format
 * [imesh*HERMITE_KNOT + iknot].
 *
 * @param sim pointer to simulation data struct
 * @param mesh pointer to mesh struct
 * @param h time step for how long this push-matrix pushes particles [s]
 * @param time current time instant (for evaluating background quantities) [s]
 * @param start the first mesh index
 * @param stop the final mesh index
 * @param r marker final R-coordinates [m]
 * @param phi marker final (periodic, not cumulative) phi-coordinates [rad]
 * @param z marker final z-coordinates [m]
 * @param ppara marker final parallel momentum coordinate [kg*m/s]
 * @param pperp marker final perpendicular momentum coordinate [kg*m/s]
 * @param fate flag indicating whether the marker terminated with error [-1],
 *        hit wall [1], hit FILD [2], or finished normally [0]
 */
void simulate_bmc_gc(
    sim_data* sim, bmc_mesh* mesh, real h, real time, size_t start, size_t stop,
    real* r, real* phi, real* z, real* ppara, real* pperp, int* fate) {

    real h_orb[NSIMD];
    real h_coll[NSIMD];
    real hermite_k[HERMITE_KNOTS] = HERMITE_K;
    real hermite_k_nsimd[5 * NSIMD * HERMITE_KNOTS];
    for(int i=0; i<NSIMD; i++) {
        h_orb[i]  = h / sim->bmc_orbit_subcycles;
        h_coll[i] = h;
        for(int k=0; k<HERMITE_KNOTS; k++) {
            hermite_k_nsimd[k*5*NSIMD + 0*NSIMD + i] = 0.0;
            hermite_k_nsimd[k*5*NSIMD + 1*NSIMD + i] = 0.0;
            hermite_k_nsimd[k*5*NSIMD + 2*NSIMD + i] = 0.0;
            hermite_k_nsimd[k*5*NSIMD + 3*NSIMD + i] = 0.0;
            hermite_k_nsimd[k*5*NSIMD + 4*NSIMD + i] = hermite_k[k];
        }
    }

    /* Go through the mesh in chunks of NSIMD */
    #pragma omp parallel for \
        shared(sim, mesh, h_orb, h_coll, hermite_k_nsimd, start, stop, \
        r, phi, z, ppara, pperp, fate)
    for(size_t iprt=start; iprt < stop + NSIMD; iprt += NSIMD) {
        /* Initialize markers from mesh */
        int lost[NSIMD];
        particle_simd_gc p;
        for(int i=0; i<NSIMD; i++) {

            /* Find the position of the node (i.e. initial marker position) */
            real origin[5];
            bmc_mesh_index2pos(mesh, iprt+i, origin);

            real ppara0 = origin[3];
            real pperp0 = origin[4];
            real pnorm = sqrt(ppara0 * ppara0 + pperp0 * pperp0);
            real xi    = ppara0 / pnorm;
            real ekin  = physlib_Ekin_pnorm(sim->bmc_mass, pnorm);
            //if(xi == -1) xi = -0.9999; //Not sure if these are needed
            //if(xi == 1) xi = 0.9999;

            particle_gc gc;
            gc.r      = origin[0];
            gc.phi    = origin[1];
            gc.z      = origin[2];
            gc.anum   = sim->bmc_anum;
            gc.znum   = sim->bmc_znum;
            gc.mass   = sim->bmc_mass;
            gc.charge = sim->bmc_charge;
            gc.time   = time;
            gc.weight = 1.0;
            gc.id     = 1;
            gc.zeta   = 0.0;
            gc.pitch  = xi;
            gc.energy = ekin;

            particle_state ps;
            particle_input_gc_to_state(&gc, &ps, &sim->B_data);
            particle_state_to_gc(&ps, 0, &p, i, &sim->B_data);
            lost[i] = 0;
            if(!wall_2d_inside(origin[0], origin[2], &sim->wall_data.w2d)) {
                lost[i] = -1;
            }

            if(iprt + i >= stop || p.err[i]) {
                /* No more mesh points to initialize; fill rest of the array
                 * with dummy markers */
                p.id[i] = -1;
                p.running[i] = 0;
            }
        }

        /* Take a number of orbit-following steps. This is then followed by
         * a single collisional step but with a larger time-step that
         * corresponds to the whole orbit-following part */
        for(int j=0; j<sim->bmc_orbit_subcycles; j++) {
            real r0[NSIMD], phi0[NSIMD], z0[NSIMD];
            for(int i=0; i<NSIMD; i++) {
                if(!lost[i]) {
                    r0[i]   = p.r[i];
                    phi0[i] = p.phi[i];
                    z0[i]   = p.z[i];
                }
            }
            step_gc_rk4(&p, h_orb, &sim->B_data, &sim->E_data);
            for(int i=0; i<NSIMD; i++) {
                if(!p.err[i] && !lost[i]) {
                    real w_coll;
                    int tile = wall_hit_wall(
                        r0[i], phi0[i], z0[i], p.r[i], p.phi[i], p.z[i],
                        &sim->wall_data, &w_coll);
                    if(tile > 0) {
                        lost[i] = 1.0;
                        if( wall_get_flag(&sim->wall_data, tile-1) == 1 ) {
                            lost[i] = 2.0;
                        }
                    }
                }
                else if(!p.err[i] && lost[i]) {

                }
                else {
                    lost[i] = -1.0;
                }
            }
        }

        /* The collisional part is evaluated with the Gauss-Hermite quadrule */
        particle_simd_gc p_knot;
        for(int i_knot=0; i_knot<HERMITE_KNOTS; i_knot++) {
            for(int i=0; i<NSIMD; i++) {
                /* p contains the marker position after the orbit step and
                 * p_knot will contain the final position. This copying is done
                 * since p_knot is computed several times from the same
                 * initial p */
                particle_copy_gc(&p, i, &p_knot, i);
            }
            mccc_gc_euler(&p_knot, h_coll, &sim->B_data, &sim->plasma_data,
                          &sim->mccc_data, &(hermite_k_nsimd[i_knot*5*NSIMD]));
            /* Store the final locations*/
            for(int i=0; i<NSIMD; i++) {
                if(iprt + i >= stop) break;
                r[(iprt+i) * HERMITE_KNOTS + i_knot]     = p_knot.r[i];
                phi[(iprt+i) * HERMITE_KNOTS + i_knot]   =
                    fmod(fmod(p_knot.phi[i], CONST_2PI) + CONST_2PI, CONST_2PI);
                z[(iprt+i) * HERMITE_KNOTS + i_knot]     = p_knot.z[i];
                ppara[(iprt+i) * HERMITE_KNOTS + i_knot] = p_knot.ppar[i];

                real Bnorm = sqrt( p_knot.B_r[i]   * p_knot.B_r[i]
                                 + p_knot.B_phi[i] * p_knot.B_phi[i]
                                 + p_knot.B_z[i]   * p_knot.B_z[i]);
                real pnorm = physlib_gc_p(p_knot.mass[i], p_knot.mu[i],
                                          p_knot.ppar[i], Bnorm);
                real pperp2 = fabs(pnorm * pnorm - p_knot.ppar[i] * p_knot.ppar[i]);
                pperp[(iprt+i) * HERMITE_KNOTS + i_knot] = sqrt(pperp2);

                if(lost[i] != 0) {
                    fate[(iprt+i) * HERMITE_KNOTS + i_knot] = lost[i];
                }
            }
        }
    }
}