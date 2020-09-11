#include "bmc_diag.h"

void diag_move_distribution(sim_offload_data* sim, diag_data* diag_dest, diag_data* diag_src, int dist_length) {
    // copy into diag_dest
    #ifdef MPI
        if (sim->diag_offload_data.dist5D_collect) {
            MPI_Allreduce(diag_src->dist5D.histogram, diag_dest->dist5D.histogram, dist_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            memset(diag_src->dist5D.histogram, 0, dist_length);
        } else if (sim->diag_offload_data.dist6D_collect) {
            MPI_Allreduce(diag_src->dist6D.histogram, diag_dest->dist6D.histogram, dist_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            memset(diag_src->dist6D.histogram, 0, dist_length);
        }
    #else
        if (sim->diag_offload_data.dist5D_collect) {
            memcpy(diag_dest->dist5D.histogram, diag_src->dist5D.histogram, dist_length * sizeof(real));
            memset(diag_src->dist5D.histogram, 0, dist_length * sizeof(real));
        } else if (sim->diag_offload_data.dist6D_collect) {
            memcpy(diag_dest->dist6D.histogram, diag_src->dist6D.histogram, dist_length * sizeof(real));
            memset(diag_src->dist6D.histogram, 0, dist_length * sizeof(real));
        }
    #endif
}

void bmc_update_particles_diag(
    int n_mpi_particles,
    particle_state* ps0,
    particle_state* ps1,
    int* ps_indexes,
    diag_data* diag0,
    diag_data* diag1,
    sim_data* sim,
    int n_montecarlo_steps
) {

    int n_updated = 0;
    for(int i = 0; i < n_mpi_particles; i++) {

        if (sim->sim_mode == simulate_mode_fo) {
            bmc_diag_update_fo(&ps0[i], &ps1[i], ps_indexes[i], diag0, diag1, n_montecarlo_steps);
        }
        else if (sim->sim_mode == simulate_mode_gc) {
            n_updated += bmc_diag_update_gc(&ps0[i], &ps1[i], ps_indexes[i], diag0, diag1, n_montecarlo_steps, &(sim->wall_data.w2d));
        }

    }

    printf("Updated %d particles\n", n_updated);
}

int bmc_diag_update_gc(
    particle_state* ps0,
    particle_state* ps1,
    int p0_index,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps,
    wall_2d_data* w2d
) {
    if(diag0->dist5D_collect) {
        return bmc_diag_5D_update_gc(&diag1->dist5D, &diag0->dist5D,  p0_index, ps1, ps0, n_montecarlo_steps, w2d);
    }

    if(diag0->dist6D_collect) {
        bmc_diag_6D_update_gc(&diag1->dist6D, &diag0->dist6D,  p0_index, ps1, ps0, n_montecarlo_steps);
    }
}

void bmc_diag_update_fo(
    particle_state* ps0,
    particle_state* ps1,
    int p0_index,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps
) {
    if(diag0->dist5D_collect) {
        bmc_diag_5D_update_fo(&diag1->dist5D, &diag0->dist5D,  p0_index, ps1, ps0, n_montecarlo_steps);
    }

    if(diag0->dist6D_collect) {
        bmc_diag_6D_update_fo(&diag1->dist6D, &diag0->dist6D,  p0_index, ps1, ps0, n_montecarlo_steps);
    }
}


int bmc_diag_5D_update_gc(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int p0_index,
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps,
        wall_2d_data* w2d
    ) {

    if (ps1->err) {
        return 0;
    }

    // check if the particle escaped the velocity space
    real vperp = sqrt(2 * sqrt(ps1->B_r*ps1->B_r
        +ps1->B_phi*ps1->B_phi
        +ps1->B_z*ps1->B_z)
        * ps1->mu / ps1->mass);
                
    if ((ps1->vpar > dist1->max_vpara) || (ps1->vpar < dist1->min_vpara) ||
        (vperp > dist1->max_vperp) || (vperp < dist1->min_vperp)) {

        // outside velocity space
        return 0;
    }

    // deposit the final state with a linear interpolation on the closest spatial vertices
    // deposit in the velocity space is nearest neighbour 
    int p1_indexes[8];
    int p1_target_hit[8];
    real p1_weights[8];
    bmc_dist5D_gc_indexes(ps0, p1_indexes, p1_weights, p1_target_hit, ps1, dist0, w2d);
    for (int i=0; i<8; i++) {
        if (p1_indexes[i] >= 0) {
            if (p1_target_hit[i]) {
                // particle hit the target domain. Set the relative probabiity to 1
                dist1->histogram[p0_index] += p1_weights[i] / n_montecarlo_steps;
            } else {
                // particle didn't hit the target domain. Weight the probability with the last probability matrix
                dist1->histogram[p0_index] += dist0->histogram[p1_indexes[i]] * p1_weights[i] / n_montecarlo_steps;
            }
        }
    }

    return 1;
}

void bmc_diag_5D_update_fo(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int p0_index,
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {

    // check if the particle hit the wall
    if (ps1->endcond & endcond_wall) {
        if (bmc_walltile_in_target(ps1->walltile)) {
            // particle ended in target domain, set the weight to 1
            dist1->histogram[p0_index] += 1. / n_montecarlo_steps;
        }
    } else {
        int p1_index = bmc_dist5D_fo_index(ps1, dist0);
        dist1->histogram[p0_index] += dist0->histogram[p1_index] / n_montecarlo_steps;
    }
}

void bmc_diag_6D_update_gc(
        dist_6D_data* dist1,
        dist_6D_data* dist0,
        int p0_index,
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {

    // check if the particle hit the wall
    if (ps1->endcond & endcond_wall) {
        if (bmc_walltile_in_target(ps1->walltile)) {
            // particle ended in target domain, set the weight to 1
            dist1->histogram[p0_index] += 1. / n_montecarlo_steps;
        }
    } else {
        int p1_index = bmc_dist6D_gc_index(ps1, dist0);
        dist1->histogram[p0_index] += dist0->histogram[p1_index] / n_montecarlo_steps;
    }
}

void bmc_diag_6D_update_fo(
        dist_6D_data* dist1,
        dist_6D_data* dist0,
        int p0_index,
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {

    // check if the particle hit the wall
    if (ps1->endcond & endcond_wall) {
        if (bmc_walltile_in_target(ps1->walltile)) {
            // particle ended in target domain, set the weight to 1
            dist1->histogram[p0_index] += 1. / n_montecarlo_steps;
        }
    } else {
        int p1_index = bmc_dist6D_fo_index(ps1, dist0);
        dist1->histogram[p0_index] += dist0->histogram[p1_index] / n_montecarlo_steps;
    }
}

/**
 * Deposit a particle with a linear interpolation in the 3 spatial dimensions
 * and to the closest neighbor in the 2 velocity dimensions.
 * 
 * Returns a list of indexes where the particle is deposited, the relative weights and
 * flags indicating whether the particle hit the target domain or not
 * 
 * @param ps0 Initial state of the particle. Used to compute the wall hit
 * @param indexes Pointer to the output indexes array. Length must be 8
 * @param weights Pointer to the output weights array. Length must be 8
 * @param target_hit Pointer to the output target hit flag array. Length must be 8
 * @param ps Pointer to the final state of the particle
 * @param dist Pointer to the 5d distribution
 * @param w2d Pointer to the 2D wall struct
 * 
 **/
 void bmc_dist5D_gc_indexes(particle_state* ps0, int* indexes, real* weights, int* target_hit, particle_state* ps, dist_5D_data* dist, wall_2d_data* w2d) {
    real phi;
    real vpara;
    real vperp;
    real weights_rphiz[3];
    int i_r;
    int i_phi;
    int i_z;
    int i_vpara;
    int i_vperp;
    int i_time;
    int i_q;

    weights_rphiz[0] = ((ps->r - dist->min_r)
                / ((dist->max_r - dist->min_r)/dist->n_r));
    i_r = floor(weights_rphiz[0]);
    weights_rphiz[0] = 1. - weights_rphiz[0] + i_r;

    phi = fmod(ps->phi, 2*CONST_PI);
    if(phi < 0) {
        phi = phi + 2*CONST_PI;
    }
    weights_rphiz[1] = ((phi - dist->min_phi)
                / ((dist->max_phi - dist->min_phi)/dist->n_phi));
    i_phi = floor(weights_rphiz[1]);
    weights_rphiz[1] = 1. - weights_rphiz[1] + i_phi;

    weights_rphiz[2] = ((ps->z - dist->min_z)
            / ((dist->max_z - dist->min_z) / dist->n_z));
    i_z = floor(weights_rphiz[2]);
    weights_rphiz[2] = 1. - weights_rphiz[2] + i_z;

    i_vpara = floor((ps->vpar - dist->min_vpara)
                / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));

    vperp = sqrt(2 * sqrt(ps->B_r*ps->B_r
                                +ps->B_phi*ps->B_phi
                                +ps->B_z*ps->B_z)
                    * ps->mu / ps->mass);
    i_vperp = floor((vperp - dist->min_vperp)
                / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));

    i_time = 0;

    i_q = 0;

    real dr = (dist->max_r - dist->min_r)/dist->n_r;
    real dphi = (dist->max_phi - dist->min_phi)/dist->n_phi;
    real dz = (dist->max_z - dist->min_z)/dist->n_z;

    int i = 0;
    real r1, phi1, z1;
    int j_phimod;
    for (int j_r=i_r; j_r<i_r + 2; j_r++)
    for (int j_phi=i_phi; j_phi<i_phi + 2; j_phi++)
    for (int j_z=i_z; j_z<i_z + 2; j_z++) {
        j_phimod = j_phi;
        if (j_phimod>=dist->n_phi) {
            j_phimod = 0;
        }
        r1 = j_r*dr + dist->min_r;
        phi1 = j_phimod*dphi + dist->min_phi;
        z1 = j_z*dz + dist->min_z;
        indexes[i] = dist_5D_index(j_r, j_phimod, j_z, i_vpara, i_vperp, i_time, i_q,
                    dist->n_phi, dist->n_z, dist->n_vpara, dist->n_vperp, 1, 1);
        weights[i] = fabs(weights_rphiz[0] - j_r + i_r) * fabs(weights_rphiz[1] - j_phi + i_phi) * fabs(weights_rphiz[2] - j_z + i_z);
        target_hit[i] = bmc_wall_2d_hit_target(ps0->r, (j_r)*dr + dist->min_r,
                        ps0->phi, (j_phimod)*dphi + dist->min_phi, ps0->z, (j_z)*dz + dist->min_z, w2d);
        i++;
    }
}
int bmc_dist5D_fo_index(particle_state* ps, dist_5D_data* dist) {
    real phi;
    real vpara;
    real vperp;
    int i_r;
    int i_phi;
    int i_z;
    int i_vpara;
    int i_vperp;
    int i_time;
    int i_q;

    i_r = floor((ps->r - dist->min_r)
                / ((dist->max_r - dist->min_r)/dist->n_r));

    phi = fmod(ps->phi, 2*CONST_PI);
    if(phi < 0) {
        phi = phi + 2*CONST_PI;
    }
    i_phi = floor((phi - dist->min_phi)
                / ((dist->max_phi - dist->min_phi)/dist->n_phi));

    i_z = floor((ps->z - dist->min_z)
                / ((dist->max_z - dist->min_z) / dist->n_z));

    vpara = (ps->rdot * ps->B_r +
                (ps->phidot * ps->r)
                * ps->B_phi + ps->zdot * ps->B_z)
                / sqrt(ps->B_r*ps->B_r
                        +ps->B_phi*ps->B_phi
                        + ps->B_z*ps->B_z);
    i_vpara = floor((vpara - dist->min_vpara)
                / ((dist->max_vpara - dist->min_vpara) / dist->n_vpara));

    vperp = sqrt(ps->rdot*ps->rdot + (ps->phidot
                                    *ps->phidot*ps->r*ps->r)
                    + ps->zdot*ps->zdot - vpara*vpara);
    i_vperp = floor((vperp - dist->min_vperp)
                / ((dist->max_vperp - dist->min_vperp) / dist->n_vperp));

    i_time = 0;
    i_q = 0;

    return dist_5D_index(i_r, i_phi, i_z,
                        i_vpara, i_vperp,
                        i_time, i_q,
                        dist->n_phi, dist->n_z,
                        dist->n_vpara, dist->n_vperp,
                        1, 1);
}


int bmc_dist6D_gc_index(particle_state* ps, dist_6D_data* dist) {
    real phi;
    int i_r;
    int i_phi;
    int i_z;
    int i_vr;
    int i_vphi;
    int i_vz;
    int i_time;
    int i_q;

    i_r = floor((ps->r - dist->min_r)
                / ((dist->max_r - dist->min_r)/dist->n_r));

    phi = fmod(ps->phi, 2*CONST_PI);
    if(phi < 0) {
        phi = phi + 2*CONST_PI;
    }
    i_phi = floor((phi - dist->min_phi)
                / ((dist->max_phi - dist->min_phi)/dist->n_phi));

    i_z = floor((ps->z - dist->min_z)
                / ((dist->max_z - dist->min_z) / dist->n_z));

    i_vr = floor((ps->rdot - dist->min_vr)
                / ((dist->max_vr - dist->min_vr) / dist->n_vr));

    i_vphi = floor((ps->phidot*ps->r - dist->min_vphi)
                / ((dist->max_vphi - dist->min_vphi) / dist->n_vphi));

    i_vz = floor((ps->zdot - dist->min_vz)
                / ((dist->max_vz - dist->min_vz) / dist->n_vz));

    i_time = 0;
    i_q = 0;

    return dist_6D_index(i_r, i_phi, i_z,
                        i_vr, i_vphi, i_vz,
                        i_time, i_q,
                        dist->n_phi, dist->n_z,
                        dist->n_vr, dist->n_vphi,
                        dist->n_vz, 1, 1);
}

int bmc_dist6D_fo_index(particle_state* ps, dist_6D_data* dist) {
    real phi;
    int i_r;
    int i_phi;
    int i_z;
    int i_vr;
    int i_vphi;
    int i_vz;
    int i_time;
    int i_q;

    i_r = floor((ps->r - dist->min_r)
                / ((dist->max_r - dist->min_r)/dist->n_r));

    phi = fmod(ps->phi, 2*CONST_PI);
    if(phi < 0) {
        phi = phi + 2*CONST_PI;
    }
    i_phi = floor((phi - dist->min_phi)
                / ((dist->max_phi - dist->min_phi)/dist->n_phi));

    i_z = floor((ps->z - dist->min_z)
                / ((dist->max_z - dist->min_z) / dist->n_z));

    i_vr = floor((ps->rdot - dist->min_vr)
                / ((dist->max_vr - dist->min_vr) / dist->n_vr));

    i_vphi = floor((ps->phidot*ps->r - dist->min_vphi)
                / ((dist->max_vphi - dist->min_vphi) / dist->n_vphi));

    i_vz = floor((ps->zdot - dist->min_vz)
                / ((dist->max_vz - dist->min_vz) / dist->n_vz));

    i_time = 0;
    i_q = 0;

    return dist_6D_index(i_r, i_phi, i_z,
                        i_vr, i_vphi, i_vz,
                        i_time, i_q,
                        dist->n_phi, dist->n_z,
                        dist->n_vr, dist->n_vphi,
                        dist->n_vz, 1, 1);
}