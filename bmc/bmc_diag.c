#include "bmc_diag.h"

int bmc_walltile_in_target(integer walltile) {
    if ((walltile > 10) && (walltile < 30)) {
        return 1;
    }

    return 0;
}

void bmc_update_particles_diag(
    int n_mpi_particles,
    particle_state* ps0,
    particle_state* ps1,
    diag_data* diag0,
    diag_data* diag1,
    sim_data* sim,
    int n_montecarlo_steps
) {
    for(int i = 0; i < n_mpi_particles; i++) {

        if (sim->sim_mode == simulate_mode_fo) {
            bmc_diag_update_fo(&ps0[i], &ps1[i], diag0, diag1, n_montecarlo_steps);
        }
        else if (sim->sim_mode == simulate_mode_gc) {
            bmc_diag_update_gc(&ps0[i], &ps1[i], diag0, diag1, n_montecarlo_steps);
        }

    }
}

void bmc_diag_update_gc(
    particle_state* ps0,
    particle_state* ps1,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps
) {
    if(diag0->dist5D_collect) {
        bmc_diag_5D_update_gc(&diag1->dist5D, &diag0->dist5D,  ps1, ps0, n_montecarlo_steps);
    }

    if(diag0->dist6D_collect) {
        bmc_diag_6D_update_gc(&diag1->dist6D, &diag0->dist6D,  ps1, ps0, n_montecarlo_steps);
    }
}

void bmc_diag_update_fo(
    particle_state* ps0,
    particle_state* ps1,
    diag_data* diag0,
    diag_data* diag1,
    int n_montecarlo_steps
) {
    if(diag0->dist5D_collect) {
        bmc_diag_5D_update_fo(&diag1->dist5D, &diag0->dist5D,  ps1, ps0, n_montecarlo_steps);
    }

    if(diag0->dist6D_collect) {
        bmc_diag_6D_update_fo(&diag1->dist6D, &diag0->dist6D,  ps1, ps0, n_montecarlo_steps);
    }
}


void bmc_diag_5D_update_gc(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {
    int p0_index = bmc_dist5D_gc_index(ps0, dist0);

    // check if the particle hit the wall
    if (ps1->endcond & endcond_wall) {
        if (bmc_walltile_in_target(ps1->walltile)) {
            // particle ended in target domain, set the weight to 1
            dist1->histogram[p0_index] += 1. / n_montecarlo_steps;
        }
    } else {
        int p1_index = bmc_dist5D_gc_index(ps1, dist0);
        dist1->histogram[p0_index] += dist0->histogram[p1_index] / n_montecarlo_steps;
    }
}

void bmc_diag_5D_update_fo(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {
    int p0_index = bmc_dist5D_fo_index(ps0, dist0);

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
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {
    int p0_index = bmc_dist6D_gc_index(ps0, dist0);

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
        particle_state* ps1,
        particle_state* ps0,
        int n_montecarlo_steps
    ) {
    int p0_index = bmc_dist6D_fo_index(ps0, dist0);

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

int bmc_dist5D_gc_index(particle_state* ps, dist_5D_data* dist) {
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

    return dist_5D_index(i_r, i_phi, i_z,
                        i_vpara, i_vperp,
                        i_time, i_q,
                        dist->n_phi, dist->n_z,
                        dist->n_vpara, dist->n_vperp,
                        1, 1);
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