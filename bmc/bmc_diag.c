#include "bmc_diag.h"

void diag_move_distribution(sim_offload_data* sim, diag_data* diag_dest, diag_data* diag_src) {
    // copy into diag_dest
    int dist_length = sim->diag_offload_data.offload_array_length;
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

int bmc_update_distr5D_from_states(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int* p0_indexes,
        particle_state* p1,
        particle_state* p0,
        int n_particles,
        wall_2d_data* w2d
    ) {

    int n_updated = 0;

    for(int i = 0; i < n_particles; i++) {

        if (p1[i].err) {
            continue;
        }


        // check if the particle escaped the velocity space
        real vperp = sqrt(2 * sqrt(p1[i].B_r*p1[i].B_r
            +p1[i].B_phi*p1[i].B_phi
            +p1[i].B_z*p1[i].B_z)
            * p1[i].mu / p1[i].mass) * p1[i].mass;
                    
        if ((p1[i].ppar > dist1->max_ppara) || (p1[i].ppar < dist1->min_ppara) ||
            (vperp > dist1->max_pperp) || (vperp < dist1->min_pperp)) {

            // outside velocity space
            continue;
        }

        // deposit the final state with a linear interpolation on the closest spatial vertices
        // deposit in the velocity space is nearest neighbour 
        int p1_indexes[8];
        int p1_target_hit[8];
        real p1_weights[8];
        // bmc_dist5D_state_indexes(p0 + i, p1_indexes, p1_weights, p1_target_hit, p1 + i, j, dist0, w2d);
        bmc_dist5D_state_indexes(p0 + i, p1_indexes, p1_weights, p1_target_hit, p1 + i, dist1, w2d);
        for (int i_nodes=0; i_nodes<8; i_nodes++) {
            if (p1_indexes[i_nodes] >= 0) {
                if (p1_target_hit[i_nodes]) {
                    // particle hit the target domain. Set the relative probabiity to 1
                    dist1->histogram[p0_indexes[i]] += p1_weights[i_nodes] * p1[i].hermite_weights;
                } else {
                    // particle didn't hit the target domain. Weight the probability with the last probability matrix
                    dist1->histogram[p0_indexes[i]] += dist0->histogram[p1_indexes[i_nodes]] * p1_weights[i_nodes] * p1[i].hermite_weights;
                }
            }
        }

        n_updated++;
    }
    return n_updated;
}

int bmc_update_distr5D(
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        int* p0_indexes,
        particle_simd_gc* p1,
        particle_simd_gc* p0,
        int n_simd_particles,
        wall_2d_data* w2d
    ) {

    int n_updated = 0;

    for(int i = 0; i < n_simd_particles; i++) {

        #pragma omp simd
        for (int j=0; j<NSIMD; j++) {


            if ((p1[i].err[j]) || (p1[i].id[j] < 0)) {
                continue;
            }

            // check if the particle escaped the velocity space
            real pperp = sqrt(2 * sqrt(p1[i].B_r[j]*p1[i].B_r[j]
                +p1[i].B_phi[j]*p1[i].B_phi[j]
                +p1[i].B_z[j]*p1[i].B_z[j])
                * p1[i].mu[j] / p1[i].mass[j]) * p1[i].mass[j];
                        
            if ((p1[i].ppar[j] > dist1->max_ppara) || (p1[i].ppar[j] < dist1->min_ppara) ||
                (pperp > dist1->max_pperp) || (pperp < dist1->min_pperp)) {

                // printf("ppar %e pperp %e\n", p1[i].ppar[j], pperp);

                // outside velocity space
                continue;
            }

            // deposit the final state with a linear interpolation on the closest spatial vertices
            // deposit in the velocity space is nearest neighbour 
            int p1_indexes[32];
            int p1_target_hit[32];
            real p1_weights[32];
            bmc_dist5D_gc_indexes(p0 + i, p1_indexes, p1_weights, p1_target_hit, p1 + i, j, dist0, w2d);
            for (int i_nodes=0; i_nodes<32; i_nodes++) {
                if (p1_indexes[i_nodes] >= 0) {
                    if (p1_target_hit[i_nodes]) {
                        // particle hit the target domain. Set the relative probabiity to 1
                        dist1->histogram[p0_indexes[i*NSIMD + j]] += p1_weights[i_nodes] * p1[i].hermite_weights[j];
                    } else {
                        // particle didn't hit the target domain. Weight the probability with the last probability matrix
                        dist1->histogram[p0_indexes[i*NSIMD + j]] += dist0->histogram[p1_indexes[i_nodes]] * p1_weights[i_nodes] * p1[i].hermite_weights[j];
                    }
                }
            }

            n_updated++;
        }

    }
    return n_updated;
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
 void bmc_dist5D_gc_indexes(particle_simd_gc* p0, int* indexes, real* weights, int* target_hit, particle_simd_gc* p, int i_simd, dist_5D_data* dist, wall_2d_data* w2d) {
    real phi;
    real ppara;
    real pperp;
    real weights_dim[5];
    int i_r;
    int i_phi;
    int i_z;
    int i_ppara;
    int i_pperp;
    int i_time;
    int i_q;

    weights_dim[0] = ((p->r[i_simd] - dist->min_r)
                / ((dist->max_r - dist->min_r)/dist->n_r));
    i_r = floor(weights_dim[0]);
    weights_dim[0] = 1. - weights_dim[0] + i_r;

    phi = fmod(p->phi[i_simd], 2*CONST_PI);
    if(phi < 0) {
        phi = phi + 2*CONST_PI;
    }
    weights_dim[1] = ((phi - dist->min_phi)
                / ((dist->max_phi - dist->min_phi)/dist->n_phi));
    i_phi = floor(weights_dim[1]);
    weights_dim[1] = 1. - weights_dim[1] + i_phi;

    weights_dim[2] = ((p->z[i_simd] - dist->min_z)
            / ((dist->max_z - dist->min_z) / dist->n_z));
    i_z = floor(weights_dim[2]);
    weights_dim[2] = 1. - weights_dim[2] + i_z;

    weights_dim[3] = (p->ppar[i_simd] - dist->min_ppara)
                / ((dist->max_ppara - dist->min_ppara) / dist->n_ppara);
    i_ppara = floor(weights_dim[3]);
    weights_dim[3] = 1. - weights_dim[3] + i_ppara;

    pperp = sqrt(2 * sqrt(p->B_r[i_simd]*p->B_r[i_simd]
                                +p->B_phi[i_simd]*p->B_phi[i_simd]
                                +p->B_z[i_simd]*p->B_z[i_simd])
                    * p->mu[i_simd] / p->mass[i_simd]) * p->mass[i_simd];
    weights_dim[4] = (pperp - dist->min_pperp)
                / ((dist->max_pperp - dist->min_pperp) / dist->n_pperp);
    i_pperp = floor(weights_dim[4]);
    weights_dim[4] = 1. - weights_dim[4] + i_pperp;
    

    i_time = 0;

    i_q = 0;

    real dr = (dist->max_r - dist->min_r)/dist->n_r;
    real dphi = (dist->max_phi - dist->min_phi)/dist->n_phi;
    real dz = (dist->max_z - dist->min_z)/dist->n_z;
    real dppara = (dist->max_ppara - dist->min_ppara)/dist->n_ppara;
    real dpperp = (dist->max_pperp - dist->min_pperp)/dist->n_pperp;

    int i = 0;
    real r1, phi1, z1, ppara1, pperp1;
    int j_phimod;
    for (int j_r=i_r; j_r<i_r + 2; j_r++)
    for (int j_phi=i_phi; j_phi<i_phi + 2; j_phi++)
    for (int j_z=i_z; j_z<i_z + 2; j_z++)
    for (int j_ppara=i_ppara; j_ppara<i_ppara + 2; j_ppara++)
    for (int j_pperp=i_pperp; j_pperp<i_pperp + 2; j_pperp++) {

        if ((j_r < 0) || (j_phi < 0) || (j_z < 0) || (j_ppara < 0) || (j_pperp < 0)) {
            indexes[i] = 0;
            weights[i] = 0;
            target_hit[i] = 0;
            i++;
            continue;
        }

        j_phimod = j_phi;
        if (j_phimod>=dist->n_phi) {
            j_phimod = 0;
        }
        r1 = j_r*dr + dist->min_r;
        phi1 = j_phimod*dphi + dist->min_phi;
        z1 = j_z*dz + dist->min_z;
        ppara1 = j_ppara*dppara + dist->min_ppara;
        pperp1 = j_pperp*dpperp + dist->min_pperp;
        indexes[i] = dist_5D_index(j_r, j_phimod, j_z, j_ppara, j_pperp, i_time, i_q,
                    dist->n_phi, dist->n_z, dist->n_ppara, dist->n_pperp, 1, 1);
        weights[i] = fabs(weights_dim[0] - j_r + i_r) * fabs(weights_dim[1] - j_phi + i_phi)
                    * fabs(weights_dim[2] - j_z + i_z) * fabs(weights_dim[3] - j_ppara + i_ppara)
                    * fabs(weights_dim[4] - j_pperp + i_pperp);
        target_hit[i] = bmc_wall_2d_hit_target(p0->r[i_simd], (j_r)*dr + dist->min_r,
                        p0->phi[i_simd], (j_phimod)*dphi + dist->min_phi, p0->z[i_simd], (j_z)*dz + dist->min_z, w2d);
        i++;
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
 void bmc_dist5D_state_indexes(particle_state* ps0, int* indexes, real* weights, int* target_hit, particle_state* ps, dist_5D_data* dist, wall_2d_data* w2d) {
    real phi;
    real ppara;
    real pperp;
    real weights_rphiz[3];
    int i_r;
    int i_phi;
    int i_z;
    int i_ppara;
    int i_pperp;
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

    i_ppara = floor((ps->ppar - dist->min_ppara)
                / ((dist->max_ppara - dist->min_ppara) / dist->n_ppara));

    pperp = sqrt(2 * sqrt(ps->B_r*ps->B_r
                                +ps->B_phi*ps->B_phi
                                +ps->B_z*ps->B_z)
                    * ps->mu / ps->mass) * ps->mass;
    i_pperp = floor((pperp - dist->min_pperp)
                / ((dist->max_pperp - dist->min_pperp) / dist->n_pperp));

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
        indexes[i] = dist_5D_index(j_r, j_phimod, j_z, i_ppara, i_pperp, i_time, i_q,
                    dist->n_phi, dist->n_z, dist->n_ppara, dist->n_pperp, 1, 1);
        weights[i] = fabs(weights_rphiz[0] - j_r + i_r) * fabs(weights_rphiz[1] - j_phi + i_phi) * fabs(weights_rphiz[2] - j_z + i_z);
        target_hit[i] = bmc_wall_2d_hit_target(ps0->r, (j_r)*dr + dist->min_r,
                        ps0->phi, (j_phimod)*dphi + dist->min_phi, ps0->z, (j_z)*dz + dist->min_z, w2d);
        i++;
    }
}

int bmc_dist6D_fo_index(particle_state* ps, dist_6D_data* dist) {
    // real phi;
    // int i_r;
    // int i_phi;
    // int i_z;
    // int i_vr;
    // int i_vphi;
    // int i_vz;
    // int i_time;
    // int i_q;

    // i_r = floor((ps->r - dist->min_r)
    //             / ((dist->max_r - dist->min_r)/dist->n_r));

    // phi = fmod(ps->phi, 2*CONST_PI);
    // if(phi < 0) {
    //     phi = phi + 2*CONST_PI;
    // }
    // i_phi = floor((phi - dist->min_phi)
    //             / ((dist->max_phi - dist->min_phi)/dist->n_phi));

    // i_z = floor((ps->z - dist->min_z)
    //             / ((dist->max_z - dist->min_z) / dist->n_z));

    // i_vr = floor((ps->rdot - dist->min_vr)
    //             / ((dist->max_vr - dist->min_vr) / dist->n_vr));

    // i_vphi = floor((ps->phidot*ps->r - dist->min_vphi)
    //             / ((dist->max_vphi - dist->min_vphi) / dist->n_vphi));

    // i_vz = floor((ps->zdot - dist->min_vz)
    //             / ((dist->max_vz - dist->min_vz) / dist->n_vz));

    // i_time = 0;
    // i_q = 0;

    // return dist_6D_index(i_r, i_phi, i_z,
    //                     i_vr, i_vphi, i_vz,
    //                     i_time, i_q,
    //                     dist->n_phi, dist->n_z,
    //                     dist->n_vr, dist->n_vphi,
    //                     dist->n_vz, 1, 1);
}