#include "bmc_diag.h"

void diag_move_distribution(sim_offload_data* sim, diag_data* diag_dest, diag_data* diag_src, int* updated, int* nloss, int* n_err) {
    // copy into diag_dest
    int dist_length = sim->diag_offload_data.offload_array_length;
    int updated_tmp, nloss_tmp, n_err_tmp;
    #ifdef MPI
        if (sim->diag_offload_data.dist5D_collect) {
            MPI_Allreduce(diag_src->dist5D.histogram, diag_dest->dist5D.histogram, dist_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(updated, &updated_tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(nloss, &nloss_tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(n_err, &n_err_tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
            memset(diag_src->dist5D.histogram, 0, sizeof(real)*dist_length);

            *updated = updated_tmp;
            *nloss = nloss_tmp;
            *n_err = n_err_tmp;
        } else if (sim->diag_offload_data.dist6D_collect) {
            MPI_Allreduce(diag_src->dist6D.histogram, diag_dest->dist6D.histogram, dist_length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            memset(diag_src->dist6D.histogram, 0, sizeof(real)*dist_length);
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

real fmc_compute_signal_from_states(
    int n_particles,
    particle_state* ps,
    int* n_updated,
    int* n_loss,
    int* n_err
) {
    *n_updated = 0;
    *n_loss = 0;
    *n_err = 0;
    real signal = 0;

    for(int i = 0; i < n_particles; i++) {

        if (ps[i].err) {
            *n_err = *n_err + 1;
            continue;
        }

        if ((ps[i].walltile > 0) && (bmc_walltile_in_target(ps[i].walltile))) {
            *n_updated = *n_updated + 1;
            signal += ps[i].weight;

        }
        else if (ps[i].walltile > 0) {
            *n_loss = *n_loss + 1;
        }
    }
    return signal;
}

int fmc_update_distr5D_from_states(
        dist_5D_data* dist1,
        int* p0_indexes,
        particle_state* p1,
        int n_particles,
        wall_2d_data* w2d,
        int *n_loss,
        int *n_err
    ) {

    int n_updated = 0;
    *n_loss = 0;
    *n_err = 0;

    for(int i = 0; i < n_particles; i++) {

        if (p1[i].err) {
            *n_err = *n_err + 1;
            continue;
        }

        if ((p1[i].walltile > 0) && (bmc_walltile_in_target(p1[i].walltile))) {
            n_updated++;
            dist1->histogram[p0_indexes[i]] += p1[i].hermite_weights;
            continue;
        }
        else if (p1[i].walltile > 0) {
            *n_loss = *n_loss + 1;
        }
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
        wall_data* wallData 
    ) {

    int n_updated = 0;

    for(int i = 0; i < n_simd_particles; i++) {

        #pragma omp simd
        for (int j=0; j<NSIMD; j++) {


            if ((p1[i].err[j]) || (p1[i].id[j] < 0)) {
                continue;
            }

            // int tile = wall_2d_hit_wall(p0[i].r[j], p0[i].phi[j], p0[i].z[j], p1[i].r[j], p1[i].phi[j], p1[i].z[j], w2d);
            if (p1[i].walltile[j] > 0) {
                n_updated++;
                if (bmc_walltile_in_target(p1[i].walltile[j])) {
                    dist1->histogram[p0_indexes[i*NSIMD + j]] += p1[i].hermite_weights[j];
                }
                continue;
            }

            // check if the particle escaped the velocity space
            real pperp = sqrt(2 * sqrt(p1[i].B_r[j]*p1[i].B_r[j]
                +p1[i].B_phi[j]*p1[i].B_phi[j]
                +p1[i].B_z[j]*p1[i].B_z[j])
                * p1[i].mu[j] / p1[i].mass[j]) * p1[i].mass[j];
                        
            if ((p1[i].ppar[j] > dist1->max_ppara) || (p1[i].ppar[j] < dist1->min_ppara) ||
                (pperp > dist1->max_pperp) || (pperp < dist1->min_pperp)) {

                // outside velocity space
                continue;
            }

            // deposit the final state with a linear interpolation on the closest spatial vertices
            // deposit in the velocity space is nearest neighbour 
            int p1_indexes[32];
            int p1_target_hit[32];
            real p1_weights[32];
            bmc_dist5D_gc_indexes(p1_indexes, p1_weights, p1_target_hit, p1 + i, j, dist0, wallData);
            for (int i_nodes=0; i_nodes<32; i_nodes++) {
                if (p1_indexes[i_nodes] >= 0) {
                    if (p1_target_hit[i_nodes] == 1) {
                        // particle hit the target domain. Set the relative probabiity to 1
                        dist1->histogram[p0_indexes[i*NSIMD + j]] += p1_weights[i_nodes] * p1[i].hermite_weights[j];
                    } else if (p1_target_hit[i_nodes] == 0) {
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
 int bmc_dist5D_gc_indexes(int* indexes, real* weights, int* target_hit, particle_simd_gc* p, int i_simd, dist_5D_data* dist, wall_data* wallData) {
    real phi;
    real pperp = sqrt(2 * sqrt(p->B_r[i_simd]*p->B_r[i_simd]
            +p->B_phi[i_simd]*p->B_phi[i_simd]
            +p->B_z[i_simd]*p->B_z[i_simd])
            * p->mu[i_simd] / p->mass[i_simd]) * p->mass[i_simd];

    return bmc_dist5D_gc_indexes_from_coordinates(indexes, weights, target_hit, p->r[i_simd], p->phi[i_simd], p->z[i_simd], p->ppar[i_simd], pperp, dist, wallData);

}

 int bmc_dist5D_gc_indexes_from_coordinates(int* indexes, real* weights, int* target_hit, real r, real phi, real z, real ppar, real pperp, dist_5D_data* dist, wall_data* wallData) {

    real weights_dim[5];
    int i_r;
    int i_phi;
    int i_z;
    int i_ppara;
    int i_pperp;
    int i_time;
    int i_q;

    weights_dim[0] = ((r - dist->min_r)
                / ((dist->max_r - dist->min_r)/(dist->n_r-1)));
    i_r = floor(weights_dim[0]);
    weights_dim[0] = 1. - weights_dim[0] + i_r;

    phi = fmod(phi, 2*CONST_PI);
    if(phi < 0) {
        phi = phi + 2*CONST_PI;
    }
    weights_dim[1] = ((phi - dist->min_phi)
                / ((dist->max_phi - dist->min_phi)/fmax(1, dist->n_phi-1)));
    i_phi = floor(weights_dim[1]);
    weights_dim[1] = 1. - weights_dim[1] + i_phi;

    weights_dim[2] = ((z - dist->min_z)
            / ((dist->max_z - dist->min_z) / (dist->n_z-1)));
    i_z = floor(weights_dim[2]);
    weights_dim[2] = 1. - weights_dim[2] + i_z;

    weights_dim[3] = (ppar - dist->min_ppara)
                / ((dist->max_ppara - dist->min_ppara) / (dist->n_ppara-1));
    i_ppara = floor(weights_dim[3]);
    weights_dim[3] = 1. - weights_dim[3] + i_ppara;

    weights_dim[4] = (pperp - dist->min_pperp)
                / ((dist->max_pperp - dist->min_pperp) / (dist->n_pperp-1));
    i_pperp = floor(weights_dim[4]);
    weights_dim[4] = 1. - weights_dim[4] + i_pperp;

    i_time = 0;

    i_q = 0;

    real dr = (dist->max_r - dist->min_r)/(dist->n_r-1);
    real dphi = (dist->max_phi - dist->min_phi)/fmax(1, dist->n_phi-1);
    real dz = (dist->max_z - dist->min_z)/(dist->n_z - 1);
    real dppara = (dist->max_ppara - dist->min_ppara)/(dist->n_ppara - 1);
    real dpperp = (dist->max_pperp - dist->min_pperp)/(dist->n_pperp - 1);

    int err = 0;
    int i = 0;
    int j_phimod;
    for (int j_r=i_r; j_r<i_r + 2; j_r++)
    for (int j_phi=i_phi; j_phi<i_phi + 2; j_phi++)
    for (int j_z=i_z; j_z<i_z + 2; j_z++)
    for (int j_ppara=i_ppara; j_ppara<i_ppara + 2; j_ppara++)
    for (int j_pperp=i_pperp; j_pperp<i_pperp + 2; j_pperp++) {

        if ((j_r < 0) || (j_phi < 0) || (j_z < 0) || (j_ppara < 0) || (j_pperp < 0)) {
            err = 1;
        }
        if ((j_r >= dist->n_r) || (j_z >= dist->n_z) || (j_ppara >= dist->n_ppara) || (j_pperp >= dist->n_pperp)) {
            err = 1;
        }

        if (!err) {
            j_phimod = j_phi;
            if (j_phimod>=dist->n_phi) {
                j_phimod = 0;
            }

            indexes[i] = dist_5D_index(j_r, j_phimod, j_z, j_ppara, j_pperp, i_time, i_q, dist->n_phi, dist->n_z, dist->n_ppara, dist->n_pperp, 1, 1);
            weights[i] = fabs(weights_dim[0] - j_r + i_r) * fabs(weights_dim[1] - j_phi + i_phi)
                        * fabs(weights_dim[2] - j_z + i_z) * fabs(weights_dim[3] - j_ppara + i_ppara)
                        * fabs(weights_dim[4] - j_pperp + i_pperp);
            target_hit[i] = bmc_wall_hit_target(r, (j_r)*dr + dist->min_r,
                            phi, (j_phimod)*dphi + dist->min_phi, z, (j_z)*dz + dist->min_z, wallData);
            i++;
        }
    }
    return err;
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
int bmc_dist5D_state_indexes(int* indexes, real* weights, int* target_hit, particle_state* ps, dist_5D_data* dist, wall_data* wdata) {

    real pperp = sqrt(2 * sqrt(ps->B_r*ps->B_r
            +ps->B_phi*ps->B_phi
            +ps->B_z*ps->B_z)
            * ps->mu / ps->mass) * ps->mass;

    return bmc_dist5D_gc_indexes_from_coordinates(indexes, weights, target_hit, ps->r, ps->phi, ps->z, ps->ppar, pperp, dist, wdata);

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

void compute_5d_indexes_from_hist_index(int i, int* i_x, dist_5D_offload_data* dist5D) {
    int i1 = i;
    i_x[4] = i1 % dist5D->n_pperp;
    i1 = i1 / dist5D->n_pperp;
    i_x[3] = i1 % dist5D->n_ppara;
    i1 = i1 / dist5D->n_ppara;
    i_x[2] = i1 % dist5D->n_z;
    i1 = i1 / dist5D->n_z;
    i_x[1] = i1 % dist5D->n_phi;
    i_x[0] = i1 / dist5D->n_phi;
}

void compute_5d_coordinates_from_hist_index(int i, int* i_x, real* r, real* phi, real* z, real* ppara, real* pperp, dist_5D_offload_data* dist5D) {
    compute_5d_indexes_from_hist_index(i, i_x, dist5D);
    compute_element_5d_coordinates(i_x, r, phi, z, ppara, pperp, dist5D);
}

void compute_element_5d_coordinates(int* i_x_new, real* r, real* phi, real* z, real* ppara, real* pperp, dist_5D_offload_data* dist) {
    *r = dist->min_r + (dist->max_r - dist->min_r) / (dist->n_r-1) * i_x_new[0];
    *phi = dist->min_phi + (dist->max_phi - dist->min_phi) / fmax(1, dist->n_phi - 1) * i_x_new[1];
    *z = dist->min_z + (dist->max_z - dist->min_z) / (dist->n_z-1) * i_x_new[2];
    *ppara = dist->min_ppara + (dist->max_ppara - dist->min_ppara) / (dist->n_ppara-1) * i_x_new[3];
    *pperp = dist->min_pperp + (dist->max_pperp - dist->min_pperp) / (dist->n_pperp - 1) * i_x_new[4];
}

void bmc_5D_to_particle_state(
        B_field_data* Bdata,
        real r, real phi, real z,
        real ppara, real pperp,
        real t,
        int id,
        particle_state* ps,
        real m,
        real q,
        int rk4_subcycles
    ) {

    a5err err = 0;

    real psi[1], rho[1], B_dB[15];
    err = B_field_eval_B_dB(B_dB, r, phi, z, t, Bdata);
    if (!err) {
        err = B_field_eval_psi(psi, r, phi, z,
                                t, Bdata);
    }
    if(!err) {
        err = B_field_eval_rho(rho, psi[0], Bdata);
    }

    if(!err) {
        ps->rho        = rho[0];
        ps->B_r        = B_dB[0];
        ps->B_phi      = B_dB[4];
        ps->B_z        = B_dB[8];
        ps->B_r_dr     = B_dB[1];
        ps->B_phi_dr   = B_dB[5];
        ps->B_z_dr     = B_dB[9];
        ps->B_r_dphi   = B_dB[2];
        ps->B_phi_dphi = B_dB[6];
        ps->B_z_dphi   = B_dB[10];
        ps->B_r_dz     = B_dB[3];
        ps->B_phi_dz   = B_dB[7];
        ps->B_z_dz     = B_dB[11];
    }

    /* Input is in (Ekin,xi) coordinates but state needs (mu,vpar) so we need to do that
        * transformation first. */
    real Bnorm = math_normc(B_dB[0], B_dB[4], B_dB[8]);
    real mu = 0.5 * pperp * pperp / m / Bnorm;
    if(!err && mu < 0)          {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        ps->n_t_subcycles = rk4_subcycles;
        ps->r        = r;
        ps->phi      = phi;
        ps->z        = z;
        ps->mu       = mu;
        ps->ppar     = ppara;
        ps->zeta     = 0;
        ps->mass     = m;
        ps->charge   = q;
        ps->anum     = 0;
        ps->znum     = 1;
        ps->weight   = 1;
        ps->time     = t;
        ps->mileage  = 0;
        ps->theta    = atan2(ps->z-B_field_get_axis_z(Bdata, ps->phi),
                                ps->r-B_field_get_axis_r(Bdata, ps->phi));
        ps->id       = id;
        ps->endcond  = 0;
        ps->walltile = 0;
        ps->cputime  = 0;
    }

    /* Guiding center transformation to get particle coordinates */
    real rprt, phiprt, zprt, pr, pphi, pz;
    if(!err) {
        real pparprt, muprt, zetaprt;
        gctransform_guidingcenter2particle(
            ps->mass, ps->charge, B_dB,
            ps->r, ps->phi, ps->z, ps->ppar, ps->mu, ps->zeta,
            &rprt, &phiprt, &zprt, &pparprt, &muprt, &zetaprt);

        B_field_eval_B_dB(B_dB, rprt, phiprt, zprt, ps->time, Bdata);

        gctransform_pparmuzeta2prpphipz(
            ps->mass, ps->charge, B_dB,
            phiprt, pparprt, muprt, zetaprt,
            &pr, &pphi, &pz);
    }
    if(!err && rprt <= 0) {err = error_raise(ERR_MARKER_UNPHYSICAL, __LINE__, EF_PARTICLE);}

    if(!err) {
        ps->rprt   = rprt;
        ps->phiprt = phiprt;
        ps->zprt   = zprt;
        ps->p_r    = pr;
        ps->p_phi  = pphi;
        ps->p_z    = pz;

        ps->err = 0;
    } else {
        ps->err = err;
    }

}

/**
 * Write the probability distribution to the output HDF file
 * 
 * @param sim_offload pointer to simulation offload data
 * @param distr pointer to distribution struct
 * @param distr_array pointer to offload array of the distribution
 * @param dist_length length of the distribution data
 * @param mpi_rank MPI node rank id
 * 
 **/
void write_probability_distribution(
    sim_offload_data* sim_offload,
    diag_data* distr,
    real* distr_array,
    int mpi_rank,
    int write_for_importance_sampling
) {

    FILE *fp;

    if (write_for_importance_sampling)
        fp = fopen("distr_prob", "wb");

    int dist_length = sim_offload->diag_offload_data.offload_array_length;

    if (write_for_importance_sampling)
        fwrite(distr->dist5D.histogram, sizeof(distr->dist5D.histogram), dist_length, fp);

    for (int i = 0; i<dist_length; i++) {
        // printf("%f\n", distr->dist5D.histogram[i]);
        if (sim_offload->diag_offload_data.dist5D_collect) {
            if (distr->dist5D.histogram[i] > 1.0001) {
                printf("Warning: unpysical probability: %f\n", distr->dist5D.histogram[i]);
            }
            if (distr->dist5D.histogram[i] < 0) {
               distr->dist5D.histogram[i] = 0; 
            }
        }
        if (sim_offload->diag_offload_data.dist6D_collect) {
            if (distr->dist6D.histogram[i] > 1.0001) {
                printf("Warning: unpysical probability: %f\n", distr->dist6D.histogram[i]);
            }
        }
    }
    if (write_for_importance_sampling)
        fclose(fp);

    /* Combine distribution and write it to HDF5 file */
    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nWriting BMC probability distribution.\n");
    if (mpi_rank == 0) {
        int err_writediag = 0;
        err_writediag = hdf5_interface_write_diagnostics(sim_offload, distr_array, sim_offload->hdf5_out);
        if(err_writediag) {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                    "\nWriting diagnostics failed.\n"
                    "See stderr for details.\n");
        }
        else {
            print_out0(VERBOSE_MINIMAL, mpi_rank,
                    "BMC distributions written.\n");
        }

    }
}

void bmc_update_distr5D_from_weights(
        particle_deposit_weights *p1_weights,
        dist_5D_data* dist1,
        dist_5D_data* dist0,
        particle_simd_gc* p1,
        int n_simd_particles,
        int* p0_indexes
    ) {

    int err;
    for(int i = 0; i < n_simd_particles; i++) {

        // #pragma omp simd
        for (int j=0; j<NSIMD; j++) {
            if (p1[i].id[j] < 0) {
                continue;
            }

            err = 0;
            int index;
            real weight;
            int target_hit;
            for (int k=0; k<32; k++) {
                if (err) break;

                index = p1_weights[i].index[j*32 + k];
                target_hit = p1_weights[i].target_hit[j*32 + k];
                weight = p1_weights[i].weight[j*32 + k];
                if (weight < 0) {
                    // if the weight is negative, the prob of the node couldn't be computed
                    // likely because it exited the velocity space
                    err = 1;

                    // just continue to the next marker.
                    // It might be that another marker in the same node will evolve to a valid location
                    continue;

                    // // mark the node as negative so the it can be caught for debug
                    // dist1->histogram[p0_indexes[i*NSIMD + j]] = weight;
                }
                else if (weight == 0) {
                    continue;
                }
                else if (target_hit == 2) {
                    // wall hit, no target
                    continue;
                }
                else if (target_hit == 1) {
                    // particle hit the target domain. Set the relative probabiity to 1
                    dist1->histogram[p0_indexes[i*NSIMD + j]] += weight * p1[i].hermite_weights[j];
                } else if (dist0->histogram[index] < 0) {
                    continue;
                    // dist1->histogram[p0_indexes[i * NSIMD + j]] = weight;
                    // err = 1;
                } else {
                    // particle didn't hit the target domain. Weight the probability with the last probability matrix
                    dist1->histogram[p0_indexes[i*NSIMD + j]] += dist0->histogram[index] * weight * p1[i].hermite_weights[j];
                }
            }
        }

    }
}

void bmc_compute_prob_weights(particle_deposit_weights *p1_weightsIndexes,
    int n_simd_particles, particle_simd_gc *p1, particle_simd_gc *p0,
    dist_5D_data* dist1, dist_5D_data* dist0, wall_data* wallData,
    int* p0_indexes, int debugExitVelocitySpace
    ) {

    for(int i = 0; i < n_simd_particles; i++) {

        #pragma omp simd
        for (int j=0; j<NSIMD; j++) {

            for (int k=0; k<32; k++) {
                p1_weightsIndexes[i].weight[32*j + k] = 0;
            }

            // check if the marker hit the wall
            // int tile = wall_2d_hit_wall(p0[i].r[j], p0[i].phi[j], p0[i].z[j], p1[i].r[j], p1[i].phi[j], p1[i].z[j], w2d);
            int tile = p1[i].walltile[j];
            if (tile > 0) {
                if (bmc_walltile_in_target(tile)) {
                    p1_weightsIndexes[i].weight[32*j] = 1.;
                    p1_weightsIndexes[i].target_hit[32*j] = 1;
                }
                continue;
            }

            if ((p1[i].err[j]) || (p1[i].id[j] < 0))
               continue; 

            // check if the particle escaped the velocity space
            real pperp = sqrt(2 * sqrt(p1[i].B_r[j]*p1[i].B_r[j]
                +p1[i].B_phi[j]*p1[i].B_phi[j]
                +p1[i].B_z[j]*p1[i].B_z[j])
                * p1[i].mu[j] / p1[i].mass[j]) * p1[i].mass[j];
                        
            if ((p1[i].ppar[j] > dist1->max_ppara) || (p1[i].ppar[j] < dist1->min_ppara) ||
                (pperp > dist1->max_pperp) || (pperp < dist1->min_pperp)) {

                // outside velocity space. Mark the weight as negative so it can be caught for debug purpose later
                if (debugExitVelocitySpace) {
                    printf("Warning: outside V: index: %d id: %d, ppar: %e, pperp %e\n", p0_indexes[i*NSIMD+j], p1[i].id[j], p1[i].ppar[j], pperp);
                }
                p1_weightsIndexes[i].weight[32 * j] = -1;
                continue; 
            }



            // deposit the final state j
            int p1_indexes[32];
            int p1_target_hit[32];
            real p1_weights[32];
            if (!bmc_dist5D_gc_indexes(p1_indexes, p1_weights, p1_target_hit, p1 + i, j, dist0, wallData)) {
                for (int k=0; k<32; k++) {
                    p1_weightsIndexes[i].weight[32*j + k] = p1_weights[k];
                    p1_weightsIndexes[i].index[32*j + k] = p1_indexes[k];
                    p1_weightsIndexes[i].target_hit[32*j + k] = p1_target_hit[k];
                }
            }
        }
    }
}

void particle_state_to_particle_gc(
    particle_state* ps,
    particle_gc* p_gc
) {
    real Brpz[3] = {ps->B_r, ps->B_phi, ps->B_z};
    real Bnorm   = math_norm(Brpz);
    real pin = physlib_gc_p( ps->mass, ps->mu, ps->ppar, Bnorm);
    p_gc->pitch = physlib_gc_xi(ps->mass, ps->mu, ps->ppar, Bnorm);
    p_gc->energy = physlib_Ekin_pnorm(ps->mass, pin);

    p_gc->r = ps->r;
    p_gc->phi = ps->phi;
    p_gc->z = ps->z;
    p_gc->zeta = ps->zeta;
    p_gc->mass = ps->mass;
    p_gc->charge = ps->charge;
    p_gc-> anum = ps->anum;
    p_gc-> znum = ps->znum;
    p_gc->weight = ps->weight;
    p_gc->time = ps->time;
    p_gc->id = ps->id;
}