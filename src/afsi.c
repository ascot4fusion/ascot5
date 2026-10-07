/**
 * @file afsi.c
 * @brief ASCOT Fusion Source Integrator AFSI
 */
#include <string.h>
#include <math.h>
#include <hdf5_hl.h>
#include "ascot5.h"
#include "print.h"
#include "gitver.h"
#include "math.h"
#include "physlib.h"
#include "consts.h"
#include "offload.h"
#include "random.h"
#include "simulate.h"
#include "boschhale.h"
#include "diag/hist.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_dist.h"
#include "hdf5io/hdf5_histogram.h"
#include "afsi.h"

/** Random number generator used by AFSI */
random_data rdata;
void B_field_cyl_to_cartesian(real* B_cart, real* B_cyl, real r, real phi);
void afsi_compute_product_velocities_3d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* B_cyl, int not_transformed_vel, real* vprod1, real* vprod2);
// void afsi_compute_product_momenta_2d(
//     int i, real mprod1, real mprod2, int prodmomspace,
//     real* vprod1, real* vprod2, real* prod1_p1, real* prod1_p2,
//     real* prod2_p1, real* prod2_p2);
void afsi_compute_product_momenta_2d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q, int prodmomspace,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* prod1_p1, real* prod1_p2, real* prod2_p1, real* prod2_p2);
void afsi_store_particle_data(int i, real r, real phi, real z, 
    real* vprod2, real mprod2, real* prod2, int position_space_id, int velocity_space_id, real weight);
void afsi_sample_reactant_momenta_2d(
    sim_data* sim, afsi_data* afsi, real mass1, real mass2,
    real vol, int nsample, size_t i0, size_t i1, size_t i2,
    real r, real phi, real z, real time, real rho,
    real* density1, real* ppara1, real* pperp1,
    real* density2, real* ppara2, real* pperp2);
void afsi_sample_beam_2d(
    histogram* hist, real mass, real vol, int nsample,
    size_t i0, size_t i1, size_t i2, real* density,
    real* ppara, real* pperp);
void afsi_sample_thermal_2d(
    sim_data* sim, int ispecies, real mass, int nsample,
    real r, real phi, real z, real time, real rho,
    real* density, real* pppara, real* ppperp);
void afsi_sample_beam_2d(histogram* hist, real mass, real vol, int nsample,
                         size_t i0, size_t i1, size_t i2, real* density, real* ppara,
                         real* pperp);
void afsi_sample_thermal_2d(sim_data* sim, int ispecies, real mass, int nsample,
                            real r, real phi, real z, real time, real rho,
                            real* density, real* pppara, real* ppperp);
void afsi_sample_reactant_momenta_2d_alt(
    sim_data* sim, afsi_data* afsi, real mass1, real mass2, real vol,
    int nsample, size_t i0, size_t i1, size_t i2,
    real r, real phi, real z, real time, real rho, real* cumdist,
    real* density1, real* ppara1, real* pperp1,
    real* density2, real* ppara2, real* pperp2);
void afsi_sample_beam_2d_alt(histogram* hist, real mass, real vol, int nsample,
                         size_t i0, size_t i1, size_t i2, real* cumdist,
                         real* density, real* ppara, real* pperp);

/**
 * @brief Calculate fusion source from two arbitrary ion distributions
 *
 * Inputs and outputs are expected to have same physical (R, phi, z)
 * dimensions.
 *
 * @param sim pointer to simulation data
 * @param reaction fusion reaction type, see the description
 * @param n number of Monte Carlo samples to be used
 * @param react1 reactant 1 distribution data
 * @param react2 reactant 2 distribution data
 * @param mult factor by which the output is scaled
 * @param prod1_data distribution data for product 1 output
 * @param prod2_data distribution data for product 2 output
 */
void afsi_run(sim_data* sim, afsi_data* afsi, int n,
              histogram* prod1, histogram* prod2) {
    /* QID for this run */
    char qid[11];
    hdf5_generate_qid(qid);
    strcpy(sim->qid, qid);

    int mpi_rank = 0, mpi_root = 0; /* AFSI does not support MPI */
    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root, "AFSI5\n");
    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);

    random_init(&rdata, time((NULL)));
    strcpy(sim->hdf5_out, sim->hdf5_in);
    simulate_init(sim);

    if( hdf5_interface_init_results(sim, qid, "afsi") ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root,
                   "\nInitializing output failed.\n"
                   "See stderr for details.\n");
        /* Free data and terminate */
        abort();
    }

    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;
    boschhale_reaction(
        afsi->reaction, &m1, &q1, &m2, &q2,
        &mprod1, &qprod1, &mprod2, &qprod2, &Q);

    int prod_mom_space;
    int i0coord, i1coord, i2coord, p1coord, p2coord;
    if(prod1->axes[0].n) {
        i0coord = 0;
        i1coord = 1;
        i2coord = 2;
    }
    else if(prod1->axes[3].n) {
        i0coord = 3;
        i1coord = 4;
        i2coord = 1;
    }
    else {
        return;
    }
    if(prod1->axes[5].n) {
        p1coord = 5;
        p2coord = 6;
        prod_mom_space = PPARPPERP;
    }
    else if(prod1->axes[10].n) {
        p1coord = 10;
        p2coord = 11;
        prod_mom_space = EKINXI;
    }
    else {
        return;
    }

    real time = 0.0;
    int fail_count = 0;
    /*
    for (int ax = 0; ax < HIST_ALLDIM; ax++) {
        printf("Axis %d: n=%ld min=%e max=%e\n",
           ax,
           prod2->axes[ax].n,
           prod2->axes[ax].min,
           prod2->axes[ax].max);
    } 
    
    printf("Mult: %e\n", afsi->mult);*/
    #pragma omp parallel for
    for(size_t i0 = 0; i0 < afsi->volshape[0]; i0++) {
        real* ppara1 = (real*) malloc(n*sizeof(real));
        real* pperp1 = (real*) malloc(n*sizeof(real));
        real* ppara2 = (real*) malloc(n*sizeof(real));
        real* pperp2 = (real*) malloc(n*sizeof(real));
        real* prod1_p1 = (real*) malloc(n*sizeof(real));
        real* prod1_p2 = (real*) malloc(n*sizeof(real));
        real* prod2_p1 = (real*) malloc(n*sizeof(real));
        real* prod2_p2 = (real*) malloc(n*sizeof(real));

        for(size_t i1 = 0; i1 < afsi->volshape[1]; i1++) {
            for(size_t i2 = 0; i2 < afsi->volshape[2]; i2++) {
                size_t spatial_index = i0*afsi->volshape[1]*afsi->volshape[2]
                                     + i1*afsi->volshape[2] + i2;
                real r = afsi->r[spatial_index];
                real z = afsi->z[spatial_index];
                real phi = afsi->phi[spatial_index];
                real vol = afsi->vol[spatial_index];

                real psi, rho[2];
                if(B_field_eval_psi(&psi, r, phi, z, time, &sim->B_data) ||
                   B_field_eval_rho(rho, psi, &sim->B_data) ) {
                    continue;
                }

                real density1, density2;
                afsi_sample_reactant_momenta_2d(
                    sim, afsi, m1, m2, vol, n, i0, i1, i2,
                    r, phi, z, time, rho[0],
                    &density1, ppara1, pperp1, &density2, ppara2, pperp2);
                if(density1 == 0 || density2 == 0) {
                    continue;
                }
                for(size_t i = 0; i < n; i++) {
                    real vcom2;
                    afsi_compute_product_momenta_2d(
                        i, m1, m2, mprod1, mprod2, Q, prod_mom_space,
                        ppara1, pperp1, ppara2, pperp2, &vcom2,
                        prod1_p1, prod1_p2, prod2_p1, prod2_p2);
                    real E = 0.5 * ( m1 * m2 ) / ( m1 + m2 ) * vcom2;

                    real weight = density1 * density2 * sqrt(vcom2)
                        * boschhale_sigma(afsi->reaction, E)/n*vol;

                    size_t ip1 = math_bin_index(
                        prod1_p1[i], prod1->axes[p1coord].n,
                        prod1->axes[p1coord].min, prod1->axes[p1coord].max);
                    size_t ip2 = math_bin_index(
                        prod1_p2[i], prod1->axes[p2coord].n,
                        prod1->axes[p2coord].min, prod1->axes[p2coord].max);
                    if( 0 <= ip1 && ip1 < prod1->axes[p1coord].n &&
                        0 <= ip2 && ip2 < prod1->axes[p2coord].n) {
                        size_t index =    i0*prod1->strides[i0coord]
                                        + i1*prod1->strides[i1coord]
                                        + i2*prod1->strides[i2coord]
                                        + ip1*prod1->strides[p1coord]
                                        + ip2*prod1->strides[p2coord];
                        prod1->bins[index] += weight * afsi->mult;
                    }

                    ip1 = math_bin_index(
                        prod2_p1[i], prod2->axes[p1coord].n,
                        prod2->axes[p1coord].min, prod2->axes[p1coord].max);
                    ip2 = math_bin_index(
                        prod2_p2[i], prod2->axes[p2coord].n,
                        prod2->axes[p2coord].min, prod2->axes[p2coord].max);
                    if( 0 <= ip1 && ip1 < prod2->axes[p1coord].n &&
                        0 <= ip2 && ip2 < prod2->axes[p2coord].n) {
                        size_t index =    i0*prod2->strides[i0coord]
                                        + i1*prod2->strides[i1coord]
                                        + i2*prod2->strides[i2coord]
                                        + ip1*prod2->strides[p1coord]
                                        + ip2*prod2->strides[p2coord];
                        prod2->bins[index] += weight * afsi->mult;
                    }
                    else {
                        printf("BIN FAIL: prod2_p1=%f prod2_p2=%f | ip1=%ld ip2=%ld\n",
                        prod2_p1[i], prod2_p2[i], ip1, ip2);
                        fail_count++;
                    }
                }
            }
        }
        free(ppara1);
        free(ppara2);
        free(pperp1);
        free(pperp2);
        free(prod1_p1);
        free(prod1_p2);
        free(prod2_p1);
        free(prod2_p2);
    }
    printf("Number of samples that did not fit in the histogram: %d\n", fail_count);

    m1     = m1 / CONST_U;
    m2     = m2 / CONST_U;
    mprod1 = mprod1 / CONST_U;
    mprod2 = mprod2 / CONST_U;
    int c1     = (int)rint(q1 / CONST_E);
    int c2     = (int)rint(q2 / CONST_E);
    int cprod1 = (int)rint(qprod1 / CONST_E);
    int cprod2 = (int)rint(qprod2 / CONST_E);

    hid_t f = hdf5_open(sim->hdf5_out);
    if(f < 0) {
        print_err("Error: File not found.\n");
        abort();
    }
    char path[300];
    sprintf(path, "/results/afsi_%s/reaction", sim->qid);
    hid_t reactiondata = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT,
                                    H5P_DEFAULT);
    if(reactiondata < 0) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    hsize_t size = 1;
    real q = Q / CONST_E;
    if(H5LTmake_dataset_double(reactiondata, "m1", 1, &size, &m1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "m2", 1, &size, &m2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "mprod1", 1, &size, &mprod1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "mprod2", 1, &size, &mprod2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "q", 1, &size, &q)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "q1", 1, &size, &c1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "q2", 1, &size, &c2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "qprod1", 1, &size, &cprod1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "qprod2", 1, &size, &cprod2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5Gclose(reactiondata)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }

    sprintf(path, "/results/afsi_%s/prod1dist5d", sim->qid);
    if( hdf5_hist_write(f, path, prod1) ) {
        print_err("Warning: 5D distribution could not be written.\n");
    }
    sprintf(path, "/results/afsi_%s/prod2dist5d", sim->qid);
    if( hdf5_hist_write(f, path, prod2) ) {
        print_err("Warning: 5D distribution could not be written.\n");
    }
    if(hdf5_close(f)) {
        print_err("Failed to close the file.\n");
        abort();
    }

    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root, "\nDone\n");
}

size_t hist_index(const histogram *hist, const size_t bins[])
{
    size_t index = 0;

    for (int ax = 0; ax < HIST_ALLDIM; ax++) {
        if (hist->axes[ax].n)
            index += bins[ax] * hist->strides[ax];
    }

    return index;
}

void afsi_run_6d(sim_data* sim, afsi_data* afsi, int n,
                 histogram* prod2){
    /* QID for this run */
    char qid[11];
    hdf5_generate_qid(qid);
    strcpy(sim->qid, qid);

    int mpi_rank = 0, mpi_root = 0; /* AFSI does not support MPI */
    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root, "AFSI5\n");
    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);

    random_init(&rdata, time((NULL)));
    strcpy(sim->hdf5_out, sim->hdf5_in);
    simulate_init(sim);

    if( hdf5_interface_init_results(sim, qid, "afsi") ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root,
                   "\nInitializing output failed.\n"
                   "See stderr for details.\n");
        /* Free data and terminate */
        abort();
    }
    
    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;

    boschhale_reaction(
        afsi->reaction,
        &m1, &q1, &m2, &q2,
        &mprod1, &qprod1, &mprod2, &qprod2, &Q
    );

    real time = 0.0;
    // printf("Mult: %e\n", afsi->mult);
    for (int ax = 0; ax < HIST_ALLDIM; ax++) {
        if (ax < HIST_ALLDIM - 1) {
            printf("Axis %2d: n=%3zu min=% .6e max=% .6e stride=%zu\n",
               ax,
               prod2->axes[ax].n,
               prod2->axes[ax].min,
               prod2->axes[ax].max,
               prod2->strides[ax]);
        }
        else {
            printf("Axis %2d: n=%3zu min=% .6e max=% .6e stride=LAST\n",
               ax,
               prod2->axes[ax].n,
               prod2->axes[ax].min,
               prod2->axes[ax].max);
        }
    }
    /* --- Ekin,theta_v,phi_v distribution ---*/
    int prod_mom_space;
    int pax0, pax1, pax2;
    if(prod2->axes[10].n) {
        pax0 = 10;
        pax1 = 16;
        pax2 = 17;
        prod_mom_space = ESPH;
        printf("Using ESPH\n");
    }
    /* --- pr,pphi,pz distribution ---*/
    else if(prod2->axes[7].n) {
        pax0 = 7;
        pax1 = 8;
        pax2 = 9;
        prod_mom_space = PCYL;
        printf("Using PCYL\n");
    }
    else {
        printf("Wrong momentum axes\n");
        return;
    }
    int fail_count = 0;
    double total_weight = 0;
    #pragma omp parallel for
    for (size_t i0 = 0; i0 < afsi->volshape[0]; i0++)
    {
        real* ppara1 = malloc(n*sizeof(real));
        real* pperp1 = malloc(n*sizeof(real));
        real* ppara2 = malloc(n*sizeof(real));
        real* pperp2 = malloc(n*sizeof(real));

        for (size_t i1 = 0; i1 < afsi->volshape[1]; i1++)
        {
            for (size_t i2 = 0; i2 < afsi->volshape[2]; i2++)
            {
                size_t spatial_index =
                    i0 * afsi->volshape[1] * afsi->volshape[2] +
                    i1 * afsi->volshape[2] + i2;

                real r   = afsi->r[spatial_index];
                real phi = afsi->phi[spatial_index];
                real z   = afsi->z[spatial_index];
                real vol = afsi->vol[spatial_index];

                /* --- Magnetic field --- */
                real psi, rho[2], B_cyl[3], B_cart[3];

                if (B_field_eval_psi(&psi, r, phi, z, time, &sim->B_data) ||
                    B_field_eval_rho(rho, psi, &sim->B_data) ||
                    B_field_eval_B(B_cyl, r, phi, z, time, &sim->B_data))
                    continue;

                B_field_cyl_to_cartesian(B_cart, B_cyl, r, phi);

                /* --- Sample reactants --- */
                real density1, density2;

                afsi_sample_reactant_momenta_2d(
                    sim, afsi, m1, m2, vol, n,
                    i0, i1, i2,
                    r, phi, z, time, rho[0],
                    &density1, ppara1, pperp1,
                    &density2, ppara2, pperp2
                );

                if (density1 == 0 || density2 == 0)
                    continue;
                for (size_t i = 0; i < n; i++)
                {
                    real vprod1[3], vprod2[3], vcom2;

                    /* --- Compute product velocities --- */
                    afsi_compute_product_velocities_3d(
                        i,
                        m1, m2, mprod1, mprod2, Q,
                        ppara1, pperp1,
                        ppara2, pperp2,
                        &vcom2,
                        B_cart,
                        0,
                        vprod1, vprod2
                    );

                    /* --- Convert to (pr, pphi, pz) --- */
                    real vx = vprod2[0];
                    real vy = vprod2[1];
                    real vz = vprod2[2];
                    
                    real v2 = vx*vx + vy*vy + vz*vz;
                    real v  = sqrt(v2);

                    /* --- Convert to (E, theta, phi) --- */
                    if (v == 0.0)
                        continue;
                    
                    real p0, p1, p2;
                    if (prod_mom_space == ESPH){
                        // Ekin (classical)
                        p0 = 0.5 * mprod2 * v2;

                        real cos_theta = vz / v;

                        /* Clamp to valid domain of acos */
                        if (cos_theta > 1.0)  cos_theta = 1.0;
                        if (cos_theta < -1.0) cos_theta = -1.0;
                        
                        // Theta_v
                        p1 = acos(cos_theta);
                        
                        // Phi_v
                        p2 = atan2(vy, vx); // Could be alternative real phi_v = atan2(vphi,vr);
                    }
                    else if (prod_mom_space == PCYL){
                        /* Cylindrical velocity */
                        real vr   =  vx*cos(phi) + vy*sin(phi);
                        real vphi = -vx*sin(phi) + vy*cos(phi);

                        /* Momentum pr, pphi, pz */ 
                        p0   = mprod2 * vr;
                        p1 = mprod2 * vphi;
                        p2   = mprod2 * vz;
                    }
                    else {
                        printf("Wrong momentum axes\n");
                        continue;
                    }

                    /* --- Reaction weight --- */
                    real Ecom = 0.5 * (m1 * m2) / (m1 + m2) * vcom2;

                    real weight =
                        density1 * density2 * sqrt(vcom2) *
                        boschhale_sigma(afsi->reaction, Ecom) / n * vol;

                    /* --- Bin indices --- */
                    size_t ip0 = math_bin_index(
                        p0,
                        prod2->axes[pax0].n,
                        prod2->axes[pax0].min,
                        prod2->axes[pax0].max
                    );

                    size_t ip1 = math_bin_index(
                        p1,
                        prod2->axes[pax1].n,
                        prod2->axes[pax1].min,
                        prod2->axes[pax1].max
                    );

                    size_t ip2 = math_bin_index(
                        p2,
                        prod2->axes[pax2].n,
                        prod2->axes[pax2].min,
                        prod2->axes[pax2].max
                    );
                    
                    if (ip0 < prod2->axes[pax0].n &&
                        ip1 < prod2->axes[pax1].n &&
                        ip2 < prod2->axes[pax2].n)
                    {
                        
                        size_t index =
                              i0 * prod2->strides[0]
                            + i1 * prod2->strides[1]
                            + i2 * prod2->strides[2]
                            + ip0 * prod2->strides[pax0]
                            + ip1 * prod2->strides[pax1]
                            + ip2;
                        #pragma omp atomic
                        prod2->bins[index] += weight * afsi->mult;
                        #pragma omp atomic
                        total_weight += weight* afsi->mult;
                    }
                    else {
                        printf("BIN FAIL: p0=%e p1=%e p2=%e | ip0=%ld ip1=%ld ip2=%ld\n",
                        p0, p1, p2,
                        ip0, ip1, ip2);
                        fail_count++;
                    }
                }
            }
        }
        free(ppara1);
        free(ppara2);
        free(pperp1);
        free(pperp2);
    }
    printf("Fail count: %d", fail_count);
    double hist_sum = 0.0;

    for (size_t i = 0; i < prod2->nbin; i++)
        hist_sum += prod2->bins[i];

    printf("Total sampled weight = %.15e\n", total_weight);
    printf("Histogram sum        = %.15e\n", hist_sum);
    printf("Difference           = %.15e\n",
       hist_sum - total_weight);
    m1     = m1 / CONST_U;
    m2     = m2 / CONST_U;
    mprod1 = mprod1 / CONST_U;
    mprod2 = mprod2 / CONST_U;
    int c1     = (int)rint(q1 / CONST_E);
    int c2     = (int)rint(q2 / CONST_E);
    int cprod1 = (int)rint(qprod1 / CONST_E);
    int cprod2 = (int)rint(qprod2 / CONST_E);

    hid_t f = hdf5_open(sim->hdf5_out);
    if(f < 0) {
        print_err("Error: File not found.\n");
        abort();
    }
    char path[300];
    sprintf(path, "/results/afsi_%s/reaction", sim->qid);
    hid_t reactiondata = H5Gcreate2(f, path, H5P_DEFAULT, H5P_DEFAULT,
                                    H5P_DEFAULT);
    if(reactiondata < 0) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    hsize_t size = 1;
    real q = Q / CONST_E;
    if(H5LTmake_dataset_double(reactiondata, "m1", 1, &size, &m1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "m2", 1, &size, &m2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "mprod1", 1, &size, &mprod1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "mprod2", 1, &size, &mprod2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_double(reactiondata, "q", 1, &size, &q)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "q1", 1, &size, &c1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "q2", 1, &size, &c2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "qprod1", 1, &size, &cprod1)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5LTmake_dataset_int(reactiondata, "qprod2", 1, &size, &cprod2)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }
    if(H5Gclose(reactiondata)) {
        print_err("Failed to write reaction data.\n");
        abort();
    }

    // sprintf(path, "/results/afsi_%s/prod1dist5d", sim->qid);
    // if( hdf5_hist_write(f, path, prod1) ) {
    //     print_err("Warning: 5D distribution could not be written.\n");
    // }
    sprintf(path, "/results/afsi_%s/prod2dist5d", sim->qid);
    if( hdf5_hist_write(f, path, prod2) ) {
        print_err("Warning: 5D distribution could not be written.\n");
    }
    if(hdf5_close(f)) {
        print_err("Failed to close the file.\n");
        abort();
    }

    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root, "\nDone\n");
}


/**
 * @brief Calculate fusion source from two arbitrary ion distributions using rejection sampling.
 *
 * Stores in prod2 the sampled neutrons in a form suitable for the use in neutronic Monte Carlo codes.
 *
 * @param sim pointer to simulation data
 * @param reaction fusion reaction type, see the description
 * @param n Total number of Monte Carlo samples required for output
 * @param react1 reactant 1 distribution data
 * @param react2 reactant 2 distribution data
 * @param mult factor by which the output is scaled
 * @param prod2 (n*9) array to store the sampled neutrons.
 */

void afsi_run_rejection(sim_data* sim, afsi_data* afsi, int n,
                int position_space_id, int velocity_space_id,
                real Smax, real* cumdist_all, real* prod2) {

    random_init(&rdata, time((NULL)));
    simulate_init(sim);

    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;

    boschhale_reaction(
        afsi->reaction, 
        &m1, &q1, &m2, &q2,
        &mprod1, &qprod1, &mprod2, &qprod2, &Q
    );
    
    real time = 0.0;

    real rmin = afsi->r[0], rmax = afsi->r[1];
    real phimin = afsi->phi[0], phimax = afsi->phi[1];
    real zmin = afsi->z[0], zmax = afsi->z[1];

    int n_accepted = 0 ;
    int n_samples = 1; 
    int not_transformed_vel = 0; // 0 -> to use B-field aligned transformation, 1 -> vz always aligned with pparall, not physical for neutronics sim.
    //int cartesian = 1; // 1 -> spatial coordinate in cartesian ref frame, good for serpent. 0 -> cylindrical coord, good to visualize results.
    int trials = 0;

    while (n_accepted < n){
        trials++;
        /* ---------------------------
         * 1. Sample position uniformly
         * --------------------------- */
        real r   = rmin + (rmax - rmin) * random_uniform(rdata);
        real phi = phimin + (phimax - phimin) * random_uniform(rdata);
        real z   = zmin + (zmax - zmin) * random_uniform(rdata);
        if (phi < 0.0) phi += 360.0;
        real phirad = phi * CONST_PI / 180.0;

        size_t i0 = 0, i1 = 0, i2 = 0;
        real vol = 1.0;

        if (afsi->type1 == 1)
        {
            i0 = math_bin_index(r, afsi->beam1->axes[0].n,
                                afsi->beam1->axes[0].min,
                                afsi->beam1->axes[0].max);

            i1 = math_bin_index(phirad, afsi->beam1->axes[1].n,
                                afsi->beam1->axes[1].min,
                                afsi->beam1->axes[1].max);

            i2 = math_bin_index(z, afsi->beam1->axes[2].n,
                                afsi->beam1->axes[2].min,
                                afsi->beam1->axes[2].max);

            size_t spatial_index =
                i0 * afsi->volshape[1] * afsi->volshape[2] +
                i1 * afsi->volshape[2] +
                i2;

            vol = afsi->vol[spatial_index];
        }

        /* ---------------------------
         * 3. Magnetic field evaluation
         * --------------------------- */
        real psi, rho[2], B_cyl[3];
        if (B_field_eval_psi(&psi, r, phirad, z, time, &sim->B_data) ||
            B_field_eval_rho(rho, psi, &sim->B_data) || B_field_eval_B(B_cyl, r, phi, z, time, &sim->B_data)) {
            continue;
        }   

        real B_cart[3];
        B_field_cyl_to_cartesian(B_cart, B_cyl, r, phirad);

        /* ---------------------------
         * 4. Sample reactant momenta
         * --------------------------- */
        real density1, density2;    
        real ppara1, pperp1, ppara2, pperp2;

        afsi_sample_reactant_momenta_2d_alt(
            sim, afsi, m1, m2, vol, n_samples,
            i0, i1, i2,
            r, phirad, z, time, rho[0],
            cumdist_all,
            &density1, &ppara1, &pperp1,
            &density2, &ppara2, &pperp2
        );
        if(density1 == 0 || density2 == 0) {
            continue;
        }

        real vcom2;
        real vprod1[3], vprod2[3];
        
        /* ---------------------------
         * 5. Product kinematics
         * --------------------------- */
        afsi_compute_product_velocities_3d(
            0, m1, m2, mprod1, mprod2, Q,
            &ppara1, &pperp1,
            &ppara2, &pperp2,
            &vcom2,
            B_cart,
            not_transformed_vel,
            vprod1, vprod2
        );

        real E = 0.5 * (m1 * m2) / (m1 + m2) * vcom2;

        real sigma = boschhale_sigma(afsi->reaction, E);
        
        real S = density1 * density2 * sqrt(vcom2)*sigma;
        /* ---------------------------
         * 6. Rejection sampling
         * --------------------------- */
        real u  =random_uniform(rdata);

        if (u < S / Smax) {
            real weight = S * vol;
            afsi_store_particle_data(
                n_accepted,
                r, phirad, z, 
                vprod2, 
                mprod2, 
                prod2, 
                position_space_id,
                velocity_space_id, 
                weight
            );
            n_accepted++;
        }
    }
    printf("Trials: %d\n", trials);
    printf("Accepted: %d\n", n_accepted);
    printf("Acceptance ratio: %.3e\n", (double)n_accepted / (double)trials);
    printf("Smax: %.3e\n", Smax);
}

/**
 * @brief Compute velocities of reaction products.
 *
 * @param i marker index on input velocity and output momentum arrays.
 * @param m1 mass of reactant 1 [kg].
 * @param m2 mass of reactant 2 [kg].
 * @param Q energy released in the reaction [eV].
 * @param ppara1 the parallel momentum of react1
 * @param pperp1 the perpendicular momentum of react1
 * @param ppara2 the parallel momentum of react2
 * @param pperp2 the perpendicular momentum of react2
 * @param vcom2 pointer for storing relative velocity of i'th reactants
 * @param B magnetic field in cylindrical coordinates.
 * @param not_transformed_vel flag to indicate if the velocities are not transformed to B-field aligned coordinates.
 * @param vprod1 array where velocity of product 1 is stored.
 * @param vprod2 array where velocity of product 2 is stored.
 */

void afsi_compute_product_velocities_3d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* B, int not_transformed_vel, real* vprod1, real* vprod2) {
    
    const real EPS = 1e-20;
    // --- Checking masses ---
    if (m1 <= 0 || m2 <= 0 || mprod1 <= 0 || mprod2 <= 0) {
        printf("Warning: Non-positive mass encountered"
               " (marker %d, m1=%e, m2=%e, mprod1=%e, mprod2=%e).\n",
               i, m1, m2, mprod1, mprod2);
        vprod1[0]=vprod1[1]=vprod1[2]=0.0;
        vprod2[0]=vprod2[1]=vprod2[2]=0.0;
        *vcom2 = 0.0;
        return;
    }
    
    real rn1 = CONST_2PI * random_uniform(rdata);
    real rn2 = CONST_2PI * random_uniform(rdata);
    real v1x, v1y, v1z, v2x, v2y, v2z;

    if (not_transformed_vel){
        v1x = cos(rn1) * pperp1[i] / m1;
        v1y = sin(rn1) * pperp1[i] / m1;
        v1z = ppara1[i] / m1; 

        v2x = cos(rn2) * pperp2[i] / m2;
        v2y = sin(rn2) * pperp2[i] / m2;
        v2z = ppara2[i] / m2;
    } else {
        real Bmag = sqrt(B[0]*B[0] + B[1]*B[1] + B[2]*B[2]);

        real bx, by, bz;
        if (Bmag < EPS || !isfinite(Bmag)){
            printf("Warning: Bmag near zero or non-finite. (marker %d)\n", i);
            bx = by = 0.0;
            bz = 1.0; // Default to z-direction if B is negligible
        }
        else {
            bx = B[0] / Bmag;
            by = B[1] / Bmag;
            bz = B[2] / Bmag;
        }
        
        // Vector e1 perpendicular to B
        real ex1, ey1, ez1;
        if (fabs(bx) > fabs(bz)){
            ex1 = -by;
            ey1 = bx;
            ez1 = 0.0;
        }
        else {
            //printf("Note: bx > bz. (marker %d)\n", i);
            ex1 = 0.0;
            ey1 = -bz;
            ez1 = by;
        }
        
        real e1mag = sqrt(ex1*ex1 + ey1*ey1 + ez1*ez1);
       if (e1mag < EPS || !isfinite(e1mag)) {
            // fallback orthogonal vector
            printf("Warning: e1mag near zero or non-finite. Using fallback orthogonal vector. (marker %d)\n", i);
            if (fabs(bz) < 0.9) {
                ex1 = -by; ey1 = bx; ez1 = 0.0;
            } else {
                ex1 = 1.0; ey1 = 0.0; ez1 = 0.0;
            }
            e1mag = sqrt(ex1*ex1 + ey1*ey1 + ez1*ez1);
        }
        ex1 /= e1mag;
        ey1 /= e1mag;
        ez1 /= e1mag;

        real ex2 = by * ez1 - bz * ey1;
        real ey2 = bz * ex1 - bx * ez1;
        real ez2 = bx * ey1 - by * ex1;

        real e2mag = sqrt(ex2*ex2 + ey2*ey2 + ez2*ez2);
        if (e2mag < EPS || !isfinite(e2mag)) {
            // fallback perpendicular
            printf("Warning: e2mag near zero or non-finite. Using fallback perpendicular vector. (marker %d)\n", i);
            ex2 = 0.0; ey2 = 1.0; ez2 = 0.0;
            e2mag = 1.0;
        }
        ex2 /= e2mag;
        ey2 /= e2mag;
        ez2 /= e2mag;

        v1x = ppara1[i] / m1 * bx + pperp1[i] / m1 * (cos(rn1) * ex1 + sin(rn1) * ex2);
        v1y = ppara1[i] / m1 * by + pperp1[i] / m1 * (cos(rn1) * ey1 + sin(rn1) * ey2);
        v1z = ppara1[i] / m1 * bz + pperp1[i] / m1 * (cos(rn1) * ez1 + sin(rn1) * ez2);

        v2x = ppara2[i] / m2 * bx + pperp2[i] / m2 * (cos(rn2) * ex1 + sin(rn2) * ex2);
        v2y = ppara2[i] / m2 * by + pperp2[i] / m2 * (cos(rn2) * ey1 + sin(rn2) * ey2);
        v2z = ppara2[i] / m2 * bz + pperp2[i] / m2 * (cos(rn2) * ez1 + sin(rn2) * ez2);
    }

    if (!isfinite(v1x) || !isfinite(v1y) || !isfinite(v1z) ||
        !isfinite(v2x) || !isfinite(v2y) || !isfinite(v2z)) {
        printf("Warning: Non-finite velocity component computed for reactants (marker %d). pppara1=%e, pperp1=%e, ppara2=%e, pperp2=%e\n", i, ppara1[i], pperp1[i], ppara2[i], pperp2[i]);

    }
    real dvx = v1x - v2x;
    real dvy = v1y - v2y;
    real dvz = v1z - v2z;
    *vcom2 = dvx*dvx + dvy*dvy + dvz*dvz;
    if (!isfinite(*vcom2)) {
        printf("DEBUG: vcom2 non-finite (marker %d)\n", i);
    }

    if ((*vcom2) > 1e200) {
        printf("WARNING: vcom2 huge (marker %d): %e\n", i, (*vcom2));
    }
    // reduced mass
    real mu = (m1 * m2) / (m1 + m2);

    real ekin = Q + 0.5 * mu * (*vcom2);

    if (!isfinite(ekin) || ekin < 0.0) {
        printf("Warning: Non-physical kinetic energy computed in CM frame (marker %d, ekin=%e). Setting ekin to zero.\n", i, ekin);
        ekin = 0.0;
    }

    // Velocity of the system's center of mass
    real v_cm[3];
    v_cm[0]  = ( m1 * v1x + m2 * v2x ) / ( m1 + m2 );
    v_cm[1]  = ( m1 * v1y + m2 * v2y ) / ( m1 + m2 );
    v_cm[2]  = ( m1 * v1z + m2 * v2z ) / ( m1 + m2 );

    // Speed and velocity of product 2 in CM frame
    rn1 = random_uniform(&rdata);
    rn2 = random_uniform(&rdata);
    real phi = CONST_2PI * rn1;

    real arg = 2.0*(rn2 - 0.5);
    arg = fmax(-1.0, fmin(1.0, arg));
    real theta = acos(arg);
    real vnorm = sqrt( 2.0 * ekin / ( mprod2 * ( 1.0 + mprod2 / mprod1 ) ) );

    real v2_cm[3];
    v2_cm[0] = vnorm * sin(theta) * cos(phi);
    v2_cm[1] = vnorm * sin(theta) * sin(phi);
    v2_cm[2] = vnorm * cos(theta);

    // Products' velocities in lab frame
    vprod1[0] = -(mprod2/mprod1) * v2_cm[0] + v_cm[0];
    vprod1[1] = -(mprod2/mprod1) * v2_cm[1] + v_cm[1];
    vprod1[2] = -(mprod2/mprod1) * v2_cm[2] + v_cm[2];
    vprod2[0] = v2_cm[0] + v_cm[0];
    vprod2[1] = v2_cm[1] + v_cm[1];
    vprod2[2] = v2_cm[2] + v_cm[2];
    if (!isfinite(vprod2[0]) || !isfinite(vprod2[1]) || !isfinite(vprod2[2])){
        printf("Warning: Non-finite velocity computed for product 2 (marker %d).\n", i);
    }
}

/**
 * @brief Store sampled particle data in the prod2 array.
 *
 * @param i index of the particle being stored.
 * @param rvec array with radial coordinates [m].
 * @param phivec array with azimuthal coordinates [rad].
 * @param zvec array with vertical coordinates [m].
 * @param vprod2 array with velocity of product 2 [m/s].
 * @param mprod2 mass of product 2 [kg].
 * @param prod2 array where the sampled particle data is stored.
 */

void afsi_store_particle_data(int i, real r, real phi, real z, real* vprod2, real mprod2, real* prod2, int position_space_id, int velocity_space_id, real weight){
    /* --------------------------
     * Velocity
     * -------------------------- */
    real vx = vprod2[0];
    real vy = vprod2[1];
    real vz = vprod2[2];
    real v2 = vx*vx + vy*vy + vz*vz;
    real v  = sqrt(v2);
    if (v <= 0.0) v = 1e-30;

    // Compute the kinetic energy of prod2
    // real energy_prod2 = physlib_Ekin_gamma(mprod2, physlib_gamma_vnorm(vprod2_magnitude));
    real energy_prod2 = 0.5 * mprod2 * v2 *6241506479963.2;// Convert to MeV

    // Direction vector
    real ux = vx/v;
    real uy = vy/v;
    real uz = vz/v;
    
    /* --------------------------
     * Spherical velocity angles
     * -------------------------- */
    real theta = acos(uz);                 // polar angle
    real phi_v = atan2(vy, vx);            // azimuthal angle
    
    /* ==========================================================
     * POSITION SPACE
     * ========================================================== */
    if (position_space_id == 0){
        /* cylindrical */
        r = r * 100; // convert to cm
        z = z * 100; // convert to cm
        phi = phi * 180.0 / CONST_PI; // convert to degrees
        prod2[(i) * 9 + 0] = r; 
        prod2[(i) * 9 + 1] = phi * 180.0 / CONST_PI; 
        prod2[(i) * 9 + 2] = z;
    }
    else if (position_space_id == 1){
        /* cartesian */
        real x = r * cos(phi) * 100; // convert to cm
        real y = r * sin(phi) * 100; // convert to cm
        z = z * 100; // convert to cm
        prod2[(i) * 9 + 0] = x;
        prod2[(i) * 9 + 1] = y;
        prod2[(i) * 9 + 2] = z;
    }
    else {
        print_err("Invalid position space ID.\n");
        abort();
    }
    
    /* ==========================================================
     * VELOCITY SPACE
     * ========================================================== */
    if (velocity_space_id == 0) {
        prod2[(i) * 9 + 3] = ux; 
        prod2[(i) * 9 + 4] = uy;
        prod2[(i) * 9 + 5] = uz;
    }
    else if (velocity_space_id == 1) {
        /* momentum */
        prod2[(i) * 9 + 3] = mprod2 * vx; 
        prod2[(i) * 9 + 4] = mprod2 * vy;
        prod2[(i) * 9 + 5] = mprod2 * vz;

    }
    else if (velocity_space_id == 2)
    {
        /* energy + angles */
        prod2[i*9 + 3] = energy_prod2;
        prod2[i*9 + 4] = theta;
        prod2[i*9 + 5] = phi_v;
    }
    else {
        print_err("Invalid velocity space ID.\n");
        abort();
    }
    /* ==========================================================
     * COMMON FIELDS
     * ========================================================== */
    prod2[i*9 + 6] = energy_prod2;
    prod2[i*9 + 7] = weight;
    prod2[i*9 + 8] = 0.0;   /* time */
}


/**
 * @brief Compute momenta of reaction products.
 *
 * @param i marker index on input velocity and output momentum arrays.
 * @param m1 mass of reactant 1 [kg].
 * @param m2 mass of reactant 2 [kg].
 * @param mprod1 mass of product 1 [kg].
 * @param mprod2 mass of product 2 [kg].
 * @param Q energy released in the reaction [eV].
 * @param ppara1 the parallel momentum of react1
 * @param pperp1 the perpendicular momentum of react1
 * @param ppara2 the parallel momentum of react2
 * @param pperp2 the perpendicular momentum of react2
 * @param vcom2 pointer for storing relative velocity of i'th reactants
 * @param pparaprod1 array where parallel momentum of product 1 is stored.
 * @param pperpprod1 array where perpendicular momentum of product 1 is stored.
 * @param pparaprod2 array where parallel momentum of product 2 is stored.
 * @param pperpprod2 array where perpendicular momentum of product 2 is stored.
 */
void afsi_compute_product_momenta_2d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q, int prodmomspace,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* prod1_p1, real* prod1_p2, real* prod2_p1, real* prod2_p2) {

    real rn1 = CONST_2PI * random_uniform(&rdata);
    real rn2 = CONST_2PI * random_uniform(&rdata);

    real v1x = cos(rn1) * pperp1[i] / m1;
    real v1y = sin(rn1) * pperp1[i] / m1;
    real v1z = ppara1[i] / m1;

    real v2x = cos(rn2) * pperp2[i] / m2;
    real v2y = sin(rn2) * pperp2[i] / m2;
    real v2z = ppara2[i] / m2;

    *vcom2 =   (v1x - v2x) * (v1x - v2x)
             + (v1y - v2y) * (v1y - v2y)
             + (v1z - v2z) * (v1z - v2z);

    // Velocity of the system's center of mass
    real v_cm[3];
    v_cm[0]  = ( m1 * v1x + m2 * v2x ) / ( m1 + m2 );
    v_cm[1]  = ( m1 * v1y + m2 * v2y ) / ( m1 + m2 );
    v_cm[2]  = ( m1 * v1z + m2 * v2z ) / ( m1 + m2 );

    // Total kinetic energy after the reaction in CM frame
    real ekin = Q
        + 0.5 * m1 * (   (v1x - v_cm[0])*(v1x - v_cm[0])
                       + (v1y - v_cm[1])*(v1y - v_cm[1])
                       + (v1z - v_cm[2])*(v1z - v_cm[2]) )
        + 0.5 * m2 * (   (v2x - v_cm[0])*(v2x - v_cm[0])
                       + (v2y - v_cm[1])*(v2y - v_cm[1])
                       + (v2z - v_cm[2])*(v2z - v_cm[2]) );

    // Speed and velocity of product 2 in CM frame
    rn1 = random_uniform(&rdata);
    rn2 = random_uniform(&rdata);
    real phi   = CONST_2PI * rn1;
    real theta = acos( 2 * ( rn2 - 0.5 ) );
    real vnorm = sqrt( 2.0 * ekin / ( mprod2 * ( 1.0 + mprod2 / mprod1 ) ) );

    real v2_cm[3];
    v2_cm[0] = vnorm * sin(theta) * cos(phi);
    v2_cm[1] = vnorm * sin(theta) * sin(phi);
    v2_cm[2] = vnorm * cos(theta);

    // Products' velocities in lab frame
    real vprod1[3], vprod2[3];
    vprod1[0] = -(mprod2/mprod1) * v2_cm[0] + v_cm[0];
    vprod1[1] = -(mprod2/mprod1) * v2_cm[1] + v_cm[1];
    vprod1[2] = -(mprod2/mprod1) * v2_cm[2] + v_cm[2];
    vprod2[0] = v2_cm[0] + v_cm[0];
    vprod2[1] = v2_cm[1] + v_cm[1];
    vprod2[2] = v2_cm[2] + v_cm[2];

    if(prodmomspace == PPARPPERP) {
        prod1_p1[i] = vprod1[2] * mprod1;
        prod1_p2[i] = sqrt(vprod1[0]*vprod1[0] + vprod1[1]*vprod1[1]) * mprod1;
        prod2_p1[i] = vprod2[2] * mprod2;
        prod2_p2[i] = sqrt(vprod2[0]*vprod2[0] + vprod2[1]*vprod2[1]) * mprod2;
    }
    else {
        real vnorm1 = math_norm(vprod1);
        prod1_p2[i] = vprod1[2] / vnorm1;
        prod1_p1[i] = physlib_Ekin_gamma(mprod1, physlib_gamma_vnorm(vnorm1));

        real vnorm2 = math_norm(vprod2);
        prod2_p2[i] = vprod2[2] / vnorm2;
        // Using classical kinetic energy
        prod2_p2[i] = 0.5*mprod2*vnorm2*vnorm2;
        //prod2_p1[i] = physlib_Ekin_gamma(mprod2, physlib_gamma_vnorm(vnorm2));
    }
}

/**
 * @brief Sample momenta from reactant distributions.
 *
 * @param react1 reactant 1 distribution data.
 * @param react2 reactant 2 distribution data.
 * @param m1 mass of reactant 1 [kg].
 * @param m2 mass of reactant 2 [kg].
 * @param n number of samples.
 * @param iR R index of the distribution cell where sampling is done.
 * @param iphi phi index of the distribution cell where sampling is done.
 * @param iz z index of the distribution cell where sampling is done.
 * @param ppara1 array where the parallel momentum of react1 will be stored
 * @param pperp1 array where the perpendicular momentum of react1 will be stored
 * @param ppara2 array where the parallel momentum of react2 will be stored
 * @param pperp2 array where the perpendicular momentum of react2 will be stored
 */
void afsi_sample_reactant_momenta_2d(
    sim_data* sim, afsi_data* afsi, real mass1, real mass2, real vol,
    int nsample, size_t i0, size_t i1, size_t i2,
    real r, real phi, real z, real time, real rho,
    real* density1, real* ppara1, real* pperp1,
    real* density2, real* ppara2, real* pperp2) {
    if(afsi->type1 == 1) {
        afsi_sample_beam_2d(afsi->beam1, mass1, vol, nsample, i0, i1, i2,
                            density1, ppara1, pperp1);
    }
    else if(afsi->type1 == 2) {
        afsi_sample_thermal_2d(sim, afsi->thermal1, mass1, nsample, r, phi, z,
                               time, rho, density1, ppara1, pperp1);
    }
    if(afsi->type2 == 1) {
        afsi_sample_beam_2d(afsi->beam2, mass2, vol, nsample, i0, i1, i2,
                            density2, ppara2, pperp2);
    }
    else if(afsi->type2 == 2) {
        afsi_sample_thermal_2d(sim, afsi->thermal2, mass2, nsample, r, phi, z,
                               time, rho, density2, ppara2, pperp2);
    }
}

/**
 * @brief Sample ppara and pperp from a 5D distribution.
 *
 * @param dist pointer to the distribution data.
 * @param nsample number of values to be sampled.
 * @param spatial_index spatial index where sampling is done.
 * @param ppara pointer to array where sampled parallel momenta are stored.
 * @param pperp pointer to array where sampled perpedicular momenta are stored.
 */
void afsi_sample_beam_2d(histogram* hist, real mass, real vol, int nsample,
                         size_t i0, size_t i1, size_t i2,
                         real* density, real* ppara, real* pperp) {
    int mom_space;
    size_t p1coord, p2coord;
    if(hist->axes[5].n) {
        p1coord = 5;
        p2coord = 6;
        mom_space = PPARPPERP;
    }
    else if(hist->axes[10].n) {
        p1coord = 10;
        p2coord = 11;
        mom_space = EKINXI;
    }
    else {
        return;
    }

    real* cumdist = (real*) malloc(
        hist->axes[p1coord].n*hist->axes[p2coord].n*sizeof(real));

    *density = 0.0;
    for(size_t ip1 = 0; ip1 < hist->axes[p1coord].n; ip1++) {
        for(size_t ip2 = 0; ip2 < hist->axes[p2coord].n; ip2++) {
            size_t index = i0*hist->strides[0]
                         + i1*hist->strides[1]
                         + i2*hist->strides[2]
                         + ip1*hist->strides[p1coord]
                         + ip2*hist->strides[p2coord];
            size_t j = ip1*hist->axes[p2coord].n+ip2;
            if(j == 0) {
                cumdist[0] = hist->bins[index];
            } else {
                cumdist[j] = cumdist[j - 1] + hist->bins[index];
            }
            *density += hist->bins[index] / vol;
        }
    }
    if(*density == 0) {
        free(cumdist);
        return;
    }

    for(size_t i = 0; i < nsample; i++) {
        real r = random_uniform(&rdata);
        r *= cumdist[hist->axes[p1coord].n*hist->axes[p2coord].n-1];
        for(size_t j = 0; j < hist->axes[p1coord].n*hist->axes[p2coord].n; j++) {
            if(cumdist[j] > r) {
                if(mom_space == PPARPPERP) {
                    ppara[i] = hist->axes[5].min + (j / hist->axes[6].n + 0.5)
                        * (hist->axes[5].max - hist->axes[5].min) / hist->axes[5].n;
                    pperp[i] = hist->axes[6].min + (j % hist->axes[6].n + 0.5)
                        * (hist->axes[6].max - hist->axes[6].min) / hist->axes[6].n;
                }
                else {
                    real ekin = hist->axes[10].min + (j / hist->axes[11].n + 0.5)
                        * (hist->axes[10].max - hist->axes[10].min) / hist->axes[10].n;
                    real pitch = hist->axes[11].min + (j % hist->axes[11].n + 0.5)
                        * (hist->axes[11].max - hist->axes[11].min) / hist->axes[11].n;
                    real gamma = physlib_gamma_Ekin(mass, ekin);
                    real pnorm = sqrt(gamma * gamma - 1.0) * mass * CONST_C;
                    ppara[i] = pitch * pnorm;
                    pperp[i] = sqrt( 1.0 - pitch*pitch ) * pnorm;
                }
                break;
            }
        }
    }
    free(cumdist);
}

/**
 * @brief Sample ppara and pperp from a thermal (Maxwellian) population.
 *
 * Sampling from Maxwellian distribution is done using Eq. (7) in
 * "Efficient Algorithm for Generating Maxwell Random Variables", N. Mohamed,
 * DOI 10.1007/s10955-011-0364-y
 *
 * @param data pointer to the thermal data.
 * @param mass mass of the particle species.
 * @param nsample number of values to be sampled.
 * @param spatial_index spatial index where sampling is done.
 * @param ppara pointer to array where sampled parallel momenta are stored.
 * @param pperp pointer to array where sampled perpedicular momenta are stored.
 */
void afsi_sample_thermal_2d(sim_data* sim, int ispecies, real mass, int nsample,
                            real r, real phi, real z, real time, real rho,
                            real* density, real* ppara, real* pperp) {
    real ni, ti;
    if( plasma_eval_dens(&ni, rho, r, phi, z, time, ispecies,
        &sim->plasma_data) ||
        plasma_eval_temp(&ti, rho, r, phi, z, time, ispecies,
        &sim->plasma_data) ) {
        *density = 0.0;
        return;
    }
    *density = ni;
    for(int i = 0; i < nsample; i++) {
        real r1, r2, r3, r4, E;

        r1 = random_uniform(&rdata);
        r2 = random_uniform(&rdata);
        r3 = cos( 0.5 * random_uniform(&rdata) * CONST_PI );
        E  = -ti * ( log(r1) + log(r2) * r3 * r3 );

        r4 = 1.0 - 2 * random_uniform(&rdata);
        pperp[i] = sqrt( ( 1 - r4*r4 ) * 2 * E * mass);
        ppara[i] = r4 * sqrt(2 * E * mass);
    }
}

/**
 * @brief Sample momenta from reactant distributions. Alternative version for rejection sampling
 *
 * @param react1 reactant 1 distribution data.
 * @param react2 reactant 2 distribution data.
 * @param m1 mass of reactant 1 [kg].
 * @param m2 mass of reactant 2 [kg].
 * @param n number of samples.
 * @param iR R index of the distribution cell where sampling is done.
 * @param iphi phi index of the distribution cell where sampling is done.
 * @param iz z index of the distribution cell where sampling is done.
 * @param r radial coordinate [m].
 * @param phi azimuthal coordinate [rad].
 * @param z vertical coordinate [m].
 * @param time time [s].
 * @param rho plasma density [m^-3].
 * @param cumdist_all array with cumulative distribution function values for all cells.
 * 
 * @param density1 array where the density of react1 will be stored
 * @param density2 array where the density of react2 will be stored
 * @param ppara1 array where the parallel momentum of react1 will be stored
 * @param pperp1 array where the perpendicular momentum of react1 will be stored
 * @param ppara2 array where the parallel momentum of react2 will be stored
 * @param pperp2 array where the perpendicular momentum of react2 will be stored
 */

void afsi_sample_reactant_momenta_2d_alt(
    sim_data* sim, afsi_data* afsi, real mass1, real mass2, real vol,
    int nsample, size_t i0, size_t i1, size_t i2,
    real r, real phi, real z, real time, real rho, real* cumdist_all,
    real* density1, real* ppara1, real* pperp1,
    real* density2, real* ppara2, real* pperp2) {
    if(afsi->type1 == 1) {
        size_t offset = (i0*afsi->volshape[1]*afsi->volshape[2]+ i1*afsi->volshape[2] + i2) * (afsi->beam1->axes[5].n*afsi->beam1->axes[6].n);
        real* cumdist = &cumdist_all[offset];
        afsi_sample_beam_2d_alt(afsi->beam1, mass1, vol, nsample, i0, i1, i2, cumdist,
                            density1, ppara1, pperp1);
    }
    else if(afsi->type1 == 2) {
        afsi_sample_thermal_2d(sim, afsi->thermal1, mass1, nsample, r, phi, z,
                               time, rho, density1, ppara1, pperp1);
    }
    if(afsi->type2 == 1) {
        size_t offset = (i0*afsi->volshape[1]*afsi->volshape[2]+ i1*afsi->volshape[2] + i2) * (afsi->beam1->axes[5].n*afsi->beam1->axes[6].n);
        real* cumdist = &cumdist_all[offset];
        afsi_sample_beam_2d_alt(afsi->beam2, mass2, vol, nsample, i0, i1, i2, cumdist,
                            density2, ppara2, pperp2);
    }
    else if(afsi->type2 == 2) {
        afsi_sample_thermal_2d(sim, afsi->thermal2, mass2, nsample, r, phi, z,
                               time, rho, density2, ppara2, pperp2);
    }
}

/**
 * @brief Sample ppara and pperp from a 5D distribution. Alternative version for rejection sampling.
 *
 * @param dist pointer to the distribution data.
 * @param nsample number of values to be sampled.
 * @param spatial_index spatial index where sampling is done.
 * @param cumdist pointer to the starting point of cumulative distribution function values for the current spatial index.
 * @param ppara pointer to array where sampled parallel momenta are stored.
 * @param pperp pointer to array where sampled perpedicular momenta are stored.
 */

void afsi_sample_beam_2d_alt(histogram* hist, real mass, real vol, int nsample,
                         size_t i0, size_t i1, size_t i2, real* cumdist,
                         real* density, real* ppara, real* pperp) {
    int mom_space;
    size_t p1coord, p2coord;
    if(hist->axes[5].n) {
        p1coord = 5;
        p2coord = 6;
        mom_space = PPARPPERP;
    }
    else if(hist->axes[10].n) {
        p1coord = 10;
        p2coord = 11;
        mom_space = EKINXI;
    }
    else {
        return;
    }

    

    *density = 0.0;
    for(size_t ip1 = 0; ip1 < hist->axes[p1coord].n; ip1++) {
        for(size_t ip2 = 0; ip2 < hist->axes[p2coord].n; ip2++) {
            size_t index = i0*hist->strides[0]
                         + i1*hist->strides[1]
                         + i2*hist->strides[2]
                         + ip1*hist->strides[p1coord]
                         + ip2*hist->strides[p2coord];
            *density += hist->bins[index] / vol;
        }
    }
    if(*density == 0) {
        return;
    }

    for(size_t i = 0; i < nsample; i++) {
        real r = random_uniform(rdata);
        r *= cumdist[hist->axes[p1coord].n*hist->axes[p2coord].n-1];
        for(size_t j = 0; j < hist->axes[p1coord].n*hist->axes[p2coord].n; j++) {
            if(cumdist[j] > r) {
                if(mom_space == PPARPPERP) {
                    ppara[i] = hist->axes[5].min + (j / hist->axes[6].n + 0.5)  
                        * (hist->axes[5].max - hist->axes[5].min) / hist->axes[5].n;
                    pperp[i] = hist->axes[6].min + (j % hist->axes[6].n + 0.5)
                        * (hist->axes[6].max - hist->axes[6].min) / hist->axes[6].n;
                }
                else {
                    real ekin = hist->axes[10].min + (j / hist->axes[10].n + 0.5)
                        * (hist->axes[10].max - hist->axes[10].min) / hist->axes[10].n;
                    real pitch = hist->axes[11].min + (j / hist->axes[11].n + 0.5) 
                        * (hist->axes[11].max - hist->axes[11].min) / hist->axes[11].n;
                    real gamma = physlib_gamma_Ekin(mass, ekin);
                    real pnorm = sqrt(gamma * gamma - 1.0) * mass * CONST_C;
                    ppara[i] = pitch * pnorm;
                    pperp[i] = sqrt( 1.0 - pitch*pitch ) * pnorm;
                }
                break;
            }
        }
    }
}

/**
 * @brief Convert magnetic field from cylindrical to Cartesian coordinates.
 *
 * @param B_cart pointer to the output Cartesian magnetic field vector.
 * @param B_cyl pointer to the input cylindrical magnetic field vector.
 * @param r radial coordinate in cylindrical coordinates [m].
 * @param phi azimuthal coordinate in cylindrical coordinates [rad].
 */
void B_field_cyl_to_cartesian(real* B_cart, real* B_cyl, real r, real phi){
    B_cart[0] = B_cyl[0] * cos(phi) - B_cyl[1] * sin(phi);
    B_cart[1] = B_cyl[0] * sin(phi) + B_cyl[1] * cos(phi);
    B_cart[2] = B_cyl[2];
}