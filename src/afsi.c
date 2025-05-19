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

void afsi_compute_product_velocities_3d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* vprod1, real* vprod2);
void afsi_compute_product_momenta_2d(
    int i, real mprod1, real mprod2, int prodmomspace,
    real* vprod1, real* vprod2, real* prod1_p1, real* prod1_p2,
    real* prod2_p1, real* prod2_p2);
void afsi_store_particle_data(int i, real* rvec, real* phivec,
     real* zvec, real* vprod2, real mprod2, real* prod2);
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
                    real vprod1[3];
                    real vprod2[3];

                    afsi_compute_product_velocities_3d(
                        i, m1, m2, mprod1, mprod2, Q,
                        ppara1, pperp1, ppara2, pperp2, &vcom2,
                        vprod1, vprod2);
                    
                    afsi_compute_product_momenta_2d(i, mprod1, mprod2, prod_mom_space,
                        vprod1, vprod2,  prod1_p1, prod1_p2, prod2_p1, prod2_p2);
                    
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

void afsi_run_new(sim_data* sim, afsi_data* afsi, int n,
                real* rvec, real* phivec, real* zvec, real* prod2, 
                int i0, int i1, int i2) {

    random_init(&rdata, time((NULL)));
    simulate_init(sim);

    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;
    boschhale_reaction(
    afsi->reaction, &m1, &q1, &m2, &q2,
    &mprod1, &qprod1, &mprod2, &qprod2, &Q);
    
    real time = 0.0;
    real* ppara1 = (real*) malloc(n*sizeof(real));
    real* pperp1 = (real*) malloc(n*sizeof(real));
    real* ppara2 = (real*) malloc(n*sizeof(real));
    real* pperp2 = (real*) malloc(n*sizeof(real));

        
    real r = afsi->r[0];
    real z = afsi->z[0];
    real phi = afsi->phi[0];
    real vol = afsi->vol[0];
    
    real psi, rho[2];
    if(B_field_eval_psi(&psi, r, phi, z, time, &sim->B_data) ||
        B_field_eval_rho(rho, psi, &sim->B_data) ) {
        return;
    }

    real density1, density2;
    afsi_sample_reactant_momenta_2d(
        sim, afsi, m1, m2, vol, n, i0, i1, i2,
        r, phi, z, time, rho[0],
        &density1, ppara1, pperp1, &density2, ppara2, pperp2);
    if(density1 == 0 || density2 == 0) {
        return;
    }
    for(size_t i = 0; i < n; i++) {
        real vcom2;
        real vprod1[3];
        real vprod2[3];
        
        afsi_compute_product_velocities_3d(
            i, m1, m2, mprod1, mprod2, Q,
            ppara1, pperp1, ppara2, pperp2, &vcom2,
            vprod1, vprod2);
        
        afsi_store_particle_data(i, rvec, phivec, zvec, vprod2, mprod2, prod2);
    }
    free(ppara1);
    free(ppara2);
    free(pperp1);
    free(pperp2);
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

void afsi_compute_product_velocities_3d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* vprod1, real* vprod2) {
    
    real rn1 = CONST_2PI * random_uniform(rdata);
    real rn2 = CONST_2PI * random_uniform(rdata);

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
    vprod1[0] = -(mprod2/mprod1) * v2_cm[0] + v_cm[0];
    vprod1[1] = -(mprod2/mprod1) * v2_cm[1] + v_cm[1];
    vprod1[2] = -(mprod2/mprod1) * v2_cm[2] + v_cm[2];
    vprod2[0] = v2_cm[0] + v_cm[0];
    vprod2[1] = v2_cm[1] + v_cm[1];
    vprod2[2] = v2_cm[2] + v_cm[2];
}

void afsi_store_particle_data(int i, real* rvec, real* phivec, real* zvec, real* vprod2, real mprod2, real* prod2){

    // Compute the magnitude of the velocity of prod2
    real vprod2_magnitude = sqrt(vprod2[0] * vprod2[0] +
    vprod2[1] * vprod2[1] +
    vprod2[2] * vprod2[2]);

    // Compute the kinetic energy of prod2
    real energy_prod2 = 0.5 * mprod2 * vprod2_magnitude * vprod2_magnitude;

    // Normalize vprod2 to get the direction vector
    real u = vprod2[0] / vprod2_magnitude;
    real v = vprod2[1] / vprod2_magnitude;
    real w = vprod2[2] / vprod2_magnitude;

    // Sample the coordinates in the bin
    real r = rvec[0] + (rvec[1] - rvec[0]) * random_uniform(rdata);
    real phi = phivec[0] + (phivec[1] - phivec[0]) * random_uniform(rdata);
    real z = zvec[0] + (zvec[1] - zvec[0]) * random_uniform(rdata);

    // Store the sampled coordinates, directions and energy in the prod2 array
    prod2[i * 7 + 0] = r; 
    prod2[i * 7 + 1] = phi; 
    prod2[i * 7 + 2] = z; 
    prod2[i * 7 + 3] = u; 
    prod2[i * 7 + 4] = v;
    prod2[i * 7 + 5] = w;
    prod2[i * 7 + 6] = energy_prod2*6241506479963.2; // Convert to MeV
}

void afsi_compute_product_momenta_2d(
    int i, real mprod1, real mprod2, int prodmomspace,
    real* vprod1, real* vprod2, real* prod1_p1, real* prod1_p2,
    real* prod2_p1, real* prod2_p2) {

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
        prod2_p1[i] = physlib_Ekin_gamma(mprod2, physlib_gamma_vnorm(vnorm2));
    }
}

/**
 * @brief Sample velocities from reactant distributions.
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
            if(ip1 == 0 && ip2 == 0) {
                cumdist[0] = hist->bins[index];
            } else {
                cumdist[ip1*hist->axes[p2coord].n+ip2] =
                      cumdist[ip1*hist->axes[p2coord].n+ip2-1]
                    + hist->bins[index];
            }
            *density += hist->bins[index] / vol;
        }
    }
    if(*density == 0) {
        return;
    }

    for(size_t i = 0; i < nsample; i++) {
        real r = random_uniform(&rdata);
        r *= cumdist[hist->axes[p1coord].n*hist->axes[p2coord].n-1];
        for(size_t j = 0; j < hist->axes[p1coord].n*hist->axes[p2coord].n; j++) {
            if(cumdist[j] > r) {
                if(mom_space == PPARPPERP) {
                    ppara[i] = hist->axes[5].min + (j / hist->axes[5].n + 0.5)
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
