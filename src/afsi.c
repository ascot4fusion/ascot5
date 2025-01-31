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

void afsi_sample_reactant_momenta_2d(
    afsi_data* react, real mass, int nsample, size_t spatial_index,
    real* ppara, real* pperp);
void afsi_compute_product_momenta_2d(
    int i, real m1, real m2, real mprod1, real mprod2, real Q, int prodmomspace,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* prod1_p1, real* prod1_p2, real* prod2_p1, real* prod2_p2);
void afsi_sample_beam_2d(histogram* hist, int nsample, size_t spatial_index,
                         real mass, real* p1, real* p2);
void afsi_sample_thermal_2d(afsi_thermal_data* data, real mass, int nsample,
                            size_t spatial_index, real* p1, real* p2);
real afsi_get_density(afsi_data* dist, size_t spatial_index, real r);
real afsi_get_volume(afsi_data* dist, real vol);

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
void afsi_run(sim_data* sim, Reaction reaction, int n,
              afsi_data* react1, afsi_data* react2, real mult,
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
        reaction, &m1, &q1, &m2, &q2, &mprod1, &qprod1, &mprod2, &qprod2, &Q);

    int n_r=0, n_phi=0, n_z=0;
    if(react1->type == 1) {
        n_r   = react1->beam->axes[0].n;
        n_phi = react1->beam->axes[1].n;
        n_z   = react1->beam->axes[2].n;
    }
    else if(react1->type == 2) {
        n_r   = react1->thermal->n_r;
        n_phi = react1->thermal->n_phi;
        n_z   = react1->thermal->n_z;
    }

    int prod_mom_space;
    int p1coord, p2coord;
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

    #pragma omp parallel for
    for(size_t ir = 0; ir < n_r; ir++) {
        real* ppara1 = (real*) malloc(n*sizeof(real));
        real* pperp1 = (real*) malloc(n*sizeof(real));
        real* ppara2 = (real*) malloc(n*sizeof(real));
        real* pperp2 = (real*) malloc(n*sizeof(real));
        real* prod1_p1 = (real*) malloc(n*sizeof(real));
        real* prod1_p2 = (real*) malloc(n*sizeof(real));
        real* prod2_p1 = (real*) malloc(n*sizeof(real));
        real* prod2_p2 = (real*) malloc(n*sizeof(real));

        real dr = (prod2->axes[0].max-prod2->axes[0].min) / prod2->axes[0].n;
        real r =  prod2->axes[0].min + 0.5 * dr + ir * dr;
        real vol = afsi_get_volume(react1, r);
        for(size_t iphi = 0; iphi < n_phi; iphi++) {
            for(size_t iz = 0; iz < n_z; iz++) {
                // size_t spatial_index = ir*prod1->strides[0]
                //              + iphi*prod1->strides[1]
                //              + iz*prod1->strides[2];
                size_t spatial_index = ir*n_phi*n_z
                             + iphi*n_z
                             + iz;
                real density1 = afsi_get_density(react1, spatial_index, vol);
                real density2 = afsi_get_density(react2, spatial_index, vol);
                if(density1 > 0 && density2 > 0) {
                    afsi_sample_reactant_momenta_2d(
                        react1, m1, n, spatial_index, ppara1, pperp1);
                    afsi_sample_reactant_momenta_2d(
                        react2, m2, n, spatial_index, ppara2, pperp2);
                    for(size_t i = 0; i < n; i++) {
                        real vcom2;
                        afsi_compute_product_momenta_2d(
                            i, m1, m2, mprod1, mprod2, Q, prod_mom_space,
                            ppara1, pperp1, ppara2, pperp2, &vcom2,
                            prod1_p1, prod1_p2, prod2_p1, prod2_p2);
                        real E = 0.5 * ( m1 * m2 ) / ( m1 + m2 ) * vcom2;

                        real weight = density1 * density2 * sqrt(vcom2)
                            * boschhale_sigma(reaction, E)/n*vol;

                        size_t ip1 = math_bin_index(
                            prod1_p1[i], prod1->axes[p1coord].n,
                            prod1->axes[p1coord].min, prod1->axes[p1coord].max);
                        size_t ip2 = math_bin_index(
                            prod1_p2[i], prod1->axes[p2coord].n,
                            prod1->axes[p2coord].min, prod1->axes[p2coord].max);
                        if( 0 <= ip1 && ip1 < prod1->axes[p1coord].n &&
                            0 <= ip2 && ip2 < prod1->axes[p2coord].n) {
                            size_t index = ir*prod1->strides[0]
                                         + iphi*prod1->strides[1]
                                         + iz*prod1->strides[2]
                                         + ip1*prod1->strides[p1coord]
                                         + ip2*prod1->strides[p2coord];
                            GPU_ATOMIC
                            prod1->bins[index] += weight * mult;
                        }

                        ip1 = math_bin_index(
                            prod2_p1[i], prod2->axes[p1coord].n,
                            prod2->axes[p1coord].min, prod2->axes[p1coord].max);
                        ip2 = math_bin_index(
                            prod2_p2[i], prod2->axes[p2coord].n,
                            prod2->axes[p2coord].min, prod2->axes[p2coord].max);
                        if( 0 <= ip1 && ip1 < prod2->axes[p1coord].n &&
                            0 <= ip2 && ip2 < prod2->axes[p2coord].n) {
                            size_t index = ir*prod2->strides[0]
                                         + iphi*prod2->strides[1]
                                         + iz*prod2->strides[2]
                                         + ip1*prod2->strides[p1coord]
                                         + ip2*prod2->strides[p2coord];
                            GPU_ATOMIC
                            prod2->bins[index] += weight * mult;
                        }
                    }
                }
                else {
                    /* Do nothing */
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
    rn1 = random_uniform(rdata);
    rn2 = random_uniform(rdata);
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
    afsi_data* react, real mass, int nsample, size_t spatial_index,
    real* ppara, real* pperp) {
    if(react->type == 1) {
        afsi_sample_beam_2d(react->beam, nsample, spatial_index, mass,
                            ppara, pperp);
    }
    else if(react->type == 2) {
        afsi_sample_thermal_2d(
            react->thermal, mass, nsample, spatial_index, ppara, pperp);
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
void afsi_sample_beam_2d(histogram* hist, int nsample, size_t spatial_index,
                         real mass, real* ppara, real* pperp) {
    int mom_space;
    int p1coord, p2coord;
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
        hist->axes[5].n*hist->axes[6].n*sizeof(real));

    for(size_t ippara = 0; ippara < hist->axes[p1coord].n; ippara++) {
        for(size_t ipperp = 0; ipperp < hist->axes[p2coord].n; ipperp++) {
            size_t index = spatial_index
                         + ippara*hist->strides[p1coord]
                         + ipperp*hist->strides[p2coord];
            if(ippara == 0 && ipperp == 0) {
                cumdist[0] = hist->bins[index];
            } else {
                cumdist[ippara*hist->axes[p2coord].n+ipperp] =
                      cumdist[ippara*hist->axes[p2coord].n+ipperp-1]
                    + hist->bins[index];
            }
        }
    }
    for(size_t ippara = 0; ippara < hist->axes[p1coord].n; ippara++) {
        for(size_t ipperp = 0; ipperp < hist->axes[p2coord].n; ipperp++) {
            cumdist[ippara*hist->axes[p2coord].n+ipperp] /=
                cumdist[hist->axes[p1coord].n*hist->axes[p2coord].n-1];
        }
    }

    for(size_t i = 0; i < nsample; i++) {
        real r = random_uniform(rdata);
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
void afsi_sample_thermal_2d(afsi_thermal_data* data, real mass, int nsample,
                            size_t spatial_index, real* ppara, real* pperp) {
    real temp = data->temperature[spatial_index];
    for(int i = 0; i < nsample; i++) {
        real r1, r2, r3, r4, E;

        r1 = random_uniform(rdata);
        r2 = random_uniform(rdata);
        r3 = cos( 0.5 * random_uniform(rdata) * CONST_PI );
        E  = -temp * ( log(r1) + log(r2) * r3 * r3 );

        r4 = 1.0 - 2 * random_uniform(rdata);
        pperp[i] = sqrt( ( 1 - r4*r4 ) * 2 * E * mass);
        ppara[i] = r4 * sqrt(2 * E * mass);
    }
}

/**
 * @brief Get particle density.
 *
 * @param dist pointer to afsi data.
 * @param spatial_index spatial index where density is queried.
 * @param vol cell volume [m^3].
 *
 * @return density [m^-3] or zero if type was unknown.
 */
real afsi_get_density(afsi_data* dist, size_t spatial_index, real vol) {
    if(dist->type == 1) {
        real density = 0.0;
        for(size_t ippara = 0; ippara < dist->beam->axes[5].n; ippara++) {
            for(size_t ipperp = 0; ipperp < dist->beam->axes[6].n; ipperp++) {
                size_t index = spatial_index
                             + ippara*dist->beam->strides[5]
                             + ipperp*dist->beam->strides[6];

                density += dist->beam->bins[index] / vol;
            }
        }
        return density;
    } else if(dist->type == 2) {
        return dist->thermal->density[spatial_index];
    } else {
        return 0;
    }
}

/**
 * @brief Get physical volume of a distribution cell.
 *
 * @param dist distribution.
 * @param r radial coordinate at the center of the cell.
 * @return cell volume [m^3].
 */
real afsi_get_volume(afsi_data* dist, real r) {
    real dr, dphi, dz;
    if(dist->type == 1) {
        dr =  (dist->beam->axes[0].max - dist->beam->axes[0].min)
             / dist->beam->axes[0].n;
        dz =  (dist->beam->axes[2].max - dist->beam->axes[2].min)
             / dist->beam->axes[2].n;
        dphi =  (dist->beam->axes[1].max - dist->beam->axes[1].min)
              / dist->beam->axes[1].n;
    }
    else {
        dr = (dist->thermal->max_r - dist->thermal->min_r) / dist->thermal->n_r;
        dz = (dist->thermal->max_z - dist->thermal->min_z) / dist->thermal->n_z;
        dphi = (dist->thermal->max_phi - dist->thermal->min_phi)
              / dist->thermal->n_phi;
    }
    return dphi*dz*dr*r;
}
