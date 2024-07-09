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
#include "consts.h"
#include "random.h"
#include "simulate.h"
#include "boschhale.h"
#include "diag/dist_5D.h"
#include "hdf5_interface.h"
#include "hdf5io/hdf5_helpers.h"
#include "hdf5io/hdf5_dist.h"
#include "afsi.h"

/** Random number generator used by AFSI */
random_data rdata;

void afsi_sample_reactant_momenta(
    afsi_data* react1, afsi_data* react2, real m1, real m2, int n, int iR,
    int iphi, int iz, real* ppara1, real* pperp1, real* ppara2, real* pperp2);
void afsi_compute_product_momenta(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* pparaprod1, real* pperpprod1, real* pparaprod2, real* pperpprod2);
void afsi_sample_5D(dist_5D_data* dist, int n, int iR, int iphi, int iz,
                    real* ppara, real* pperp);
void afsi_sample_thermal(afsi_thermal_data* data, real mass, int n, int iR,
                         int iphi, int iz, real* ppara, real* pperp);
real afsi_get_density(afsi_data* dist, int iR, int iphi, int iz);
real afsi_get_volume(afsi_data* dist, int iR);

/**
 * @brief Calculate fusion source from two arbitrary ion distributions
 *
 * Inputs and outputs are expected to have same physical (R, phi, z)
 * dimensions.
 *
 * @param sim pointer to simulation offload data
 * @param reaction fusion reaction type, see the description
 * @param n number of Monte Carlo samples to be used
 * @param react1 reactant 1 distribution data
 * @param react2 reactant 2 distribution data
 * @param mult factor by which the output is scaled
 * @param prod1_offload_data distribution data for product 1 output
 * @param prod2_offload_data distribution data for product 2 output
 * @param prod1_offload_array array where product 1 distribution is stored
 * @param prod2_offload_array array where product 2 distribution is stored
 */
void afsi_run(sim_offload_data* sim, Reaction reaction, int n,
              afsi_data* react1, afsi_data* react2, real mult,
              dist_5D_offload_data* prod1_offload_data,
              dist_5D_offload_data* prod2_offload_data,
              real* prod1_offload_array, real* prod2_offload_array) {
    /* QID for this run */
    char qid[11];
    hdf5_generate_qid(qid);
    strcpy(sim->qid, qid);

    int mpi_rank = 0, mpi_root = 0; /* AFSI does not support MPI */
    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root, "AFSI5\n");
    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);

    dist_5D_data prod1, prod2;
    dist_5D_init(&prod1, prod1_offload_data, prod1_offload_array);
    dist_5D_init(&prod2, prod2_offload_data, prod2_offload_array);

    simulate_init_offload(sim);
    random_init(&rdata, time((NULL)));
    sim_data sim_data;
    strcpy(sim->hdf5_out, sim->hdf5_in);
    sim_init(&sim_data, sim);

    if( hdf5_interface_init_results(sim, qid, "afsi") ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root,
                   "\nInitializing output failed.\n"
                   "See stderr for details.\n");
        /* Free offload data and terminate */
        abort();
    }

    real m1, q1, m2, q2, mprod1, qprod1, mprod2, qprod2, Q;
    boschhale_reaction(
        reaction, &m1, &q1, &m2, &q2, &mprod1, &qprod1, &mprod2, &qprod2, &Q);

    int n_r=0, n_phi=0, n_z=0;
    if(react1->type == 1) {
        n_r   = react1->dist_5D->n_r;
        n_phi = react1->dist_5D->n_phi;
        n_z   = react1->dist_5D->n_z;
    }
    else if(react1->type == 2) {
        n_r   = react1->dist_thermal->n_r;
        n_phi = react1->dist_thermal->n_phi;
        n_z   = react1->dist_thermal->n_z;
    }

    #pragma omp parallel for
    for(int iR = 0; iR < n_r; iR++) {
        real* ppara1     = (real*) malloc(n*sizeof(real));
        real* pperp1     = (real*) malloc(n*sizeof(real));
        real* ppara2     = (real*) malloc(n*sizeof(real));
        real* pperp2     = (real*) malloc(n*sizeof(real));
        real* pparaprod1 = (real*) malloc(n*sizeof(real));
        real* pperpprod1 = (real*) malloc(n*sizeof(real));
        real* pparaprod2 = (real*) malloc(n*sizeof(real));
        real* pperpprod2 = (real*) malloc(n*sizeof(real));

        real vol = afsi_get_volume(react1, iR);
        for(int iphi = 0; iphi < n_phi; iphi++) {
            for(int iz = 0; iz < n_z; iz++) {
                real density1 = afsi_get_density(react1, iR, iphi, iz);
                real density2 = afsi_get_density(react2, iR, iphi, iz);
                if(density1 > 0 && density2 > 0) {
                    afsi_sample_reactant_momenta(
                        react1, react2, m1, m2, n, iR, iphi, iz,
                        ppara1, pperp1, ppara2, pperp2);
                    for(int i = 0; i < n; i++) {
                        real vcom2;
                        afsi_compute_product_momenta(
                            i, m1, m2, mprod1, mprod2, Q,
                            ppara1, pperp1, ppara2, pperp2, &vcom2,
                            pparaprod1, pperpprod1, pparaprod2, pperpprod2);
                        real E = 0.5 * ( m1 * m2 ) / ( m1 + m2 ) * vcom2;

                        real weight = density1 * density2 * sqrt(vcom2)
                            * boschhale_sigma(reaction, E)/n*vol;

                        int ippara = floor(
                            (pparaprod1[i] - prod1.min_ppara) * prod1.n_ppara
                            / ( prod1.max_ppara - prod1.min_ppara ) );
                        int ipperp = floor(
                            (pperpprod1[i] - prod1.min_pperp) * prod1.n_pperp
                            / ( prod1.max_pperp - prod1.min_pperp ) );
                        if( 0 <= ippara && ippara < prod1.n_ppara &&
                            0 <= ipperp && ipperp < prod1.n_pperp) {
                            prod1.histogram[dist_5D_index(
                                    iR, iphi, iz, ippara, ipperp, 0, 0,
                                    prod1.step_6, prod1.step_5, prod1.step_4,
                                    prod1.step_3, prod1.step_2, prod1.step_1)]
                                += weight * mult;
                        }

                        ippara = floor(
                            (pparaprod2[i] - prod2.min_ppara) * prod2.n_ppara
                            / ( prod2.max_ppara - prod2.min_ppara ) );
                        ipperp = floor(
                            (pperpprod2[i] - prod2.min_pperp) * prod2.n_pperp
                            / ( prod2.max_pperp - prod2.min_pperp ) );
                        if( 0 <= ippara && ippara < prod2.n_ppara &&
                            0 <= ipperp && ipperp < prod2.n_pperp) {
                            prod2.histogram[dist_5D_index(
                                    iR, iphi, iz, ippara, ipperp, 0, 0,
                                    prod2.step_6, prod2.step_5, prod2.step_4,
                                    prod2.step_3, prod2.step_2, prod2.step_1)]
                                += weight * mult;
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
        free(pparaprod1);
        free(pparaprod2);
        free(pperpprod1);
        free(pperpprod2);
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
    if( hdf5_dist_write_5D(f, path, prod1_offload_data, prod1_offload_array) ) {
        print_err("Warning: 5D distribution could not be written.\n");
    }
    sprintf(path, "/results/afsi_%s/prod2dist5d", sim->qid);
    if( hdf5_dist_write_5D(f, path, prod2_offload_data, prod2_offload_array) ) {
        print_err("Warning: 5D distribution could not be written.\n");
    }
    if(hdf5_close(f)) {
        print_err("Failed to close the file.\n");
        abort();
    }

    print_out0(VERBOSE_MINIMAL, mpi_rank, mpi_root, "\nDone\n");
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
void afsi_sample_reactant_momenta(
    afsi_data* react1, afsi_data* react2, real m1, real m2, int n, int iR,
    int iphi, int iz, real* ppara1, real* pperp1, real* ppara2, real* pperp2) {

    if(react1->type == 1) {
        afsi_sample_5D(react1->dist_5D, n, iR, iphi, iz, ppara1, pperp1);
    }
    else if(react1->type == 2) {
        afsi_sample_thermal(
            react1->dist_thermal, m1, n, iR, iphi, iz, ppara1, pperp1);
    }

    if(react2->type == 1) {
        afsi_sample_5D(react2->dist_5D, n, iR, iphi, iz, ppara2, pperp2);
    }
    else if(react2->type == 2) {
        afsi_sample_thermal(
            react2->dist_thermal, m2, n, iR, iphi, iz, ppara2, pperp2);
    }
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
void afsi_compute_product_momenta(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2, real* vcom2,
    real* pparaprod1, real* pperpprod1, real* pparaprod2, real* pperpprod2) {

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

    // ppara and pperp
    pparaprod1[i] = vprod1[2] * mprod1;
    pperpprod1[i] = sqrt( vprod1[0]*vprod1[0] + vprod1[1]*vprod1[1] ) * mprod1;
    pparaprod2[i] = vprod2[2] * mprod2;
    pperpprod2[i] = sqrt( vprod2[0]*vprod2[0] + vprod2[1]*vprod2[1] ) * mprod2;
}

/**
 * @brief Sample ppara and pperp from a 5D distribution.
 *
 * @param dist pointer to the distribution data.
 * @param n number of values to be sampled.
 * @param iR R index where sampling is done.
 * @param iphi phi index where sampling is done.
 * @param iz z index where sampling is done.
 * @param ppara pointer to array where sampled parallel momenta are stored.
 * @param pperp pointer to array where sampled perpedicular momenta are stored.
 */
void afsi_sample_5D(dist_5D_data* dist, int n, int iR, int iphi, int iz,
                    real* ppara, real* pperp) {
    real* cumdist = (real*) malloc(dist->n_ppara*dist->n_pperp*sizeof(real));

    for(int ippara = 0; ippara < dist->n_ppara; ippara++) {
        for(int ipperp = 0; ipperp < dist->n_pperp; ipperp++) {
            if(ippara == 0 && ipperp == 0) {
                cumdist[0] = dist->histogram[dist_5D_index(
                        iR, iphi, iz, 0, 0, 0, 0, dist->step_6, dist->step_5,
                        dist->step_4, dist->step_3, dist->step_2,
                        dist->step_1)];
            } else {
                cumdist[ippara*dist->n_pperp+ipperp] =
                    cumdist[ippara*dist->n_pperp+ipperp-1]
                    + dist->histogram[dist_5D_index(
                        iR, iphi, iz, ippara, ipperp, 0, 0, dist->step_6,
                        dist->step_5, dist->step_4, dist->step_3, dist->step_2,
                        dist->step_1)];
            }
        }
    }
    for(int ippara = 0; ippara < dist->n_ppara; ippara++) {
        for(int ipperp = 0; ipperp < dist->n_pperp; ipperp++) {
            cumdist[ippara*dist->n_pperp+ipperp] /=
                cumdist[dist->n_ppara*dist->n_pperp-1];
        }
    }

    for(int i = 0; i < n; i++) {
        real r = random_uniform(rdata);
        for(int j = 0; j < dist->n_ppara*dist->n_pperp; j++) {
            if(cumdist[j] > r) {
                pperp[i] = dist->min_pperp + (j % dist->n_pperp + 0.5)
                    * (dist->max_pperp - dist->min_pperp) / dist->n_pperp;
                ppara[i] = dist->min_ppara + (j / dist->n_pperp + 0.5)
                    * (dist->max_ppara - dist->min_ppara) / dist->n_ppara;
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
 * @param n number of values to be sampled.
 * @param iR R index where sampling is done.
 * @param iphi phi index where sampling is done.
 * @param iz z index where sampling is done.
 * @param ppara pointer to array where sampled parallel momenta are stored.
 * @param pperp pointer to array where sampled perpedicular momenta are stored.
 */
void afsi_sample_thermal(afsi_thermal_data* data, real mass, int n, int iR,
                         int iphi, int iz, real* ppara, real* pperp) {
    int ind = iR * (data->n_phi * data->n_z) + iphi * data->n_z + iz;
    real temp = data->temperature[ind];

    for(int i = 0; i < n; i++) {
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
 * @param iR R index at which density is queried.
 * @param iphi phi index at which density is queried.
 * @param iz z index at which density is queried.
 * @return density [m^-3] or zero if type was unknown.
 */
real afsi_get_density(afsi_data* dist, int iR, int iphi, int iz) {
    if(dist->type == 1) {
        real vol = afsi_get_volume(dist, iR);

        real density = 0.0;
        for(int ippara = 0; ippara < dist->dist_5D->n_ppara; ippara++) {
            for(int ipperp = 0; ipperp < dist->dist_5D->n_pperp; ipperp++) {
                density += dist->dist_5D->histogram[dist_5D_index(
                        iR, iphi, iz, ippara, ipperp, 0, 0,
                        dist->dist_5D->step_6, dist->dist_5D->step_5,
                        dist->dist_5D->step_4, dist->dist_5D->step_3,
                        dist->dist_5D->step_2, dist->dist_5D->step_1)] / vol;
            }
        }
        return density;
    }

    else if(dist->type == 2) {
        return dist->dist_thermal->density[
              iR * dist->dist_thermal->n_phi * dist->dist_thermal->n_z
            + iphi * dist->dist_thermal->n_z + iz];
    }

    else {
        return 0;
    }
}

/**
 * @brief Get physical volume of a distribution cell.
 *
 * @param dist distribution.
 * @param iR radial index of the cell.
 * @return cell volume [m^3].
 */
real afsi_get_volume(afsi_data* dist, int iR) {
    real dR, dz;

    if(dist->type == 1) {
        dR = (dist->dist_5D->max_r - dist->dist_5D->min_r) / dist->dist_5D->n_r;
        dz = (dist->dist_5D->max_z - dist->dist_5D->min_z) / dist->dist_5D->n_z;
        return CONST_2PI*(dist->dist_5D->min_r + iR * dR + 0.5*dR)*dR*dz;
    }

    else if(dist->type == 2) {
        dR = (dist->dist_thermal->max_r - dist->dist_thermal->min_r)
            / dist->dist_thermal->n_r;
        dz = (dist->dist_thermal->max_z - dist->dist_thermal->min_z)
            / dist->dist_thermal->n_z;
        return CONST_2PI*(dist->dist_thermal->min_r + iR * dR + 0.5*dR)*dR*dz;
    }
    return 0;
}

/**
 * @brief Test distribution.
 *
 * @param dist1 distribution to be tested.
 */
void afsi_test_dist(dist_5D_data* dist1) {
    printf("%d %le %le\n", dist1->n_r, dist1->min_r, dist1->max_r);
    printf("%d %le %le\n", dist1->n_phi, dist1->min_phi, dist1->max_phi);
    printf("%d %le %le\n", dist1->n_z, dist1->min_z, dist1->max_z);
    printf("%d %le %le\n", dist1->n_ppara, dist1->min_ppara, dist1->max_ppara);
    printf("%d %le %le\n", dist1->n_pperp, dist1->min_pperp, dist1->max_pperp);
    printf("%d %le %le\n", dist1->n_time, dist1->min_time, dist1->max_time);
    printf("%d %le %le\n", dist1->n_q, dist1->min_q, dist1->max_q);

    real sum = 0.0;

    int ncell = dist1->n_r * dist1->n_phi * dist1->n_z * dist1->n_ppara
        * dist1->n_pperp * dist1->n_time * dist1->n_q;
    for(int i = 0; i < ncell; i++) {
        sum += dist1->histogram[i];
    }

    printf("%le\n", sum);
}

/**
 * @brief Test thermal source.
 */
void afsi_test_thermal() {
    afsi_thermal_data data;

    data.n_r     = 1;
    data.min_r   = 0.1;
    data.max_r   = 1;
    data.n_phi   = 1;
    data.min_phi = 0;
    data.max_phi = 360;
    data.n_z     = 1;
    data.min_z   = -1;
    data.max_z   = 1;

    real temperature = 1e3;
    real density = 1e19;

    data.temperature = &temperature;
    data.density = &density;

    int n = 100000;
    real* ppara = (real*) malloc(n*sizeof(real));
    real* pperp = (real*) malloc(n*sizeof(real));

    afsi_sample_thermal(&data, 3.343e-27, n, 0, 0, 0, ppara, pperp);

    for(int i = 0; i < n; i++) {
        printf("%le %le\n", ppara[i], pperp[i]);
    }
}
