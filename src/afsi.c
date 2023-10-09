/**
 * @file afsi.c
 * @brief ASCOT Fusion Source Integrator AFSI
 */
#define _XOPEN_SOURCE 500 /**< rand48 requires POSIX 1995 standard */

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

void afsi_sample_reactant_velocities(
    afsi_data* react1, afsi_data* react2, real m1, real m2, int n,
    int iR, int iphi, int iz, real* v1x, real* v1y,
    real* v1z, real* v2x, real* v2y, real* v2z);
void afsi_compute_product_momenta(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* v1x, real* v1y, real* v1z, real* v2x, real* v2y, real* v2z,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2);
void afsi_sample_5D(dist_5D_data* dist, int n, int iR, int iphi, int iz,
                    real* ppara, real* pperp);
void afsi_sample_thermal(afsi_thermal_data* data, real mass, int n, int iR,
                         int iphi, int iz, real* ppara, real* pperp);
real afsi_get_density(afsi_data* dist, int iR, int iphi, int iz);
real afsi_get_volume(afsi_data* dist, int iR);

/**
 * @brief Calculate fusion source from two arbitrary ion distributions.
 *
 * Possible fusion reactions are listed below in a format
 * reactant 1 + reactant 2 -> product1 + product2:
 *
 * 1: D + T   -> He4 + n
 * 2: D + He3 -> He4 + p
 * 3: D + D   -> T + p
 * 4: D + D   -> He3 + n
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
void afsi_run(sim_offload_data* sim, int reaction, int n, afsi_data* react1,
              afsi_data* react2, real mult,
              dist_5D_offload_data* prod1_offload_data,
              dist_5D_offload_data* prod2_offload_data,
              real* prod1_offload_array, real* prod2_offload_array) {
    /* QID for this run */
    char qid[11];
    hdf5_generate_qid(qid);
    strcpy(sim->qid, qid);

    int mpi_rank = 0; /* AFSI does not support MPI */
    print_out0(VERBOSE_MINIMAL, mpi_rank, "AFSI5\n");

#ifdef GIT_VERSION
    print_out0(VERBOSE_MINIMAL, mpi_rank,
               "Tag %s\nBranch %s\n\n", GIT_VERSION, GIT_BRANCH);
#else
    print_out0(VERBOSE_MINIMAL, mpi_rank, "Not under version control\n\n");
#endif

    dist_5D_data prod1, prod2;
    dist_5D_init(&prod1, prod1_offload_data, prod1_offload_array);
    dist_5D_init(&prod2, prod2_offload_data, prod2_offload_array);

    simulate_init_offload(sim);
    random_init(&rdata, time((NULL)));
    sim_data sim_data;
    strcpy(sim->hdf5_out, sim->hdf5_in);
    sim_init(&sim_data, sim);

    if( hdf5_interface_init_results(sim, qid, "afsi") ) {
        print_out0(VERBOSE_MINIMAL, mpi_rank,
                   "\nInitializing output failed.\n"
                   "See stderr for details.\n");
        /* Free offload data and terminate */
        abort();
    }

    real m1=0, m2=0, mprod1=0, mprod2=0, Q=0;
    switch(reaction) {
        case 1: /* DT */
            m1     = 3.344e-27; // D
            m2     = 5.008e-27; // T
            mprod1 = 6.645e-27; // He4
            mprod2 = 1.675e-27; // n
            Q      = 17.6e6*CONST_E;
            break;
        case 2: /* D-He3 */
            m1     = 3.344e-27; // D
            m2     = 5.008e-27; // He3
            mprod1 = 6.645e-27; // He4
            mprod2 = 1.673e-27; // p
            Q      = 18.3e6*CONST_E;
            break;
        case 3: /* DDp */
            m1     = 3.344e-27; // D
            m2     = 3.344e-27; // D
            mprod1 = 5.008e-27; // T
            mprod2 = 1.673e-27; // p
            Q      = 4.03e6*CONST_E;
            break;
        case 4: /* DDn */
            m1     = 3.344e-27; // D
            m2     = 3.344e-27; // D
            mprod1 = 5.008e-27; // He3
            mprod2 = 1.675e-27; // n
            Q      = 3.27e6*CONST_E;
            break;
    }

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
        real* v1x = (real*) malloc(n*sizeof(real));
        real* v1y = (real*) malloc(n*sizeof(real));
        real* v1z = (real*) malloc(n*sizeof(real));
        real* v2x = (real*) malloc(n*sizeof(real));
        real* v2y = (real*) malloc(n*sizeof(real));
        real* v2z = (real*) malloc(n*sizeof(real));
        real* ppara1 = (real*) malloc(n*sizeof(real));
        real* pperp1 = (real*) malloc(n*sizeof(real));
        real* ppara2 = (real*) malloc(n*sizeof(real));
        real* pperp2 = (real*) malloc(n*sizeof(real));

        real vol = afsi_get_volume(react1, iR);
        for(int iphi = 0; iphi < n_phi; iphi++) {
            for(int iz = 0; iz < n_z; iz++) {
                real density1 = afsi_get_density(react1, iR, iphi, iz);
                real density2 = afsi_get_density(react2, iR, iphi, iz);
                if(density1 > 0 && density2 > 0) {
                    afsi_sample_reactant_velocities(
                        react1, react2, m1, m2, n, iR, iphi, iz,
                        v1x, v1y, v1z, v2x, v2y, v2z);
                    for(int i = 0; i < n; i++) {
                        real vcom2 =   (v1x[i] - v2x[i]) * (v1x[i] - v2x[i])
                                     + (v1y[i] - v2y[i]) * (v1y[i] - v2y[i])
                                     + (v1z[i] - v2z[i]) * (v1z[i] - v2z[i]);

                        real E_keV = 0.5 * ( m1 * m2 ) / ( m1 + m2 ) * vcom2
                            / (1.e3 * CONST_E);

                        real weight = density1 * density2 * sqrt(vcom2)
                            * boschhale_sigma(reaction, E_keV)/n*vol;

                        afsi_compute_product_momenta(
                            i, m1, m2, mprod1, mprod2, Q,
                            v1x, v1y, v1z, v2x, v2y, v2z,
                            ppara1, pperp1, ppara2, pperp2);

                        int ippara = floor(
                            (ppara1[i] - prod1.min_ppara) * prod1.n_ppara
                            / ( prod1.max_ppara - prod1.min_ppara ) );
                        int ipperp = floor(
                            (pperp1[i] - prod1.min_pperp) * prod1.n_pperp
                            / ( prod1.max_pperp - prod1.min_pperp ) );
                        prod1.histogram[dist_5D_index(
                                iR, iphi, iz, ippara, ipperp, 0, 0,
                                prod1.n_phi, prod1.n_z, prod1.n_ppara,
                                prod1.n_pperp, 1, 1)] += weight * mult;

                        ippara = floor(
                            (ppara2[i] - prod2.min_ppara) * prod2.n_ppara
                            / ( prod2.max_ppara - prod2.min_ppara ) );
                        ipperp = floor(
                            (pperp2[i] - prod2.min_pperp) * prod2.n_pperp
                            / ( prod2.max_pperp - prod2.min_pperp ) );
                        prod2.histogram[dist_5D_index(
                                iR, iphi, iz, ippara, ipperp, 0, 0,
                                prod2.n_phi, prod2.n_z, prod2.n_ppara,
                                prod2.n_pperp, 1, 1)] += weight * mult;
                    }
                }
                else {
                    /* Do nothing */
                }
            }
        }
        free(v1x);
        free(v1y);
        free(v1z);
        free(v2x);
        free(v2y);
        free(v2z);
        free(ppara1);
        free(ppara2);
        free(pperp1);
        free(pperp2);
    }

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
    if(H5LTmake_dataset_double(reactiondata, "q",      1, &size, &q)) {
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

    print_out0(VERBOSE_MINIMAL, mpi_rank, "\nDone\n");
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
 * @param v1x array where first velocity components of react1 will be stored.
 * @param v1y array where second velocity components of react1 will be stored.
 * @param v1z array where third velocity components of react1 will be stored.
 * @param v2x array where first velocity components of react2 will be stored.
 * @param v2y array where second velocity components of react2 will be stored.
 * @param v2z array where third velocity components of react2 will be stored.
 */
void afsi_sample_reactant_velocities(
    afsi_data* react1, afsi_data* react2, real m1, real m2, int n, int iR,
    int iphi, int iz, real* v1x, real* v1y, real* v1z, real* v2x, real* v2y,
    real* v2z) {
    real* ppara1 = (real*) malloc(n*sizeof(real));
    real* pperp1 = (real*) malloc(n*sizeof(real));
    real* ppara2 = (real*) malloc(n*sizeof(real));
    real* pperp2 = (real*) malloc(n*sizeof(real));

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
    for(int i = 0; i < n; i++) {
        real rx = 2*round(random_uniform(rdata))-1;
        real ry = 2*round(random_uniform(rdata))-1;
        real rz = random_uniform(rdata);
        v1x[i] = rx * pperp1[i]/m1 * sqrt(rz);
        v1y[i] = ry * pperp1[i]/m1 * sqrt(1-rz);
        v1z[i] = ppara1[i]/m1;

        rx = 2*round(random_uniform(rdata))-1;
        ry = 2*round(random_uniform(rdata))-1;
        rz = random_uniform(rdata);
        v2x[i] = rx * pperp2[i]/m2 * sqrt(rz);
        v2y[i] = ry * pperp2[i]/m2 * sqrt(1-rz);
        v2z[i] = ppara2[i]/m2;
    }
    free(ppara1);
    free(pperp1);
    free(ppara2);
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
 * @param v1x reactant 1 velocity x components.
 * @param v1y reactant 1 velocity y components.
 * @param v1z reactant 1 velocity z components.
 * @param v2x reactant 2 velocity x components.
 * @param v2y reactant 2 velocity y components.
 * @param v2z reactant 2 velocity z components.
 * @param ppara1 array where parallel momentum of product 1 is stored.
 * @param pperp1 array where perpendicular momentum of product 1 is stored.
 * @param ppara2 array where parallel momentum of product 2 is stored.
 * @param pperp2 array where perpendicular momentum of product 2 is stored.
 */
void afsi_compute_product_momenta(
    int i, real m1, real m2, real mprod1, real mprod2, real Q,
    real* v1x, real* v1y, real* v1z, real* v2x, real* v2y, real* v2z,
    real* ppara1, real* pperp1, real* ppara2, real* pperp2) {

    // Velocity of the system's center of mass
    real v_cm[3];
    v_cm[0]  = ( m1 * v1x[i] + m2 * v2x[i] ) / ( m1 + m2 );
    v_cm[1]  = ( m1 * v1y[i] + m2 * v2y[i] ) / ( m1 + m2 );
    v_cm[2]  = ( m1 * v1z[i] + m2 * v2z[i] ) / ( m1 + m2 );

    // Total kinetic energy after the reaction in CM frame
    real ekin = Q
        + 0.5 * m1 * (   (v1x[i] - v_cm[0])*(v1x[i] - v_cm[0])
                       + (v1y[i] - v_cm[1])*(v1y[i] - v_cm[1])
                       + (v1z[i] - v_cm[2])*(v1z[i] - v_cm[2]) )
        + 0.5 * m2 * (   (v2x[i] - v_cm[0])*(v2x[i] - v_cm[0])
                       + (v2y[i] - v_cm[1])*(v2y[i] - v_cm[1])
                       + (v2z[i] - v_cm[2])*(v2z[i] - v_cm[2]) );

    // Speed and velocity of product 2 in CM frame
    real rn1 = random_uniform(rdata);
    real rn2 = random_uniform(rdata);
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
    ppara1[i] = vprod1[2] * mprod1;
    pperp1[i] = sqrt( vprod1[0]*vprod1[0] + vprod1[1]*vprod1[1] ) * mprod1;
    ppara2[i] = vprod2[2] * mprod2;
    pperp2[i] = sqrt( vprod2[0]*vprod2[0] + vprod2[1]*vprod2[1] ) * mprod2;
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
                cumdist[0] = dist->histogram[dist_5D_index(iR, iphi, iz,
                    0, 0, 0, 0, dist->n_phi, dist->n_z, dist->n_ppara,
                    dist->n_pperp, 1, 1)];
            } else {
                cumdist[ippara*dist->n_pperp+ipperp] =
                    cumdist[ippara*dist->n_pperp+ipperp-1]
                    + dist->histogram[dist_5D_index(iR, iphi, iz,
                        ippara, ipperp, 0, 0, dist->n_phi, dist->n_z,
                        dist->n_ppara, dist->n_pperp, 1, 1)];
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
        real w = 0;
        real r1;
        real r2;

        while(w <= 1.0) {
            r1 = random_uniform(rdata);
            r2 = random_uniform(rdata);
            w = sqrt(r1*r1 + r2*r2);
        }

        real r3 = random_uniform(rdata);
        real r4 = random_uniform(rdata);

        real E = -temp * (r1*r1 * log(r3) / w + log(r4));
        real absv = sqrt(2*E*CONST_E/mass);

        r1 = random_uniform(rdata);
        r2 = random_uniform(rdata);

        real theta = CONST_2PI * r1;
        real phi = acos(1 - 2*r2);
        real vx = absv * sin(phi) * cos(theta);
        real vy = absv * sin(phi) * sin(theta);
        real vz = absv * cos(phi);

        ppara[i] = vz * mass;
        pperp[i] = sqrt(vx*vx + vy*vy) * mass;
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
                density += dist->dist_5D->histogram[dist_5D_index(iR, iphi, iz,
                           ippara, ipperp, 0, 0, dist->dist_5D->n_phi,
                           dist->dist_5D->n_z, dist->dist_5D->n_ppara,
                           dist->dist_5D->n_pperp, 1, 1)] / vol;
            }
        }
        return density;
    }

    else if(dist->type == 2) {
        return dist->dist_thermal->density[
            iR*(  dist->dist_thermal->n_phi*dist->dist_thermal->n_z
                + iphi*dist->dist_thermal->n_z)+iz];
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
