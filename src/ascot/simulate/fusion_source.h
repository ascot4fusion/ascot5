/**
 * @file fusion_source.h
 * ASCOT Fusion Source Integrator AFSI.
 */
#ifndef FUSION_SOURCE_H
#define FUSION_SOURCE_H

#include "consts.h"
#include "data/diag.h"
#include "datatypes.h"
#include "defines.h"
#include "utils/boschhale.h"
#include "utils/random.h"
#include <math.h>

/**
 * Valid momentum space basis.
 */
typedef enum
{
    PPARPPERP,
    EKINXI
} mom_space_basis;

/**
 * Compute momenta of reaction products.
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
    size_t i, real m1, real m2, real mprod1, real mprod2, real Q,
    int prodmomspace, real *ppara1, real *pperp1, real *ppara2, real *pperp2,
    real *vcom2, real *prod1_p1, real *prod1_p2, real *prod2_p1,
    real *prod2_p2);

/**
 * Sample velocities from reactant distributions.
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
    Simulation *sim, FusionSource *afsi, real mass1, real mass2, real vol,
    size_t nsample, size_t i0, size_t i1, size_t i2, real r, real phi, real z,
    real time, real rho, real *density1, real *ppara1, real *pperp1,
    real *density2, real *ppara2, real *pperp2);

/**
 * Sample ppara and pperp from a 5D distribution.
 *
 * @param dist pointer to the distribution data.
 * @param nsample number of values to be sampled.
 * @param spatial_index spatial index where sampling is done.
 * @param ppara pointer to array where sampled parallel momenta are stored.
 * @param pperp pointer to array where sampled perpedicular momenta are stored.
 */
void afsi_sample_beam_2d(
    DiagHist *hist, real mass, real vol, size_t nsample, size_t i0, size_t i1,
    size_t i2, real *density, real *ppara, real *pperp);

/**
 * Sample ppara and pperp from a thermal (Maxwellian) population.
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
void afsi_sample_thermal_2d(
    Simulation *sim, int ispecies, real mass, size_t nsample, real r, real phi,
    real z, real time, real rho, real *density, real *pppara, real *ppperp);

#endif
