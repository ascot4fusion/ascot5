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
#include "utils/random.h"
#include "utils/ratecoefficients.h"

/**
 * Compute reaction products.
 *
 * @param reaction Reaction that takes place.
 * @param vol Physical volume of the cell where reactions take place [m^3].
 * @param react1 First reactant.
 * @param react2 Second reactant.
 * @param prod1 First product.
 * @param prod2 Second product.
 */
void fusion_source_compute_product(
    Reaction reaction, real vol, MarkerGyroOrbit *react1,
    MarkerGyroOrbit *react2, MarkerGyroOrbit *prod1, MarkerGyroOrbit *prod2);

/**
 * Sample markers from an arbitrary distribution.
 *
 * @param hist Histogram to sample from.
 * @param spatial_index Spatial index within the histogram where sampling is
 *        done.
 * @param vol Physical volume of the cell where sampling takes place [m^3].
 * @param sampled Sampled markers.
 * @param density Density of the volume where sampling takes place [m^-3].
 */
void fusion_source_sample_beam(
    DiagHist *hist, size_t spatial_index, real vol, MarkerGyroOrbit *sampled,
    real density[1]);

/**
 * Sample markers from a thermal (Maxwellian) distribution.
 *
 * @param sim Simulation data.
 * @param species_index Index of the plasma species to sample.
 * @param sampled Sampled markers.
 * @param density Density of the volume where sampling takes place [m^-3].
 */
void fusion_source_sample_thermal(
    Simulation *sim, size_t species_index, MarkerGyroOrbit *sampled,
    real density[1]);

#endif
