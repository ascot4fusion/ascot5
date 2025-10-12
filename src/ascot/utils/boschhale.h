/**
 * @file boschhale.h
 * Formulas for fusion cross-sections and thermal reactivities.
 *
 * The model is adapted from here:
 * https://www.doi.org/10.1088/0029-5515/32/4/I07
 */
#ifndef BOSCHHALE_H
#define BOSCHHALE_H

#include "defines.h"
#include "datatypes.h"

/**
 * Get masses and charges of particles participating in the reaction and
 * the released energy.
 *
 * @param reaction reaction enum.
 * @param m1 mass of the first reactant [kg].
 * @param q1 charge of the first reactant [C].
 * @param m2 mass of the second reactant [kg].
 * @param q2 charge of the second reactant [C].
 * @param mprod1 mass of the first product [kg].
 * @param qprod1 charge of the first product [C].
 * @param mprod2 mass of the second product [kg].
 * @param qprod2 charge of the second product [C].
 * @param Q energy released [J].
 */
void boschhale_reaction(
    Reaction reaction, real *m1, real *q1, real *m2, real *q2, real *mprod1,
    real *qprod1, real *mprod2, real *qprod2, real *Q);

/**
 * Estimate cross-section for a given fusion reaction.
 *
 * @param reaction Reaction for which the cross-section is estimated.
 * @param E Ion energy [J].
 *
 * @return Cross-section [m^2].
 */
real boschhale_sigma(Reaction reaction, real E);

/**
 * Estimate reactivity for a given fusion reaction.
 *
 * @param reaction Reaction for which the reactivity is estimated.
 * @param Ti Ion temperature [keV].
 *
 * @return Reactivity.
 */
real boschhale_sigmav(Reaction reaction, real Ti);

#endif
