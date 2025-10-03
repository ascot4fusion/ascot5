/**
 * @file gctransform.h
 * Module for performing guiding center transformations
 *
 * The guiding center transformation is done both ways between particle
 * phase-space [r, phi, z, pr, pphi, pz] and guiding center phase-space
 * [R, Phi, Z, ppar, mu, zeta] to first order (both in spatial and momentum
 * space coordinates).
 *
 * The guiding genter motion is defined from a basis {bhat, e1, e2}, where bhat
 * is magnetic field unit vector. e1 and e2 are chosen so that
 * e1 = (bhat x z) / |bhat x z| and bhat x e1 = e2. Thus we are assuming that
 * magnetic field is *never* parallel to z axis.
 *
 * This module works as follows:
 *
 * - particle to guiding center transformation is accomplished by calling
 *   gctransform_particle2guidingcenter()
 *
 * - guiding center to particle transformation is accomplished by calling
 *   gctransform_guidingcenter2particle() which gives [r, phi, z, ppar, mu,
 *   zeta] in particle coordinates. To obtain [r, phi, z, pr, pphi, pz],
 *   first evaluate magnetic field at particle position and then call
 *   gctransform_pparmuzeta2pRpphipz().
 *
 * The transformation is relativistic.
 *
 * Reference: "Guiding-center transformation of the radiation-reaction force in
 * a nonuniform magnetic field", E. Hirvijoki et. al.
 * https://arxiv.org/pdf/1412.1966.pdf
 */
#ifndef GCTRANSFORM_H
#define GCTRANSFORM_H

#include "ascot5.h"
#include "offload.h"

/**
 * Set the order of the transformation.
 *
 * @param order Order which is either zero or one.
 */
void gctransform_setorder(int order);

/**
 * Transform particle to guiding center  phase space.
 *
 * The transformation is done from coordinates [r, phi, z, pr, pphi, pz] to
 * [R, Phi, Z, ppar, mu, zeta]. If particle is neutral, the transformation
 * will be in zeroth order both in real and momentum space as that
 * transformation is valid when q = 0.
 *
 * @param mass Mass [kg].
 * @param charge Charge [C].
 * @param b_db Gradient of magnetic field vector at particle position.
 * @param r Particle R coordinate [m].
 * @param phi Particle phi coordinate [rad].
 * @param z Particle z coordinate [m].
 * @param pr Particle momentum R component [kg m/s].
 * @param pphi Particle momentum phi component [kg m/s].
 * @param pz Particle momentum z component [kg m/s].
 * @param R Pointer to guiding center R coordinate [m].
 * @param Phi Pointer to guiding center phi coordinate [rad].
 * @param Z Pointer to guiding center z coordinate [m].
 * @param ppar Pointer to guiding center parallel momentum [kg m/s].
 * @param mu Pointer to guiding center magnetic moment [J/T].
 * @param zeta Pointer to guiding center gyroangle [rad].
 */
DECLARE_TARGET_SIMD
void gctransform_particle2guidingcenter(
    real mass, real charge, real *b_db, real r, real phi, real z, real pr,
    real pphi, real pz, real *R, real *Phi, real *Z, real *ppar, real *mu,
    real *zeta);

/**
 * Transform guiding center to particle phase space.
 *
 * The transformation is done from coordinates [R, Phi, Z, ppar, mu] to
 * [r, phi, z, ppar_prt, mu_prt, zeta_prt]. If particle is neutral,
 * the transformation will be in zeroth order both in real and momentum space
 * as that transformation is valid when q = 0.
 *
 * @param mass Mass [kg].
 * @param charge Charge [C].
 * @param b_db Gradient of magnetic field vector at guiding center position.
 * @param R Guiding center R coordinate [m].
 * @param Phi Guiding center phi coordinate [rad].
 * @param Z Guiding center z coordinate [m].
 * @param ppar Guiding center parallel momentum [kg m/s].
 * @param mu Guiding center magnetic moment [J/T].
 * @param zeta Guiding center gyroangle [rad].
 * @param r Pointer to particle R coordinate [m].
 * @param phi Pointer to particle phi coordinate [rad].
 * @param z Pointer to particle z coordinate [m].
 * @param pparprt Pointer to particle parallel momentum [kg m/s].
 * @param muprt Pointer to particle magnetic moment [J/T].
 * @param zetaprt Pointer to particle gyroangle [rad].
 */
DECLARE_TARGET_SIMD
void gctransform_guidingcenter2particle(
    real mass, real charge, real *b_db, real R, real Phi, real Z, real ppar,
    real mu, real zeta, real *r, real *phi, real *z, real *pparprt, real *muprt,
    real *zetaprt);

/**
 * Transform particle ppar, mu, and zeta to momentum vector.
 *
 * The transformation is done from coordinates [R, Phi, Z, ppar, mu] to
 * [r, phi, z, pr, pphi, pz]. The transformation is done to first order.
 *
 * @param mass Mass [kg].
 * @param charge Charge [C].
 * @param b_db Gradient of magnetic field vector at particle position.
 * @param phi Particle phi coordinate [rad].
 * @param ppar Particle parallel momentum [kg m/s].
 * @param mu Particle magnetic moment [J/T].
 * @param zeta Particle gyroangle [rad].
 * @param pr Pointer to particle momentum R-component [kg m/s].
 * @param pphi Pointer to particle momentum phi-component [kg m/s].
 * @param pz Pointer to particle momentum z-component [kg m/s].
 */
DECLARE_TARGET_SIMD
void gctransform_pparmuzeta2prpphipz(
    real mass, real charge, real *b_db, real phi, real ppar, real mu, real zeta,
    real *pr, real *pphi, real *pz);

#endif
