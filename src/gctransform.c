/**
 * Implements gctransform.h.
 */
#include "gctransform.h"
#include "ascot5.h"
#include "consts.h"
#include "math.h"
#include "offload.h"
#include "physlib.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Order to which guiding center transformation is done in momentum space. */
// static int GCTRANSFORM_ORDER = 1;
#ifdef _OPENACC
static int GCTRANSFORM_ORDER = 1;
#pragma acc declare copyin(GCTRANSFORM_ORDER)
#elif defined(_OPENMP)
DECLARE_TARGET
static int GCTRANSFORM_ORDER = 1;
DECLARE_TARGET_END
#else
static int GCTRANSFORM_ORDER = 1;
#endif


void gctransform_setorder(int order) { GCTRANSFORM_ORDER = order; }


void gctransform_particle2guidingcenter(
    real mass, real charge, real *b_db, real r, real phi, real z, real pr,
    real pphi, real pz, real *R, real *Phi, real *Z, real *ppar, real *mu,
    real *zeta)
{

    /* |B| */
    real Bnorm =
        sqrt(b_db[0] * b_db[0] + b_db[4] * b_db[4] + b_db[8] * b_db[8]);

    /* Guiding center transformation is more easily done in cartesian
     * coordinates so we switch to using those */
    real rpz[3] = {r, phi, z};
    real prpz[3] = {pr, pphi, pz};
    real xyz[3];
    real pxyz[3];
    real b_dbxyz[12];
    math_rpz2xyz(rpz, xyz);
    math_vec_rpz2xyz(prpz, pxyz, phi);
    math_jac_rpz2xyz(b_db, b_dbxyz, r, phi);

    /* bhat = Unit vector of B */
    real bhat[3] = {b_dbxyz[0] / Bnorm, b_dbxyz[4] / Bnorm, b_dbxyz[8] / Bnorm};

    /* Magnetic field norm gradient */
    real gradB[3];
    gradB[0] =
        (bhat[0] * b_dbxyz[1] + bhat[1] * b_dbxyz[5] + bhat[2] * b_dbxyz[9]);
    gradB[1] =
        (bhat[0] * b_dbxyz[2] + bhat[1] * b_dbxyz[6] + bhat[2] * b_dbxyz[10]);
    gradB[2] =
        (bhat[0] * b_dbxyz[3] + bhat[1] * b_dbxyz[7] + bhat[2] * b_dbxyz[11]);

    /* nabla x |B| */
    real curlB[3] = {
        b_dbxyz[10] - b_dbxyz[7], b_dbxyz[3] - b_dbxyz[9],
        b_dbxyz[5] - b_dbxyz[2]};

    /* Magnetic field torsion = bhat dot ( nabla X bhat ) which is equivalent to
       bhat dot ( nabla x |B| ) / |B|
       because nabla x |B| = |B| nabla x bhat + (nabla |B|) x bhat */
    real tau = math_dot(bhat, curlB) / Bnorm;

    /* Gradient of magnetic field unit vector */
    real nablabhat[9];
    nablabhat[0] = (b_dbxyz[1] - gradB[0] * bhat[0]) / Bnorm;
    nablabhat[1] = (b_dbxyz[2] - gradB[0] * bhat[1]) / Bnorm;
    nablabhat[2] = (b_dbxyz[3] - gradB[0] * bhat[2]) / Bnorm;
    nablabhat[3] = (b_dbxyz[5] - gradB[1] * bhat[0]) / Bnorm;
    nablabhat[4] = (b_dbxyz[6] - gradB[1] * bhat[1]) / Bnorm;
    nablabhat[5] = (b_dbxyz[7] - gradB[1] * bhat[2]) / Bnorm;
    nablabhat[6] = (b_dbxyz[9] - gradB[2] * bhat[0]) / Bnorm;
    nablabhat[7] = (b_dbxyz[10] - gradB[2] * bhat[1]) / Bnorm;
    nablabhat[8] = (b_dbxyz[11] - gradB[2] * bhat[2]) / Bnorm;

    /* Magnetic field curvature vector = bhat dot nablabhat.
       Note that nabla x bhat = tau bhat + bhat x kappa */
    real kappa[3];
    kappa[0] = (bhat[0] * b_dbxyz[1] + bhat[1] * b_dbxyz[5] +
                bhat[2] * b_dbxyz[9] - math_dot(bhat, gradB) * bhat[0]) /
               Bnorm;
    kappa[1] = (bhat[0] * b_dbxyz[2] + bhat[1] * b_dbxyz[6] +
                bhat[2] * b_dbxyz[10] - math_dot(bhat, gradB) * bhat[1]) /
               Bnorm;
    kappa[2] = (bhat[0] * b_dbxyz[3] + bhat[1] * b_dbxyz[7] +
                bhat[2] * b_dbxyz[11] - math_dot(bhat, gradB) * bhat[2]) /
               Bnorm;

    /* Zeroth order mu and ppar */
    real ppar0 = pxyz[0] * bhat[0] + pxyz[1] * bhat[1] + pxyz[2] * bhat[2];

    /* Gyrovector rhohat = bhat X perphat */
    real pperp[3] = {
        pxyz[0] - ppar0 * bhat[0], pxyz[1] - ppar0 * bhat[1],
        pxyz[2] - ppar0 * bhat[2]};
    real mu0 = math_dot(pperp, pperp) / (2 * mass * Bnorm);
    real perphat[3];
    math_unit(pperp, perphat);
    if (charge < 0)
    {
        perphat[0] = -perphat[0];
        perphat[1] = -perphat[1];
        perphat[2] = -perphat[2];
    }
    real rhohat[3];
    math_cross(bhat, perphat, rhohat);

    /* Double product of dyadic a1 = -(1/2) * (rhohat perphat + perphat rhohat)
       and gradb (gradient of magnetic field unit vector) */
    real a1ddotgradb =
        -0.5 * (2 * (rhohat[0] * perphat[0] * nablabhat[0] +
                     rhohat[1] * perphat[1] * nablabhat[4] +
                     rhohat[2] * perphat[2] * nablabhat[8]) +
                (rhohat[0] * perphat[1] + rhohat[1] * perphat[0]) *
                    (nablabhat[1] + nablabhat[3]) +
                (rhohat[0] * perphat[2] + rhohat[2] * perphat[0]) *
                    (nablabhat[2] + nablabhat[6]) +
                (rhohat[1] * perphat[2] + rhohat[2] * perphat[1]) *
                    (nablabhat[5] + nablabhat[7]));

    /* Double product of dyadic a2 = (1/4) * (perphat perphat - rhohat rhohat)
       and gradb (gradient of magnetic field unit vector) */
    real a2ddotgradb = 0.25 * (perphat[0] * perphat[0] * nablabhat[0] +
                               perphat[1] * perphat[0] * nablabhat[3] +
                               perphat[2] * perphat[0] * nablabhat[6] +
                               perphat[0] * perphat[1] * nablabhat[1] +
                               perphat[1] * perphat[1] * nablabhat[4] +
                               perphat[2] * perphat[1] * nablabhat[7] +
                               perphat[0] * perphat[2] * nablabhat[2] +
                               perphat[1] * perphat[2] * nablabhat[5] +
                               perphat[2] * perphat[2] * nablabhat[8]) -
                       0.25 * (rhohat[0] * rhohat[0] * nablabhat[0] +
                               rhohat[0] * rhohat[1] * nablabhat[3] +
                               rhohat[0] * rhohat[2] * nablabhat[6] +
                               rhohat[1] * rhohat[0] * nablabhat[1] +
                               rhohat[1] * rhohat[1] * nablabhat[4] +
                               rhohat[1] * rhohat[2] * nablabhat[7] +
                               rhohat[2] * rhohat[0] * nablabhat[2] +
                               rhohat[2] * rhohat[1] * nablabhat[5] +
                               rhohat[2] * rhohat[2] * nablabhat[8]);

    /* Choose a fixed basis so that e2 = bhat x e1. We choose e1 = bhat x z */
    real e1[3] = {0, 0, 1};
    real e2[3];
    math_cross(bhat, e1, e2);
    math_unit(e2, e1);
    math_cross(bhat, e1, e2);

    /* Gyrolength */
    real rho0 = 0.0;
    if (charge != 0)
    {
        rho0 = sqrt((2 * mass * mu0) / Bnorm) / fabs(charge);
    }

    /* First order position */
    real XYZ[3];
    XYZ[0] = xyz[0] - rho0 * rhohat[0];
    XYZ[1] = xyz[1] - rho0 * rhohat[1];
    XYZ[2] = xyz[2] - rho0 * rhohat[2];

    real RPZ[3];
    math_xyz2rpz(XYZ, RPZ);
    *R = RPZ[0];
    *Phi = RPZ[1];
    *Z = RPZ[2];

    /* Zeroth order gyroangle is directly defined from basis {e1, e2} and
       gyrovector as tan(zeta) = -rhohat dot e2 / rhohat dot e1 */
    real zeta0 = atan2(-math_dot(rhohat, e2), math_dot(rhohat, e1));

    /* First order momentum terms ppar, mu1, and zeta1 */
    real ppar1 = 0.0;
    if (charge != 0)
    {
        ppar1 = -ppar0 * rho0 * math_dot(rhohat, kappa) +
                (mass * mu0 / charge) * (tau + a1ddotgradb);
    }

    real ppar2 = ppar0 * ppar0; // Power, not second order component ;)
    real temp[3] = {
        mu0 * gradB[0] + ppar2 * kappa[0] / mass,
        mu0 * gradB[1] + ppar2 * kappa[1] / mass,
        mu0 * gradB[2] + ppar2 * kappa[2] / mass};
    real mu1 = 0.0;
    if (charge != 0)
    {
        mu1 = (rho0 / Bnorm) * math_dot(rhohat, temp) -
              (ppar0 * mu0 / (charge * Bnorm)) * (tau + a1ddotgradb);
    }

    /* This monster is Littlejohn's gyro gauge vector (nabla e1) dot e2 */
    real bx2by2 = bhat[0] * bhat[0] + bhat[1] * bhat[1];
    real b1x = -bhat[0] * bhat[2] / bx2by2;
    real b2x = -bhat[1] / (bx2by2 * bx2by2);
    real b1y = -bhat[1] * bhat[2] / bx2by2;
    real b2y = -bhat[0] / (bx2by2 * bx2by2);
    real Rvec[3];
    Rvec[0] = b1x * (nablabhat[1] + b2x * (nablabhat[0] + nablabhat[1])) +
              b1y * (nablabhat[0] + b2y * (nablabhat[0] + nablabhat[1]));
    Rvec[1] = b1x * (nablabhat[4] + b2x * (nablabhat[3] + nablabhat[4])) +
              b1y * (nablabhat[3] + b2y * (nablabhat[3] + nablabhat[4]));
    Rvec[2] = b1x * (nablabhat[7] + b2x * (nablabhat[6] + nablabhat[7])) +
              b1y * (nablabhat[6] + b2y * (nablabhat[6] + nablabhat[7]));

    temp[0] = gradB[0] + kappa[0] * (ppar2 / mass) / (2 * mu0);
    temp[1] = gradB[1] + kappa[1] * (ppar2 / mass) / (2 * mu0);
    temp[2] = gradB[2] + kappa[2] * (ppar2 / mass) / (2 * mu0);
    real zeta1 = 0.0;
    if (charge != 0)
    {
        zeta1 = -rho0 * math_dot(rhohat, Rvec) +
                (ppar0 / (charge * Bnorm)) * a2ddotgradb +
                (rho0 / Bnorm) * math_dot(perphat, temp);
    }

    /* Choose whether to use first order transformation in momentum space */
    if (GCTRANSFORM_ORDER)
    {
        *ppar = ppar0 + ppar1;
        *mu = fabs(mu0 + mu1);
        *zeta = zeta0 + zeta1;
    }
    else
    {
        *ppar = ppar0;
        *mu = mu0;
        *zeta = zeta0;
    }

    /* zeta is defined to be in interval [0, 2pi] */
    *zeta = fmod(CONST_2PI + (*zeta), CONST_2PI);
}


void gctransform_guidingcenter2particle(
    real mass, real charge, real *b_db, real R, real Phi, real Z, real ppar,
    real mu, real zeta, real *r, real *phi, real *z, real *pparprt, real *muprt,
    real *zetaprt)
{

    /* |B| */
    real Bnorm =
        sqrt(b_db[0] * b_db[0] + b_db[4] * b_db[4] + b_db[8] * b_db[8]);

    /* Guiding center transformation is more easily done in cartesian
     * coordinates so we switch to using those */
    real RPZ[3] = {R, Phi, Z};
    real XYZ[3];
    real b_dbxyz[12];
    math_rpz2xyz(RPZ, XYZ);
    math_jac_rpz2xyz(b_db, b_dbxyz, R, Phi);

    /* bhat = Unit vector of B */
    real bhat[3] = {b_dbxyz[0] / Bnorm, b_dbxyz[4] / Bnorm, b_dbxyz[8] / Bnorm};

    /* Magnetic field norm gradient */
    real gradB[3];
    gradB[0] =
        (bhat[0] * b_dbxyz[1] + bhat[1] * b_dbxyz[5] + bhat[2] * b_dbxyz[9]);
    gradB[1] =
        (bhat[0] * b_dbxyz[2] + bhat[1] * b_dbxyz[6] + bhat[2] * b_dbxyz[10]);
    gradB[2] =
        (bhat[0] * b_dbxyz[3] + bhat[1] * b_dbxyz[7] + bhat[2] * b_dbxyz[11]);

    /* nabla x |B| */
    real curlB[3] = {
        b_dbxyz[10] - b_dbxyz[7], b_dbxyz[3] - b_dbxyz[9],
        b_dbxyz[5] - b_dbxyz[2]};

    /* Magnetic field torsion = bhat dot ( nabla X bhat ) which is equivalent to
       bhat dot ( nabla x |B| ) / |B|
       because nabla x |B| = |B| nabla x bhat + (nabla |B|) x bhat */
    real tau = math_dot(bhat, curlB) / Bnorm;

    /* Gradient of magnetic field unit vector */
    real nablabhat[9];
    nablabhat[0] = (b_dbxyz[1] - gradB[0] * bhat[0]) / Bnorm;
    nablabhat[1] = (b_dbxyz[2] - gradB[0] * bhat[1]) / Bnorm;
    nablabhat[2] = (b_dbxyz[3] - gradB[0] * bhat[2]) / Bnorm;
    nablabhat[3] = (b_dbxyz[5] - gradB[1] * bhat[0]) / Bnorm;
    nablabhat[4] = (b_dbxyz[6] - gradB[1] * bhat[1]) / Bnorm;
    nablabhat[5] = (b_dbxyz[7] - gradB[1] * bhat[2]) / Bnorm;
    nablabhat[6] = (b_dbxyz[9] - gradB[2] * bhat[0]) / Bnorm;
    nablabhat[7] = (b_dbxyz[10] - gradB[2] * bhat[1]) / Bnorm;
    nablabhat[8] = (b_dbxyz[11] - gradB[2] * bhat[2]) / Bnorm;

    /* Magnetic field curvature vector = bhat dot nablabhat.
       Note that nabla x bhat = tau bhat + bhat x kappa */
    real kappa[3];
    kappa[0] = (bhat[0] * b_dbxyz[1] + bhat[1] * b_dbxyz[5] +
                bhat[2] * b_dbxyz[9] - math_dot(bhat, gradB) * bhat[0]) /
               Bnorm;
    kappa[1] = (bhat[0] * b_dbxyz[2] + bhat[1] * b_dbxyz[6] +
                bhat[2] * b_dbxyz[10] - math_dot(bhat, gradB) * bhat[1]) /
               Bnorm;
    kappa[2] = (bhat[0] * b_dbxyz[3] + bhat[1] * b_dbxyz[7] +
                bhat[2] * b_dbxyz[11] - math_dot(bhat, gradB) * bhat[2]) /
               Bnorm;

    /* Choose a fixed basis so that e2 = bhat x e1. We choose e1 = bhat x z */
    real e1[3] = {0, 0, 1};
    real e2[3];
    math_cross(bhat, e1, e2);
    math_unit(e2, e1);
    math_cross(bhat, e1, e2);

    /* Gyrolength */
    real rho0 = 0.0;
    if (charge != 0)
    {
        rho0 = sqrt((2 * mass * mu) / Bnorm) / fabs(charge);
    }

    /* Gyrovector rhohat and pperphat*/
    real rhohat[3];
    real perphat[3];
    real c = cos(zeta);
    real s = sin(zeta);

    rhohat[0] = c * e1[0] - s * e2[0];
    rhohat[1] = c * e1[1] - s * e2[1];
    rhohat[2] = c * e1[2] - s * e2[2];

    perphat[0] = -s * e1[0] - c * e2[0];
    perphat[1] = -s * e1[1] - c * e2[1];
    perphat[2] = -s * e1[2] - c * e2[2];

    /* Double product of dyadic a1 = -(1/2) * (rhohat perphat + perphat rhohat)
       and gradb (gradient of magnetic field unit vector) */
    real a1ddotgradb =
        -0.5 * (2 * (rhohat[0] * perphat[0] * nablabhat[0] +
                     rhohat[1] * perphat[1] * nablabhat[4] +
                     rhohat[2] * perphat[2] * nablabhat[8]) +
                (rhohat[0] * perphat[1] + rhohat[1] * perphat[0]) *
                    (nablabhat[1] + nablabhat[3]) +
                (rhohat[0] * perphat[2] + rhohat[2] * perphat[0]) *
                    (nablabhat[2] + nablabhat[6]) +
                (rhohat[1] * perphat[2] + rhohat[2] * perphat[1]) *
                    (nablabhat[5] + nablabhat[7]));

    /* Double product of dyadic a2 = (1/4) * (perphat perphat - rhohat rhohat)
       and gradb (gradient of magnetic field unit vector) */
    real a2ddotgradb = 0.25 * (perphat[0] * perphat[0] * nablabhat[0] +
                               perphat[1] * perphat[0] * nablabhat[3] +
                               perphat[2] * perphat[0] * nablabhat[6] +
                               perphat[0] * perphat[1] * nablabhat[1] +
                               perphat[1] * perphat[1] * nablabhat[4] +
                               perphat[2] * perphat[1] * nablabhat[7] +
                               perphat[0] * perphat[2] * nablabhat[2] +
                               perphat[1] * perphat[2] * nablabhat[5] +
                               perphat[2] * perphat[2] * nablabhat[8]) -
                       0.25 * (rhohat[0] * rhohat[0] * nablabhat[0] +
                               rhohat[0] * rhohat[1] * nablabhat[3] +
                               rhohat[0] * rhohat[2] * nablabhat[6] +
                               rhohat[1] * rhohat[0] * nablabhat[1] +
                               rhohat[1] * rhohat[1] * nablabhat[4] +
                               rhohat[1] * rhohat[2] * nablabhat[7] +
                               rhohat[2] * rhohat[0] * nablabhat[2] +
                               rhohat[2] * rhohat[1] * nablabhat[5] +
                               rhohat[2] * rhohat[2] * nablabhat[8]);

    /* First order terms */
    real ppar1 = 0.0;
    if (charge != 0)
    {
        ppar1 = -ppar * rho0 * math_dot(rhohat, kappa) +
                (mass * mu / charge) * (tau + a1ddotgradb);
    }

    real ppar2 = ppar * ppar; // Second power, not second order term
    real temp[3] = {
        mu * gradB[0] + ppar2 * kappa[0] / mass,
        mu * gradB[1] + ppar2 * kappa[1] / mass,
        mu * gradB[2] + ppar2 * kappa[2] / mass};
    real mu1 = 0.0;
    if (charge != 0)
    {
        mu1 = (rho0 / Bnorm) * math_dot(rhohat, temp) -
              (ppar * mu / (charge * Bnorm)) * (tau + a1ddotgradb);
    }

    temp[0] = gradB[0] + kappa[0] * ppar2 / (2 * mass * mu);
    temp[1] = gradB[1] + kappa[1] * ppar2 / (2 * mass * mu);
    temp[2] = gradB[2] + kappa[2] * ppar2 / (2 * mass * mu);

    /* This monster is Littlejohn's gyro gauge vector (nabla e1) dot e2 */
    real bx2by2 = bhat[0] * bhat[0] + bhat[1] * bhat[1];
    real b1x = -bhat[0] * bhat[2] / bx2by2;
    real b2x = -bhat[1] / (bx2by2 * bx2by2);
    real b1y = -bhat[1] * bhat[2] / bx2by2;
    real b2y = -bhat[0] / (bx2by2 * bx2by2);
    real Rvec[3];
    Rvec[0] = b1x * (nablabhat[1] + b2x * (nablabhat[0] + nablabhat[1])) +
              b1y * (nablabhat[0] + b2y * (nablabhat[0] + nablabhat[1]));
    Rvec[1] = b1x * (nablabhat[4] + b2x * (nablabhat[3] + nablabhat[4])) +
              b1y * (nablabhat[3] + b2y * (nablabhat[3] + nablabhat[4]));
    Rvec[2] = b1x * (nablabhat[7] + b2x * (nablabhat[6] + nablabhat[7])) +
              b1y * (nablabhat[6] + b2y * (nablabhat[6] + nablabhat[7]));
    real zeta1 = 0.0;
    if (charge != 0)
    {
        zeta1 = -rho0 * math_dot(rhohat, Rvec) +
                (ppar / (charge * Bnorm)) * a2ddotgradb +
                (rho0 / Bnorm) * math_dot(perphat, temp);
    }

    /* Choose whether to use first or zeroth order momentum transform */
    if (GCTRANSFORM_ORDER)
    {
        mu = fabs(mu - mu1);
        ppar -= ppar1;
        zeta -= zeta1;
        mu = fabs(mu);

        /* Calculate new unit vector for position */
        c = cos(zeta);
        s = sin(zeta);

        rhohat[0] = c * e1[0] - s * e2[0];
        rhohat[1] = c * e1[1] - s * e2[1];
        rhohat[2] = c * e1[2] - s * e2[2];
    }

    /* First order position */
    real xyz[3];
    xyz[0] = XYZ[0] + rho0 * rhohat[0];
    xyz[1] = XYZ[1] + rho0 * rhohat[1];
    xyz[2] = XYZ[2] + rho0 * rhohat[2];

    real rpz[3];
    math_xyz2rpz(xyz, rpz);

    *r = rpz[0];
    *phi = rpz[1];
    *z = rpz[2];
    *pparprt = ppar;
    *muprt = mu;
    *zetaprt = zeta;
}


void gctransform_pparmuzeta2prpphipz(
    real mass, real charge, real *b_db, real phi, real ppar, real mu, real zeta,
    real *pr, real *pphi, real *pz)
{
    /* Find magnetic field norm and unit vector */
    real Brpz[3] = {b_db[0], b_db[4], b_db[8]};
    real Bxyz[3];
    math_vec_rpz2xyz(Brpz, Bxyz, phi);

    real bhat[3];
    math_unit(Bxyz, bhat);
    real Bnorm = math_norm(Bxyz);

    /* Find the basis vectors e1 and e2 */
    real e1[3] = {0, 0, 1};
    real e2[3];
    math_cross(bhat, e1, e2);
    math_unit(e2, e1);
    math_cross(bhat, e1, e2);

    /* Perpendicular basis vector */
    real c = cos(zeta);
    real s = sin(zeta);
    real perphat[3];
    perphat[0] = -s * e1[0] - c * e2[0];
    perphat[1] = -s * e1[1] - c * e2[1];
    perphat[2] = -s * e1[2] - c * e2[2];

    /* Perpendicular momentum, negative particles travel opposite to perphat */
    real pperp = sqrt(2.0 * mass * Bnorm * mu);
    if (charge < 0)
    {
        pperp = -pperp;
    }

    /* Evaluate the momentum vector from ppar and pperp */
    real pxyz[3];
    pxyz[0] = ppar * bhat[0] + pperp * perphat[0];
    pxyz[1] = ppar * bhat[1] + pperp * perphat[1];
    pxyz[2] = ppar * bhat[2] + pperp * perphat[2];

    /* Back to cylindrical coordinates */
    real prpz[3];
    math_vec_xyz2rpz(pxyz, prpz, phi);

    *pr = prpz[0];
    *pphi = prpz[1];
    *pz = prpz[2];
}
