/**
 * Implements bfield_analytical.h.
 */
#include "bfield_analytical.h"
#include "defines.h"
#include "consts.h"
#include "errors.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int BfieldAnalytical_init(
    BfieldAnalytical *bfield, int nripple, real bphi, real rmajor, real rminor,
    real psiscaling, real ripplescaling, real rippledamping, real axisrz[2],
    real psilimits[2], real coefficients[13])
{

    bfield->bphi = bphi;
    bfield->rmajor = rmajor;
    bfield->rminor = rminor;
    bfield->nripple = nripple;
    bfield->axisrz[0] = axisrz[0];
    bfield->axisrz[1] = axisrz[1];
    bfield->psiscaling = psiscaling;
    bfield->psilimits[0] = psilimits[0];
    bfield->psilimits[1] = psilimits[1];
    bfield->ripplescaling = ripplescaling;
    bfield->rippledamping = rippledamping;

    for (int i = 0; i < 13; i++)
    {
        bfield->coefficients[i] = coefficients[i];
    }

    return 0;
}

void BfieldAnalytical_free(BfieldAnalytical *bfield)
{
    (void)bfield;
    // No resources were dynamically allocated
}

void BfieldAnalytical_offload(BfieldAnalytical *bfield)
{
    (void)bfield;
    GPU_MAP_TO_DEVICE(bfield->coefficients [0:14])
}

a5err BfieldAnalytical_eval_psi(
    real psi[1], real r, real z, BfieldAnalytical *bfield)
{
    z -= bfield->axisrz[1];
    r /= bfield->rmajor;
    z /= bfield->rmajor;

    real logr = log(r);
    real r2 = r * r, r3 = r2 * r, r4 = r3 * r, r5 = r4 * r, r6 = r5 * r;
    real z2 = z * z, z3 = z2 * z, z4 = z3 * z, z5 = z4 * z, z6 = z5 * z;

    real *C = bfield->coefficients;
    psi[0] =
        bfield->psiscaling *
        ((1 - C[12]) * (r4 / 8) + C[12] * (r2 * logr / 2) + C[0] * (1) +
         C[1] * (r2) + C[2] * (r2 * logr - z2) + C[3] * (r4 - 4 * r2 * z2) +
         C[4] * (3 * r4 * logr - 9 * r2 * z2 - 12 * r2 * logr * z2 + 2 * z4) +
         C[5] * (r6 - 12 * r4 * z2 + 8 * r2 * z4) +
         C[6] * (8 * z6 - 140 * r2 * z4 - 120 * r2 * logr * z4 +
                 180 * r4 * logr * z2 + 75 * r4 * z2 - 15 * r6 * logr) +
         C[7] * (z) + C[8] * (z * r2) + C[9] * (z3 - 3 * z * r2 * logr) +
         C[10] * (3 * z * r4 - 4 * z3 * r2) +
         C[11] *
             (8 * z5 - 45 * z * r4 - 80 * z3 * r2 * logr + 60 * z * r4 * logr));

    return 0;
}

a5err BfieldAnalytical_eval_psi_dpsi(
    real psi_dpsi[4], real r, real z, BfieldAnalytical *bfield)
{

    z -= bfield->axisrz[1];
    r /= bfield->rmajor;
    z /= bfield->rmajor;

    real logr = log(r);
    real r2 = r * r, r3 = r2 * r, r4 = r3 * r, r5 = r4 * r, r6 = r5 * r;
    real z2 = z * z, z3 = z2 * z, z4 = z3 * z, z5 = z4 * z, z6 = z5 * z;

    real *C = bfield->coefficients;
    psi_dpsi[0] =
        bfield->psiscaling *
        ((1 - C[12]) * (r4 / 8) + C[12] * (r2 * logr / 2) + C[0] * (1) +
         C[1] * (r2) + C[2] * (r2 * logr - z2) + C[3] * (r4 - 4 * r2 * z2) +
         C[4] * (3 * r4 * logr - 9 * r2 * z2 - 12 * r2 * logr * z2 + 2 * z4) +
         C[5] * (r6 - 12 * r4 * z2 + 8 * r2 * z4) +
         C[6] * (8 * z6 - 140 * r2 * z4 - 120 * r2 * logr * z4 +
                 180 * r4 * logr * z2 + 75 * r4 * z2 - 15 * r6 * logr) +
         C[7] * (z) + C[8] * (z * r2) + C[9] * (z3 - 3 * z * r2 * logr) +
         C[10] * (3 * z * r4 - 4 * z3 * r2) +
         C[11] *
             (8 * z5 - 45 * z * r4 - 80 * z3 * r2 * logr + 60 * z * r4 * logr));

    psi_dpsi[1] =
        (bfield->psiscaling / (r * bfield->rmajor)) *
        ((1 - C[12]) * (r3 / 2) + C[12] * (r / 2 + r * logr) + C[1] * (2 * r) +
         C[2] * (2 * r * logr + r) + C[3] * (4 * r3 - 8 * r * z2) +
         C[4] * (12 * r3 * logr + 3 * r3 - 30 * r * z2 - 24 * r * logr * z2) +
         C[5] * (6 * r5 - 48 * r3 * z2 + 16 * r * z4) +
         C[6] * (-400 * r * z4 - 240 * r * logr * z4 + 720 * r3 * logr * z2 +
                 480 * r3 * z2 - 90 * r5 * logr - 15 * r5) +
         C[8] * (2 * z * r) + C[9] * (-6 * z * r * logr - 3 * z * r) +
         C[10] * (12 * z * r3 - 8 * z3 * r) +
         C[11] * (-120 * z * r3 - 160 * z3 * r * logr - 80 * z3 * r +
                  240 * z * r3 * logr));

    psi_dpsi[2] = 0;

    psi_dpsi[3] =
        (bfield->psiscaling / (r * bfield->rmajor)) *
        (C[2] * (-2 * z) + C[3] * (-8 * r2 * z) +
         C[4] * (-18 * r2 * z - 24 * r2 * logr * z + 8 * z3) +
         C[5] * (-24 * r4 * z + 32 * r2 * z3) +
         C[6] * (48 * z5 - 560 * r2 * z3 - 480 * r2 * logr * z3 +
                 360 * r4 * logr * z + 150 * r4 * z) +
         C[7] * (1) + C[8] * (r2) + C[9] * (3 * z2 - 3 * r2 * logr) +
         C[10] * (3 * r4 - 12 * z2 * r2) +
         C[11] * (40 * z4 - 45 * r4 - 240 * z2 * r2 * logr + 60 * r4 * logr));

    return 0;
}

a5err BfieldAnalytical_eval_rho_drho(
    real rho_drho[4], real r, real z, BfieldAnalytical *bfield)
{
    real psi_dpsi[4];

    BfieldAnalytical_eval_psi_dpsi(psi_dpsi, r, z, bfield);

    real delta = bfield->psilimits[1] - bfield->psilimits[0];
    if ((psi_dpsi[0] - bfield->psilimits[0]) / delta < 0)
    {
        return error_raise(ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_GS);
    }

    rho_drho[0] = sqrt((psi_dpsi[0] - bfield->psilimits[0]) / delta);
    rho_drho[1] = psi_dpsi[1] / (2 * delta * rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[3] / (2 * delta * rho_drho[0]);

    return 0;
}

a5err BfieldAnalytical_eval_b(
    real b[3], real r, real phi, real z, BfieldAnalytical *bfield)
{
    z -= bfield->axisrz[1];
    r /= bfield->rmajor;
    z /= bfield->rmajor;

    real *C = bfield->coefficients;

    real logr = log(r);
    real r2 = r * r, r3 = r2 * r, r4 = r3 * r, r5 = r4 * r;
    real z2 = z * z, z3 = z2 * z, z4 = z3 * z, z5 = z4 * z;

    b[0] = C[2] * (-2 * z) + C[3] * (-8 * r2 * z) +
           C[4] * (-18 * r2 * z - 24 * r2 * logr * z + 8 * z3) +
           C[5] * (-24 * r4 * z + 32 * r2 * z3) +
           C[6] * (48 * z5 - 560 * r2 * z3 - 480 * r2 * logr * z3 +
                   360 * r4 * logr * z + 150 * r4 * z) +
           C[7] * (1) + C[8] * (r2) + C[9] * (3 * z2 - 3 * r2 * logr) +
           C[10] * (3 * r4 - 12 * z2 * r2) +
           C[11] * (40 * z4 - 45 * r4 - 240 * z2 * r2 * logr + 60 * r4 * logr);
    b[0] *= -bfield->psiscaling / (r * bfield->rmajor * bfield->rmajor);

    b[1] = bfield->bphi / r;

    b[2] = (1 - C[12]) * (r3 / 2) + C[12] * (r / 2 + r * logr) +
           C[1] * (2 * r) + C[2] * (2 * r * logr + r) +
           C[3] * (4 * r3 - 8 * r * z2) +
           C[4] * (12 * r3 * logr + 3 * r3 - 30 * r * z2 - 24 * r * logr * z2) +
           C[5] * (6 * r5 - 48 * r3 * z2 + 16 * r * z4) +
           C[6] * (-400 * r * z4 - 240 * r * logr * z4 + 720 * r3 * logr * z2 +
                   480 * r3 * z2 - 90 * r5 * logr - 15 * r5) +
           C[8] * (2 * z * r) + C[9] * (-6 * z * r * logr - 3 * z * r) +
           C[10] * (12 * z * r3 - 8 * z3 * r) +
           C[11] * (-120 * z * r3 - 160 * z3 * r * logr - 80 * z3 * r +
                    240 * z * r3 * logr);
    b[2] *= bfield->psiscaling / (r * bfield->rmajor * bfield->rmajor);

    if (bfield->nripple > 0)
    {
        r *= bfield->rmajor;
        z *= bfield->rmajor;
        real radius = sqrt(
            (r - bfield->rmajor) * (r - bfield->rmajor) +
            (z - bfield->axisrz[1]) * (z - bfield->axisrz[1]));
        real theta = atan2(z - bfield->axisrz[1], r - bfield->rmajor);
        real delta = bfield->ripplescaling * exp(-0.5 * theta * theta) *
                     pow(radius / bfield->rminor, bfield->rippledamping);
        b[1] = b[1] * (1 + delta * cos(bfield->nripple * phi));
    }

    return 0;
}

a5err BfieldAnalytical_eval_b_db(
    real b_db[12], real r, real phi, real z, BfieldAnalytical *bfield)
{

    real *C = bfield->coefficients;

    z -= bfield->axisrz[1];
    real R0 = bfield->rmajor;
    real z0 = bfield->axisrz[1];
    real B_phi0 = bfield->bphi;
    real psi_mult = bfield->psiscaling;

    /* Precalculated terms used in the components */
    real logr = log(r);
    real r2 = r * r, r3 = r2 * r, r4 = r3 * r, r5 = r4 * r;
    real z2 = z * z, z3 = z2 * z, z4 = z3 * z, z5 = z4 * z;
    real R02 = R0 * R0;
    real R03 = R02 * R0;
    real R04 = R03 * R0;
    real R05 = R04 * R0;
    real R06 = R05 * R0;

    real B0 = C[2] * (-2 * z) / R02 + C[3] * (-8 * r2 * z) / R04 +
              C[4] * (-18 * r2 * z - 24 * r2 * logr * z + 8 * z3) / R04 +
              C[5] * (-24 * r4 * z + 32 * r2 * z3) / R06 +
              C[6] *
                  (48 * z5 - 560 * r2 * z3 - 480 * r2 * logr * z3 +
                   360 * r4 * logr * z + 150 * r4 * z) /
                  R06 +
              C[7] * (1) / R0 + C[8] * (r2) / R03 +
              C[9] * (3 * z2 - 3 * r2 * logr) / R03 +
              C[10] * (3 * r4 - 12 * z2 * r2) / R05 +
              C[11] *
                  (40 * z4 - 45 * r4 - 240 * z2 * r2 * logr + 60 * r4 * logr) /
                  R05;
    B0 *= -psi_mult / r;
    real B1 =
        C[3] * (-16 * r * z) / R04 +
        C[4] * (-60 * r * z - 48 * r * logr * z) / R04 +
        C[5] * (-96 * r3 * z + 64 * r * z3) / R06 +
        C[6] *
            (-1600 * r * z3 - 960 * r * logr * z3 + 1440 * r3 * logr * z +
             960 * r3 * z) /
            R06 +
        C[8] * (2 * r) / R03 + C[9] * (-6 * r * logr - 3 * r) / R03 +
        C[10] * (12 * r3 - 24 * z2 * r) / R05 +
        C[11] *
            (-120 * r3 - 480 * z2 * r * logr - 240 * z2 * r + 240 * r3 * logr) /
            R05;
    B1 = -B0 / r - B1 * psi_mult / r;
    real B2 = 0;
    real B3 = C[2] * (-2) / R02 + C[3] * (-8 * r2) / R04 +
              C[4] * (-18 * r2 - 24 * r2 * logr + 24 * z2) / R04 +
              C[5] * (-24 * r4 + 96 * r2 * z2) / R06 +
              C[6] *
                  (240 * z4 - 1680 * r2 * z2 - 1440 * r2 * logr * z2 +
                   360 * r4 * logr + 150 * r4) /
                  R06 +
              C[9] * (6 * z) / R03 + C[10] * (-24 * z * r2) / R05 +
              C[11] * (160 * z3 - 480 * z * r2 * logr) / R05;
    B3 *= -psi_mult / r;

    real B4 = B_phi0 * R0 / r;
    real B5 = -B_phi0 * R0 / r2;
    real B6 = 0;
    real B7 = 0;

    real B8 = (1 - C[12]) * (r3 / 2) / R04 + C[12] * (r / 2 + r * logr) / R02 +
              C[1] * (2 * r) / R02 + C[2] * (2 * r * logr + r) / R02 +
              C[3] * (4 * r3 - 8 * r * z2) / R04 +
              C[4] *
                  (12 * r3 * logr + 3 * r3 - 30 * r * z2 - 24 * r * logr * z2) /
                  R04 +
              C[5] * (6 * r5 - 48 * r3 * z2 + 16 * r * z4) / R06 +
              C[6] *
                  (-400 * r * z4 - 240 * r * logr * z4 + 720 * r3 * logr * z2 +
                   480 * r3 * z2 - 90 * r5 * logr - 15 * r5) /
                  R06 +
              C[8] * (2 * z * r) / R03 +
              C[9] * (-6 * z * r * logr - 3 * z * r) / R03 +
              C[10] * (12 * z * r3 - 8 * z3 * r) / R05 +
              C[11] *
                  (-120 * z * r3 - 160 * z3 * r * logr - 80 * z3 * r +
                   240 * z * r3 * logr) /
                  R05;
    B8 *= psi_mult / r;
    real B9 =
        (1 - C[12]) * (3 * r2 / 2) / R04 + C[12] * (1.5 + logr) / R02 +
        C[1] * (2) / R02 + C[2] * (2 * logr + 3) / R02 +
        C[3] * (12 * r2 - 8 * z2) / R04 +
        C[4] * (36 * r2 * logr + 21 * r2 - 54 * z2 - 24 * logr * z2) / R04 +
        C[5] * (30 * r4 - 144 * r2 * z2 + 16 * z4) / R06 +
        C[6] *
            (-640 * z4 - 240 * logr * z4 + 2160 * r2 * logr * z2 +
             2160 * r2 * z2 - 450 * r4 * logr - 165 * r4) /
            R06 +
        C[8] * (2 * z) / R03 + C[9] * (-6 * z * logr - 9 * z) / R03 +
        C[10] * (36 * z * r2 - 8 * z3) / R05 +
        C[11] *
            (-120 * z * r2 - 160 * z3 * logr - 240 * z3 + 720 * z * r2 * logr) /
            R05;
    B9 = B9 * psi_mult / r - B8 / r;
    real B10 = 0;
    real B11 = -B1 - B0 / r;

    /* Ripple */
    if (bfield->nripple > 0)
    {
        real radius = sqrt((r - R0) * (r - R0) + (z - z0) * (z - z0));
        real theta = atan2(z - z0, r - R0);
        real delta = bfield->ripplescaling * exp(-0.5 * theta * theta) *
                     pow(radius / bfield->rminor, bfield->rippledamping);

        real Bphi = B4;
        real Bpert = Bphi * delta * cos(bfield->nripple * phi);
        B4 += Bpert;
        B6 += -Bphi * delta * bfield->nripple * sin(bfield->nripple * phi);

        real dBpertdR =
            Bpert * (((r - R0) / radius) * (bfield->rippledamping / radius) +
                     ((z - z0) / (radius * radius)) * theta);

        real dBpertdz =
            Bpert * (((z - z0) / radius) * (bfield->rippledamping / radius) -
                     ((r - R0) / (radius * radius)) * theta);

        B5 += B5 * Bpert / Bphi + dBpertdR;
        B7 += dBpertdz;
    }

    b_db[0] = B0;
    b_db[1] = B1;
    b_db[2] = B2;
    b_db[3] = B3;
    b_db[4] = B4;
    b_db[5] = B5;
    b_db[6] = B6;
    b_db[7] = B7;
    b_db[8] = B8;
    b_db[9] = B9;
    b_db[10] = B10;
    b_db[11] = B11;

    return 0;
}

a5err BfieldAnalytical_eval_axisrz(real axisrz[2], BfieldAnalytical *bfield)
{
    axisrz[0] = bfield->axisrz[0];
    axisrz[1] = bfield->axisrz[1];

    return 0;
}
