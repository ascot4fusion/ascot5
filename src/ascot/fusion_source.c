/**
 * Implements fusion_source.h.
 */
#include "fusion_source.h"
#include "defines.h"
#include "boschhale.h"
#include "consts.h"
#include "hist.h"
#include "mathlib.h"
#include "parallel.h"
#include "physlib.h"
#include "random.h"
#include "simulate.h"
#include <math.h>
#include <string.h>

/** Random number generator used by AFSI */
random_data rdata;

void afsi_compute_product_momenta_2d(
    size_t i, real m1, real m2, real mprod1, real mprod2, real Q,
    int prodmomspace, real *ppara1, real *pperp1, real *ppara2, real *pperp2,
    real *vcom2, real *prod1_p1, real *prod1_p2, real *prod2_p1, real *prod2_p2)
{

    real rn1 = CONST_2PI * random_uniform(rdata);
    real rn2 = CONST_2PI * random_uniform(rdata);

    real v1x = cos(rn1) * pperp1[i] / m1;
    real v1y = sin(rn1) * pperp1[i] / m1;
    real v1z = ppara1[i] / m1;

    real v2x = cos(rn2) * pperp2[i] / m2;
    real v2y = sin(rn2) * pperp2[i] / m2;
    real v2z = ppara2[i] / m2;

    *vcom2 = (v1x - v2x) * (v1x - v2x) + (v1y - v2y) * (v1y - v2y) +
             (v1z - v2z) * (v1z - v2z);

    // Velocity of the system's center of mass
    real v_cm[3];
    v_cm[0] = (m1 * v1x + m2 * v2x) / (m1 + m2);
    v_cm[1] = (m1 * v1y + m2 * v2y) / (m1 + m2);
    v_cm[2] = (m1 * v1z + m2 * v2z) / (m1 + m2);

    // Total kinetic energy after the reaction in CM frame
    real ekin = Q +
                0.5 * m1 *
                    ((v1x - v_cm[0]) * (v1x - v_cm[0]) +
                     (v1y - v_cm[1]) * (v1y - v_cm[1]) +
                     (v1z - v_cm[2]) * (v1z - v_cm[2])) +
                0.5 * m2 *
                    ((v2x - v_cm[0]) * (v2x - v_cm[0]) +
                     (v2y - v_cm[1]) * (v2y - v_cm[1]) +
                     (v2z - v_cm[2]) * (v2z - v_cm[2]));

    // Speed and velocity of product 2 in CM frame
    rn1 = random_uniform(rdata);
    rn2 = random_uniform(rdata);
    real phi = CONST_2PI * rn1;
    real theta = acos(2 * (rn2 - 0.5));
    real vnorm = sqrt(2.0 * ekin / (mprod2 * (1.0 + mprod2 / mprod1)));

    real v2_cm[3];
    v2_cm[0] = vnorm * sin(theta) * cos(phi);
    v2_cm[1] = vnorm * sin(theta) * sin(phi);
    v2_cm[2] = vnorm * cos(theta);

    // Products' velocities in lab frame
    real vprod1[3], vprod2[3];
    vprod1[0] = -(mprod2 / mprod1) * v2_cm[0] + v_cm[0];
    vprod1[1] = -(mprod2 / mprod1) * v2_cm[1] + v_cm[1];
    vprod1[2] = -(mprod2 / mprod1) * v2_cm[2] + v_cm[2];
    vprod2[0] = v2_cm[0] + v_cm[0];
    vprod2[1] = v2_cm[1] + v_cm[1];
    vprod2[2] = v2_cm[2] + v_cm[2];

    if (prodmomspace == PPARPPERP)
    {
        prod1_p1[i] = vprod1[2] * mprod1;
        prod1_p2[i] =
            sqrt(vprod1[0] * vprod1[0] + vprod1[1] * vprod1[1]) * mprod1;
        prod2_p1[i] = vprod2[2] * mprod2;
        prod2_p2[i] =
            sqrt(vprod2[0] * vprod2[0] + vprod2[1] * vprod2[1]) * mprod2;
    }
    else
    {
        real vnorm1 = math_norm(vprod1);
        prod1_p2[i] = vprod1[2] / vnorm1;
        prod1_p1[i] = physlib_Ekin_gamma(mprod1, physlib_gamma_vnorm(vnorm1));

        real vnorm2 = math_norm(vprod2);
        prod2_p2[i] = vprod2[2] / vnorm2;
        prod2_p1[i] = physlib_Ekin_gamma(mprod2, physlib_gamma_vnorm(vnorm2));
    }
}

void afsi_sample_reactant_momenta_2d(
    sim_data *sim, afsi_data *afsi, real mass1, real mass2, real vol,
    size_t nsample, size_t i0, size_t i1, size_t i2, real r, real phi, real z,
    real time, real rho, real *density1, real *ppara1, real *pperp1,
    real *density2, real *ppara2, real *pperp2)
{
    if (afsi->type1 == 1)
    {
        afsi_sample_beam_2d(
            afsi->beam1, mass1, vol, nsample, i0, i1, i2, density1, ppara1,
            pperp1);
    }
    else if (afsi->type1 == 2)
    {
        afsi_sample_thermal_2d(
            sim, afsi->thermal1, mass1, nsample, r, phi, z, time, rho, density1,
            ppara1, pperp1);
    }
    if (afsi->type2 == 1)
    {
        afsi_sample_beam_2d(
            afsi->beam2, mass2, vol, nsample, i0, i1, i2, density2, ppara2,
            pperp2);
    }
    else if (afsi->type2 == 2)
    {
        afsi_sample_thermal_2d(
            sim, afsi->thermal2, mass2, nsample, r, phi, z, time, rho, density2,
            ppara2, pperp2);
    }
}

void afsi_sample_beam_2d(
    histogram *hist, real mass, real vol, size_t nsample, size_t i0, size_t i1,
    size_t i2, real *density, real *ppara, real *pperp)
{
    int mom_space;
    size_t p1coord, p2coord;
    if (hist->axes[5].n)
    {
        p1coord = 5;
        p2coord = 6;
        mom_space = PPARPPERP;
    }
    else if (hist->axes[10].n)
    {
        p1coord = 10;
        p2coord = 11;
        mom_space = EKINXI;
    }
    else
    {
        return;
    }

    real *cumdist = (real *)malloc(
        hist->axes[p1coord].n * hist->axes[p2coord].n * sizeof(real));

    *density = 0.0;
    for (size_t ip1 = 0; ip1 < hist->axes[p1coord].n; ip1++)
    {
        for (size_t ip2 = 0; ip2 < hist->axes[p2coord].n; ip2++)
        {
            size_t index = i0 * hist->strides[0] + i1 * hist->strides[1] +
                           i2 * hist->strides[2] +
                           ip1 * hist->strides[p1coord] +
                           ip2 * hist->strides[p2coord];
            if (ip1 == 0 && ip2 == 0)
            {
                cumdist[0] = hist->bins[index];
            }
            else
            {
                cumdist[ip1 * hist->axes[p2coord].n + ip2] =
                    cumdist[ip1 * hist->axes[p2coord].n + ip2 - 1] +
                    hist->bins[index];
            }
            *density += hist->bins[index] / vol;
        }
    }
    if (*density == 0)
    {
        return;
    }

    for (size_t i = 0; i < nsample; i++)
    {
        real r = random_uniform(rdata);
        r *= cumdist[hist->axes[p1coord].n * hist->axes[p2coord].n - 1];
        for (size_t j = 0; j < hist->axes[p1coord].n * hist->axes[p2coord].n;
             j++)
        {
            if (cumdist[j] > r)
            {
                if (mom_space == PPARPPERP)
                {
                    ppara[i] = hist->axes[5].min +
                               (j / hist->axes[6].n + 0.5) *
                                   (hist->axes[5].max - hist->axes[5].min) /
                                   hist->axes[5].n;
                    pperp[i] = hist->axes[6].min +
                               (j % hist->axes[6].n + 0.5) *
                                   (hist->axes[6].max - hist->axes[6].min) /
                                   hist->axes[6].n;
                }
                else
                {
                    real ekin = hist->axes[10].min +
                                (j / hist->axes[10].n + 0.5) *
                                    (hist->axes[10].max - hist->axes[10].min) /
                                    hist->axes[10].n;
                    real pitch = hist->axes[11].min +
                                 (j / hist->axes[11].n + 0.5) *
                                     (hist->axes[11].max - hist->axes[11].min) /
                                     hist->axes[11].n;
                    real gamma = physlib_gamma_Ekin(mass, ekin);
                    real pnorm = sqrt(gamma * gamma - 1.0) * mass * CONST_C;
                    ppara[i] = pitch * pnorm;
                    pperp[i] = sqrt(1.0 - pitch * pitch) * pnorm;
                }
                break;
            }
        }
    }
    free(cumdist);
}

void afsi_sample_thermal_2d(
    sim_data *sim, int ispecies, real mass, size_t nsample, real r, real phi,
    real z, real time, real rho, real *density, real *ppara, real *pperp)
{
    real ni, ti;
    if (plasma_eval_dens(
            &ni, rho, r, phi, z, time, ispecies, &sim->plasma) ||
        plasma_eval_temp(
            &ti, rho, r, phi, z, time, ispecies, &sim->plasma))
    {
        *density = 0.0;
        return;
    }
    *density = ni;
    for (size_t i = 0; i < nsample; i++)
    {
        real r1, r2, r3, r4, E;

        r1 = random_uniform(rdata);
        r2 = random_uniform(rdata);
        r3 = cos(0.5 * random_uniform(rdata) * CONST_PI);
        E = -ti * (log(r1) + log(r2) * r3 * r3);

        r4 = 1.0 - 2 * random_uniform(rdata);
        pperp[i] = sqrt((1 - r4 * r4) * 2 * E * mass);
        ppara[i] = r4 * sqrt(2 * E * mass);
    }
}
