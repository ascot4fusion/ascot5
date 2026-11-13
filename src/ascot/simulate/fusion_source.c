/**
 * Implements "fusion_source.h".
 */
#include "fusion_source.h"
#include "consts.h"
#include "data/diag_hist.h"
#include "data/marker.h"
#include "data/plasma.h"
#include "datatypes.h"
#include "defines.h"
#include "parallel.h"
#include "utils/ratecoefficients.h"
#include "utils/gctransform.h"
#include "utils/mathlib.h"
#include "utils/physlib.h"
#include "utils/random.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

/**
 * Helper function to sample coordinates from 2D histogram.
 */
void fusion_source_sample_2d(
    real out[2], real rnd, size_t spatial_index, HistCoordinate p1coord,
    HistCoordinate p2coord, DiagHist *hist)
{
    for (size_t ip1 = 0; ip1 < hist->axes[p1coord].n; ip1++)
    {
        for (size_t ip2 = 0; ip2 < hist->axes[p2coord].n; ip2++)
        {
            size_t index = spatial_index + ip1 * hist->strides[p1coord] +
                           ip2 * hist->strides[p2coord];
            rnd -= hist->values[index];
            if (rnd < 0)
            {
                out[0] =
                    hist->axes[p1coord].min +
                    (0.5 + (double)ip1) *
                        (hist->axes[p1coord].max - hist->axes[p1coord].min) /
                        hist->axes[p1coord].n;
                out[1] =
                    hist->axes[p2coord].min +
                    (0.5 + (double)ip2) *
                        (hist->axes[p2coord].max - hist->axes[p2coord].min) /
                        hist->axes[p2coord].n;
                return;
            }
        }
    }
}

/**
 * Helper function to sample coordinates from 3D histogram.
 */
void fusion_source_sample_3d(
    real out[3], real rnd, size_t spatial_index, HistCoordinate p1coord,
    HistCoordinate p2coord, HistCoordinate p3coord, DiagHist *hist)
{
    for (size_t ip1 = 0; ip1 < hist->axes[p1coord].n; ip1++)
    {
        for (size_t ip2 = 0; ip2 < hist->axes[p2coord].n; ip2++)
        {
            for (size_t ip3 = 0; ip3 < hist->axes[p3coord].n; ip3++)
            {
                size_t index = spatial_index + ip1 * hist->strides[p1coord] +
                               ip2 * hist->strides[p2coord] +
                               ip3 * hist->strides[p3coord];
                rnd -= hist->values[index];
                if (rnd < 0)
                {
                    out[0] = hist->axes[p1coord].min +
                             (0.5 + (double)ip1) *
                                 (hist->axes[p1coord].max -
                                  hist->axes[p1coord].min) /
                                 hist->axes[p1coord].n;
                    out[1] = hist->axes[p2coord].min +
                             (0.5 + (double)ip2) *
                                 (hist->axes[p2coord].max -
                                  hist->axes[p2coord].min) /
                                 hist->axes[p2coord].n;
                    out[2] = hist->axes[p3coord].min +
                             (0.5 + (double)ip3) *
                                 (hist->axes[p3coord].max -
                                  hist->axes[p3coord].min) /
                                 hist->axes[p3coord].n;
                    return;
                }
            }
        }
    }
}

void fusion_source_compute_product(
    Reaction reaction, real vol, MarkerGyroOrbit *react1,
    MarkerGyroOrbit *react2, MarkerGyroOrbit *prod1, MarkerGyroOrbit *prod2)
{
    size_t nsample = react1->size;
    for (size_t i = 0; i < nsample; i++)
    {
        real v1[3] = {
            react1->p_r[i] / react1->mass,
            react1->p_phi[i] / react1->mass,
            react1->p_z[i] / react1->mass,
        };
        real v2[3] = {
            react2->p_r[i] / react2->mass,
            react2->p_phi[i] / react2->mass,
            react2->p_z[i] / react2->mass,
        };

        real vrel = sqrt(
            (v1[0] - v2[0]) * (v1[0] - v2[0]) +
            (v1[1] - v2[1]) * (v1[1] - v2[1]) +
            (v1[2] - v2[2]) * (v1[2] - v2[2]));

        real v_com[3] = {
            (react1->mass * v1[0] + react2->mass * v2[0]) /
                (react1->mass + react2->mass),
            (react1->mass * v1[1] + react2->mass * v2[1]) /
                (react1->mass + react2->mass),
            (react1->mass * v1[2] + react2->mass * v2[2]) /
                (react1->mass + react2->mass),
        };

        real ekin_before = 0.5 * react1->mass *
                               ((v1[0] - v_com[0]) * (v1[0] - v_com[0]) +
                                (v1[1] - v_com[1]) * (v1[1] - v_com[1]) +
                                (v1[2] - v_com[2]) * (v1[2] - v_com[2])) +
                           0.5 * react2->mass *
                               ((v2[0] - v_com[0]) * (v2[0] - v_com[0]) +
                                (v2[1] - v_com[1]) * (v2[1] - v_com[1]) +
                                (v2[2] - v_com[2]) * (v2[2] - v_com[2]));
        real weight =
            vrel * boschhale_sigma(reaction, ekin_before) / nsample * vol;
        real ekin_after = ekin_before + Q_value[reaction] * CONST_E;

        real rnd[2];
        random_uniform(rnd, 2, react1->id, react1->rng_count, i);
        real phi = CONST_2PI * rnd[0];
        real theta = acos(2 * (rnd[1] - 0.5));
        real vnorm = sqrt(
            2 * ekin_after / (prod2->mass * (1. + prod2->mass / prod1->mass)));

        real vprod2_com[3] = {
            vnorm * sin(theta) * cos(phi),
            vnorm * sin(theta) * sin(phi),
            vnorm * cos(theta),
        };

        real vprod1[3] = {
            -(prod2->mass / prod1->mass) * vprod2_com[0] + v_com[0],
            -(prod2->mass / prod1->mass) * vprod2_com[1] + v_com[1],
            -(prod2->mass / prod1->mass) * vprod2_com[2] + v_com[2],
        };
        real vprod2[3] = {
            vprod2_com[0] + v_com[0],
            vprod2_com[1] + v_com[1],
            vprod2_com[2] + v_com[2],
        };

        prod1->p_r[i] = vprod1[0] * prod1->mass;
        prod1->p_phi[i] = vprod1[1] * prod1->mass;
        prod1->p_z[i] = vprod1[2] * prod1->mass;
        prod1->weight[i] = weight;
        prod2->p_r[i] = vprod2[0] * prod2->mass;
        prod2->p_phi[i] = vprod2[1] * prod2->mass;
        prod2->p_z[i] = vprod2[2] * prod2->mass;
        prod2->weight[i] = weight;
    }
}

void fusion_source_sample_beam(
    DiagHist *hist, size_t spatial_index, real vol, MarkerGyroOrbit *sampled,
    real density[1])
{
    HistCoordinate p1coord = HIST_PR, p2coord = HIST_PPHI, p3coord = HIST_PZ;
    if (hist->axes[HIST_PPAR].n && hist->axes[HIST_PPERP].n)
    {
        p1coord = HIST_PPAR;
        p2coord = HIST_PPERP;
    }
    else if (hist->axes[HIST_XI].n && hist->axes[HIST_EKIN].n)
    {
        p1coord = HIST_XI;
        p2coord = HIST_EKIN;
    }
    int is3d = p1coord == HIST_PR;

    real cumsum = 0;
    for (size_t ip1 = 0; ip1 < hist->axes[p1coord].n; ip1++)
    {
        for (size_t ip2 = 0; ip2 < hist->axes[p2coord].n; ip2++)
        {
            if (is3d)
            {
                for (size_t ip3 = 0; ip3 < hist->axes[p3coord].n; ip3++)
                {
                    size_t index = spatial_index +
                                   ip1 * hist->strides[p1coord] +
                                   ip2 * hist->strides[p2coord] +
                                   ip3 * hist->strides[p3coord];
                    cumsum += hist->values[index];
                }
            }
            else
            {
                size_t index = spatial_index + ip1 * hist->strides[p1coord] +
                               ip2 * hist->strides[p2coord];
                cumsum += hist->values[index];
            }
        }
    }
    density[0] = cumsum / vol;
    if (density[0] == 0)
        return;

    for (size_t i = 0; i < sampled->size; i++)
    {
        real rnd[1];
        random_uniform(rnd, 1, sampled->id, sampled->rng_count, i);
        rnd[0] *= cumsum;
        if (is3d)
        {
            real out[3];
            fusion_source_sample_3d(
                out, rnd[0], spatial_index, p1coord, p2coord, p3coord, hist);
            sampled->p_r[i] = out[0];
            sampled->p_z[i] = out[1];
            sampled->p_phi[i] = out[2];
        }
        else
        {
            real out[2];
            fusion_source_sample_2d(
                out, rnd[0], spatial_index, p1coord, p2coord, hist);
        }
    }
}

void fusion_source_sample_thermal(
    Simulation *sim, size_t species_index, MarkerGyroOrbit *sampled,
    real density[1])
{
    real ni, Ti;
    if (Plasma_eval_density(
            &ni, sampled->rho[0], sampled->r[0], sampled->phi[0], sampled->z[0],
            sampled->time[0], species_index, &sim->plasma) ||
        Plasma_eval_temperature(
            &Ti, sampled->rho[0], sampled->r[0], sampled->phi[0], sampled->z[0],
            sampled->time[0], species_index, &sim->plasma))
    {
        density[0] = 0;
        return;
    }
    density[0] = ni;

    for (size_t i = 0; i < sampled->size; i++)
    {
        real rnd[3];
        real sigma = sqrt(Ti * sampled->mass);
        random_norm(rnd, 3, sampled->id, sampled->rng_count, i);
        sampled->p_r[i] = rnd[0] * sigma;
        sampled->p_z[i] = rnd[1] * sigma;
        sampled->p_phi[i] = rnd[2] * sigma;
    }
}
