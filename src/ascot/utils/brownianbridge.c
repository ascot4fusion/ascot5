/**
 * Implements "brownianbridge.h".
 */
#include "brownianbridge.h"
#include "defines.h"
#include "parallel.h"
#include <math.h>
#include <stdio.h>

int BrownianBridge_init(
    BrownianBridge *bbridge, size_t marker_slots, size_t time_instances,
    size_t wiener_dimensions)
{
    int err = 0;
    bbridge->time =
        (real *)xmalloc(&err, marker_slots * time_instances * sizeof(real));
    bbridge->wiener = (real *)xmalloc(
        &err, marker_slots * time_instances * wiener_dimensions * sizeof(real));

    if (err)
    {
        BrownianBridge_free(bbridge);
        return 1;
    }

    bbridge->nmrk = marker_slots;
    bbridge->ndim = wiener_dimensions;
    bbridge->ntime = time_instances;
    return 0;
}

void BrownianBridge_free(BrownianBridge *bbridge)
{
    free(bbridge->time);
    free(bbridge->wiener);
}

void BrownianBridge_offload(BrownianBridge *bbridge)
{
    SUPPRESS_UNUSED_WARNING(bbridge);
    GPU_MAP_TO_DEVICE(
        bbridge->time [0:bbridge->nmrk * bbridge->ntime],
        bbridge->wiener [0:bbridge->nmrk * bbridge->ntime * bbridge->ndim])
}

void BrownianBridge_generate0th(BrownianBridge *bbridge, size_t imrk, real t)
{
    size_t idx = imrk * bbridge->ntime;
    bbridge->time[idx] = t;
    for (size_t i = 0; i < bbridge->ndim; i++)
        bbridge->wiener[idx * bbridge->ndim + i] = 0.0;
    for (size_t i = 1; i < bbridge->ntime; i++)
        bbridge->time[idx + i] = -1.0;
}

int BrownianBridge_generate5(
    real wiener[5], real t, size_t imrk, BrownianBridge *bbridge)
{
    size_t idx = imrk * bbridge->ntime;
    real t0 = bbridge->time[idx], t1 = -1;
    real *W0 = &bbridge->wiener[idx * bbridge->ndim], *W1 = W0;

    // Find values for t0, t1, W0, and W1. Also find index for a free slot
    // where the new process can be inserted (zero if there is no free slot).
    // Check if there is an existing process at time t (ntime if not).
    size_t free_slot = 0, existing_process = bbridge->ntime;
    for (size_t i = 1; i < bbridge->ntime; i++)
    {
        real t_ = bbridge->time[idx + i];
        int update_t0 = t_ > t0 && t_ < t;
        int update_t1 = (t_ < t1 || t1 == -1) && t_ >= t;

        t0 = update_t0 ? t_ : t0;
        t1 = update_t1 ? t_ : t1;
        W0 = update_t0
                 ? &bbridge->wiener[idx * bbridge->ndim + i * bbridge->ndim]
                 : W0;
        W1 = update_t1
                 ? &bbridge->wiener[idx * bbridge->ndim + i * bbridge->ndim]
                 : W1;

        free_slot = t_ == -1 ? i : free_slot;
        existing_process = (t_ == t) ? i : existing_process;
    }

    // Set free slot to existing process if we found an existing process
    free_slot =
        existing_process < bbridge->ntime ? existing_process : free_slot;
    idx = imrk * bbridge->ntime * bbridge->ndim + free_slot * bbridge->ndim;
    int bridge = t1 != -1;
    real var = (t - t0) * (1 + bridge * ((t1 - t) / (t1 - t0) - 1));
    real temp = bridge * (t - t0) / (t1 - t0);
    for (size_t i = 0; i < 5; i++)
    {
        wiener[i] = W0[i] + temp * (W1[i] - W0[i]) + sqrt(var) * wiener[i];
        wiener[i] = existing_process < bbridge->ntime ? bbridge->wiener[idx + i]
                                                      : wiener[i];
        wiener[i] = free_slot ? wiener[i] : W0[i];
        bbridge->wiener[idx + i] = wiener[i];
    }
    bbridge->time[imrk * bbridge->ntime + free_slot] = free_slot ? t : t0;

    return free_slot == 0 && existing_process == bbridge->ntime;
}

void BrownianBridge_clear(BrownianBridge *bbridge, real t, size_t imrk)
{

    for (size_t i = 0; i < bbridge->ntime; i++)
    {
        real t_ = bbridge->time[imrk * bbridge->ntime + i];
        bbridge->time[imrk * bbridge->ntime + i] = t_ <= t ? -1 : t_;

        for (size_t j = 0; j < bbridge->ndim; j++)
        {
            size_t idx =
                imrk * bbridge->ntime * bbridge->ndim + i * bbridge->ndim + j;
            bbridge->wiener[0] =
                t_ == t ? bbridge->wiener[idx] : bbridge->wiener[0];
        }
    }
    bbridge->time[imrk * bbridge->ntime] = t;
}
