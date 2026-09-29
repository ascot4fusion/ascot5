/**
 * @file brownianbridge.h
 * Module for storing and generating Wiener processes taking into account
 * possible existing processes.
 *
 * When sampling Wiener processes it is important to take into account the
 * existing Wiener processes when time steps are rejected as these bias the¨
 * distribution from where a new process is sampled. Assuming that there are
 * no realizations of the Wiener process at time t > t0, the Wiener process is
 * drawn from normal distribution N(W(t0), t - t0). However, if there is a
 * realization of the Wiener process at time t1 > t0, then new process is
 * sampled using the Brownian bridge; the Wiener process is drawn from normal
 * distribution with mean
 * E[W(t)] = W(t0) + (t - t0) * (W(t1) - W(t0)) / (t1 - t0),
 * and variance
 * Var[W(t)] = (t - t0) * (t1 - t) / (t1 - t0).
 *
 * For these reason any existing process must be kept stored until the
 * simulation time passes the time of the process.
 */
#ifndef BROWNIANBRIDGE_H
#define BROWNIANBRIDGE_H

#include "defines.h"
#include <stdlib.h>

/**
 * Struct for storing Wiener processes.
 */
typedef struct
{
    /** Number of marker slots. */
    size_t nmrk;

    /** Number of dimensions in the Wiener process. */
    size_t ndim;

    /** Number of time instances per marker. */
    size_t ntime;

    /**
     * Time instances for different Wiener processes.
     *
     * The earliest existing process is always at the first time slot. Remaining
     * slots are unordered. Free slots have value -1.
     *
     * This array has format i_marker * ntime + i_time.
     */
    real *time;

    /**
     * Wiener processes.
     *
     * This array has format i_marker * ntime * ndim + i_time * ndim + i_dim.
     */
    real *wiener;
} BrownianBridge;

/**
 * Initialize the brownian bridge data.
 *
 * @param bbridge The brownian bridge data.
 * @param marker_slots The number of marker slots to be stored.
 * @param time_instances The number of time instances per marker.
 * @param wiener_dimensions The number of dimensions in the Wiener process.
 * @return Non-zero value if there was not enough memory to allocate the data.
 */
int BrownianBridge_init(
    BrownianBridge *bbridge, size_t marker_slots, size_t time_instances,
    size_t wiener_dimensions);

/**
 * Free the Brownian bridge data.
 *
 * @param bbridge The Brownian bridge data.
 */
void BrownianBridge_free(BrownianBridge *bbridge);

/**
 * Offload the Brownian bridge data to the GPU.
 *
 * @param bbridge The Brownian bridge data.
 */
void BrownianBridge_offload(BrownianBridge *bbridge);

/**
 * Generate the initial Wiener process for a marker.
 *
 * This function should be called whenever a simulation for a new marker is
 * started.
 *
 * @param bbridge The Brownian bridge data.
 * @param imrk The index of the marker in the Brownian bridge data.
 * @param t The current simulation time.
 */
void BrownianBridge_generate0th(BrownianBridge *bbridge, size_t imrk, real t);

/**
 * Generate the 5-dimensional Wiener processes for a marker at new
 * time-instance.
 *
 * @param wiener The generated Wiener processes.
 * @param t The time instance for which the processes are generated.
 * @param imrk The index of the marker in the Brownian bridge data.
 * @param bbridge The Brownian bridge data (must have `ndim == 5`).
 * @return Non-zero value if the process could not be generated due to the
 *         Wiener array being full.
 */
int BrownianBridge_generate5(
    real wiener[5], real t, size_t imrk, BrownianBridge *bbridge);

/**
 * Clear the Wiener processes up to a given time point.
 *
 * Processes W(t') are redundant if t' <  t, where t is the current simulation
 * time. Note that W(t) should exist before W(t') are removed. W(t) is moved to
 * the beginning of the array.
 *
 * @param bbridge The Brownian bridge data.
 * @param t The current simulation time.
 * @param imrk The index of the marker in the Brownian bridge data.
 */
void BrownianBridge_clear(BrownianBridge *bbridge, real t, size_t imrk);

#endif
