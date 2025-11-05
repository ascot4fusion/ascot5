/**
 * @file simulate.h
 * Contains declarations of simulation_offload_data and simulation_data structs.
 * Also simulation mode enums are declared here.
 */
#ifndef SIMULATE_H
#define SIMULATE_H

#include "ascot.h"
#include "data/atomic.h"
#include "data/bfield.h"
#include "data/boozer.h"
#include "defines.h"
#include "data/diag.h"
#include "data/efield.h"
#include "data/mhd.h"
#include "data/nbi.h"
#include "data/neutral.h"
#include "options.h"
#include "data/plasma.h"
#include "utils/random.h"
#include "data/rfof.h"
#include "data/wall.h"
#include "datatypes.h"

/**
 * Simulate gyro orbits and neutrals using fixed time-step.
 *
 * The simulation includes:
 * - orbit-following with Volume-Preserving Algorithm
 * - Coulomb collisions with Euler-Maruyama method
 *
 * The simulation is carried until all markers have met some
 * end condition or are aborted/rejected. The final state of the
 * markers is stored in the given marker array. Other output
 * is stored in the diagnostic array.
 *
 * The time-step is user-defined: either a directly given fixed value
 * or a given fraction of gyrotime.
 *
 * @param queue Simulated markers in a queue.
 *        Once a marker has finished simulation, it is returned to the queue
 *        with the state it had at the end of the simulation.
 * @param sim Simulation data.
 * @param vector_size Size of the marker vector array, i.e. how many markers are
 *        processed in parallel with SIMD instructions.
 */
int simulate_go_fixed(Simulation *sim, MarkerQueue *queue, size_t vector_size);

/**
 * Trace guiding centers using fixed time-step.
 *
 * The simulation loop proceeds as follows:
 *
 * 1. The current state of the marker is stored.
 * 2. Marker motion due to EM-fields, collisions, and ICRH kick is solved.
 * 4. End conditions are checked.
 * 5. Diagnostics are updated if end conditions are not met.
 * 5. Return to the beginning of the loop.
 *
 * @param queue Simulated markers in a queue.
 *        Once a marker has finished simulation, it is returned to the queue
 *        with the state it had at the end of the simulation.
 * @param sim Simulation data.
 * @param vector_size Size of the marker vector array, i.e. how many markers are
 *        processed in parallel with SIMD instructions.
 */
int simulate_gc_fixed(Simulation *sim, MarkerQueue *queue, size_t vector_size);

/**
 * Simulate guiding centers using adaptive time-step.
 *
 * The simulation loop proceeds as follows:
 *
 * 1. The current state of the marker is stored.
 * 2. Marker motion due to EM-fields, collisions, and ICRH kick is solved.
 * 3. Time step is either accepted or rejected.
 * 4. End conditions are checked.
 * 5. Diagnostics are updated if the time step was accepted and end conditions
 *    are not met.
 * 5. Return to the beginning of the loop.
 *
 * @param queue Simulated markers in a queue.
 *        Once a marker has finished simulation, it is returned to the queue
 *        with the state it had at the end of the simulation.
 * @param sim Simulation data.
 * @param vector_size Size of the marker vector array, i.e. how many markers are
 *        processed in parallel with SIMD instructions.
 */
int simulate_gc_adaptive(
    Simulation *sim, MarkerQueue *queue, size_t vector_size);

/**
 * Trace magnetic field lines using adaptive time-step.
 *
 * The simulation loop proceeds as follows:
 *
 * 1. The current state of the marker is stored.
 * 2. Marker motion along the field line is integrated for a single time-step.
 * 3. Time step is either accepted or rejected based on the error tolerance.
 * 4. End conditions are checked.
 * 5. Diagnostics are updated if the time step was accepted and end conditions
 *    are not met.
 * 5. Return to the beginning of the loop.
 *
 * Note that even though the integration time-step has "time" in it, in reality
 * the step is measured in distance. Marker time does not advance in the
 * simulation as the background is assumed to be frozen.
 *
 * @param queue Simulated markers in a queue.
 *        Once a marker has finished simulation, it is returned to the queue
 *        with the state it had at the end of the simulation.
 * @param sim Simulation data.
 * @param vector_size Size of the marker vector array, i.e. how many markers are
 *        processed in parallel with SIMD instructions.
 */
int simulate_fl_adaptive(
    Simulation *sim, MarkerQueue *queue, size_t vector_size);

#endif
