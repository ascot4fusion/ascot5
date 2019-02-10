/**
 * @file mccc.c
 * @brief Interface for using mccc package within ascot5
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../../ascot5.h"
#include "../../error.h"
#include "../../particle.h"
#include "../../B_field.h"
#include "../../plasma.h"
#include "../../random.h"
#include "../../math.h"
#include "../../consts.h"
#include "../../physlib.h"
#include "mccc.h"
#include "mccc_wiener.h"
#include "mccc_push.h"
#include "mccc_coefs.h"

#pragma omp declare target
/** Let collisions change energy */
static int MCCC_INCLUDE_ENERGY = 1;

/** Let collisions change pitch */
static int MCCC_INCLUDE_PITCH = 1;

/** Let collisions change guiding center position */
static int MCCC_INCLUDE_GCDIFF = 1;
#pragma omp end declare target

/**
 * @brief Set which quantities are affected by the collisions.
 *
 * @param include_energy can collisions change marker energy, either 0 or 1
 * @param include_pitch  can collisions change marker pitch, either 0 or 1
 * @param include_gcdiff can collisions change GC position, either 0 or 1
 */
void mccc_setoperator(int include_energy, int include_pitch,
                      int include_gcdiff) {
    MCCC_INCLUDE_ENERGY = include_energy;
    MCCC_INCLUDE_PITCH  = include_pitch;
    MCCC_INCLUDE_GCDIFF = include_gcdiff;
}

/**
 * @brief Evaluate collision coefficients
 *
 * This function is not called during the simulation but is used as a way
 * to get easy access to collision coefficients. Coefficients are evaluated
 * for plasma parameters found on given coordinates for an array of given
 * velocities.
 *
 * Evaluated coefficients are stored in the given arrays as:
 * [nv*i_species + i_v].
 *
 * @param m particle mass [kg]
 * @param q particle charge [C]
 * @param r R-coordinate where plasma is evaluated [m]
 * @param phi phi-coordinate where plasma is evaluated [rad]
 * @param z z-coordinate where plasma is evaluated [m]
 * @param t time-coordinate where plasma is evaluated [s]
 * @param v array of velocities [m/s]
 * @param nv number of velocities
 * @param pdata pointer to plasma data
 * @param F
 * @param Dpara
 * @param Dperp
 * @param K
 * @param nu
 * @param err error flags (nv length array) which are zero if evaluation
 *        succeeded for that velocity
 */
void mccc_eval_coefs(real m, real q, real r, real phi, real z, real t, real v,
                     int nv, plasma_data* pdata, real* F, real* Dpara,
                     real* Dperp, real* K, real* nu, int* err) {

}
