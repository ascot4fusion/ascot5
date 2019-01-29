/**
 * @file boozer.c
 * @brief Module for transforming between cylindrical and Boozer coordinates.
 */
#include "ascot5.h"
#include "error.h"
#include "boozer.h"

/**
 * @brief Load Boozer data and prepare parameters for offload.
 *
 * This function fills the boozer offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 *
 * @todo Konsta will write this.
 */
int boozer_init_offload(boozer_offload_data* offload_data,
                        real** offload_array) {
    return 0;
}

/**
 * @brief Initialize boozer data struct on target
 *
 * @param boozerdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 *
 * @todo Konsta will write this.
 */
void boozer_init(boozer_data* boozerdata, boozer_offload_data* offload_data,
                     real* offload_array) {

}

/**
 * @brief Transform cylindrical coordinates to Boozer coordinates.
 *
 * @todo This is just a dummy.
 *
 * @param ptz Boozer coordinates as [psi, theta, zeta].
 *
 * @return zero on success
 */
a5err boozer_cyl2booz(real ptz[3], real r, real phi, real z,
                      boozer_data* boozerdata) {
    return 0;
}

/**
 * @brief Transform Boozer coordinates to cylindrical coordinates.
 *
 * @todo This is just a dummy.
 *
 * @param rz cylindrical coordinates as [R, z].
 *
 * @return zero on success
 */
a5err boozer_booz2cyl(real rz[2], real psi, real theta, real zeta,
                      boozer_data* boozerdata) {
    return 0;
}

/**
 * @brief Evaluate Boozer coordinates and gradients on a given location.
 *
 * @todo This is just a dummy.
 *
 * The values are stored in the given array as:
 * - ptz_dptz[0]  = psi
 * - ptz_dptz[1]  = dpsi/dR
 * - ptz_dptz[2]  = dpsi/dphi
 * - ptz_dptz[3]  = dpsi/dz
 * - ptz_dptz[4]  = theta
 * - ptz_dptz[5]  = dtheta/dR
 * - ptz_dptz[6]  = dtheta/dphi
 * - ptz_dptz[7]  = dtheta/dz
 * - ptz_dptz[8]  = zeta
 * - ptz_dptz[9]  = dzeta/dR
 * - ptz_dptz[10] = dzeta/dphi
 * - ptz_dptz[11] = dzeta/dz
 *
 * @param ptz_dptz evaluated Boozer coordinates and their gradients.
 *
 * @return zero on success
 */
a5err boozer_eval_gradients(real ptz_dptz[12], real r, real phi, real z,
                            boozer_data* boozerdata) {
    return 0;
}
