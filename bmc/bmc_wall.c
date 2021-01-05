#include "bmc_wall.h"

/**
 * Check if a given tile index is withing the target domain
 * For now, this is just a dummy functions that always return true
 * 
 * @param walltile walltile index
 * 
 * @todo read the target domain definition from the input HDF file
 * @todo make sure that the wall2d ascot functions return the walltile (or something similar)
 * @todo implement this check for 3d walls
 **/
int bmc_walltile_in_target(integer walltile) {
    return walltile == 141;
    // return 1;
}

/**
 * check if a particle trajectory hit the target domain
 * @param r0 coordinate r of the initial state
 * @param r1 coordinate r of the final state
 * @param phi0 coordinate phi of the initial state
 * @param phi1 coordinate phi of the final state
 * @param z0 coordinate z of the initial state
 * @param z1 coordinate z of the final state
 * @param w2d pointer to the 2d wall struct
 * 
 * @todo make sure that the wall2d ascot functions return the walltile (or something similar)
 **/
int bmc_wall_2d_hit_target(real r0, real r1, real phi0, real phi1, real z0, real z1, wall_2d_data* w2d) {
    int tile = wall_2d_hit_wall(r0, phi0, z0, r1, phi1, z1, w2d);
    if (tile > 0) {
        return bmc_walltile_in_target(tile);
    }
    return 0;
}