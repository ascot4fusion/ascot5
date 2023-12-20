#include "../particle.h"
#include "../wall.h"
#include "../endcond.h"

int bmc_walltile_in_target(integer walltile);
int bmc_wall_hit_target(real r0, real r1, real phi0, real phi1, real z0, real z1, wall_data* wdata);
void bmc_check_simd_particle_wallhit(
    particle_simd_gc *p,
    particle_simd_gc *p0,
    wall_data *wdata);