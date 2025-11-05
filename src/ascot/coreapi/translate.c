#include "datatypes.h"
#include "defines.h"
#include <string.h>

const char *file_names[C_FILE_COUNT] = {
    "./coreapi/interpolate.c",
    "./coreapi/solve.c",
    "./coreapi/tools.c",
    "./data/atomic.c",
    "./data/bfield.c",
    "./data/bfield_analytical.c",
    "./data/bfield_cartesian.c",
    "./data/bfield_spline2d.c",
    "./data/bfield_spline3d.c",
    "./data/bfield_stellarator.c",
    "./data/boozer.c",
    "./data/diag.c",
    "./data/diag_orb.c",
    "./data/diag_transcoef.c",
    "./data/dist_5D.c",
    "./data/dist_6D.c",
    "./data/dist_com.c",
    "./data/dist_rho5D.c",
    "./data/dist_rho6D.c",
    "./data/efield.c",
    "./data/efield_cartesian.c",
    "./data/efield_potential1d.c",
    "./data/hist.c",
    "./data/marker.c",
    "./data/mhd.c",
    "./data/mhd_dynamic.c",
    "./data/mhd_stationary.c",
    "./data/nbi.c",
    "./data/neutral.c",
    "./data/neutral_arbitrary.c",
    "./data/neutral_radial.c",
    "./data/plasma.c",
    "./data/plasma_dynamic1d.c",
    "./data/plasma_linear1d.c",
    "./data/rfof.c",
    "./data/wall.c",
    "./data/wall_contour2d.c",
    "./data/wall_triangular3d.c",
    "./endcond.c",
    "./simulate/atomic_reactions.c",
    "./simulate/ccoll_gc_euler.c",
    "./simulate/ccoll_gc_milstein.c",
    "./simulate/ccoll_go_euler.c",
    "./simulate/ccoll_wiener.c",
    "./simulate/coulomb_collisions.c",
    "./simulate/field_line_adaptive.c",
    "./simulate/fusion_source.c",
    "./simulate/guiding_center_adaptive.c",
    "./simulate/guiding_center_fixed.c",
    "./simulate/gyro_orbit_fixed.c",
    "./simulate/nbi_source.c",
    "./simulate/orbit_fl_cashkarp.c",
    "./simulate/orbit_gc_cashkarp.c",
    "./simulate/orbit_gc_rk4.c",
    "./simulate/orbit_go_vpa.c",
    "./utils/boschhale.c",
    "./utils/gctransform.c",
    "./utils/interp1Dcomp.c",
    "./utils/interp2Dcomp.c",
    "./utils/interp3Dcomp.c",
    "./utils/linint.c",
    "./utils/list.c",
    "./utils/math.c",
    "./utils/octree.c",
    "./utils/random.c",
    "./utils/splinecomp.c",
    "./utils/suzuki.c"};

int get_endcond(const char *name)
{
    if (strcmp(name, "reached_time_limit") == 0)
        return ENDCOND_TLIM;
    if (strcmp(name, "below_min_energy") == 0)
        return ENDCOND_EMIN;
    if (strcmp(name, "thermalized") == 0)
        return ENDCOND_THERM;
    if (strcmp(name, "hit_wall") == 0)
        return ENDCOND_WALL;
    if (strcmp(name, "above_rho_limit") == 0)
        return ENDCOND_RHOMAX;
    if (strcmp(name, "below_rho_limit") == 0)
        return ENDCOND_RHOMIN;
    if (strcmp(name, "completed_poloidal_orbits") == 0)
        return ENDCOND_POLMAX;
    if (strcmp(name, "completed_toroidal_orbits") == 0)
        return ENDCOND_TORMAX;
    if (strcmp(name, "simulation_not_finished") == 0)
        return ENDCOND_CPUMAX;
    if (strcmp(name, "finished_gc_in_hybrid_mode") == 0)
        return ENDCOND_HYBRID;
    if (strcmp(name, "neutralized") == 0)
        return ENDCOND_NEUTR;
    if (strcmp(name, "ionized") == 0)
        return ENDCOND_IONIZ;
    return -1;
}