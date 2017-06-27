#include "simulate.h"

void sim_init(sim_data* sim, sim_offload_data* offload_data) {
    sim->sim_mode             = offload_data->sim_mode;
    sim->enable_ada           = offload_data->enable_ada;
    sim->record_GOasGC        = offload_data->record_GOasGC;

    sim->fix_usrdef_use       = offload_data->fix_usrdef_use;
    sim->fix_usrdef_val       = offload_data->fix_usrdef_val;
    sim->fix_stepsPerGO       = offload_data->fix_stepsPerGO;

    sim->ada_tol_orbfol       = offload_data->ada_tol_orbfol;
    sim->ada_tol_clmbcol      = offload_data->ada_tol_clmbcol;
    sim->ada_max_drho         = offload_data->ada_max_drho;
    sim->ada_max_dphi         = offload_data->ada_max_dphi;
    sim->ada_max_acc          = offload_data->ada_max_acc;

    sim->enable_orbfol        = offload_data->enable_orbfol;
    sim->enable_clmbcol       = offload_data->enable_clmbcol;

    sim->endcond_active       = offload_data->endcond_active;
    sim->endcond_maxSimTime   = offload_data->endcond_maxSimTime;
    sim->endcond_maxCpuTime   = offload_data->endcond_maxCpuTime;
    sim->endcond_minRho       = offload_data->endcond_minRho;
    sim->endcond_maxRho       = offload_data->endcond_maxRho;
    sim->endcond_minEkin      = offload_data->endcond_minEkin;
    sim->endcond_minEkinPerTe = offload_data->endcond_minEkinPerTe;
    sim->endcond_maxTorOrb    = offload_data->endcond_maxTorOrb;
    sim->endcond_maxPolOrb    = offload_data->endcond_maxPolOrb;
}
