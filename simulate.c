#include "simulate.h"

void sim_init(sim_data* sim, sim_offload_data* offload_data) {
    sim->t0 = offload_data->t0;
    sim->tstep = offload_data->tstep;
    sim->tcollstep = offload_data->tcollstep;
    sim->tmax = offload_data->tmax;
    sim->trstep = offload_data->trstep;
    sim->active_endcond = offload_data->active_endcond;
    sim->emin = offload_data->emin;
}
