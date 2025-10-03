#include "diag.h"
#include "B_field.h"
#include "ascot5.h"
#include "diag/diag_orb.h"
#include "diag/diag_transcoef.h"
#include "diag/dist_5D.h"
#include "diag/dist_6D.h"
#include "diag/dist_com.h"
#include "diag/dist_rho5D.h"
#include "diag/dist_rho6D.h"
#include "options.h"
#include "particle.h"
#include "simulate.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void diag_free(Diagnostics *diag, sim_parameters *params)
{
    if (params->collect_dist5d)
    {
        dist_5D_free(diag->dist5d);
    }
    if (params->collect_dist6d)
    {
        dist_6D_free(diag->dist6d);
    }
    if (params->collect_dist5drho)
    {
        dist_rho5D_free(diag->dist5drho);
    }
    if (params->collect_dist6drho)
    {
        dist_rho6D_free(diag->dist6drho);
    }
    if (params->collect_distcom)
    {
        dist_COM_free(diag->distcom);
    }
    if (params->collect_orbit)
    {
        diag_orb_free(diag->orbit, params);
    }
    if (params->collect_transport_coefficient)
    {
        diag_transcoef_free(diag->transport_coefficient);
    }
}

void diag_offload(Diagnostics *diag, sim_parameters *params)
{
    if (params->collect_dist5d)
    {
        dist_5D_offload(diag->dist5d);
    }
    if (params->collect_dist6d)
    {
        dist_6D_offload(diag->dist6d);
    }
    if (params->collect_dist5drho)
    {
        dist_rho5D_offload(diag->dist5drho);
    }
    if (params->collect_dist6drho)
    {
        dist_rho6D_offload(diag->dist6drho);
    }
    if (params->collect_distcom)
    {
        dist_COM_offload(diag->distcom);
    }
}

void diag_update_fo(
    Diagnostics *diag, sim_parameters *params, B_field_data *Bdata,
    particle_simd_fo *p_f, particle_simd_fo *p_i)
{
    if (params->collect_orbit)
    {
        diag_orb_update_fo(diag->orbit, params, p_f, p_i);
    }
    if (params->collect_dist5d)
    {
        dist_5D_update_fo(diag->dist5d, params, p_f, p_i);
    }
    if (params->collect_dist6d)
    {
        dist_6D_update_fo(diag->dist6d, p_f, p_i);
    }
    if (params->collect_dist5drho)
    {
        dist_rho5D_update_fo(diag->dist5drho, p_f, p_i);
    }
    if (params->collect_dist6drho)
    {
        dist_rho6D_update_fo(diag->dist6drho, p_f, p_i);
    }
    if (params->collect_distcom)
    {
        dist_COM_update_fo(diag->distcom, Bdata, p_f, p_i);
    }
    if (params->collect_transport_coefficient)
    {
        diag_transcoef_update_fo(diag->transport_coefficient, params, p_f, p_i);
    }
}

void diag_update_gc(
    Diagnostics *diag, sim_parameters *params, B_field_data *Bdata,
    particle_simd_gc *p_f, particle_simd_gc *p_i)
{
    if (params->collect_orbit)
    {
        diag_orb_update_gc(diag->orbit, params, p_f, p_i);
    }
    if (params->collect_dist5d)
    {
        dist_5D_update_gc(diag->dist5d, params, p_f, p_i);
    }
    if (params->collect_dist6d)
    {
        dist_6D_update_gc(diag->dist6d, p_f, p_i);
    }
    if (params->collect_dist5drho)
    {
        dist_rho5D_update_gc(diag->dist5drho, p_f, p_i);
    }
    if (params->collect_dist6drho)
    {
        dist_rho6D_update_gc(diag->dist6drho, p_f, p_i);
    }
    if (params->collect_distcom)
    {
        dist_COM_update_gc(diag->distcom, Bdata, p_f, p_i);
    }
    if (params->collect_transport_coefficient)
    {
        diag_transcoef_update_gc(diag->transport_coefficient, params, p_f, p_i);
    }
}

void diag_update_ml(
    Diagnostics *diag, sim_parameters *params, particle_simd_ml *p_f,
    particle_simd_ml *p_i)
{
    if (params->collect_orbit)
    {
        diag_orb_update_ml(diag->orbit, params, p_f, p_i);
    }
    if (params->collect_transport_coefficient)
    {
        diag_transcoef_update_ml(diag->transport_coefficient, params, p_f, p_i);
    }
}
