/**
 * Implements diag.h.
 */
#include "diag.h"
#include "bfield.h"
#include "defines.h"
#include "diag_orb.h"
#include "diag_transcoef.h"
#include "dist_5D.h"
#include "dist_6D.h"
#include "dist_com.h"
#include "dist_rho5D.h"
#include "dist_rho6D.h"
#include "marker.h"
#include "options.h"
#include "datatypes.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void diag_free(Diagnostics *diag, Options *options)
{
    if (options->collect_dist5d)
    {
        dist_5D_free(diag->dist5d);
    }
    if (options->collect_dist6d)
    {
        dist_6D_free(diag->dist6d);
    }
    if (options->collect_dist5drho)
    {
        dist_rho5D_free(diag->dist5drho);
    }
    if (options->collect_dist6drho)
    {
        dist_rho6D_free(diag->dist6drho);
    }
    if (options->collect_distcom)
    {
        dist_COM_free(diag->distcom);
    }
    if (options->collect_orbit)
    {
        diag_orb_free(diag->orbit, options);
    }
    if (options->collect_transport_coefficient)
    {
        diag_transcoef_free(diag->transport_coefficient);
    }
}

void diag_offload(Diagnostics *diag, Options *options)
{
    if (options->collect_dist5d)
    {
        dist_5D_offload(diag->dist5d);
    }
    if (options->collect_dist6d)
    {
        dist_6D_offload(diag->dist6d);
    }
    if (options->collect_dist5drho)
    {
        dist_rho5D_offload(diag->dist5drho);
    }
    if (options->collect_dist6drho)
    {
        dist_rho6D_offload(diag->dist6drho);
    }
    if (options->collect_distcom)
    {
        dist_COM_offload(diag->distcom);
    }
}

void diag_update_go(
    Diagnostics *diag, Options *options, Bfield *bfield, MarkerGyroOrbit *p_f,
    MarkerGyroOrbit *p_i)
{
    if (options->collect_orbit)
    {
        diag_orb_update_fo(diag->orbit, options, p_f, p_i);
    }
    if (options->collect_dist5d)
    {
        dist_5D_update_fo(diag->dist5d, options, p_f, p_i);
    }
    if (options->collect_dist6d)
    {
        dist_6D_update_fo(diag->dist6d, p_f, p_i);
    }
    if (options->collect_dist5drho)
    {
        dist_rho5D_update_fo(diag->dist5drho, p_f, p_i);
    }
    if (options->collect_dist6drho)
    {
        dist_rho6D_update_fo(diag->dist6drho, p_f, p_i);
    }
    if (options->collect_distcom)
    {
        dist_COM_update_fo(diag->distcom, bfield, p_f, p_i);
    }
    if (options->collect_transport_coefficient)
    {
        diag_transcoef_update_fo(
            diag->transport_coefficient, options, p_f, p_i);
    }
}

void diag_update_gc(
    Diagnostics *diag, Options *options, Bfield *bfield,
    MarkerGuidingCenter *p_f, MarkerGuidingCenter *p_i)
{
    if (options->collect_orbit)
    {
        diag_orb_update_gc(diag->orbit, options, p_f, p_i);
    }
    if (options->collect_dist5d)
    {
        dist_5D_update_gc(diag->dist5d, options, p_f, p_i);
    }
    if (options->collect_dist6d)
    {
        dist_6D_update_gc(diag->dist6d, p_f, p_i);
    }
    if (options->collect_dist5drho)
    {
        dist_rho5D_update_gc(diag->dist5drho, p_f, p_i);
    }
    if (options->collect_dist6drho)
    {
        dist_rho6D_update_gc(diag->dist6drho, p_f, p_i);
    }
    if (options->collect_distcom)
    {
        dist_COM_update_gc(diag->distcom, bfield, p_f, p_i);
    }
    if (options->collect_transport_coefficient)
    {
        diag_transcoef_update_gc(
            diag->transport_coefficient, options, p_f, p_i);
    }
}

void diag_update_fl(
    Diagnostics *diag, Options *options, MarkerFieldLine *p_f,
    MarkerFieldLine *p_i)
{
    if (options->collect_orbit)
    {
        diag_orb_update_ml(diag->orbit, options, p_f, p_i);
    }
    if (options->collect_transport_coefficient)
    {
        diag_transcoef_update_ml(
            diag->transport_coefficient, options, p_f, p_i);
    }
}
