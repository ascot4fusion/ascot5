/**
 * Implements diag.h.
 */
#include "diag.h"
#include "bfield.h"
#include "datatypes.h"
#include "defines.h"
#include "diag_hist.h"
#include "diag_orbit.h"
#include "marker.h"

void Diag_offload(Diagnostics *diag)
{
    for (size_t i = 0; i < diag->nhist; i++)
    {
        DiagHist_offload(&diag->hist[i]);
    }
    if (diag->orbit != NULL)
    {
        DiagOrbit_offload(diag->orbit);
    }
}

void Diag_update_go(
    Diagnostics *diag, Bfield *bfield, MarkerGyroOrbit *mrk_f,
    MarkerGyroOrbit *mrk_i)
{
    for (size_t i = 0; i < diag->nhist; i++)
    {
        DiagHist_update_go(&diag->hist[i], mrk_f, mrk_i);
    }
    if (diag->orbit != NULL)
    {
        DiagOrbit_update_go(diag->orbit, bfield, mrk_f, mrk_i);
    }
}

void Diag_update_gc(
    Diagnostics *diag, Bfield *bfield, MarkerGuidingCenter *mrk_f,
    MarkerGuidingCenter *mrk_i)
{
    for (size_t i = 0; i < diag->nhist; i++)
    {
        DiagHist_update_gc(&diag->hist[i], mrk_f, mrk_i);
    }
    if (diag->orbit != NULL)
    {
        DiagOrbit_update_gc(diag->orbit, bfield, mrk_f, mrk_i);
    }
}

void Diag_update_fl(
    Diagnostics *diag, Bfield *bfield, MarkerFieldLine *mrk_f,
    MarkerFieldLine *mrk_i)
{
    if (diag->orbit != NULL)
    {
        DiagOrbit_update_fl(diag->orbit, bfield, mrk_f, mrk_i);
    }
}
