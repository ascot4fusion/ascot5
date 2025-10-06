/**
 * Implements errors.h.
 */
#include "errors.h"
#include <math.h>
#include <stdio.h>

void error_parse(a5err err, int *msg, int *line, int *file)
{
    *file = (int)(floor(err / (256 * 1024)));
    err -= *file * (256 * 1024);
    *line = (int)(floor(err / 256));
    err -= *line * (256);
    *msg = err;
}

void error_parse2str(a5err err, char *msg, char *line, char *file)
{
    int msg_i, line_i, file_i;
    error_parse(err, &msg_i, &line_i, &file_i);

    /* Find the error message */
    switch (msg_i)
    {
    case ERR_INPUT_EVALUATION:
        sprintf(
            msg, "Input evaluation failed "
                 "(marker could be outside input data grid)");
        break;
    case ERR_UNKNOWN_INPUT:
        sprintf(
            msg, "Input was not regonized "
                 "(offload or target data could be uninitialized)");
        break;
    case ERR_INPUT_UNPHYSICAL:
        sprintf(
            msg, "Input evaluation yields unphysical results "
                 "(something could be wrong with the input)");
        break;
    case ERR_MARKER_UNPHYSICAL:
        sprintf(
            msg, "One or more of marker's fields are unphysical or "
                 "inconsistent (marker input could be corrupted)");
        break;
    case ERR_INVALID_TIMESTEP:
        sprintf(
            msg, "Time step is zero, NaN, or smaller than MIN_ALLOWED"
                 "_TIME_STEP in ascot5.h (time step limits could be too "
                 "conservative)");
        break;
    case ERR_WIENER_ARRAY:
        sprintf(
            msg, "Wiener array is full of rejected steps or  "
                 "limits could be too conservative or initial step too "
                 "large)");
        break;
    case ERR_ATOMIC_EVALUATION:
        sprintf(
            msg, "Atomic reaction evaluation failed "
                 "(atomic reactions for the current charge state might"
                 "not have been implemented)");
        break;
    default:
        sprintf(msg, "Unknown error");
        break;
    }

    /* Line number to string */
    sprintf(line, "%4d", line_i);

    /* Find the file where error originated from */
    switch (file_i)
    {
    case EF_MCCC_WIENER:
        sprintf(file, "mccc_wiener.c");
        break;
    case EF_MCCC_PUSH:
        sprintf(file, "mccc_push.c");
        break;
    case EF_MCCC_COEFS:
        sprintf(file, "mccc_coefs.c");
        break;
    case EF_MCCC:
        sprintf(file, "mccc.c");
        break;
    case EF_STEP_FO_VPA:
        sprintf(file, "step_fo_vpa.c");
        break;
    case EF_STEP_GC_RK4:
        sprintf(file, "step_gc_rk4.c");
        break;
    case EF_STEP_GC_CASHKARP:
        sprintf(file, "step_gc_cashkarp.c");
        break;
    case EF_PLASMA:
        sprintf(file, "plasma.c");
        break;
    case EF_PLASMA_1D:
        sprintf(file, "plasma_1D.c");
        break;
    case EF_PLASMA_1DS:
        sprintf(file, "plasma_1DS.c");
        break;
    case EF_E_FIELD:
        sprintf(file, "E_field.c");
        break;
    case EF_E_1DS:
        sprintf(file, "E_1DS.c");
        break;
    case EF_NEUTRAL:
        sprintf(file, "neutral.c");
        break;
    case EF_N0_1D:
        sprintf(file, "N0_1D.c");
        break;
    case EF_N0_3D:
        sprintf(file, "N0_3D.c");
        break;
    case EF_N0_ST:
        sprintf(file, "N0_ST.c");
        break;
    case EF_B_FIELD:
        sprintf(file, "B_field.c");
        break;
    case EF_B_GS:
        sprintf(file, "B_GS.c");
        break;
    case EF_B_STS:
        sprintf(file, "B_STS.c");
        break;
    case EF_B_2DS:
        sprintf(file, "B_2DS.c");
        break;
    case EF_B_3DS:
        sprintf(file, "B_3DS.c");
        break;
    case EF_PARTICLE:
        sprintf(file, "particle.c");
        break;
    case EF_BOOZER:
        sprintf(file, "boozer.c");
        break;
    case EF_MHD:
        sprintf(file, "mhd.c");
        break;
    case EF_ATOMIC:
        sprintf(file, "atomic.c");
        break;
    case EF_ASIGMA:
        sprintf(file, "asigma.c");
        break;
    case EF_ASIGMA_LOC:
        sprintf(file, "asigma_loc.c");
        break;
    default:
        sprintf(file, "unknown file");
        break;
    }
}
