/**
 * Implements asigma.h.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ascot5.h"
#include "error.h"
#include "math.h"
#include "asigma.h"
#include "asigma/asigma_loc.h"
#include "consts.h"

#pragma omp declare target
/** Set values outside abscissae to zero instead of raising an error. */
static int ASIGMA_EXTRAPOLATE = 0;
#pragma omp end declare target


void asigma_free(asigma_data* atomic) {
    switch(atomic->type) {
        case asigma_type_loc:
            asigma_loc_free(atomic->asigma_loc);
            break;
    }
}


void asigma_offload(asigma_data* atomic) {
    switch(atomic->type) {
        case asigma_type_loc:
            asigma_loc_offload(atomic->asigma_loc);
            break;
    }
}


void asigma_extrapolate(int extrapolate) {
    ASIGMA_EXTRAPOLATE = extrapolate;
}


a5err asigma_eval_sigma(
    real* sigma, int z_1, int a_1, int z_2, int a_2, real E_coll_per_amu,
    asigma_reac_type reac_type, asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_sigma(
                sigma, z_1, a_1, z_2, a_2, E_coll_per_amu, reac_type,
                ASIGMA_EXTRAPOLATE, asigma_data->asigma_loc);
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || sigma[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        sigma[0] = 0.0;
    }

    return err;
}


a5err asigma_eval_sigmav(
    real* sigmav, int z_1, int a_1, real m_1, int z_2, int a_2,
    real E, real T_e, real T_0, real n_i, asigma_reac_type reac_type,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_sigmav(
                sigmav, z_1, a_1, m_1, z_2, a_2, E, T_e, T_0, n_i,
                reac_type, ASIGMA_EXTRAPOLATE, asigma_data->asigma_loc);
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || sigmav[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        sigmav[0] = 0.0;
    }

    return err;
}


a5err asigma_eval_cx(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nspec,
    const int* znum, const int* anum, real T_0, real* n_0,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_cx(
                    ratecoeff, z_1, a_1, E, mass, nspec, znum, anum, T_0, n_0,
                    ASIGMA_EXTRAPOLATE, asigma_data->asigma_loc);
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || ratecoeff[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        ratecoeff[0] = 0.0;
    }

    return err;
}


a5err asigma_eval_bms(
    real* ratecoeff, int z_1, int a_1, real E, real mass, int nion,
    const int* znum, const int* anum, real T_e, real* n_i,
    asigma_data* asigma_data) {
    a5err err = 0;

    switch(asigma_data->type) {
        case asigma_type_loc:
            err = asigma_loc_eval_bms(
                    ratecoeff, z_1, a_1, E, mass, nion, znum, anum, T_e, n_i,
                    ASIGMA_EXTRAPOLATE, asigma_data->asigma_loc);
            break;

        default:
            /* Unrecognized input. Produce error. */
            err = error_raise( ERR_UNKNOWN_INPUT, __LINE__, EF_ASIGMA );
            break;
    }
    if(err || ratecoeff[0] < 0.0) {
        /* In case of error or unphysical negative value, return zero value
           to avoid further complications */
        ratecoeff[0] = 0.0;
    }

    return err;
}
