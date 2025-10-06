/**
 * Implements N0_1D.h.
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "mathlib.h"
#include "defines.h"
#include "errors.h"
#include "neutral.h"
#include "linint.h"


int N0_1D_init(N0_1D_data* data, int n_rho, real rho_min, real rho_max,
               int n_species, int* anum, int* znum,
               real* density, real* temperature) {

    data->n_species = n_species;
    data->anum = (int*) malloc(n_species * sizeof(int));
    data->znum = (int*) malloc(n_species * sizeof(int));
    data->n0 = (linint1D_data*) malloc( n_species * sizeof(linint1D_data) );
    data->t0 = (linint1D_data*) malloc( n_species * sizeof(linint1D_data) );
    int err = 0;
    for(int i = 0; i < data->n_species; i++) {
        data->anum[i] = anum[i];
        data->znum[i] = znum[i];

        real* c = (real*) malloc(n_rho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for(int i = 0; i < n_rho; i++) {
            c[i] = density[i];
        }
        linint1D_init(&data->n0[i], c, n_rho, NATURALBC, rho_min, rho_max);
        c = (real*) malloc(n_rho * sizeof(real));
        err += c == NULL ? 1 : 0;
        for(int i = 0; i < n_rho; i++) {
            c[i] = temperature[i];
        }
        linint1D_init(&data->t0[i], c, n_rho, NATURALBC, rho_min, rho_max);
    }
    return err;
}


void N0_1D_free(N0_1D_data* data) {
    free(data->anum);
    free(data->znum);
    for(int i = 0; i < data->n_species; i++) {
        free(data->n0[i].c);
        free(data->t0[i].c);
    }
    free(data->n0);
    free(data->t0);
}


void N0_1D_offload(N0_1D_data* data) {
    SUPPRESS_UNUSED_WARNING(data);
}


a5err N0_1D_eval_n0(real* n0, real rho, N0_1D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for(int i=0; i<ndata->n_species; i++) {
        interperr += linint1D_eval_f(&n0[i], &ndata->n0[i], rho);
    }

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}


a5err N0_1D_eval_t0(real* t0, real rho, N0_1D_data* ndata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    for(int i=0; i<ndata->n_species; i++) {
        interperr += linint1D_eval_f(&t0[i], &ndata->t0[i], rho);
    }

    if(interperr) {
        return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_N0_1D);
    }

    return err;
}


int N0_1D_get_n_species(N0_1D_data* ndata) {
    return ndata->n_species;
}
