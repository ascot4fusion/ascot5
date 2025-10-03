#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "B_field.h"
#include "interp.h"
#include "E_1DS.h"


int EfieldPotential1D_init(
    EfieldPotential1D* data, int nrho, real reff, real rholim[2],
    real dvdrho[nrho]
) {
    data->reff = reff;
    int err = interp1Dcomp_setup(
        &data->dv, dvdrho, nrho, NATURALBC, rholim[0], rholim[1]
    );
    return err;
}


void EfieldPotential1D_free(EfieldPotential1D* efield) {
    free(efield->dv.c);
}


void EfieldPotential1D_offload(EfieldPotential1D* efield) {
    (void)efield;
    GPU_MAP_TO_DEVICE( efield->dV, efield->dV.c[0:efield->dv.n_x*NSIZE_COMP1D] )
}


a5err EfieldPotential1D_eval_e(
    real e[3], real r, real phi, real z, EfieldPotential1D* efield,
    B_field_data* bfield
) {
    a5err err = 0;
    int interperr = 0;
    real rho_drho[4];
    err = B_field_eval_rho_drho(rho_drho, r, phi, z, bfield);
    if(err) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_1DS );
    }
    /* Convert partial derivative to gradient */
    rho_drho[2] = rho_drho[2]/r;
    /* We set the field to zero if outside the profile. */
    if (rho_drho[0] < efield->dv.x_min || rho_drho[0] > efield->dv.x_max ) {
        e[0] = 0;
        e[1] = 0;
        e[2] = 0;
    } else {
        real dV;
        interperr += interp1Dcomp_eval_f(&dV, &efield->dv, rho_drho[0]);
        e[0] = -dV * rho_drho[1] * efield->reff;
        e[1] = -dV * rho_drho[2] * efield->reff;
        e[2] = -dV * rho_drho[3] * efield->reff;

        if(interperr) {
            err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_E_1DS );
        }
    }
    return err;
}
