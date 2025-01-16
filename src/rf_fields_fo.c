/**
 * @file rf_fields_fo.c
 * @brief 2D electromagnetic fields from TORIC RF wave evaluation.
 *
 * Evaluates the electromagnetic wave components from a typical TORIC
 * solution. The TORIC solution provides with the 
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "spline/interp.h"
#include "rf_fields_fo.h"

int rffield_init(RF2D_fields* rffield_data, real rmin, real rmax, int nr, \
              real zmin, real zmax, int nz, \
              int ntor, real omega, \
              real* Er_real, real* Er_imag, real* Ephi_real, \
              real* Ephi_imag, real* Ez_real, real* Ez_imag, \
              real* Br_real, real* Br_imag, real* Bphi_real, \
              real* Bphi_imag, real* Bz_real, real* Bz_imag){
    a5err err = 0; // Error flag

    // We will simplify the way to initialize all the structures,
    // by using a for loop and a pointer to the structure.
    real* data[12];

    rffield_data->introbj[0] = &rffield_data->Er_real;
    data[0] = Er_real;
    rffield_data->introbj[1] = &rffield_data->Er_imag;
    data[1] = Er_imag;
    rffield_data->introbj[2] = &rffield_data->Ephi_real;
    data[2] = Ephi_real;
    rffield_data->introbj[3] = &rffield_data->Ephi_imag;
    data[3] = Ephi_imag;
    rffield_data->introbj[4] = &rffield_data->Ez_real;
    data[4] = Ez_real;
    rffield_data->introbj[5] = &rffield_data->Ez_imag;
    data[5] = Ez_imag;
    rffield_data->introbj[6] = &rffield_data->Br_real;
    data[6] = Br_real;
    rffield_data->introbj[7] = &rffield_data->Br_imag;
    data[7] = Br_imag;
    rffield_data->introbj[8] = &rffield_data->Bphi_real;
    data[8] = Bphi_real;
    rffield_data->introbj[9] = &rffield_data->Bphi_imag;
    data[9] = Bphi_imag;
    rffield_data->introbj[10] = &rffield_data->Bz_real;
    data[10] = Bz_real;
    rffield_data->introbj[11] = &rffield_data->Bz_imag;
    data[11] = Bz_imag;

    // Setting some data.
    rffield_data->ntor = ntor;
    rffield_data->omega = omega;

    // Initialize the structures
    for(int i = 0; i < 12; i++){
        err += interp2Dcomp_setup(rffield_data->introbj[i], \
                                  data[i], nr, nz, \
                                  NATURALBC, NATURALBC, \
                                  rmin, rmax, zmin, zmax);
        if(err) break;
    }

    rffield_data->initialized = 1;
    return err;

}

void rffield_free(RF2D_fields* rffield_data){
    if(!rffield_data->initialized) return; // Nothing to do here.

    // Free the allocated memory for the splines.
    if(rffield_data->Er_real.c != NULL) free(rffield_data->Er_real.c);
    if(rffield_data->Er_imag.c != NULL) free(rffield_data->Er_imag.c);
    if(rffield_data->Ez_real.c != NULL) free(rffield_data->Ez_real.c);
    if(rffield_data->Ez_imag.c != NULL) free(rffield_data->Ez_imag.c);
    if(rffield_data->Ephi_real.c != NULL) free(rffield_data->Ephi_real.c);
    if(rffield_data->Ephi_imag.c != NULL) free(rffield_data->Ephi_imag.c);
    if(rffield_data->Br_real.c != NULL) free(rffield_data->Br_real.c);
    if(rffield_data->Br_imag.c != NULL) free(rffield_data->Br_imag.c);
    if(rffield_data->Bz_real.c != NULL) free(rffield_data->Bz_real.c);
    if(rffield_data->Bz_imag.c != NULL) free(rffield_data->Bz_imag.c);
    if(rffield_data->Ephi_real.c != NULL) free(rffield_data->Bphi_real.c);
    if(rffield_data->Ephi_imag.c != NULL) free(rffield_data->Bphi_imag.c);
    rffield_data->ntor = 0;
    rffield_data->omega = 0.0;
    rffield_data->initialized = 0;

    // Clearing the pointers.
    for(int i = 0; i < 12; i++){
        rffield_data->introbj[i] = NULL;
    }
}

void RF2D_offload(RF2D_fields* rffield_data){
    if(!rffield_data->initialized) return; // Nothing to do here.

    int nsize = rffield_data->Er_real.n_x * rffield_data->Er_real.n_y;
    // Offload the data to the accelerator.
    GPU_MAP_TO_DEVICE(rffield_data->Er_real, rffield_data->Er_real.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Er_imag, rffield_data->Er_imag.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Ez_real, rffield_data->Ez_real.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Ez_imag, rffield_data->Ez_imag.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Ephi_real, rffield_data->Ephi_real.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Ephi_imag, rffield_data->Ephi_imag.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Br_real, rffield_data->Br_real.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Br_imag, rffield_data->Br_imag.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Bz_real, rffield_data->Bz_real.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Bz_imag, rffield_data->Bz_imag.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Bphi_real, rffield_data->Bphi_real.c[0:nsize*NSIZE_COMP2D], \
                      rffield_data->Bphi_imag, rffield_data->Bphi_imag.c[0:nsize*NSIZE_COMP2D],
                      rffield_data->ntor, rffield_data->omega);
}

a5err RF_field_eval(real E[3], real B[3], real r, real phi,\
                    real z, real t, RF2D_fields* rffield_data){
    a5err err = 0; // Error flag

    // Interpolation error flag
    int interperr = 0;

    // Interpolated values, ordered as follows:
    // 0:Er_real, 1:Er_imag, 2:Ephi_real, 3:Ephi_imag, 4:Ez_real, 5:Ez_imag,
    // 6:Br_real, 7:Br_imag, 8:Bphi_real, 9:Bphi_imag, 10:z_real, 11:Bz_imag
    real interpolated[12];

    if(!rffield_data->initialized){
        // The RF fields have not been initialized, so
        // RF fields are not available.
        E[0] = 0.0;
        E[1] = 0.0;
        E[2] = 0.0;
        B[0] = 0.0;
        B[1] = 0.0;
        B[2] = 0.0;
        return 0;
    }


    // Interpolating the values.
    for(int k = 0; k < 12; k++){
        interperr = interp2Dcomp_eval_f(&interpolated[k], rffield_data->introbj[k], r, z);
        if(interperr){
            E[0] = 0.0; E[1] = 0.0; E[2] = 0.0;
            B[0] = 0.0; B[1] = 0.0; B[2] = 0.0;
            return 0; 
        }
    }

    // Computing the actual electric and magnetic waves.
    real phase = rffield_data->omega * t + rffield_data->ntor * phi;
    real cosphase = cos(phase);
    real sinphase = sin(phase);
    E[0] = interpolated[0]  * cosphase + interpolated[1]  * sinphase; // Er
    E[1] = interpolated[2]  * cosphase + interpolated[3]  * sinphase; // Ephi
    E[2] = interpolated[4]  * cosphase + interpolated[5]  * sinphase; // Ez
    B[0] = interpolated[6]  * cosphase + interpolated[7]  * sinphase; // Br
    B[1] = interpolated[8]  * cosphase + interpolated[9]  * sinphase; // Bphi
    B[2] = interpolated[10] * cosphase + interpolated[11] * sinphase; // Bz

    return err;
}