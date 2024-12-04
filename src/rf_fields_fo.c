/**
 * @file rf_fields_fo.c
 * @brief 2D electromagnetic fields from TORIC RF wave evaluation.
 *
 * Evaluates the electromagnetic wave components from a typical TORIC
 * solution. The TORIC solution provides with the 
 */
#include <stdio.h>
#include <stdlib.h>
#include "ascot5.h"
#include "error.h"
#include "print.h"
#include "spline/interp.h"
#include "rf_fields_fo.h"

int RF2D_init(RF2D_fields* rfdata, real rmin, real rmax, int nr, \
              real zmin, real zmax, int nz, \
              int ntor, int omega, \
              real* Er_real, real* Ez_real, real* Ephi_real, \
              real* Er_imag, real* Ez_imag, real* Ephi_imag, \
              real* Br_real, real* Bz_real, real* Bphi_real, \
              real* Br_imag, real* Bz_imag, real* Bphi_imag){
    a5err err = 0; // Error flag

    // We will simplify the way to initialize all the structures,
    // by using a for loop and a pointer to the structure.
    real* data[12];

    rfdata->introbj[0] = &rfdata->Er_real;
    data[0] = Er_real;
    rfdata->introbj[1] = &rfdata->Er_imag;
    data[1] = Er_imag;
    rfdata->introbj[2] = &rfdata->Ephi_real;
    data[2] = Ephi_real;
    rfdata->introbj[3] = &rfdata->Ephi_imag;
    data[3] = Ephi_imag;
    rfdata->introbj[4] = &rfdata->Ez_real;
    data[4] = Ez_real;
    rfdata->introbj[5] = &rfdata->Ez_imag;
    data[5] = Ez_imag;
    rfdata->introbj[6] = &rfdata->Br_real;
    data[6] = Br_real;
    rfdata->introbj[7] = &rfdata->Br_imag;
    data[7] = Br_imag;
    rfdata->introbj[8] = &rfdata->Bphi_real;
    data[8] = Bphi_real;
    rfdata->introbj[9] = &rfdata->Bphi_imag;
    data[9] = Bphi_imag;
    rfdata->introbj[10] = &rfdata->Bz_real;
    data[10] = Bz_real;
    rfdata->introbj[11] = &rfdata->Bz_imag;
    data[11] = Bz_imag;

    // Setting some data.
    rfdata->ntor = ntor;
    rfdata->omega = omega;

    // Initialize the structures
    for(int i = 0; i < 12; i++){
        err += interp2Dcomp_setup(rfdata->introbj[i], data[i], nr, nz, \
                                  NATURALBC, NATURALBC, \
                                  rmin, rmax, zmin, zmax);
        if(err){
            // Setting the error message.
            print_err("Error: Failed to initialize splines.\n");
            break;
        }
    }
    
    free(rfdata->introbj);
    free(data);
    return err;

}

void RF2D_free(RF2D_fields* rfdata){
    // Free the allocated memory for the splines.
    free(rfdata->Er_real.c);
    free(rfdata->Er_imag.c);
    free(rfdata->Ez_real.c);
    free(rfdata->Ez_imag.c);
    free(rfdata->Ephi_real.c);
    free(rfdata->Ephi_imag.c);
    free(rfdata->Br_real.c);
    free(rfdata->Br_imag.c);
    free(rfdata->Bz_real.c);
    free(rfdata->Bz_imag.c);
    free(rfdata->Bphi_real.c);
    free(rfdata->Bphi_imag.c);
    rfdata->ntor = 0;
    rfdata->omega = 0;

    // Clearing the pointers.
    for(int i = 0; i < 12; i++){
        rfdata->introbj[i] = NULL;
    }
}

void RF_offload(RF2D_fields* rfdata){
    int nsize = rfdata->Er_real.n_x * rfdata->Er_real.n_y;
    // Offload the data to the accelerator.
    GPU_MAP_TO_DEVICE(rfdata->Er_real, rfdata->Er_real.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Er_imag, rfdata->Er_imag.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Ez_real, rfdata->Ez_real.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Ez_imag, rfdata->Ez_imag.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Ephi_real, rfdata->Ephi_real.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Ephi_imag, rfdata->Ephi_imag.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Br_real, rfdata->Br_real.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Br_imag, rfdata->Br_imag.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Bz_real, rfdata->Bz_real.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Bz_imag, rfdata->Bz_imag.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Bphi_real, rfdata->Bphi_real.c[0:nsize*NSIZE_COMP2D], \
                      rfdata->Bphi_imag, rfdata->Bphi_imag.c[0:nsize*NSIZE_COMP2D],
                      rfdata->ntor, rfdata->omega);
}

a5err RF_field_eval(real E[3], real B[3], real r, real phi,\
                    real z, real t, RF2D_fields* rfdata){
    a5err err = 0; // Error flag

    // Interpolation error flag
    int interperr = 0;

    // Interpolated values, ordered as follows:
    // 0:Er_real, 1:Er_imag, 2:Ephi_real, 3:Ephi_imag, 4:Ez_real, 5:Ez_imag,
    // 6:Br_real, 7:Br_imag, 8:Bphi_real, 9:Bphi_imag, 10:z_real, 11:Bz_imag
    real interpolated[12];


    // Interpolating the values.
    for(int k = 0; k < 12; k++){
        interperr = interp2Dcomp_eval_f(&interpolated[k], rfdata->introbj[k], r, z);
        if(interperr){
            // Setting the error message.
            print_err("Error: Failed to evaluate the field.\n");
            return error_raise(ERR_INPUT_EVALUATION, __LINE__, EF_RF_FIELDS_FO);
        }
    }

    // Computing the actual electric and magnetic waves.
    real phase = rfdata->omega * t + rfdata->ntor * phi;
    real cosphase = cos(phase);
    real sinphase = sin(phase);
    E[0] = interpolated[0] * cosphase + interpolated[1] * sinphase; // Er
    E[1] = interpolated[2] * cosphase + interpolated[3] * sinphase; // Ephi
    E[2] = interpolated[4] * cosphase + interpolated[5] * sinphase; // Ez
    B[0] = interpolated[6] * cosphase + interpolated[7] * sinphase; // Br
    B[1] = interpolated[8] * cosphase + interpolated[9] * sinphase; // Bphi
    B[2] = interpolated[10] * cosphase + interpolated[11] * sinphase; // Bz

    return err;
}