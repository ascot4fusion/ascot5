/**
 * @author Pablo Oyola oyola@pppl.gov
 * @file rf_fields_fo.h
 * @brief Header file for rf_fields_fo.c
 *
 * Contains the declaration of the input structure for the evaluation
 * of the RF fields for the full-orbit cases.
 */

#ifndef RF_FIELDS_FO_H
#define RF_FIELDS_FO_H

#include "offload.h"
#include "ascot5.h"
#include "error.h"
#include "interp.h"


// Defining the object structure of the RF fields.
typedef struct RF2D_fields{
    interp2D_data Er_real;     /**< Interpolation 2D spline struct for radial electric field    */
    interp2D_data Er_imag;     /**< Interpolation 2D spline struct for radial electric field    */
    interp2D_data Ez_real;     /**< Interpolation 2D spline struct for vertical electric field  */
    interp2D_data Ez_imag;     /**< Interpolation 2D spline struct for vertical electric field  */
    interp2D_data Ephi_real;   /**< Interpolation 2D spline struct for toroidal electric field  */
    interp2D_data Ephi_imag;   /**< Interpolation 2D spline struct for toroidal electric field  */
    interp2D_data Br_real;     /**< Interpolation 2D spline struct for radial magnetic field    */
    interp2D_data Br_imag;     /**< Interpolation 2D spline struct for radial magnetic field    */
    interp2D_data Bz_imag;     /**< Interpolation 2D spline struct for vertical magnetic field  */
    interp2D_data Bz_real;     /**< Interpolation 2D spline struct for vertical magnetic field  */
    interp2D_data Bphi_imag;   /**< Interpolation 2D spline struct for toroidal magnetic field  */
    interp2D_data Bphi_real;   /**< Interpolation 2D spline struct for toroidal magnetic field  */
    interp2D_data Bphi_imag;   /**< Interpolation 2D spline struct for toroidal magnetic field  */
    interp2D_data* introbj[12];   /**< Simplified pointer to the interpolation objects.           */
    int ntor;                  /**< Toroidal mode number.                           */
    int omega;                 /**< Wave frequency                                  */
}RF2D_fields;


int RF2D_init( RF2D_fields* rfdata, real rmin, real rmax, int nr, \
               real zmin, real zmax, int nz, \
               int ntor, int omega, \
               real* Er_real, real* Ez_real, real* Ephi_real, \
               real* Er_imag, real* Ez_imag, real* Ephi_imag, \
               real* Br_real, real* Bz_real, real* Bphi_real, \
                real* Br_imag, real* Bz_imag, real* Bphi_imag);
void RF2D_free(RF2D_fields* rfdata);
void RF2D_offload(RF2D_fields* rfdata);

GPU_DECLARE_TARGET_SIMD_UNIFORM(rfdata)
a5err RF_field_eval(real E[3], real B[3], real r, real phi,\
                    real z, real t, RF2D_fields* rfdata);
DECLARE_TARGET_END

#endif