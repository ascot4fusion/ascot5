/**
 * @author Pablo Oyola oyola@pppl.gov
 * @file RF_full_orbit_3d.h
 * @brief Header file for 3D RF fields evaluation for the full-orbit case.
 *
 * Contains the declaration of the input structure for the evaluation
 * of the RF fields for the full-orbit cases.
 */

#ifndef RF_FULL_ORBIT_3D_H
#define RF_FULL_ORBIT_3D_H

#include "../offload.h"
#include "../ascot5.h"
#include "../error.h"
#include "../spline/interp.h"


// Defining the object structure of the RF fields.
typedef struct RF3D_fields{
    interp3D_data Er_real;     /**< Interpolation 3D spline struct for radial electric field    */
    interp3D_data Er_imag;     /**< Interpolation 3D spline struct for radial electric field    */
    interp3D_data Ez_real;     /**< Interpolation 3D spline struct for vertical electric field  */
    interp3D_data Ez_imag;     /**< Interpolation 3D spline struct for vertical electric field  */
    interp3D_data Ephi_real;   /**< Interpolation 3D spline struct for toroidal electric field  */
    interp3D_data Ephi_imag;   /**< Interpolation 3D spline struct for toroidal electric field  */
    interp3D_data Br_real;     /**< Interpolation 3D spline struct for radial magnetic field    */
    interp3D_data Br_imag;     /**< Interpolation 3D spline struct for radial magnetic field    */
    interp3D_data Bz_imag;     /**< Interpolation 3D spline struct for vertical magnetic field  */
    interp3D_data Bz_real;     /**< Interpolation 3D spline struct for vertical magnetic field  */
    interp3D_data Bphi_imag;   /**< Interpolation 3D spline struct for toroidal magnetic field  */
    interp3D_data Bphi_real;   /**< Interpolation 3D spline struct for toroidal magnetic field  */
    interp3D_data* introbj[12];   /**< Simplified pointer to the interpolation objects.           */
    int ntor;                  /**< Toroidal mode number.                           */
    real omega;                 /**< Wave frequency                                  */
    int initialized;       /**< Flag to check if the object has been initialized */
}RF3D_fields;

a5err RF3D_fields_init_from_file(RF3D_fields* rffield_data, hid_t f, char* qid);
int RF3D_fields_init( RF3D_fields* rffield_data, real rmin, real rmax, int nr, \
               real zmin, real zmax, int nz, int nphi, \
               int ntor, real omega, \
               real* Er_real, real* Er_imag, real* Ephi_real, \
               real* Ephi_imag, real* Ez_real, real* Ez_imag, \
               real* Br_real, real* Br_imag, real* Bphi_real, \
               real* Bphi_imag, real* Bz_real, real* Bz_imag);
void RF3D_fields_free(RF3D_fields* rffield_data);
void RF3D_fields_offload(RF3D_fields* rffield_data);

GPU_DECLARE_TARGET_SIMD_UNIFORM(rffield_data)
a5err RF3D_field_eval(real E[3], real B[3], real r, real phi,\
                       real z, real t, RF3D_fields* rffield_data);
DECLARE_TARGET_END

#endif