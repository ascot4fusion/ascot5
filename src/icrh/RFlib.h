/**
 * @author Pablo Oyola oyola@pppl.gov
 * @file RFlib.h
 * @brief Header file for RF fields evaluation - includes all necessary headers for all the
 * implemented RF models.
 */

#ifndef RFLIB_H
#define RFLIB_H
#include "RF_full_orbit_2d.h"
#include "RF_full_orbit_3d.h"
#include "RF2D_gc_stix.h"
#include "RF_stix_particle_history.h"
#include <hdf5.h>

typedef enum RF_type{
    RF_NONE,
    RF_FULL_ORBIT_2D,
    RF_FULL_ORBIT_3D,
    RF2D_GC_STIX
} RF_type;

typedef struct RF_fields{
    RF_type type; /** Type of RF field */
    union {
        RF2D_fields rf2d;       /** 2D full-orbit RF fields */
        RF3D_fields rf3d;       /** 3D full-orbit RF fields */
        RF2D_gc_stix stix;      /** 2D drift-kinetic Stix RF fields */
    };
} RF_fields;

a5err RF_fields_init(RF_fields* rf, hid_t f, char* qid,
                    int lhigh, B_field_data* bdata);
void RF_fields_free(RF_fields* rf);
void RF_fields_offload(RF_fields* rf);

GPU_DECLARE_TARGET_SIMD_UNIFORM(rf)
a5err RF_fields_eval(real E[3], real B[3], real r, real phi, real z, real t, RF_fields* rf);
DECLARE_TARGET_END

#endif