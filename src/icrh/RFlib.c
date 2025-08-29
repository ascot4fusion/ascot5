#include <string.h>
#include "../hdf5io/hdf5_helpers.h"
#include "RFlib.h"
#include "RF2D_gc_stix.h"


a5err RF_fields_init(RF_fields* rf, hid_t f, char* qid,
                    int lhigh, B_field_data* bdata){
    if(!rf) return 0;
    char path[256]; // Storage array required for hdf5_gen_path() calls
    a5err err;

    rf->type = RF_NONE; // Setting the default.

    hdf5_gen_path("/rf/RF_2D_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        rf->type = RF_FULL_ORBIT_2D;
        err = RF2D_fields_init_from_file(&rf->rf2d, f, qid);
    }

    hdf5_gen_path("/rf/RF_3D_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        rf->type = RF_FULL_ORBIT_3D;
        err = RF3D_fields_init_from_file(&rf->rf3d, f, qid);
    }

    hdf5_gen_path("/rf/RF_STIX_2D_XXXXXXXXXX", qid, path);
    if( !hdf5_find_group(f, path) ) {
        rf->type = RF2D_GC_STIX;
        err = RF2D_gc_stix_init_from_file(&rf->stix, f, qid, lhigh, bdata);
    }

    return err;
}

void RF_fields_free(RF_fields* rf){
    if(!rf) return;

    if(rf->type == RF_FULL_ORBIT_2D){
        RF2D_fields_free(&rf->rf2d);
    } else if(rf->type == RF_FULL_ORBIT_3D){
        RF3D_fields_free(&rf->rf3d);
    } else if(rf->type == RF2D_GC_STIX){
        RF2D_gc_stix_free(&rf->stix);
    }
}

void RF_fields_offload(RF_fields* rf){
    if(!rf) return;

    if(rf->type == RF_FULL_ORBIT_2D){
        RF2D_fields_offload(&rf->rf2d);
    } else if(rf->type == RF_FULL_ORBIT_3D){
        RF3D_fields_offload(&rf->rf3d);
    } else if(rf->type == RF2D_GC_STIX){
        RF2D_gc_stix_offload(&rf->stix);
    }
}

a5err RF_fields_eval(real E[3], real B[3], real r, real phi,\
                       real z, real t, RF_fields* rf){
    if(!rf) return 0;
    
    switch(rf->type){
        case RF_FULL_ORBIT_2D:
            return RF2D_field_eval(E, B, r, phi, z, t, &rf->rf2d);
        case RF_FULL_ORBIT_3D:
            return RF3D_field_eval(E, B, r, phi, z, t, &rf->rf3d);
        case RF2D_GC_STIX:
            return -1;
        default:
            E[0] = 0.0;
            E[1] = 0.0;
            E[2] = 0.0;
            B[0] = 0.0;
            B[1] = 0.0;
            B[2] = 0.0;
            return 0;
    }
    return 0;
}