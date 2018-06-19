/**
 * @file symmetry.c
 * @brief Magnetic field interface
 */
#include "ascot5.h"
#include "error.h"
#include "symmetry.h"

a5err symmetry_stellarator_apply_scalar(real r_phi_z[], real r, real phi, real z,
                                real period_length) {
    a5err err = 0;

    r_phi_z[0] = r;
    
    /* We need to use the stellarator symmetry here. */
    /* http://dx.doi.org/10.1016/S0167-2789(97)00216-9 */
    /* The data is expected to include half a period. */
    
    /* First figure out the effective phi angle */
    r_phi_z[1] = fmod(phi, period_length);
    if ( 2*r_phi_z[1] > period_length) {
        r_phi_z[1] = period_length - r_phi_z[1];
        r_phi_z[2] = -z;
    } else {
        r_phi_z[2] = z;
    }

    return err;
}

a5err symmetry_periodic_apply_scalar(real r_phi_z[], real r, real phi, real z,
                                     real period_length) {
    a5err err = 0;
    r_phi_z[0] = r;
    r_phi_z[1] = fmod(phi, period_length);
    r_phi_z[2] = z;
    return err;
}

a5err symmetry_stellarator_apply_vector(real r_phi_z[], real scaling[], real r, real phi, real z,
                                real period_length) {
    a5err err = 0;

    r_phi_z[0] = r;
    
    /* We need to use the stellarator symmetry here. */
    /* http://dx.doi.org/10.1016/S0167-2789(97)00216-9 */
    /* The data is expected to include half a period. */
    
    /* First figure out the effective phi angle */
    r_phi_z[1] = fmod(phi, period_length);
    if ( 2*r_phi_z[1] > period_length) {
        r_phi_z[1] = period_length - r_phi_z[1];
        r_phi_z[2] = -z;
        scaling[0] = -1;
        scaling[1] = 1;
        scaling[2] = 1;
    } else {
        r_phi_z[2] = z;
        scaling[0] = 1;
        scaling[1] = 1;
        scaling[2] = 1;
    }

    return err;
}

a5err symmetry_periodic_apply_vector(real r_phi_z[], real scaling[], real r, real phi, real z,
                                     real period_length) {
    a5err err = 0;

    r_phi_z[0] = r;
    r_phi_z[1] = fmod(phi, period_length);
    r_phi_z[2] = z;
    scaling[0] = 1;
    scaling[1] = 1;
    scaling[2] = 1;

    return err;
}

a5err symmetry_apply_scalar(real r_phi_z[], real r, real phi, real z,
                            symmetry_type type, real period_length) {
    a5err err = 0;

    switch(type) {
        case symmetry_type_periodic:
        err = symmetry_periodic_apply_scalar(r_phi_z, r, phi, z, period_length);
        break;
	
        case symmetry_type_stellarator:
        err = symmetry_stellarator_apply_scalar(r_phi_z, r, phi, z, period_length);
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	r_phi_z[0] = r;
	r_phi_z[1] = phi;
	r_phi_z[2] = z;
    }

    return err;
}

a5err symmetry_apply_vector(real r_phi_z[], real scaling[], real r, real phi, real z,
                            symmetry_type type, real period_length) {
    a5err err = 0;

    switch(type) {
        case symmetry_type_periodic:
        err = symmetry_periodic_apply_vector(r_phi_z, scaling, r, phi, z, period_length);
        break;
	
        case symmetry_type_stellarator:
        err = symmetry_stellarator_apply_vector(r_phi_z, scaling, r, phi, z, period_length);
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	r_phi_z[0] = r;
	r_phi_z[1] = phi;
	r_phi_z[2] = z;
	scaling[0] = 1;
	scaling[1] = 1;
	scaling[2] = 1;
    }

    return err;
}
