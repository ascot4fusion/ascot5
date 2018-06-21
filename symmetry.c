/**
 * @file symmetry.c
 * @brief Symmetry interface
 */
#include "ascot5.h"
#include "error.h"
#include "symmetry.h"

a5err symmetry_stellarator_apply_scalar(real r_phi_z[], real r0, real phi0, real z0,
                                real period_length) {
    a5err err = 0;

    r_phi_z[0] = r0;
    
    /* We need to use the stellarator symmetry here. */
    /* http://dx.doi.org/10.1016/S0167-2789(97)00216-9 */
    
    r_phi_z[1] = fmod(phi0, period_length);
    if ( 2*r_phi_z[1] > period_length) {
        r_phi_z[1] = period_length - r_phi_z[1];
        r_phi_z[2] = -z0;
    } else {
        r_phi_z[2] = z0;
    }

    return err;
}

a5err symmetry_periodic_apply_scalar(real r_phi_z[], real r0, real phi0, real z0,
                                     real period_length) {
    a5err err = 0;
    r_phi_z[0] = r0;
    r_phi_z[1] = fmod(phi0, period_length);
    if(r_phi_z[1] < 0) {
        r_phi_z[1] += period_length;
    }
    r_phi_z[2] = z0;
    return err;
}

a5err symmetry_stellarator_apply_vector(real r_phi_z[], real scaling[], real r0, real phi0, real z0,
                                real period_length) {
    a5err err = 0;

    r_phi_z[0] = r0;
    
    /* We need to use the stellarator symmetry here. */
    /* http://dx.doi.org/10.1016/S0167-2789(97)00216-9 */

    r_phi_z[1] = fmod(phi0, period_length);
    if ( 2*r_phi_z[1] > period_length) {
        r_phi_z[1] = period_length - r_phi_z[1];
        r_phi_z[2] = -z0;
        scaling[0] = -1;
        scaling[1] = 1;
        scaling[2] = 1;
    } else {
        r_phi_z[2] = z0;
        scaling[0] = 1;
        scaling[1] = 1;
        scaling[2] = 1;
    }

    return err;
}

a5err symmetry_periodic_apply_vector(real r_phi_z[], real scaling[], real r0, real phi0, real z0,
                                     real period_length) {
    a5err err = 0;

    r_phi_z[0] = r0;
    r_phi_z[1] = fmod(phi0, period_length);
    if(r_phi_z[1] < 0) {
        r_phi_z[1] += period_length;
    }
    r_phi_z[2] = z0;
    scaling[0] = 1;
    scaling[1] = 1;
    scaling[2] = 1;

    return err;
}

/**
 * @brief Symmetry transformation for scalar variables
 *
 * This function returns the symmetric (r, phi, z) coordinates for given (r0, phi0, z0)
 * so that f(r0,phi0,z0) = f(r,phi,z)
 *
 * @param r_phi_z transformed coordinates (r_phi_z[0] = r, r_phi_z[1] = phi, r_phi_z[2] = z)
 * @param r0 original value of r
 * @param phi0 original value of phi
 * @param z0 original value of z
 * @param type symmetry transformation type
 * @param period_length length of one symmetric period
 */
a5err symmetry_apply_scalar(real r_phi_z[], real r0, real phi0, real z0,
                            symmetry_type type, real period_length) {
    a5err err = 0;

    switch(type) {
    case symmetry_type_none:
	r_phi_z[0] = r0;
	r_phi_z[1] = phi0;
	r_phi_z[2] = z0;
        break;

    case symmetry_type_periodic:
        err = symmetry_periodic_apply_scalar(r_phi_z, r0, phi0, z0, period_length);
        break;
	
    case symmetry_type_stellarator:
        err = symmetry_stellarator_apply_scalar(r_phi_z, r0, phi0, z0, period_length);
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	r_phi_z[0] = r0;
	r_phi_z[1] = phi0;
	r_phi_z[2] = z0;
    }

    return err;
}

/**
 * @brief Symmetry transformation for vector valued variables
 *
 * This function returns the symmetric (r, phi, z) coordinates for given (r0, phi0, z0).
 * Also, the scaling coefficients for vector valued quantities are given, so that
 * (f_r,f_phi,f_z)_(r0,phi0,z0) = (scaling[0]*f_r,scaling[1]*f_phi,scaling[2]*f_z)_(r,phi,z)
 *
 * @param r_phi_z transformed coordinates (r_phi_z[0] = r, r_phi_z[1] = phi, r_phi_z[2] = z)
 * @param scaling scaling factor for each vector component
 * @param r0 original value of r
 * @param phi0 original value of phi
 * @param z0 original value of z
 * @param type symmetry transformation type
 * @param period_length length of one symmetric period
 */
a5err symmetry_apply_vector(real r_phi_z[], real scaling[], real r0, real phi0, real z0,
                            symmetry_type type, real period_length) {
    a5err err = 0;

    switch(type) {
    case symmetry_type_none:
	r_phi_z[0] = r0;
	r_phi_z[1] = phi0;
	r_phi_z[2] = z0;
	scaling[0] = 1;
	scaling[1] = 1;
	scaling[2] = 1;
        break;
        
    case symmetry_type_periodic:
        err = symmetry_periodic_apply_vector(r_phi_z, scaling, r0, phi0, z0, period_length);
        break;
	
    case symmetry_type_stellarator:
        err = symmetry_stellarator_apply_vector(r_phi_z, scaling, r0, phi0, z0, period_length);
        break;
    }

    if(err) {
	/* Return some reasonable values to avoid further errors */
	r_phi_z[0] = r0;
	r_phi_z[1] = phi0;
	r_phi_z[2] = z0;
	scaling[0] = 1;
	scaling[1] = 1;
	scaling[2] = 1;
    }

    return err;
}
