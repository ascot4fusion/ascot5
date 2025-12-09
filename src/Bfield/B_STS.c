/**
 * @file B_STS.c
 * @brief Stellarator magnetic field with cubic spline interpolation
 *
 * This module represents a magnetic field where data is given in \f$R\phi z\f$-
 * grid from which it is interpolated with tricubic splines.
 *
 * The magnetic field is evaluated from magnetic field strength \f$\mathbf{B}\f$
 * which may not be divergence free. The poloidal magnetic flux \f$\psi\f$ is
 * interpolated using tricubic splines as well. \f$\psi\f$ and \f$\mathbf{B}\f$
 * are given in separate grids.
 *
 * The magnetic axis location for stellarators varies with the \f$\phi\f$ angle
 * and is evaluated using linear interpolation.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown. The \f$\phi\f$-grid is assumed to be
 * periodic. Periodic boundary conditions are used but it is user's
 * responsibility to provide input whose \f$\phi\f$-grid makes sense (in that it
 * actually represents a periodic field), i.e.,
 * \f$\phi_\mathrm{max}-\phi_\mathrm{min} = 2\pi/(N+1)\f$.  However, do note
 * that in this module \f$\phi_\mathrm{max}\f$ is not the "last" grid point but
 * the second last, e.g. if \f$\phi_\mathrm{min}=0\f$ and \f$n_\phi = 360\f$,
 * then \f$\phi_\mathrm{max}=359\f$ if periodicity is \f$N=0\f$.
 *
 * @see B_field.c linint1D.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "../consts.h"
#include "B_STS.h"
#include "../linint/linint.h"
#include "../spline/interp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/**
 * @brief Initialize magnetic field data
 *
 * @param data pointer to the data struct
 * @param p_n_r number of r grid points in psi data
 * @param p_r_min minimum R coordinate in psi data grid [m]
 * @param p_r_max maximum R coordinate in psi data grid [m]
 * @param p_n_phi number of phi grid points in psi data
 * @param p_phi_min minimum phi coordinate in psi data grid [rad]
 * @param p_phi_max maximum phi coordinate in psi data grid [rad]
 * @param p_n_z number of z grid points in psi data
 * @param p_z_min minimum z coordinate in psi data grid [m]
 * @param p_z_max maximum z coordinate in psi data grid [m]
 * @param p_n_r number of r grid points in B data
 * @param b_r_min minimum R coordinate in B data grid [m]
 * @param b_r_max maximum R coordinate in B data grid [m]
 * @param b_n_phi number of phi grid points in B data
 * @param b_phi_min minimum phi coordinate in B data grid [rad]
 * @param b_phi_max maximum phi coordinate in B data grid [rad]
 * @param b_n_z number of z grid points in B data
 * @param b_z_min minimum z coordinate in B data grid [m]
 * @param b_z_max maximum z coordinate in B data grid [m]
 * @param naxis number of phi coordinates in magnetic axis data
 * @param axis_min ,inimum phi grid point in axis data [rad]
 * @param axis_max maximum phi grid point in axis data [rad]
 * @param axis_r R coordinates of the magnetic axis [m]
 * @param axis_z z coordinates of the magnetic axis [m]
 * @param psi0 poloidal flux at magnetic axis [Vs/m]
 * @param psi1 poloidal flux at separatrix [Vs/m]
 * @param psi poloidal flux
 *        psi(R_i,phi_j,z_k) = arr[k*p_n_r*p_n_phi + j*p_n_r + i] [Vs/m]
 * @param B_r Magnetic field R component
 *        B_r(R_i,phi_j,z_k) = arr[k*b_n_r*b_n_phi + j*b_n_r + i] [T]
 * @param B_phi Magnetic field phi component
 *        B_phi(R_i,phi_j,z_k) = arr[k*b_n_r*b_n_phi + j*b_n_r + i] [T]
 * @param B_z Magnetic field z component
 *        B_z(R_i,phi_j,z_k) = arr[k*b_n_r*b_n_phi + j*b_n_r + i] [T]
 * @param Nperiods Number of periods in the stellarator
 * @param stell_sym When set, it uses the stellarator symmetry
 *
 * @return zero if initialization succeeded
 */
int B_STS_init(B_STS_data* data,
               int p_n_r, real p_r_min, real p_r_max,
               int p_n_phi, real p_phi_min, real p_phi_max,
               int p_n_z, real p_z_min, real p_z_max,
               int b_n_r, real b_r_min, real b_r_max,
               int b_n_phi, real b_phi_min, real b_phi_max,
               int b_n_z, real b_z_min, real b_z_max,
               int naxis, real axis_min, real axis_max,
               real* axis_r, real* axis_z, real psi0, real psi1,
               real* psi, real* B_r, real* B_phi, real* B_z,
               int Nperiods, int stell_sym) {

    /* Spline initialization. */
    int err = 0;
    data->psi0 = psi0;
    data->psi1 = psi1;
    int PHI_INTERP_CONDITION = PERIODICBC;
    if(Nperiods > 0 && stell_sym){
        PHI_INTERP_CONDITION = NATURALBC;
    }

    err = interp3Dcomp_setup(&data->psi, psi, p_n_r, p_n_phi, p_n_z,
                             NATURALBC, PHI_INTERP_CONDITION, NATURALBC,
                             p_r_min, p_r_max, p_phi_min, p_phi_max,
                             p_z_min, p_z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }
    interp3Dcomp_setup(
        &data->B_r, B_r, b_n_r, b_n_phi, b_n_z,
        NATURALBC, PHI_INTERP_CONDITION, NATURALBC,
        b_r_min, b_r_max, b_phi_min, b_phi_max, b_z_min, b_z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }
    interp3Dcomp_setup(
        &data->B_phi, B_phi, b_n_r, b_n_phi, b_n_z,
        NATURALBC, PHI_INTERP_CONDITION, NATURALBC,
        b_r_min, b_r_max, b_phi_min, b_phi_max, b_z_min, b_z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }
    interp3Dcomp_setup(
        &data->B_z, B_z, b_n_r, b_n_phi, b_n_z,
        NATURALBC, PHI_INTERP_CONDITION, NATURALBC,
        b_r_min, b_r_max, b_phi_min, b_phi_max, b_z_min, b_z_max);
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return 1;
    }

    real* c1 = (real*)malloc(naxis*sizeof(real));
    real* c2 = (real*)malloc(naxis*sizeof(real));
    for(int i = 0; i < naxis; i++) {
        c1[i] = axis_r[i];
        c2[i] = axis_z[i];
    }
    linint1D_init(&data->axis_r, c1, naxis, PERIODICBC, axis_min, axis_max);
    linint1D_init(&data->axis_z, c2, naxis, PERIODICBC, axis_min, axis_max);

    /** Setting the symmetry properties */
    data->Nperiods = Nperiods;
    data->stell_sym = stell_sym;

    /* Evaluate psi and magnetic field on axis for checks */
    real psival[1], Bval[3], axis[2];
    err += B_STS_get_axis_rz(axis, data, 0);
    err += B_STS_eval_psi(psival, axis[0], 0, axis[1], data);
    err += B_STS_eval_B(Bval, axis[0], 0, axis[1], data);
    if(err) {
        print_err("Error: Initialization failed.\n");
        return err;
    }

    printf("\nStellarator magnetic field (B_STS)\n");
    print_out(VERBOSE_IO, "Psi-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n",
              p_n_r, p_r_min, p_r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f m zmax = %3.3f m\n",
              p_n_z, p_z_min, p_z_max);
    print_out(VERBOSE_IO, "nphi = %4.d phimin = %3.3f deg phimax = %3.3f deg\n",
              p_n_phi, math_rad2deg(p_phi_min), math_rad2deg(p_phi_max));
    print_out(VERBOSE_IO, "B-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n",
              b_n_r, b_r_min, b_r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f m zmax = %3.3f m\n",
              b_n_z, b_z_min, b_z_max);
    print_out(VERBOSE_IO, "nphi = %4.d phimin = %3.3f deg phimax = %3.3f deg\n",
              b_n_phi, math_rad2deg(b_phi_min), math_rad2deg(b_phi_max));
    print_out(VERBOSE_IO, "Psi at magnetic axis (phi=0) (%1.3f m, %1.3f m)\n",
              axis[0], axis[1]);
    print_out(VERBOSE_IO, "%3.3f (evaluated)\n%3.3f (given)\n",
              psival[0], data->psi0);
    print_out(VERBOSE_IO, "Magnetic field on axis:\n"
              "B_R = %3.3f B_phi = %3.3f B_z = %3.3f\n",
              Bval[0], Bval[1], Bval[2]);
    if(data->Nperiods > 0){
        print_out(VERBOSE_IO, "Stellarator with %d periods.\n",
                  data->Nperiods);
    }

    if(data->Nperiods > 0 && data->stell_sym){
        print_out(VERBOSE_IO,
                  "Stellarator symmetry enabled: assuming the grid is a half-period.\n");
    }

    return 0;
}

/**
 * @brief Free allocated resources
 *
 * @param data pointer to the data struct
 */
void B_STS_free(B_STS_data* data) {
    free(data->psi.c);
    free(data->B_r.c);
    free(data->B_phi.c);
    free(data->B_z.c);
}

/**
 * @brief Offload data to the accelerator.
 *
 * @param data pointer to the data struct
 */
void B_STS_offload(B_STS_data* data) {
    GPU_MAP_TO_DEVICE(
        data->axis_r, data->axis_r.c[0:data->axis_r.n_x], \
        data->axis_z, data->axis_z.c[0:data->axis_z.n_x], \
        data->psi, data->B_r, data->B_phi, data->B_z, \
        data->psi.c[0:data->psi.n_x*data->psi.n_y*data->psi.n_z*NSIZE_COMP3D], \
        data->B_r.c[0:data->B_r.n_x*data->B_r.n_y*data->B_r.n_z*NSIZE_COMP3D],\
        data->B_phi.c[0:data->B_phi.n_x*data->B_phi.n_y *data->B_phi.n_z*NSIZE_COMP3D], \
        data->B_z.c[0:data->B_z.n_x*data->B_z.n_y*data->B_z.n_z*NSIZE_COMP3D]
    )
}

/**
 * @brief Reduce coordinates to stellarator symmetry region
 * 
 * In the stellarator symmetry, the following symmetry holds:
 * \Psi(R, \phi, z) = \Psi(R, -\phi, -z)
 * B_R(R, \phi, z)  = - B_R(R, -\phi, -z)
 * B_\phi(R, \phi, z) =  B_\phi(R, -\phi, -z)
 * B_z(R, \phi, z)  = B_z(R, -\phi, -z)
 * 
 * This function will reduce the coordinates (R, phi, z) to the 
 * fundamental domain of the stellarator symmetry, i.e.,
 * phi in [0, pi/Nperiods] and z >= 0.
 *
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param r_out pointer to reduced R coordinate [m]
 * @param phi_out pointer to reduced phi coordinate [rad]
 * @param z_out pointer to reduced z coordinate [m]
 * @param flip pointer to sign flip factor for B_R and derivatives
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_reduce_symm(real r, real phi, real z,
                        real* r_out, real* phi_out, real* z_out, real* flip,
                        B_STS_data* Bdata){
    a5err err = 0;

    real phi_int = fmod(phi, 2.0*M_PI);

    if(Bdata->Nperiods <= 0 || Bdata->stell_sym == 0){
        // No stellarator symmetry to apply
        *r_out = r;
        *phi_out = phi_int;
        *z_out = z;
        *flip = 1.0;
        return err;
    }
    
    real half_period = (M_PI / Bdata->Nperiods);
    real phi_mod = fmod(phi_int, 2.0*half_period);

    if(phi_mod > half_period){
        // Map to negative phi
        *phi_out = 2.0*half_period - phi_mod;
        *z_out = -z;
        *flip = -1.0;
    } else {
        *phi_out = phi_mod;
        *z_out = z;
        *flip = 1.0;
    }

    *r_out = r;

    return err;
}

/**
 * @brief Evaluate poloidal flux psi
 *
 * @param psi pointer where psi [V*s*m^-1] value will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_eval_psi(real* psi, real r, real phi, real z,
                     B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    real r_red, phi_red, z_red, flip;
    err += B_STS_reduce_symm(r, phi, z, &r_red, &phi_red, &z_red, &flip, Bdata);

    interperr += interp3Dcomp_eval_f(&psi[0], &Bdata->psi, r_red, phi_red, z_red);

#ifdef B_STS_CLAMP_RHO_NONNEGATIVE
    if ( psi[0] < Bdata->psi0 ){
        psi[0] = Bdata->psi0;
    }
#endif

    /* Test for psi interpolation error */
    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * @param psi_dpsi pointer for storing psi [V*s*m^-1] and its derivatives
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp[10];

    real r_red, phi_red, z_red, flip;
    err += B_STS_reduce_symm(r, phi, z, &r_red, &phi_red, &z_red, &flip, Bdata);

    interperr += interp3Dcomp_eval_df(psi_dpsi_temp, &Bdata->psi, r_red, phi_red, z_red);

    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = flip * psi_dpsi_temp[2];
    psi_dpsi[3] = flip * psi_dpsi_temp[3];

#ifdef B_STS_CLAMP_RHO_NONNEGATIVE
    if ( psi_dpsi_temp[0] < Bdata->psi0 ){
        psi_dpsi[0] = Bdata->psi0;
        psi_dpsi[1] = 0.0;
        psi_dpsi[2] = 0.0;
        psi_dpsi[3] = 0.0;
    }
#endif


    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho and its derivatives
 *
 * @param rho_drho pointer where rho and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_STS_data* Bdata) {
    a5err err = 0;
    real psi_dpsi[4];

    real r_red, phi_red, z_red, flip;
    err += B_STS_reduce_symm(r, phi, z, &r_red, &phi_red, &z_red, &flip, Bdata);

    err = B_STS_eval_psi_dpsi(psi_dpsi, r_red, phi_red, z_red, Bdata);
    if(err){
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (psi_dpsi[0] - Bdata->psi0) / delta <= 0 ) {
#ifdef B_STS_CLAMP_RHO_NONNEGATIVE
        /* Set rho and the partial derivatives to zero, because otherwise one
           would need to divide by rho, which is zero. Of course, this problem
           persists when B_STS_CLAMP_RHO_NONNEGATIVE is not defined and
           the evaluation happens exactly at rho=0.0 */
        rho_drho[0] = 0.0;
        rho_drho[1] = 0.0;
        rho_drho[2] = 0.0;
        rho_drho[3] = 0.0;
        return err;
#else
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
#endif
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt( (psi_dpsi[0] - Bdata->psi0) / delta );

    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = flip * psi_dpsi[2] / (2*delta*rho_drho[0]);
    rho_drho[3] = flip * psi_dpsi[3] / (2*delta*rho_drho[0]);

    return err;
}

/**
 * @brief Evaluate magnetic field
 *
 * @param B pointer to array where magnetic field values are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_eval_B(real B[3], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    real r_red, phi_red, z_red, flip;
    err += B_STS_reduce_symm(r, phi, z, &r_red, &phi_red, &z_red, &flip, Bdata);

    interperr += interp3Dcomp_eval_f(&B[0], &Bdata->B_r, r_red, phi_red, z_red);
    interperr += interp3Dcomp_eval_f(&B[1], &Bdata->B_phi, r_red, phi_red, z_red);
    interperr += interp3Dcomp_eval_f(&B[2], &Bdata->B_z, r_red, phi_red, z_red);
    /* Test for B field interpolation error */
    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    B[0] = flip * B[0]; // B_R changes sign under stellarator symmetry

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    return err;

}

/**
 * @brief Evaluate magnetic field and its derivatives
 *
 * @param B_dB pointer to array where the field and its derivatives are stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [deg]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];

    real r_red, phi_red, z_red, flip;
    err += B_STS_reduce_symm(r, phi, z, &r_red, &phi_red, &z_red, &flip, Bdata);

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_r, r_red, phi_red, z_red);

    B_dB[0] = flip * B_dB_temp[0];
    B_dB[1] = flip * B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_phi, r_red, phi_red, z_red);

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = flip * B_dB_temp[2];
    B_dB[7] = flip * B_dB_temp[3];

    interperr += interp3Dcomp_eval_df(B_dB_temp, &Bdata->B_z, r_red, phi_red, z_red);

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = flip * B_dB_temp[2];
    B_dB[11] = flip * B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {
        return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );
    }

    return 0;
}

/**
 * @brief Return magnetic axis Rz-coordinates
 *
 * @param rz pointer where axis R and z [m] values will be stored
 * @param Bdata pointer to magnetic field data struct
 * @param phi phi coordinate [rad]
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_get_axis_rz(real rz[2], B_STS_data* Bdata, real phi) {
    a5err err = 0;

    int interperr = 0; /* If error happened during interpolation */
    interperr += linint1D_eval_f(&rz[0], &Bdata->axis_r, phi);
    interperr += linint1D_eval_f(&rz[1], &Bdata->axis_z, phi);
    if(interperr) {
        err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );
    }
    return err;
}
