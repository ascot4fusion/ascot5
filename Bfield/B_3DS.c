/**
 * @file B_3DS.c
 * @brief 3D magnetic field with tricubic spline interpolation
 *
 * This module represents a magnetic field where data is given in \f$R\phi z\f$-
 * grid from which it is interpolated with tricubic splines.
 *
 * The magnetic field is evaluated from magnetic field strength \f$\mathbf{B}\f$
 * which may not be divergence free. However, \f$B_R\f$ and \f$B_z\f$ components
 * are also evaluated from poloidal magnetic flux \f$\psi(R,z)\f$ as
 *
 * \f{align*}{
 * B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
 * B_z &= \frac{1}{R}\frac{\partial\psi}{\partial R}
 * \f}
 *
 * The total field is then a sum of components interpolated directly from
 * \f$\mathbf{B}\f$ and components calculated via interpolated \f$\psi\f$.
 * Note that \f$\psi\f$ is assumed to be axisymmetric and is interpolated with
 * bicubic splines. \f$\psi\f$ and \f$\mathbf{B}\f$ are given in separate grids.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown. For \f$\phi\f$-grid, periodic boundary
 * conditions are used but it is user's responsibility to provide input
 * whose \f$\phi\f$-grid makes sense (in that it actually represents a periodic
 * field), i.e., \f$\phi_\mathrm{max}-\phi_\mathrm{min} = 2\pi/(N+1)\f$.
 * However, do note that in this module \f$\phi_\mathrm{max}\f$ is not the
 * "last" grid point but the second last, e.g. if \f$\phi_\mathrm{min}=0\f$
 * and \f$n_\phi = 360\f$, then \f$\phi_\mathrm{max}=359\f$ if periodicity is
 * \f$N=0\f$.
 *
 * The splines may either have compact or explicit forms which is toggled by
 * INTERP_SPL_EXPL in ascot5.h. Compact forms require 1/8 th of memory (in 3D)
 * but require more floating point operations.
 *
 * @see B_field.c B_2DS.c interp3Dexpl.c interp3Dcomp.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "B_3DS.h"
#include "../spline/interp2D.h"
#include "../spline/interp3D.h"
#include "../spline/interp2Dcomp.h"
#include "../spline/interp3Dcomp.h"
#include "../spline/interp2Dexpl.h"
#include "../spline/interp3Dexpl.h"

/**
 * @brief Initialize magnetic field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct and to construct
 * bicubic splines whose coefficients are stored in re-allocated offload array.
 *
 * The offload data struct must have the following fields initialized:
 * - B_3DS_offload_data.psigrid_n_r
 * - B_3DS_offload_data.psigrid_n_z
 * - B_3DS_offload_data.psigrid_r_min
 * - B_3DS_offload_data.psigrid_r_max
 * - B_3DS_offload_data.psigrid_z_min
 * - B_3DS_offload_data.psigrid_z_max
 *
 * - B_3DS_offload_data.Bgrid_n_r
 * - B_3DS_offload_data.Bgrid_n_z
 * - B_3DS_offload_data.Bgrid_r_min
 * - B_3DS_offload_data.Bgrid_r_max
 * - B_3DS_offload_data.Bgrid_z_min
 * - B_3DS_offload_data.Bgrid_z_max
 * - B_3DS_offload_data.n_phi
 * - B_3DS_offload_data.phi_min
 * - B_3DS_offload_data.phi_max
 *
 * - B_3DS_offload_data.psi0
 * - B_3DS_offload_data.psi1
 * - B_3DS_offload_data.axis_r
 * - B_3DS_offload_data.axis_z
 *
 * B_3DS_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[                     z*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_R(R_i, phi_z, z_j)   [T]
 * - offload_array[  Bn_r*Bn_z*Bn_phi + z*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_phi(R_i, phi_z, z_j)   [T]
 * - offload_array[2*Bn_r*Bn_z*Bn_phi + z*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_z(R_i, phi_z, z_j)   [T]
 * - offload_array[3*Bn_r*Bn_z*Bn_phi + j*n_r + i]
 *   = psi(R_i, z_j)   [V*s*m^-1]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int B_3DS_init_offload(B_3DS_offload_data* offload_data, real** offload_array) {

    /* Fill rest of the offload data struct */
    offload_data->psigrid_r_grid =
        (offload_data->psigrid_r_max - offload_data->psigrid_r_min)
        / (offload_data->psigrid_n_r - 1);
    offload_data->psigrid_z_grid =
        (offload_data->psigrid_z_max - offload_data->psigrid_z_min)
        / (offload_data->psigrid_n_z - 1);

    offload_data->Bgrid_r_grid =
        (offload_data->Bgrid_r_max - offload_data->Bgrid_r_min)
        / (offload_data->Bgrid_n_r - 1);
    offload_data->Bgrid_z_grid =
        (offload_data->Bgrid_z_max - offload_data->Bgrid_z_min)
        / (offload_data->Bgrid_n_z - 1);
    offload_data->phi_grid =
        (offload_data->phi_max - offload_data->phi_min)
        / (offload_data->n_phi);// phi_max is not the last but second last point

    /* Spline initialization. Use spline structs for temporary storage */
    int err = 0;
    int psi_size = offload_data->psigrid_n_r * offload_data->psigrid_n_z;
    int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->n_phi;

    interp2D_data psi;
    interp3D_data B_r;
    interp3D_data B_phi;
    interp3D_data B_z;

#if INTERP_SPL_EXPL
    err += interp2Dexpl_init(
        &psi, *offload_array + 3*B_size,
        offload_data->psigrid_n_r, offload_data->psigrid_n_z,
        offload_data->psigrid_r_min, offload_data->psigrid_r_max,
        offload_data->psigrid_r_grid,
        offload_data->psigrid_z_min, offload_data->psigrid_z_max,
        offload_data->psigrid_z_grid);

    err += interp3Dexpl_init(
        &B_r, *offload_array + 0*B_size,
        offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    err += interp3Dexpl_init(
        &B_phi, *offload_array + 1*B_size,
        offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    err += interp3Dexpl_init(
        &B_z, *offload_array + 2*B_size,
        offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    /* The data is now presented with splines, each data point has
     * 64 spline coefficients in 3D and 16 in 2D */
    psi_size *= 16;
    B_size *= 64;

#else
    err += interp2Dcomp_init(
        &psi, *offload_array + 3*B_size,
        offload_data->psigrid_n_r, offload_data->psigrid_n_z,
        offload_data->psigrid_r_min, offload_data->psigrid_r_max,
        offload_data->psigrid_r_grid,
        offload_data->psigrid_z_min, offload_data->psigrid_z_max,
        offload_data->psigrid_z_grid);

    err += interp3Dcomp_init(
        &B_r, *offload_array + 0*B_size,
        offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    err += interp3Dcomp_init(
        &B_phi, *offload_array + 1*B_size,
        offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    err += interp3Dcomp_init(
        &B_z, *offload_array + 2*B_size,
        offload_data->Bgrid_n_r, offload_data->n_phi, offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->phi_min, offload_data->phi_max, offload_data->phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    /* The data is now presented with splines, each data point has
     * 16 spline coefficients in 3D and 4 in 2D */
    psi_size *= 4;
    B_size *= 8;
#endif
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);

    offload_data->offload_array_length = psi_size + B_size*3;
    *offload_array =
        (real*) malloc( offload_data->offload_array_length * sizeof(real) );

    for(int i = 0; i < B_size; i++) {
        (*offload_array)[0*B_size + i] = B_r.c[i];
        (*offload_array)[1*B_size + i] = B_phi.c[i];
        (*offload_array)[2*B_size + i] = B_z.c[i];
    }
    for(int i = 0; i < psi_size; i++) {
        (*offload_array)[3*B_size + i] = psi.c[i];
    }

    interp2D_free(&psi);
    interp3D_free(&B_r);
    interp3D_free(&B_phi);
    interp3D_free(&B_z);

    /* Evaluate psi and magnetic field on axis for checks */
    B_3DS_data Bdata;
    B_3DS_init(&Bdata, offload_data, *offload_array);
    real psival[1], Bval[3];
    err = B_3DS_eval_psi(psival, offload_data->axis_r, 0, offload_data->axis_z,
                         &Bdata);
    err = B_3DS_eval_B(Bval, offload_data->axis_r, 0, offload_data->axis_z,
                       &Bdata);
    if(err) {
        print_err("Error: Initialization failed.\n");
        return err;
    }

    /* Print some sanity check on data */
    printf("\n3D magnetic field (B_3DS)\n");
    print_out(VERBOSE_IO, "Psi-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n",
              offload_data->psigrid_n_r,
              offload_data->psigrid_r_min, offload_data->psigrid_r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f m zmax = %3.3f m\n",
              offload_data->psigrid_n_z,
              offload_data->psigrid_z_min, offload_data->psigrid_z_max);
    print_out(VERBOSE_IO, "B-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n",
              offload_data->Bgrid_n_r,
              offload_data->Bgrid_r_min, offload_data->Bgrid_r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f m zmax = %3.3f m\n",
              offload_data->Bgrid_n_z,
              offload_data->Bgrid_z_min, offload_data->Bgrid_z_max);
    print_out(VERBOSE_IO, "nphi = %4.d phimin = %3.3f deg phimax = %3.3f deg\n",
              offload_data->n_phi,
              math_rad2deg(offload_data->phi_min),
              math_rad2deg(offload_data->phi_max));
    print_out(VERBOSE_IO, "Psi at magnetic axis (%1.3f m, %1.3f m)\n",
              offload_data->axis_r, offload_data->axis_z);
    print_out(VERBOSE_IO, "%3.3f (evaluated)\n%3.3f (given)\n",
              psival[0], offload_data->psi0);
    print_out(VERBOSE_IO, "Magnetic field on axis:\n"
              "B_R = %3.3f B_phi = %3.3f B_z = %3.3f\n",
              Bval[0], Bval[1], Bval[2]);

    return err;
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_3DS_free_offload(B_3DS_offload_data* offload_data,
                        real** offload_array) {
    free(*offload_array);
    *offload_array = NULL;
}

/**
 * @brief Initialize magnetic field data struct on target
 *
 * @param Bdata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array offload array
 */
void B_3DS_init(B_3DS_data* Bdata, B_3DS_offload_data* offload_data,
                real* offload_array) {

    int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->n_phi ;
#if INTERP_SPL_EXPL
    B_size   *= 64;
#else
    B_size   *= 8;
#endif

    /* Initialize target data struct */
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;

    /* Copy parameters and assign pointers to offload array to initialize the
       spline structs */
    Bdata->B_r.n_r        = offload_data->Bgrid_n_r;
    Bdata->B_r.r_min      = offload_data->Bgrid_r_min;
    Bdata->B_r.r_max      = offload_data->Bgrid_r_max;
    Bdata->B_r.r_grid     = offload_data->Bgrid_r_grid;
    Bdata->B_r.n_z        = offload_data->Bgrid_n_z;
    Bdata->B_r.z_min      = offload_data->Bgrid_z_min;
    Bdata->B_r.z_max      = offload_data->Bgrid_z_max;
    Bdata->B_r.z_grid     = offload_data->Bgrid_z_grid;
    Bdata->B_r.n_phi      = offload_data->n_phi;
    Bdata->B_r.phi_min    = offload_data->phi_min;
    Bdata->B_r.phi_max    = offload_data->phi_max;
    Bdata->B_r.phi_grid   = offload_data->phi_grid;
    Bdata->B_r.c          = &(offload_array[0*B_size]);

    Bdata->B_phi.n_r      = offload_data->Bgrid_n_r;
    Bdata->B_phi.r_min    = offload_data->Bgrid_r_min;
    Bdata->B_phi.r_max    = offload_data->Bgrid_r_max;
    Bdata->B_phi.r_grid   = offload_data->Bgrid_r_grid;
    Bdata->B_phi.n_z      = offload_data->Bgrid_n_z;
    Bdata->B_phi.z_min    = offload_data->Bgrid_z_min;
    Bdata->B_phi.z_max    = offload_data->Bgrid_z_max;
    Bdata->B_phi.z_grid   = offload_data->Bgrid_z_grid;
    Bdata->B_phi.n_phi    = offload_data->n_phi;
    Bdata->B_phi.phi_min  = offload_data->phi_min;
    Bdata->B_phi.phi_max  = offload_data->phi_max;
    Bdata->B_phi.phi_grid = offload_data->phi_grid;
    Bdata->B_phi.c        = &(offload_array[1*B_size]);

    Bdata->B_z.n_r        = offload_data->Bgrid_n_r;
    Bdata->B_z.r_min      = offload_data->Bgrid_r_min;
    Bdata->B_z.r_max      = offload_data->Bgrid_r_max;
    Bdata->B_z.r_grid     = offload_data->Bgrid_r_grid;
    Bdata->B_z.n_z        = offload_data->Bgrid_n_z;
    Bdata->B_z.z_min      = offload_data->Bgrid_z_min;
    Bdata->B_z.z_max      = offload_data->Bgrid_z_max;
    Bdata->B_z.z_grid     = offload_data->Bgrid_z_grid;
    Bdata->B_z.n_phi      = offload_data->n_phi;
    Bdata->B_z.phi_min    = offload_data->phi_min;
    Bdata->B_z.phi_max    = offload_data->phi_max;
    Bdata->B_z.phi_grid   = offload_data->phi_grid;
    Bdata->B_z.c          = &(offload_array[2*B_size]);

    Bdata->psi.n_r        = offload_data->psigrid_n_r;
    Bdata->psi.r_min      = offload_data->psigrid_r_min;
    Bdata->psi.r_max      = offload_data->psigrid_r_max;
    Bdata->psi.r_grid     = offload_data->psigrid_r_grid;
    Bdata->psi.n_z        = offload_data->psigrid_n_z;
    Bdata->psi.z_min      = offload_data->psigrid_z_min;
    Bdata->psi.z_max      = offload_data->psigrid_z_max;
    Bdata->psi.z_grid     = offload_data->psigrid_z_grid;
    Bdata->psi.c          = &(offload_array[3*B_size]);
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
a5err B_3DS_eval_psi(real psi[1], real r, real phi, real z,
                   B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_B(&psi[0], &Bdata->psi, r, z);
#else
    interperr += interp2Dcomp_eval_B(&psi[0], &Bdata->psi, r, z);
#endif

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate poloidal flux psi and its derivatives
 *
 * @param psi pointer where psi [V*s*m^-1] and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                   B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp[6];
#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(psi_dpsi_temp, &Bdata->psi, r, z);
#else
    interperr += interp2Dcomp_eval_dB(psi_dpsi_temp, &Bdata->psi, r, z);
#endif
    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = 0;
    psi_dpsi[3] = psi_dpsi_temp[2];

    if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho
 *
 * @param rho pointer where rho value will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_rho(real rho[1], real psi, B_3DS_data* Bdata) {
    a5err err = 0;

    /* Check that the values seem valid */
    real delta = (Bdata->psi1 - Bdata->psi0);
    if( (psi - Bdata->psi0) / delta < 0 ) {
         err = error_raise( ERR_UNPHYSICAL_PSI, __LINE__ );
    }

    if(!err) {
        rho[0] = sqrt( (psi - Bdata->psi0) / delta );
    }
    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho and its derivatives
 *
 * @param rho pointer where rho and its derivatives will be stored
 * @param r R coordinate [m]
 * @param phi phi coordinate [rad]
 * @param z z coordinate [m]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_3DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_3DS_data* Bdata) {
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi[6];

#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#else
    interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#endif

    if(interperr) {
        return error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (rho_drho[0] - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_UNPHYSICAL_PSI, __LINE__ );
    }

    /* Normalize psi to get rho */
    rho_drho[0] = sqrt(fabs((psi_dpsi[0] - Bdata->psi0) / delta));

    rho_drho[1] = psi_dpsi[1] / (2*delta*rho_drho[0]);
    rho_drho[2] = 0;
    rho_drho[3] = psi_dpsi[2] / (2*delta*rho_drho[0]);

    return 0;
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
a5err B_3DS_eval_B(real B[3], real r, real phi, real z, B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
#if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_B(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dexpl_eval_B(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dexpl_eval_B(&B[2], &Bdata->B_z, r, phi, z);
#else
    interperr += interp3Dcomp_eval_B(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[2], &Bdata->B_z, r, phi, z);
#endif

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

#ifndef NOPSI
    if(!err) {
        real psi_dpsi[6];
#if INTERP_SPL_EXPL
        interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#else
        interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#endif
        B[0] = B[0] - psi_dpsi[2]/r;
        B[2] = B[2] + psi_dpsi[1]/r;

        /* Test for psi interpolation error */
        if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}
    }
#endif

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

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
a5err B_3DS_eval_B_dB(real B_dB[], real r, real phi, real z,
                      B_3DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];
#if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->B_r, r, phi, z);
#else
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_r, r, phi, z);
#endif

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];


#if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->B_phi, r, phi, z);
#else
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_phi, r, phi, z);
#endif

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];


#if INTERP_SPL_EXPL
    interperr += interp3Dexpl_eval_dB(B_dB_temp, &Bdata->B_z, r, phi, z);
#else
    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_z, r, phi, z);
#endif

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_OUTSIDE_BFIELD, __LINE__ );}

#ifndef NOPSI
    if(!err) {
        real psi_dpsi[6];
#if INTERP_SPL_EXPL
        interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#else
        interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#endif

        B_dB[0] = B_dB[0] - psi_dpsi[2]/r;
        B_dB[1] = B_dB[1] + psi_dpsi[2]/(r*r)-psi_dpsi[5]/r;
        B_dB[3] = B_dB[3] - psi_dpsi[4]/r;
        B_dB[8] = B_dB[8] + psi_dpsi[1]/r;
        B_dB[9] = B_dB[9] - psi_dpsi[1]/(r*r) + psi_dpsi[3]/r;
        B_dB[11] = B_dB[11] + psi_dpsi[5]/r;

        /* Test for psi interpolation error */
        if(interperr) {err = error_raise( ERR_OUTSIDE_PSIFIELD, __LINE__ );}
    }
#endif

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {err = error_raise( ERR_UNPHYSICAL_B, __LINE__ );}

    return err;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis R-coordinate [m]
 */
real B_3DS_get_axis_r(B_3DS_data* Bdata) {
    return Bdata->axis_r;
}

/**
 * @brief Return magnetic axis z-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis z-coordinate [m]
 */
real B_3DS_get_axis_z(B_3DS_data* Bdata) {
    return Bdata->axis_z;
}
