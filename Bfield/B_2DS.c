/**
 * @file B_2DS.c
 * @brief 2D magnetic field with bicubic spline interpolation
 *
 * This module represents a magnetic field where data is given in \f$Rz\f$-grid
 * from which it is interpolated with bicubic splines. The field is
 * axisymmetric.
 *
 * The magnetic field is evaluated from magnetic field strength \f$\mathbf{B}\f$
 * which may not be divergence free. However, \f$B_R\f$ and \f$B_z\f$ components
 * are also evaluated from poloidal magnetic flux \f$\psi\f$ as
 *
 * \f{align*}{
 * B_R &= -\frac{1}{R}\frac{\partial\psi}{\partial z}\\
 * B_z &= \frac{1}{R}\frac{\partial\psi}{\partial R}
 * \f}
 *
 * The total field is then a sum of components interpolated directly from
 * \f$\mathbf{B}\f$ and components calculated via interpolated \f$\psi\f$.
 *
 * This module does no extrapolation so if queried value is outside the
 * \f$Rz\f$-grid an error is thrown.
 *
 * The splines may either have compact or explicit forms which is toggled by
 * INTERP_SPL_EXPL in ascot5.h. Compact forms require 1/4 th of memory (in 2D)
 * but require more floating point operations.
 *
 * @see B_field.c interp2Dexpl.c interp2Dcomp.c
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../math.h"
#include "../ascot5.h"
#include "../error.h"
#include "../print.h"
#include "B_2DS.h"
#include "../spline/interp2D.h"
#include "../spline/interp2Dcomp.h"
#include "../spline/interp2Dexpl.h"

/**
 * @brief Initialize magnetic field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct and to construct
 * bicubic splines whose coefficients are stored in re-allocated offload array.
 *
 * The offload data struct must have the following fields initialized:
 * - B_2DS_offload_data.n_r
 * - B_2DS_offload_data.n_z
 * - B_2DS_offload_data.r_min
 * - B_2DS_offload_data.r_max
 * - B_2DS_offload_data.z_min
 * - B_2DS_offload_data.z_max
 * - B_2DS_offload_data.psi0
 * - B_2DS_offload_data.psi1
 * - B_2DS_offload_data.axis_r
 * - B_2DS_offload_data.axis_z
 *
 * B_2DS_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[            j*n_r + i] = psi(R_i, z_j)   [V*s*m^-1]
 * - offload_array[  n_r*n_z + j*n_r + i] = B_R(R_i, z_j)   [T]
 * - offload_array[2*n_r*n_z + j*n_r + i] = B_phi(R_i, z_j) [T]
 * - offload_array[3*n_r*n_z + j*n_r + i] = B_z(R_i, z_j)   [T]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int B_2DS_init_offload(B_2DS_offload_data* offload_data,
                       real** offload_array) {

    /* Fill rest of the offload data struct */
    offload_data->r_grid = (offload_data->r_max - offload_data->r_min)
        / (offload_data->n_r - 1);
    offload_data->z_grid = (offload_data->z_max - offload_data->z_min)
        / (offload_data->n_z - 1);

    /* Spline initialization. Use spline structs for temporary storage */
    int err = 0;
    int splinesize = offload_data->n_r*offload_data->n_z;
    interp2D_data psi;
    interp2D_data B_r;
    interp2D_data B_phi;
    interp2D_data B_z;
#if INTERP_SPL_EXPL
    splinesize *= 16;
    err += interp2Dexpl_init(
        &psi, *offload_array,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp2Dexpl_init(
        &B_r, *offload_array+offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp2Dexpl_init(
        &B_phi, *offload_array+2*offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp2Dexpl_init(
        &B_z, *offload_array+3*offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);
#else
    splinesize *= 4;
    err += interp2Dcomp_init(
        &psi, *offload_array,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp2Dcomp_init(
        &B_r, *offload_array+offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp2Dcomp_init(
        &B_phi, *offload_array+2*offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);

    err += interp2Dcomp_init(
        &B_z, *offload_array+3*offload_data->n_z*offload_data->n_r,
        offload_data->n_r, offload_data->n_z,
        offload_data->r_min, offload_data->r_max, offload_data->r_grid,
        offload_data->z_min, offload_data->z_max, offload_data->z_grid);
#endif
    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);

    offload_data->offload_array_length = splinesize*4;
    *offload_array =
        (real*) malloc( offload_data->offload_array_length * sizeof(real) );

    for(int i = 0; i < splinesize; i++) {
        (*offload_array)[0*splinesize + i] = psi.c[i];
        (*offload_array)[1*splinesize + i] = B_r.c[i];
        (*offload_array)[2*splinesize + i] = B_phi.c[i];
        (*offload_array)[3*splinesize + i] = B_z.c[i];
    }

    interp2D_free(&psi);
    interp2D_free(&B_r);
    interp2D_free(&B_phi);
    interp2D_free(&B_z);

    /* Evaluate psi and magnetic field on axis for checks */
    B_2DS_data Bdata;
    B_2DS_init(&Bdata, offload_data, *offload_array);
    real psival[1], Bval[3];
    err = B_2DS_eval_psi(psival, offload_data->axis_r, 0, offload_data->axis_z,
                         &Bdata);
    err = B_2DS_eval_B(Bval, offload_data->axis_r, 0, offload_data->axis_z,
                       &Bdata);
    if(err) {
        print_err("Error: Initialization failed.\n");
        return err;
    }

    /* Print some sanity check on data */
    printf("\n2D magnetic field (B_2DS)\n");
    print_out(VERBOSE_IO, "Grid: nR = %4.d Rmin = %3.3f Rmax = %3.3f\n",
              offload_data->n_r,
              offload_data->r_min, offload_data->r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f zmax = %3.3f\n",
              offload_data->n_z,
              offload_data->z_min, offload_data->z_max);
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
void B_2DS_free_offload(B_2DS_offload_data* offload_data,
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
void B_2DS_init(B_2DS_data* Bdata, B_2DS_offload_data* offload_data,
                real* offload_array) {

    int splinesize = offload_data->n_r * offload_data->n_z;
#if INTERP_SPL_EXPL
    splinesize *= 16;
#else
    splinesize *= 4;
#endif

    /* Initialize target data struct */
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;
    Bdata->axis_r = offload_data->axis_r;
    Bdata->axis_z = offload_data->axis_z;

    /* Copy parameters and assign pointers to offload array to initialize the
       spline structs */
    Bdata->psi.n_r    = offload_data->n_r;
    Bdata->psi.r_min  = offload_data->r_min;
    Bdata->psi.r_max  = offload_data->r_max;
    Bdata->psi.r_grid = offload_data->r_grid;
    Bdata->psi.n_z    = offload_data->n_z;
    Bdata->psi.z_min  = offload_data->z_min;
    Bdata->psi.z_max  = offload_data->z_max;
    Bdata->psi.z_grid = offload_data->z_grid;
    Bdata->psi.c      = &(offload_array[0*splinesize]);

    Bdata->B_r.n_r    = offload_data->n_r;
    Bdata->B_r.r_min  = offload_data->r_min;
    Bdata->B_r.r_max  = offload_data->r_max;
    Bdata->B_r.r_grid = offload_data->r_grid;
    Bdata->B_r.n_z    = offload_data->n_z;
    Bdata->B_r.z_min  = offload_data->z_min;
    Bdata->B_r.z_max  = offload_data->z_max;
    Bdata->B_r.z_grid = offload_data->z_grid;
    Bdata->B_r.c      = &(offload_array[1*splinesize]);

    Bdata->B_phi.n_r    = offload_data->n_r;
    Bdata->B_phi.r_min  = offload_data->r_min;
    Bdata->B_phi.r_max  = offload_data->r_max;
    Bdata->B_phi.r_grid = offload_data->r_grid;
    Bdata->B_phi.n_z    = offload_data->n_z;
    Bdata->B_phi.z_min  = offload_data->z_min;
    Bdata->B_phi.z_max  = offload_data->z_max;
    Bdata->B_phi.z_grid = offload_data->z_grid;
    Bdata->B_phi.c      = &(offload_array[2*splinesize]);

    Bdata->B_z.n_r    = offload_data->n_r;
    Bdata->B_z.r_min  = offload_data->r_min;
    Bdata->B_z.r_max  = offload_data->r_max;
    Bdata->B_z.r_grid = offload_data->r_grid;
    Bdata->B_z.n_z    = offload_data->n_z;
    Bdata->B_z.z_min  = offload_data->z_min;
    Bdata->B_z.z_max  = offload_data->z_max;
    Bdata->B_z.z_grid = offload_data->z_grid;
    Bdata->B_z.c      = &(offload_array[3*splinesize]);
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
a5err B_2DS_eval_psi(real psi[1], real r, real phi, real z, B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_B(&psi[0], &Bdata->psi, r, z);
#else
    interperr += interp2Dcomp_eval_B(&psi[0], &Bdata->psi, r, z);
#endif

    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );}

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
a5err B_2DS_eval_psi_dpsi(real psi_dpsi[4], real r, real phi, real z,
                          B_2DS_data* Bdata) {
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

    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );}

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
a5err B_2DS_eval_rho(real rho[1], real psi, B_2DS_data* Bdata) {

    /* Check that the values seem valid */
    real delta = (Bdata->psi1 - Bdata->psi0);
    if( (psi - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
    }

    rho[0] = sqrt( (psi - Bdata->psi0) / delta );
    return 0;
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
a5err B_2DS_eval_rho_drho(real rho_drho[4], real r, real phi, real z,
                          B_2DS_data* Bdata) {
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi[6];

#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#else
    interperr += interp2Dcomp_eval_dB(psi_dpsi, &Bdata->psi, r, z);
#endif

    if(interperr) {
        return error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );
    }

    /* Check that the values seem valid */
    real delta = Bdata->psi1 - Bdata->psi0;
    if( (rho_drho[0] - Bdata->psi0) / delta < 0 ) {
         return error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );
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
a5err B_2DS_eval_B(real B[3], real r, real phi, real z, B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_B(&B[0], &Bdata->B_r, r, z);
    interperr += interp2Dexpl_eval_B(&B[1], &Bdata->B_phi, r, z);
    interperr += interp2Dexpl_eval_B(&B[2], &Bdata->B_z, r, z);
#else
    interperr += interp2Dcomp_eval_B(&B[0], &Bdata->B_r, r, z);
    interperr += interp2Dcomp_eval_B(&B[1], &Bdata->B_phi, r, z);
    interperr += interp2Dcomp_eval_B(&B[2], &Bdata->B_z, r, z);
#endif

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );}

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
        if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );}
    }
#endif

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );}

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
a5err B_2DS_eval_B_dB(real B_dB[12], real r, real phi, real z,
                      B_2DS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[6];

#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(B_dB_temp, &Bdata->B_r, r, z);
#else
    interperr += interp2Dcomp_eval_dB(B_dB_temp, &Bdata->B_r, r, z);
#endif

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = 0;
    B_dB[3] = B_dB_temp[2];

#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(B_dB_temp, &Bdata->B_phi, r, z);
#else
    interperr += interp2Dcomp_eval_dB(B_dB_temp, &Bdata->B_phi, r, z);
#endif

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = 0;
    B_dB[7] = B_dB_temp[2];

#if INTERP_SPL_EXPL
    interperr += interp2Dexpl_eval_dB(B_dB_temp, &Bdata->B_z, r, z);
#else
    interperr += interp2Dcomp_eval_dB(B_dB_temp, &Bdata->B_z, r, z);
#endif

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = 0;
    B_dB[11] = B_dB_temp[2];

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );}


#ifndef NOPSI
    real psi_dpsi[6];

    if(!err) {
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
        if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_2DS );}
    }
#endif

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_2DS );}

    return err;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis R-coordinate [m]
 */
real B_2DS_get_axis_r(B_2DS_data* Bdata) {
    return Bdata->axis_r;
}

/**
 * @brief Return magnetic axis z-coordinate
 *
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Magnetic axis z-coordinate [m]
 */
real B_2DS_get_axis_z(B_2DS_data* Bdata) {
    return Bdata->axis_z;
}
