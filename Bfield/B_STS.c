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
 * \f$Rz\f$-grid an error is thrown. For \f$\phi\f$-grid, the input can be
 * either stellarator-symmetric or periodic. The stellarator symmetric field is
 * convertedto periodic during initialization. Periodic boundary conditions are
 * used but it is user's responsibility to provide input whose \f$\phi\f$-grid
 * makes sense (in that it actually represents a periodic field), i.e.,
 * \f$\phi_\mathrm{max}-\phi_\mathrm{min} = 2\pi/(N+1)\f$.
 * However, do note that in this module \f$\phi_\mathrm{max}\f$ is not the
 * "last" grid point but the second last, e.g. if \f$\phi_\mathrm{min}=0\f$
 * and \f$n_\phi = 360\f$, then \f$\phi_\mathrm{max}=359\f$ if periodicity is
 * \f$N=0\f$.
 *
 * @see B_field.c B_2DS.c B_3DS.c interp3Dexpl.c interp3Dcomp.c linint1D.c
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
#include "../linint/linint1D.h"
#include "../spline/interp3Dcomp.h"

/**
 * @brief Initialize magnetic field offload data
 *
 * This function takes pre-initialized offload data struct and offload array as
 * inputs. The data is used to fill rest of the offload struct and to construct
 * bicubic splines whose coefficients are stored in re-allocated offload array.
 *
 * The offload data struct must have the following fields initialized:
 * - B_STS_offload_data.psigrid_n_r
 * - B_STS_offload_data.psigrid_n_z
 * - B_STS_offload_data.psigrid_r_min
 * - B_STS_offload_data.psigrid_r_max
 * - B_STS_offload_data.psigrid_z_min
 * - B_STS_offload_data.psigrid_z_max
 *
 * - B_STS_offload_data.Bgrid_n_r
 * - B_STS_offload_data.Bgrid_n_z
 * - B_STS_offload_data.Bgrid_r_min
 * - B_STS_offload_data.Bgrid_r_max
 * - B_STS_offload_data.Bgrid_z_min
 * - B_STS_offload_data.Bgrid_z_max
 * - B_STS_offload_data.Bgrid_n_phi
 * - B_STS_offload_data.Bgrid_phi_min
 * - B_STS_offload_data.Bgrid_phi_max
 *
 * - B_STS_offload_data.n_axis
 * - B_STS_offload_data.axis_min
 * - B_STS_offload_data.axis_max
 *
 * - B_STS_offload_data.psi0
 * - B_STS_offload_data.psi1
 * - B_STS_offload_data.axis_r
 * - B_STS_offload_data.axis_z
 *
 * - B_STS_offload_data.n_periods
 * - B_STS_offload_data.symmetry_mode
 *
 * B_STS_offload_data.offload_array_length is set here.
 *
 * The offload array must contain the following data:
 * - offload_array[                     z*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_R(R_i, phi_z, z_j)   [T]
 * - offload_array[  Bn_r*Bn_z*Bn_phi + z*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_phi(R_i, phi_z, z_j)   [T]
 * - offload_array[2*Bn_r*Bn_z*Bn_phi + z*Bn_r*Bn_z + j*Bn_r + i]
 *   = B_z(R_i, phi_z, z_j)   [T]
 * - offload_array[3*Bn_r*Bn_z*Bn_phi + z*psin_r*psin_z + j*psin_r + i]
 *   = psi(R_i, phi_z, z_j)   [V*s*m^-1]
 * - offload_array[3*Bn_r*Bn_z*Bn_phi + psin_r*psin_z*psin_phi + z]
 *   = axis_R(phi_z)   [m]
 * - offload_array[3*Bn_r*Bn_z*Bn_phi + psin_r*psin_z*psin_phi + n_axis + z]
 *   = axis_z(phi_z)   [m]
 *
 * Sanity checks are printed if data was initialized succesfully.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array which is reallocated here
 *
 * @return zero if initialization succeeded
 */
int B_STS_init_offload(B_STS_offload_data* offload_data, real** offload_array) {

    /* Fill rest of the offload data struct */
    offload_data->psigrid_r_grid =
        (offload_data->psigrid_r_max - offload_data->psigrid_r_min)
        / (offload_data->psigrid_n_r - 1);
    offload_data->psigrid_z_grid =
        (offload_data->psigrid_z_max - offload_data->psigrid_z_min)
        / (offload_data->psigrid_n_z - 1);
    offload_data->psigrid_phi_grid =
        (offload_data->psigrid_phi_max - offload_data->psigrid_phi_min)
        / (offload_data->psigrid_n_phi - 1);// phi_max currently last point

    offload_data->Bgrid_r_grid =
        (offload_data->Bgrid_r_max - offload_data->Bgrid_r_min)
        / (offload_data->Bgrid_n_r - 1);
    offload_data->Bgrid_z_grid =
        (offload_data->Bgrid_z_max - offload_data->Bgrid_z_min)
        / (offload_data->Bgrid_n_z - 1);
    offload_data->Bgrid_phi_grid =
        (offload_data->Bgrid_phi_max - offload_data->Bgrid_phi_min)
        / (offload_data->Bgrid_n_phi - 1);// phi_max currently last point

    offload_data->axis_grid = (offload_data->axis_max - offload_data->axis_min)
        / (offload_data->n_axis - 1);

    offload_data->period_length = CONST_2PI/offload_data->n_periods;

    int offload_B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->Bgrid_n_phi;
    int offload_psi_size = offload_data->psigrid_n_r * offload_data->psigrid_n_z
        * offload_data->psigrid_n_phi;
    real* offload_B_r = &(*offload_array)[0*offload_B_size];
    real* offload_B_phi = &(*offload_array)[1*offload_B_size];
    real* offload_B_z = &(*offload_array)[2*offload_B_size];
    real* offload_psi = &(*offload_array)[3*offload_B_size];
    real* temp_array = NULL;

    int B_size = 0;
    int psi_size = 0;
    int axis_size = offload_data->n_axis;
    if (offload_data->symmetry_mode == symmetry_type_stellarator) {
        /* We need to use stellarator symmetry here.
         * http://dx.doi.org/10.1016/S0167-2789(97)00216-9
         * The data is expected to include half a period.
         * Stellarator symmetry: i_phi        <=> 2*(n_phi-1)-i_phi
         *                       i_z          <=> n_z-i_z-1
         * So, a point at:
         *     (i_r, i_phi, i_z)  <=> (i_r, 2*(n_phi-1)-i_phi, n_z-i_z-1)
         *
         * temp_B_x data is in the format:
         *     (i_r, i_phi, i_z) = temp_B_x(i_z*n_phi*n_r + i_phi*n_r + i_r)
         * offload_array data -"-"-"-"-  :
         *     (i_r, i_phi, i_z) =
         *                    (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ]
         * => (*offload_array)[i_phi*n_z*n_r + i_z*n_r + i_r ] =
         *                           temp_B_x(i_z*n_phi*n_r + i_phi*n_r + i_r);
         * The values are: Sym[B_r, B_phi, B_z] = [-B_r, B_phi, B_z]
         */
        B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
            * (2 * (offload_data->Bgrid_n_phi - 1));
        psi_size = offload_data->psigrid_n_r * offload_data->psigrid_n_z
            * (2 * (offload_data->psigrid_n_phi - 1));
        temp_array = (real*) malloc((3 * B_size + psi_size) * sizeof(real));

        int i_phi;
        int i_z;
        int i_r;
        int temp_ind, off_ind, sym_ind;

        /* Bfield */
        int n_r = offload_data->Bgrid_n_r;
        int n_phi = offload_data->Bgrid_n_phi;
        int n_z = offload_data->Bgrid_n_z;
        for (i_phi = 0; i_phi < n_phi; i_phi++) {
            for (i_z = 0; i_z < n_z; i_z++) {
                for (i_r = 0; i_r < n_r; i_r++) {
                    /* Index of data point in temp arrays */
                    temp_ind = i_z*n_phi*n_r + i_phi*n_r + i_r;

                    /* Index of data point in offload_array and
                       corresponding stel.-symmetric index  */
                    off_ind = i_phi*n_z*n_r + i_z*n_r + i_r;
                    sym_ind = (2*(n_phi - 1) - i_phi)*n_z*n_r
                        + (n_z - i_z - 1)*n_r + i_r;

                    /* B_r */
                    temp_array[off_ind] =
                        offload_B_r[temp_ind];
                    if (i_phi != 0) {
                        temp_array[sym_ind] =
                            - offload_B_r[temp_ind];
                    }
                    // B_phi
                    temp_array[B_size + off_ind] =
                        offload_B_phi[temp_ind];
                    if (i_phi != 0) {
                        temp_array[B_size + sym_ind] =
                            offload_B_phi[temp_ind];
                    }
                    // B_z
                    temp_array[2*B_size + off_ind] =
                        offload_B_z[temp_ind];
                    if (i_phi != 0) {
                        temp_array[2*B_size + sym_ind] =
                            offload_B_z[temp_ind];
                    }
                }
            }
        }
        /* Bfield data is now for one toroidal period */
        offload_data->Bgrid_n_phi = 2*(offload_data->Bgrid_n_phi - 1);
        offload_data->Bgrid_phi_max = 2*offload_data->Bgrid_phi_max;

        /* Psi */
        n_r = offload_data->psigrid_n_r;
        n_phi = offload_data->psigrid_n_phi;
        n_z = offload_data->psigrid_n_z;
        for (i_phi = 0; i_phi < n_phi; i_phi++) {
            for (i_z = 0; i_z < n_z; i_z++) {
                for (i_r = 0; i_r < n_r; i_r++) {
                    /* Index of data point in temp arrays */
                    temp_ind = i_z*n_phi*n_r + i_phi*n_r + i_r;

                    /* Index of data point in offload_array and
                       corresponding stel.-symmetric index  */
                    off_ind = i_phi*n_z*n_r + i_z*n_r + i_r;
                    sym_ind = (2*(n_phi - 1) - i_phi)*n_z*n_r
                        + (n_z - i_z - 1)*n_r + i_r;

                    // psi
                    temp_array[3*B_size + off_ind] =
                        offload_psi[temp_ind];
                    if (i_phi != 0) {
                        temp_array[3*B_size + sym_ind] =
                            offload_psi[temp_ind];
                    }
                }
            }
        }
        /* Psi data is now for one toroidal period */
        offload_data->psigrid_n_phi = 2*(offload_data->psigrid_n_phi - 1);
        offload_data->psigrid_phi_max = 2*offload_data->psigrid_phi_max;
        /* Since the data is now for a whole toroidal period, we can set the
           symmetry mode as pediodic */
        offload_data->symmetry_mode = symmetry_type_periodic;
    }
    else if (offload_data->symmetry_mode == symmetry_type_periodic) {
        /* Remove duplicate phi point from data */
        B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
            * (offload_data->Bgrid_n_phi - 1);
        psi_size = offload_data->psigrid_n_r * offload_data->psigrid_n_z
            * (offload_data->psigrid_n_phi - 1);
        temp_array = (real*) malloc((3 * B_size + psi_size)
                                        * sizeof(real));

        int i_phi;
        int i_z;
        int i_r;
        int temp_ind, off_ind;

        /* Bfield */
        int n_r = offload_data->Bgrid_n_r;
        int n_phi = offload_data->Bgrid_n_phi - 1;
        int n_z = offload_data->Bgrid_n_z;
        for (i_phi = 0; i_phi < n_phi - 1; i_phi++) {
            for (i_z = 0; i_z < n_z; i_z++) {
                for (i_r = 0; i_r < n_r; i_r++) {
                    /* Index of data point in temp arrays */
                    temp_ind = i_z*n_phi*n_r + i_phi*n_r + i_r;
                    /* Index of data point in offload_array */
                    off_ind = i_phi*n_z*n_r + i_z*n_r + i_r;
                    /* B_r */
                    temp_array[off_ind] =
                        offload_B_r[temp_ind];
                    // B_phi
                    temp_array[B_size + off_ind] =
                        offload_B_phi[temp_ind];
                    // B_z
                    temp_array[2*B_size + off_ind] =
                        offload_B_z[temp_ind];
                }
            }
        }
        /* Bfield data is now for one point less */
        offload_data->Bgrid_n_phi = offload_data->Bgrid_n_phi - 1;
        offload_data->Bgrid_phi_max =
            offload_data->Bgrid_phi_max - offload_data->Bgrid_phi_grid;

        /* Psi */
        n_r = offload_data->psigrid_n_r;
        n_phi = offload_data->psigrid_n_phi - 1;
        n_z = offload_data->psigrid_n_z;
        for (i_phi = 0; i_phi < n_phi; i_phi++) {
            for (i_z = 0; i_z < n_z; i_z++) {
                for (i_r = 0; i_r < n_r; i_r++) {
                    /* Index of data point in temp arrays */
                    temp_ind = i_z*n_phi*n_r + i_phi*n_r + i_r;
                    /* Index of data point in offload_array and
                       corresponding stel.-symmetric index  */
                    off_ind = i_phi*n_z*n_r + i_z*n_r + i_r;
                    // psi
                    temp_array[3*B_size + off_ind] =
                        offload_psi[temp_ind];
                }
            }
        }
        /* Psi data is now for one point less */
        offload_data->psigrid_n_phi = offload_data->psigrid_n_phi - 1;
        offload_data->psigrid_phi_max =
            offload_data->psigrid_phi_max - offload_data->psigrid_phi_grid;
    }

    /* Spline initialization. Use spline structs for temporary storage */
    int err = 0;

    interp3D_data psi;
    interp3D_data B_r;
    interp3D_data B_phi;
    interp3D_data B_z;
    linint1D_data axis_r;
    linint1D_data axis_z;

    err += interp3Dcomp_init(
        &psi, temp_array + 3*B_size,
        offload_data->psigrid_n_r, offload_data->psigrid_n_phi,
        offload_data->psigrid_n_z,
        offload_data->psigrid_r_min, offload_data->psigrid_r_max,
        offload_data->psigrid_r_grid,
        offload_data->psigrid_phi_min, offload_data->psigrid_phi_max,
        offload_data->psigrid_phi_grid,
        offload_data->psigrid_z_min, offload_data->psigrid_z_max,
        offload_data->psigrid_z_grid);

    err += interp3Dcomp_init(
        &B_r, temp_array + 0*B_size,
        offload_data->Bgrid_n_r, offload_data->Bgrid_n_phi,
        offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->Bgrid_phi_min, offload_data->Bgrid_phi_max,
        offload_data->Bgrid_phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    err += interp3Dcomp_init(
        &B_phi, temp_array + 1*B_size,
        offload_data->Bgrid_n_r, offload_data->Bgrid_n_phi,
        offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->Bgrid_phi_min, offload_data->Bgrid_phi_max,
        offload_data->Bgrid_phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    err += interp3Dcomp_init(
        &B_z, temp_array + 2*B_size,
        offload_data->Bgrid_n_r, offload_data->Bgrid_n_phi,
        offload_data->Bgrid_n_z,
        offload_data->Bgrid_r_min, offload_data->Bgrid_r_max,
        offload_data->Bgrid_r_grid,
        offload_data->Bgrid_phi_min, offload_data->Bgrid_phi_max,
        offload_data->Bgrid_phi_grid,
        offload_data->Bgrid_z_min, offload_data->Bgrid_z_max,
        offload_data->Bgrid_z_grid);

    /* Magnetic axis */
    int periodic = 1;
    err += linint1D_init(
        &axis_r,
        *offload_array + 3*offload_B_size + offload_psi_size,
        offload_data->n_axis, offload_data->axis_min, offload_data->axis_max,
        offload_data->axis_grid, periodic);

    err += linint1D_init(
        &axis_z,
        *offload_array + 3*offload_B_size + offload_psi_size + axis_size,
        offload_data->n_axis, offload_data->axis_min, offload_data->axis_max,
        offload_data->axis_grid, periodic);

    /* The data is now presented with splines, each data point has
     * 8 spline coefficients in 3D */
    psi_size *= 8;
    B_size *= 8;

    if(err) {
        print_err("Error: Failed to initialize splines.\n");
        return err;
    }

    /* Re-allocate the offload array and store spline coefficients there */
    free(*offload_array);

    offload_data->offload_array_length = psi_size + B_size*3 + 2*axis_size;
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
    for(int i = 0; i < axis_size; i++) {
        (*offload_array)[3*B_size + psi_size + i] = axis_r.c[i];
        (*offload_array)[3*B_size + psi_size + axis_size + i] = axis_z.c[i];
    }

    interp3D_free(&psi);
    interp3D_free(&B_r);
    interp3D_free(&B_phi);
    interp3D_free(&B_z);
    linint1D_free(&axis_r);
    linint1D_free(&axis_z);

    /* Evaluate psi and magnetic field on axis for checks */
    B_STS_data Bdata;
    B_STS_init(&Bdata, offload_data, *offload_array);
    real psival[1], Bval[3], axis[2];
    err = B_STS_get_axis_r(axis + 0, &Bdata, 0);
    if(err) {
        print_err("Error: Axis r evaluation failed.\n");
        return err;
    }
    err = B_STS_get_axis_z(axis + 1, &Bdata, 0);
    if(err) {
        print_err("Error: Axis z evaluation failed.\n");
        return err;
    }
    err = B_STS_eval_psi(psival, axis[0], 0, axis[1],
                         &Bdata);
    if(err) {
        print_err("Error: psi evaluation failed.\n");
        return err;
    }
    err = B_STS_eval_B(Bval, axis[0], 0, axis[1],
                       &Bdata);
    if(err) {
        print_err("Error: B evaluation failed.\n");
        return err;
    }

    printf("\nStellarator magnetic field (B_STS)\n");
    print_out(VERBOSE_IO, "Psi-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n",
              offload_data->psigrid_n_r,
              offload_data->psigrid_r_min, offload_data->psigrid_r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f m zmax = %3.3f m\n",
              offload_data->psigrid_n_z,
              offload_data->psigrid_z_min, offload_data->psigrid_z_max);
    print_out(VERBOSE_IO, "nphi = %4.d phimin = %3.3f deg phimax = %3.3f deg\n",
              offload_data->psigrid_n_phi,
              math_rad2deg(offload_data->psigrid_phi_min),
              math_rad2deg(offload_data->psigrid_phi_max));
    print_out(VERBOSE_IO, "B-grid: nR = %4.d Rmin = %3.3f m Rmax = %3.3f m\n",
              offload_data->Bgrid_n_r,
              offload_data->Bgrid_r_min, offload_data->Bgrid_r_max);
    print_out(VERBOSE_IO, "      nz = %4.d zmin = %3.3f m zmax = %3.3f m\n",
              offload_data->Bgrid_n_z,
              offload_data->Bgrid_z_min, offload_data->Bgrid_z_max);
    print_out(VERBOSE_IO, "nphi = %4.d phimin = %3.3f deg phimax = %3.3f deg\n",
              offload_data->Bgrid_n_phi,
              math_rad2deg(offload_data->Bgrid_phi_min),
              math_rad2deg(offload_data->Bgrid_phi_max));
    print_out(VERBOSE_IO, "Psi at magnetic axis (phi=0) (%1.3f m, %1.3f m)\n",
              axis[0], axis[1]);
    print_out(VERBOSE_IO, "%3.3f (evaluated)\n%3.3f (given)\n",
              psival[0], offload_data->psi0);
    print_out(VERBOSE_IO, "Magnetic field on axis:\n"
              "B_R = %3.3f B_phi = %3.3f B_z = %3.3f\n",
              Bval[0], Bval[1], Bval[2]);
    print_out(VERBOSE_IO, "period length = %1.3f rad, symmetry_mode = %4.d\n",
              offload_data->period_length, offload_data->symmetry_mode);

    return 0;
}

/**
 * @brief Free offload array
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 */
void B_STS_free_offload(B_STS_offload_data* offload_data, real** offload_array) {
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
void B_STS_init(B_STS_data* Bdata, B_STS_offload_data* offload_data,
               real* offload_array) {

    int B_size = offload_data->Bgrid_n_r * offload_data->Bgrid_n_z
        * offload_data->Bgrid_n_phi ;
    int psi_size = offload_data->psigrid_n_r * offload_data->psigrid_n_z
        * offload_data->psigrid_n_phi;
    int axis_size = offload_data->n_axis;

    B_size   *= 8;
    psi_size *= 8;

    /* Initialize target data struct */
    Bdata->period_length = offload_data->period_length;
    Bdata->symmetry_mode = offload_data->symmetry_mode;
    Bdata->psi0 = offload_data->psi0;
    Bdata->psi1 = offload_data->psi1;

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
    Bdata->B_r.n_phi      = offload_data->Bgrid_n_phi;
    Bdata->B_r.phi_min    = offload_data->Bgrid_phi_min;
    Bdata->B_r.phi_max    = offload_data->Bgrid_phi_max;
    Bdata->B_r.phi_grid   = offload_data->Bgrid_phi_grid;
    Bdata->B_r.c          = &(offload_array[0*B_size]);

    Bdata->B_phi.n_r      = offload_data->Bgrid_n_r;
    Bdata->B_phi.r_min    = offload_data->Bgrid_r_min;
    Bdata->B_phi.r_max    = offload_data->Bgrid_r_max;
    Bdata->B_phi.r_grid   = offload_data->Bgrid_r_grid;
    Bdata->B_phi.n_z      = offload_data->Bgrid_n_z;
    Bdata->B_phi.z_min    = offload_data->Bgrid_z_min;
    Bdata->B_phi.z_max    = offload_data->Bgrid_z_max;
    Bdata->B_phi.z_grid   = offload_data->Bgrid_z_grid;
    Bdata->B_phi.n_phi    = offload_data->Bgrid_n_phi;
    Bdata->B_phi.phi_min  = offload_data->Bgrid_phi_min;
    Bdata->B_phi.phi_max  = offload_data->Bgrid_phi_max;
    Bdata->B_phi.phi_grid = offload_data->Bgrid_phi_grid;
    Bdata->B_phi.c        = &(offload_array[1*B_size]);

    Bdata->B_z.n_r        = offload_data->Bgrid_n_r;
    Bdata->B_z.r_min      = offload_data->Bgrid_r_min;
    Bdata->B_z.r_max      = offload_data->Bgrid_r_max;
    Bdata->B_z.r_grid     = offload_data->Bgrid_r_grid;
    Bdata->B_z.n_z        = offload_data->Bgrid_n_z;
    Bdata->B_z.z_min      = offload_data->Bgrid_z_min;
    Bdata->B_z.z_max      = offload_data->Bgrid_z_max;
    Bdata->B_z.z_grid     = offload_data->Bgrid_z_grid;
    Bdata->B_z.n_phi      = offload_data->Bgrid_n_phi;
    Bdata->B_z.phi_min    = offload_data->Bgrid_phi_min;
    Bdata->B_z.phi_max    = offload_data->Bgrid_phi_max;
    Bdata->B_z.phi_grid   = offload_data->Bgrid_phi_grid;
    Bdata->B_z.c          = &(offload_array[2*B_size]);

    Bdata->psi.n_r        = offload_data->psigrid_n_r;
    Bdata->psi.r_min      = offload_data->psigrid_r_min;
    Bdata->psi.r_max      = offload_data->psigrid_r_max;
    Bdata->psi.r_grid     = offload_data->psigrid_r_grid;
    Bdata->psi.n_z        = offload_data->psigrid_n_z;
    Bdata->psi.z_min      = offload_data->psigrid_z_min;
    Bdata->psi.z_max      = offload_data->psigrid_z_max;
    Bdata->psi.z_grid     = offload_data->psigrid_z_grid;
    Bdata->psi.n_phi      = offload_data->psigrid_n_phi;
    Bdata->psi.phi_min    = offload_data->psigrid_phi_min;
    Bdata->psi.phi_max    = offload_data->psigrid_phi_max;
    Bdata->psi.phi_grid   = offload_data->psigrid_phi_grid;
    Bdata->psi.c          = &(offload_array[3*B_size]);

    Bdata->axis_r.n_r     = offload_data->n_axis;
    Bdata->axis_r.r_min   = offload_data->axis_min;
    Bdata->axis_r.r_max   = offload_data->axis_max;
    Bdata->axis_r.r_grid  = offload_data->axis_grid;
    Bdata->axis_r.c       = &(offload_array[3*B_size + psi_size]);

    Bdata->axis_z.n_r     = offload_data->n_axis;
    Bdata->axis_z.r_min   = offload_data->axis_min;
    Bdata->axis_z.r_max   = offload_data->axis_max;
    Bdata->axis_z.r_grid  = offload_data->axis_grid;
    Bdata->axis_z.c       = &(offload_array[3*B_size + psi_size + axis_size]);
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
a5err B_STS_eval_psi(real psi[], real r, real phi, real z,
                   B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += interp3Dcomp_eval_B(&psi[0], &Bdata->psi, r, phi, z);

    /* Test for psi interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );}

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
a5err B_STS_eval_psi_dpsi(real psi_dpsi[], real r, real phi, real z,
                          B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real psi_dpsi_temp[10];

    interperr += interp3Dcomp_eval_dB(psi_dpsi_temp, &Bdata->psi, r, phi, z);

    psi_dpsi[0] = psi_dpsi_temp[0];
    psi_dpsi[1] = psi_dpsi_temp[1];
    psi_dpsi[2] = psi_dpsi_temp[2];
    psi_dpsi[3] = psi_dpsi_temp[3];

    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );}

    return err;
}

/**
 * @brief Evaluate normalized poloidal flux rho
 *
 * @param rho pointer where rho value will be stored
 * @param psi corresponding poloidal flux psi [V*s*m^-1]
 * @param Bdata pointer to magnetic field data struct
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_eval_rho(real rho[], real psi, B_STS_data* Bdata) {
    a5err err = 0;

    /* Normalize psi to get rho */
    rho[0] = sqrt(fabs( (psi - Bdata->psi0) / (Bdata->psi1 - Bdata->psi0) ));

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
a5err B_STS_eval_rho_drho(real rho_drho[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    real rho;

    B_STS_eval_psi_dpsi(rho_drho, r, phi, z, Bdata);
    /* Convert: rho = sqrt(psi), drho = dpsi/(2 * sqrt(psi)) */
    rho = sqrt(rho_drho[0]);
    rho_drho[0] = rho;
    rho_drho[1] = rho_drho[1] / (2*rho);
    rho_drho[2] = rho_drho[2] / (2*rho);
    rho_drho[3] = rho_drho[3] / (2*rho);

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
a5err B_STS_eval_B(real B[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */

    interperr += interp3Dcomp_eval_B(&B[0], &Bdata->B_r, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[1], &Bdata->B_phi, r, phi, z);
    interperr += interp3Dcomp_eval_B(&B[2], &Bdata->B_z, r, phi, z);

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );}

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B[0]*B[0] + B[1]*B[1] + B[2]*B[2]) == 0);
    if(!err && check) {err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );}

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
a5err B_STS_eval_B_dB(real B_dB[], real r, real phi, real z, B_STS_data* Bdata) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    real B_dB_temp[10];

    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_r, r, phi, z);

    B_dB[0] = B_dB_temp[0];
    B_dB[1] = B_dB_temp[1];
    B_dB[2] = B_dB_temp[2];
    B_dB[3] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_phi, r, phi, z);

    B_dB[4] = B_dB_temp[0];
    B_dB[5] = B_dB_temp[1];
    B_dB[6] = B_dB_temp[2];
    B_dB[7] = B_dB_temp[3];

    interperr += interp3Dcomp_eval_dB(B_dB_temp, &Bdata->B_z, r, phi, z);

    B_dB[8] = B_dB_temp[0];
    B_dB[9] = B_dB_temp[1];
    B_dB[10] = B_dB_temp[2];
    B_dB[11] = B_dB_temp[3];

    /* Test for B field interpolation error */
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );}

    /* Check that magnetic field seems valid */
    int check = 0;
    check += ((B_dB[0]*B_dB[0] + B_dB[4]*B_dB[4] + B_dB[8]*B_dB[8]) == 0);
    if(!err && check) {err = error_raise( ERR_INPUT_UNPHYSICAL, __LINE__, EF_B_STS );}

    return err;
}

/**
 * @brief Return magnetic axis R-coordinate
 *
 * @param axis_r pointer where axis R [m] value will be stored
 * @param Bdata pointer to magnetic field data struct
 * @param phi phi coordinate [deg]
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_get_axis_r(real* axis_r, B_STS_data* Bdata, real phi) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    phi = fmod(phi, 2*CONST_PI);
    if(phi < 0) {
        phi += 2*CONST_PI;
    }
    interperr += linint1D_eval(axis_r, &Bdata->axis_r, phi);
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );}
    return err;
}

/**
 * @brief Return magnetic axis z-coordinate
 *
 * @param axis_z pointer where axis z [m] value will be stored
 * @param Bdata pointer to magnetic field data struct
 * @param phi phi coordinate [deg]
 *
 * @return Non-zero a5err value if evaluation failed, zero otherwise
 */
a5err B_STS_get_axis_z(real* axis_z, B_STS_data* Bdata, real phi) {
    a5err err = 0;
    int interperr = 0; /* If error happened during interpolation */
    phi = fmod(phi, 2*CONST_PI);
    if(phi < 0) {
        phi += 2*CONST_PI;
    }
    interperr += linint1D_eval(axis_z, &Bdata->axis_z, phi);
    if(interperr) {err = error_raise( ERR_INPUT_EVALUATION, __LINE__, EF_B_STS );}
    return err;
}
