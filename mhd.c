/**
 * @file mhd.c
 * @brief Module for evaluating MHD parameters.
 */
#include <stdlib.h>
#include "ascot5.h"
#include "print.h"
#include "error.h"
#include "mhd.h"
#include "spline/interp1D.h"
#include "spline/interp1Dcomp.h"

/**
 * @brief Load MHD data and prepare parameters for offload.
 *
 * This function fills the MHD offload struct with parameters and allocates
 * and fills the offload array. Sets offload array length in the offload struct.
 *
 * It is assumed that the offload_data struct is completely filled before
 * calling this function (except for the offload_array_length and psigrid).
 * Furthermore, offload array should contain following data:
 *
 * - offload_array[j*npsi + i] : alpha(mode_j, psi_i).
 * - offload_array[n_modes*npsi + j*npsi + i] : phi(mode_j, psi_i).
 *
 * 1D splines are constructed here and stored to offload array which is
 * reallocated.
 *
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to pointer to offload array
 *
 * @return zero if initialization succeeded.
 */
int mhd_init_offload(mhd_offload_data* offload_data,
                     real** offload_array) {

    offload_data->psigrid = (offload_data->psimax - offload_data->psimin)
        / (offload_data->npsi - 1);

    /* Allocate array for storing coefficients (which later replaces the
       offload array) */
    int cpernode = 2; // Number of coefficients in one node.
    real* coeff_array = (real*)malloc(2 * cpernode * offload_data->n_modes
                                      * offload_data->npsi * sizeof(real));

    /* Go through all modes, initialize splines for each, copy coefficients
       and free splines. */
    int err     = 0;
    int npsi    = offload_data->npsi;
    int n_modes = offload_data->n_modes;
    for(int j=0; j<offload_data->n_modes; j++) {
        interp1D_data spline;
        err += interp1Dcomp_init(&spline, (offload_array[j*npsi]),
                                 offload_data->npsi, offload_data->psimin,
                                 offload_data->psimax,
                                 offload_data->psigrid);
        for(int i=0; i<cpernode*npsi; i++) {
            coeff_array[j*cpernode*npsi + i] = spline.c[i];
        }
        interp1Dcomp_free(&spline);

        err += interp1Dcomp_init(&spline, (offload_array[npsi*n_modes+j*npsi]),
                                 offload_data->npsi, offload_data->psimin,
                                 offload_data->psimax,
                                 offload_data->psigrid);
        for(int i=0; i<cpernode*npsi; i++) {
            coeff_array[cpernode*npsi*n_modes + j*cpernode*npsi + i] = spline.c[i];
        }
        interp1Dcomp_free(&spline);
    }

    free(offload_array);
    *offload_array = coeff_array;
    offload_data->offload_array_length = n_modes*npsi*cpernode;

    /* Print some sanity check on data */
    print_out(VERBOSE_IO, "\nMHD input\n");
    print_out(VERBOSE_IO, "Grid: npsi = %4.d psimin = %3.3f psimax = %3.3f\n",
              offload_data->npsi,
              offload_data->psimin, offload_data->psimax);

    print_out(VERBOSE_IO, "\nModes:\n");
    for(int j=0; j<n_modes; j++) {
        print_out(VERBOSE_IO, "(n,m) = (%2.d,%2.d) Amplitude = %3.3g Frequency = %3.3g\n",
                  offload_data->nmode[j], offload_data->mmode[j],
                  offload_data->amplitude_nm[j], offload_data->omega_nm[j]);
    }

    return err;
}

/**
 * @brief Initialize MHD data struct on target
 *
 * @param mhddata pointer to data struct on target
 * @param offload_data pointer to offload data struct
 * @param offload_array pointer to offload array
 */
void mhd_init(mhd_data* MHDdata, mhd_offload_data* offload_data,
              real* offload_array) {
    MHDdata->n_modes = offload_data->n_modes;

    int npsi     = offload_data->npsi;
    int n_modes  = offload_data->n_modes;
    int cpernode = 2;
    for(int j=0; j<MHDdata->n_modes; j++) {
        MHDdata->nmode[j]        = offload_data->nmode[j];
        MHDdata->mmode[j]        = offload_data->mmode[j];
        MHDdata->amplitude_nm[j] = offload_data->amplitude_nm[j];
        MHDdata->omega_nm[j]     = offload_data->omega_nm[j];

        MHDdata->alpha_nm[j].n_r    = offload_data->npsi;
        MHDdata->alpha_nm[j].r_min  = offload_data->psimin;
        MHDdata->alpha_nm[j].r_max  = offload_data->psimax;
        MHDdata->alpha_nm[j].r_grid = offload_data->psigrid;
        MHDdata->alpha_nm[j].c      = &(offload_array[j*npsi]);

        MHDdata->phi_nm[j].n_r    = offload_data->npsi;
        MHDdata->phi_nm[j].r_min  = offload_data->psimin;
        MHDdata->phi_nm[j].r_max  = offload_data->psimax;
        MHDdata->phi_nm[j].r_grid = offload_data->psigrid;
        MHDdata->phi_nm[j].c      =
            &(offload_array[cpernode*npsi*n_modes + j*cpernode*npsi]);
    }
}

/**
 * @brief MHD magic here.
 *
 * @todo This is just a dummy.
 */
a5err mhd_eval(mhd_data* mhddata) {
    return 0;
}
