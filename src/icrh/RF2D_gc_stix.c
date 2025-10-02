/**
 * @file RF2D_gc_stix.c
 * @brief Implementation of the Stix diffusion operator for ICRH
 *
 * stix_data file contains the implementation of the Stix diffusion operator
 * for the ASCOT code, specifically for the 2D representation of the
 * radio-frequency (RF) heating in the context of gyrokinetic simulations.
 * It includes functions for initializing the Stix data, computing
 * resonances, scattering particles, and evaluating fields.
 */

#include <hdf5.h>
#include <hdf5_hl.h>
#include <string.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_airy.h>
#include <gsl/gsl_errno.h>
#include <math.h>
#include "RF2D_gc_stix.h"
#include "../error.h"
#include "../print.h"
#include "../hdf5io/hdf5_helpers.h"
#include "../ascot5.h"
#include "../B_field.h"
#include "../math.h"
#include "../spline/interp.h"
#include "../particle.h"

#define M_PI 3.14159265358979323846

/**
 * @brief Initializes the RF2D_gc_stix structure with data from an HDF5 file.
 *
 * stix_data function initializes the RF2D_gc_stix structure with the necessary
 * data for the Stix diffusion operator, including wave fields, magnetic field,
 * and other parameters.
 *
 * @param stix_data Pointer to the RF2D_gc_stix structure to be initialized.
 * @param f HDF5 file identifier containing the required data.
 * @param qid Active ID of the RF field to activate.
 * @param lhigh Flag indicating whether high-frequency data is used.
 * @param bdata Pointer to the B_field_data structure containing magnetic field data.
 * @return Error code indicating success or failure of the initialization.
 */
a5err RF2D_gc_stix_init_from_file(RF2D_gc_stix* stix_data, hid_t f, char* qid,
                        int lhigh, B_field_data* bdata){
    // Initialize the Stix data structure with data from the HDF5 file.
    a5err err = 0;
    
    // Check if the stix_data pointer is valid.
    if (stix_data == NULL) {
        print_err("RF2D_gc_stix_init: stix_data pointer is NULL.\n");
        return 1;
    }

    if(bdata == NULL) {
        print_err("RF2D_gc_stix_init: bdata pointer is NULL.\n");
        return 1;
    }

    if(lhigh < 1){
        print_err("RF2D_gc_stix_init: lhigh must be at least 1.\n");
        return 1;
    }

    if (f < 0) {
        print_err("RF2D_gc_stix_init: Invalid file_id.\n");
        return 1;
    }

    #ifdef RFPATHSTIX
    #undef RFPATHSTIX
    #endif
    #define RFPATHSTIX "/RF/RF2D_Stix_XXXXXXXXXX/" 

    // Internal parameters to store temporarily what 
    // we read from the HDF5 file
    char tmp[256];
    int nr, nz;

    /* Read data */
    if( hdf5_read_int(RFPATHSTIX "nr", &nr,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(RFPATHSTIX "nz", &nz,
                      f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(RFPATHSTIX "nwaves", &stix_data->nwaves,
                      f, qid, __FILE__, __LINE__) ) {return 1;}

    if(stix_data->nwaves <= 0) {
        print_err("RF2D_gc_stix_init: Number of waves must be positive.\n");
        return 1;
    }

    stix_data->omega = (real*) malloc(stix_data->nwaves * sizeof(real));
    stix_data->ntor = (int*) malloc(stix_data->nwaves * sizeof(real));


    if( hdf5_read_double(RFPATHSTIX "omega", stix_data->omega,
                            f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(RFPATHSTIX "ntor", stix_data->ntor,
                        f, qid, __FILE__, __LINE__) ) {return 1;}


    if( hdf5_read_double(RFPATHSTIX "rmin", &stix_data->rmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATHSTIX "rmax", &stix_data->rmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATHSTIX "zmin", &stix_data->zmin,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_double(RFPATHSTIX "zmax", &stix_data->zmax,
                         f, qid, __FILE__, __LINE__) ) {return 1;}
    if( hdf5_read_int(RFPATHSTIX "include_Eminus", &stix_data->include_Eminus,
                         f, qid, __FILE__, __LINE__) ){
        print_err("RF2D_gc_stix_init: The include_Eminus flag is missing in the HDF5 file: assuming it is 1.\n");
        stix_data->include_Eminus = 1; // Default to 1 if not found
    }
    if( hdf5_read_int(RFPATHSTIX "include_stochastic", &stix_data->include_stochastic,
                         f, qid, __FILE__, __LINE__) ){
        print_err("RF2D_gc_stix_init: The include_stochastic flag is missing in the HDF5 file: assuming it is 1.\n");
        stix_data->include_stochastic = 1; // Default to 1 if not found
    }

    if( hdf5_read_int(RFPATHSTIX "include_vpara_kick", &stix_data->include_vpara_kick,
                         f, qid, __FILE__, __LINE__) ){
        print_err("RF2D_gc_stix_init: The include_vpara_kick flag is missing in the HDF5 file: assuming it is 1.\n");
        stix_data->include_vpara_kick = 1; // Default to 1 if not found
    }
    if( hdf5_read_int(RFPATHSTIX "include_phase_factor", &stix_data->include_phase_factor,
                         f, qid, __FILE__, __LINE__) ){
        print_err("RF2D_gc_stix_init: The include_phase_factor flag is missing in the HDF5 file: assuming it is 1.\n");
        stix_data->include_phase_factor = 1; // Default to 1 if not found
    }
    
    int nsize = nr * nz * stix_data->nwaves;
    real* fields[7]; // Pointers to the fields in the buffer
    char *fieldnames[] = {
        "Eplus_re", "Eplus_im", "Eminus_re", "Eminus_im",
        "kperp", "costheta",  "sintheta"
    };

    for(int i = 0; i < 7; i++) fields[i] = NULL;

    for(int i = 0; i < 7; i++) {
        fields[i] = (real*) malloc(nsize * sizeof(real));
        snprintf(tmp, 256, "%s%s", RFPATHSTIX, fieldnames[i]);
        if( hdf5_read_double(tmp, fields[i],
                             f, qid, __FILE__, __LINE__) ){
            err = 1;
            break;
            }
    }
    if(err) {
        print_err("RF2D_gc_stix_init: Error reading fields from HDF5 file.\n");
        for(int i = 0; i < 7; i++) {
            if(fields[i] != NULL) free(fields[i]);
        }
        return 1;
    }

    // We make some aliases to the fields for easier access.
    real* Eplus_re = fields[0];
    real* Eplus_im = fields[1];
    real* Eminus_re = fields[2];
    real* Eminus_im = fields[3];
    real* kperp = fields[4];
    real* costheta = fields[5];
    real* sintheta = fields[6];

    // Initialize the Stix data structure.
    err = RF2D_gc_stix_init(stix_data, Eplus_re, Eplus_im,
                            Eminus_re, Eminus_im,
                            kperp, costheta, sintheta,
                            stix_data->rmin, stix_data->rmax,
                            stix_data->zmin, stix_data->zmax,
                            nr, nz, stix_data->nwaves,
                            stix_data->include_Eminus,
                            stix_data->include_stochastic,
                            stix_data->include_vpara_kick,
                            stix_data->include_phase_factor,
                            bdata);
    for(int i = 0; i < 7; i++) {
        if(fields[i] != NULL) free(fields[i]);
    }
    if(err) {
        print_err("RF2D_gc_stix_init: Error initializing Stix data from fields.\n");
        return err;
    }

    // Set the B_field_data pointer.
    stix_data->bdata = bdata;
    err = 0;
    return err;
}


/**
 * @brief Initializes the RF2D_gc_stix structure with data from provided arrays.
 *
 * RF2D_gc_stix_init_from_data function initializes the RF2D_gc_stix structure
 * with the necessary data for the Stix diffusion operator, including wave fields,
 * magnetic field, and other parameters. It computes the required fields and sets up
 * the interpolation objects for the wave fields.
 *
 * @param stix_data Pointer to the RF2D_gc_stix structure to be initialized.
 * @param Eplus_re Pointer to the real part of the Eplus field.
 * @param Eplus_im Pointer to the imaginary part of the Eplus field.
 * @param Eminus_re Pointer to the real part of the Eminus field.
 * @param Eminus_im Pointer to the imaginary part of the Eminus field.
 * @param kperp Pointer to the kperp field.
 * @param costerm Pointer to the cosine term for the phase factor.
 * @param sinterm Pointer to the sine term for the phase factor.
 * @param rmin Minimum radial coordinate.
 * @param rmax Maximum radial coordinate.
 * @param zmin Minimum vertical coordinate.
 * @param zmax Maximum vertical coordinate.
 * @param nr Number of radial grid points.
 * @param nz Number of vertical grid points.
 * @param nwaves Number of waves in the Stix diffusion operator.
 * @param include_Eminus Flag to include Eminus in the Stix diffusion operator.
 * @param include_stochastic Flag to include stochastic diffusion in the Stix operator.
 * @param include_vpara_kick Flag to include vpara kick in the Stix operator
 * @param include_phase_factor Flag to include phase factor in the Stix operator.
 * @param bdata Pointer to the B_field_data structure containing magnetic field data.
 * @return Error code indicating success or failure of the initialization.
 */
a5err RF2D_gc_stix_init(RF2D_gc_stix* stix_data,
                                  real* Eplus_re, real* Eplus_im,
                                  real* Eminus_re, real* Eminus_im,
                                  real* kperp, real* costerm, real* sinterm,
                                  real rmin, real rmax, real zmin, real zmax,
                                  int nr, int nz, int nwaves,
                                  int include_Eminus, int include_stochastic,
                                  int include_vpara_kick, int include_phase_factor,
                                  B_field_data* bdata){
    a5err err = 0;
    
    // Allocating memory for the Stix data structure.
    if (stix_data == NULL) {
        print_err("RF2D_gc_stix_init_from_data: stix_data pointer is NULL.\n");
        err = 1;
        return err;
    }
    if(!Eplus_re || !Eplus_im || !Eminus_re || !Eminus_im ||
       !kperp || !costerm || !sinterm) {
        print_err("RF2D_gc_stix_init_from_data: One or more field pointers are NULL.\n");
        err = 1;
        return err;
    }

    if(nr <= 0 || nz <= 0 || nwaves <= 0) {
        print_err("RF2D_gc_stix_init_from_data: nr, nz, and nwaves must be positive.\n");
        return err;
    }

    int nsize = nr * nz * nwaves;
    stix_data->rmin = rmin;
    stix_data->rmax = rmax;
    stix_data->zmin = zmin;
    stix_data->zmax = zmax;
    stix_data->nr = nr;
    stix_data->nz = nz;
    stix_data->nwaves = nwaves;
    stix_data->include_Eminus = include_Eminus;
    stix_data->include_stochastic = include_stochastic;
    stix_data->include_vpara_kick = include_vpara_kick;
    stix_data->include_phase_factor = include_phase_factor;
    stix_data->bdata = bdata;

    real* Eplus2 = (real*) malloc(nsize * sizeof(real));
    real* Eminus2 = (real*) malloc(nsize * sizeof(real));
    real* Ecross2 = (real*) malloc(nsize * sizeof(real));

    if(!Eplus2 || !Eminus2 || !Ecross2) {
        print_err("RF2D_gc_stix_init_from_data: Memory allocation failed.\n");
        return err;
    }

    // Copy the fields into the Stix data structure.
    // Compute the Eplus2, Eminus2 and the cross term.
    #pragma omp parallel for
    for(int i = 0; i < nsize; i++) {
        // 2.0d0 * ( Eplus_re[i]**2 + Eplus_im[i]**2 )
        Eplus2[i] = 2.0 * (Eplus_re[i] * Eplus_re[i] + Eplus_im[i] * Eplus_im[i]);
        // 2.0d0 * ( Eminus_re[i]**2 + Eminus_im[i]**2 )
        Eminus2[i] = 2.0 * (Eminus_re[i] * Eminus_re[i] + Eminus_im[i] * Eminus_im[i]);

        real cos2term = 1.0;
        real sin2term = 0.0;
        if(include_phase_factor) {
            // Re(Eplus2 * conj(Eminus))
            sin2term = 2.0 * costerm[i] * sinterm[i];
            cos2term = 2.0 * costerm[i] * costerm[i] - 1.0;
        }
        Ecross2[i] = + 4.0 * cos2term * ( Eplus_re[i] * Eminus_re[i]   + \
                                          Eplus_im[i] * Eminus_im[i] )   \
                     - 4.0 * sin2term * ( Eplus_re[i] * Eminus_im[i]  - \
                                          Eplus_im[i] * Eminus_re[i] );

    }
    // Set the enabled flag to 1.
    stix_data->enabled = 1;

    // Generating the interpolation objects for the fields.
    err = 0;
    stix_data->Eplus_2 = (interp2D_data*) malloc(sizeof(interp2D_data));
    stix_data->Eminus_2 = (interp2D_data*) malloc(sizeof(interp2D_data));
    stix_data->E2cross = (interp2D_data*) malloc(sizeof(interp2D_data));
    stix_data->kperp = (interp2D_data*) malloc(sizeof(interp2D_data));
    err += interp2Dcomp_setup(stix_data->Eplus_2, Eplus2, nr, nz, 
                              NATURALBC, NATURALBC,
                              rmin, rmax, zmin, zmax);
    err += interp2Dcomp_setup(stix_data->Eminus_2, Eminus2, nr, nz, 
                              NATURALBC, NATURALBC,
                              rmin, rmax, zmin, zmax);
    err += interp2Dcomp_setup(stix_data->E2cross, Ecross2, nr, nz, 
                              NATURALBC, NATURALBC,
                              rmin, rmax, zmin, zmax);
    err += interp2Dcomp_setup(stix_data->kperp, kperp, nr, nz, 
                              NATURALBC, NATURALBC,
                              rmin, rmax, zmin, zmax);

    free(Eplus2);
    free(Eminus2);
    free(Ecross2);

    print_out(VERBOSE_IO, "RF2D_gc_stix_init: Stix data initialized with %d waves.\n", nwaves);
    if(err) {
        print_err("RF2D_gc_stix_init: Error initializing interpolation objects.\n");
        return err;
    }
    print_out(VERBOSE_IO, "RF2D_gc_stix_init: Stix data initialized successfully.\n");
    print_out(VERBOSE_IO, "  Radial grid: %d points from %f to %f\n", nr, rmin, rmax);
    print_out(VERBOSE_IO, "  Vertical grid: %d points from %f to %f\n", nz, zmin, zmax);
    print_out(VERBOSE_IO, "  Number of waves: %d\n", nwaves);
    print_out(VERBOSE_IO, "  Include Eminus: %d\n", include_Eminus);
    print_out(VERBOSE_IO, "  Include stochastic diffusion: %d\n", include_stochastic);
    print_out(VERBOSE_IO, "  Include vpara kick: %d\n", include_vpara_kick);
    print_out(VERBOSE_IO, "  Include phase factor: %d\n", include_phase_factor);

    // Disabling the error raising from the GSL library.
    gsl_set_error_handler_off();
    
    err = 0;
    return err;
}

/**
 * @brief Frees the memory allocated for the RF2D_gc_stix structure.
 *
 * @param stix_data Pointer to the RF2D_gc_stix structure to be freed.
 */
void RF2D_gc_stix_free(RF2D_gc_stix* stix_data){
    // Free the memory allocated for the Stix data structure.
    if (stix_data == NULL) {
        print_err("RF2D_gc_stix_free: stix_data pointer is NULL.");
        return;
    }

    if (stix_data->Eplus_2 != NULL) {
        free(stix_data->Eplus_2->c);
        free(stix_data->Eplus_2);
        stix_data->Eplus_2 = NULL;
    }
    if (stix_data->Eminus_2 != NULL) {
        free(stix_data->Eminus_2->c);
        free(stix_data->Eminus_2);
        stix_data->Eminus_2 = NULL;
    }
    if (stix_data->E2cross != NULL) {
        free(stix_data->E2cross->c);
        free(stix_data->E2cross);
        stix_data->E2cross = NULL;
    }
    if (stix_data->kperp != NULL) {
        free(stix_data->kperp->c);
        free(stix_data->kperp);
        stix_data->kperp = NULL;
    }

    if(stix_data->omega !=NULL) free(stix_data->omega);
    if(stix_data->ntor !=NULL) free(stix_data->ntor);

    // Free resonance arrays.
    if(stix_data->R_resonances != NULL) {
        for(int i = 0; i < stix_data->nwaves; i++) {
            if(stix_data->R_resonances[i] != NULL) {
                free(stix_data->R_resonances[i]);
            }
        }
        free(stix_data->R_resonances);
    }

    if(stix_data->nres != NULL) free(stix_data->nres);
    if(stix_data->res_nums != NULL) {
        for (int i = 0; i < stix_data->nwaves; i++) {
            if(stix_data->res_nums[i] != NULL) {
                free(stix_data->res_nums[i]);
            }
        }
    }
    stix_data->res_nums = NULL;
    stix_data->nres = NULL;


    // Free the Stix data structure itself.
    stix_data->enabled = 0;

    // Setting all elements pointer to NULL
    stix_data->Eplus_2 = NULL;
    stix_data->Eminus_2 = NULL;
    stix_data->E2cross = NULL;
    stix_data->kperp = NULL;
    stix_data->omega = NULL;
    stix_data->ntor = NULL;
}


/**
 * @brief Offloads the RF2D_gc_stix structure to the GPU or other device.
 *
 * @param stix_data Pointer to the RF2D_gc_stix structure to be offloaded.
 */
void RF2D_gc_stix_offload(RF2D_gc_stix* stix_data){
    // Offload the Stix data structure to the GPU or other device.
    if (stix_data == NULL) {
        print_err("RF2D_gc_stix_offload: stix_data pointer is NULL.");
        return;
    }

    int nsize = stix_data->nr * stix_data->nz * stix_data->nwaves * NSIZE_COMP2D;
    GPU_MAP_TO_DEVICE(stix_data[0:1],
                      stix_data->Eplus_2, stix_data->Eminus_2,
                      stix_data->E2cross, stix_data->kperp,
                      stix_data->omega[0:stix_data->nwaves],
                      stix_data->ntor[0:stix_data->nwaves],
                      stix_data->Eplus_2->c[0:nsize],
                      stix_data->Eminus_2->c[0:nsize],
                      stix_data->E2cross->c[0:nsize],
                      stix_data->kperp->c[0:nsize],
                      stix_data->R_resonances[0:stix_data->nwaves],
                      stix_data->res_nums[0:stix_data->nwaves],
                      stix_data->nres[0:stix_data->nwaves]
                    );

    for (int i = 0; i < stix_data->nwaves; i++) {
        int nres = stix_data->nres[i];
        if (nres <= 0) continue; // Skip if no resonances
        GPU_MAP_TO_DEVICE(stix_data->R_resonances[i][0:stix_data->nres],
                          stix_data->res_nums[i][0:stix_data->nres]);
    }
}

/**
 * @brief Computes the cold resonance locations based on the frequencies and the
 * magnetic field strenght.
 * 
 * It will look for the cold resonance condition, i.e., the condition ignoring
 * the Doppler shift. The resonance condition is given by:
 * \[ \omega_\mathrm{RF} = l_i \Omega_\mathrm{cyclotron} \]
 * being l_i the resonance number, which is an integer.
 * 
 * The purpose of this function is to compute the radial locations of the cold
 * resonances, assuming approximately that B \propto 1/R.
 *
 * @param stix_data Pointer to the RF2D_gc_stix structure to be initialized.
 * @param bdata Pointer to the B_field_data structure containing magnetic field data.
 * @param n_max_res Maximum number of resonance locations to compute.
 * @param qm Charge-to-mass ratio of the particles.
 */
a5err RF2D_gc_stix_compute_cold_resonances(RF2D_gc_stix* stix_data,
                                           B_field_data* bdata, int n_max_res,
                                           real qm){
    if(stix_data == NULL || bdata == NULL) {
        print_err("RF2D_gc_stix_compute_cold_resonances: stix_data or bdata pointer is NULL.");
        return 1;
    }

    if(n_max_res <= 0) {
        print_err("RF2D_gc_stix_compute_cold_resonances: n_max_res must be positive.");
        return 1;
    }

    stix_data->n_max_res = n_max_res;
    stix_data->R_resonances = (real**) malloc(stix_data->nwaves * sizeof(real*));
    stix_data->nres = (int*) malloc(stix_data->nwaves * sizeof(int));
    stix_data->res_nums = (int**) malloc(stix_data->nwaves * sizeof(int*));

    // Checking that we are dealing with a 2D magnetic field, as otherwise the resonances
    // would not be only function of the major radius coordinate.
    if(bdata->type != B_field_type_2DS) {
        print_err("RF2D_gc_stix_compute_cold_resonances: B_field_data must be of type B_field_type_2DS.");
        // Creating dummy data.
        for(int i = 0; i < stix_data->nwaves; i++) {
            stix_data->R_resonances[i] = (real*) malloc(n_max_res * sizeof(real));
            stix_data->res_nums[i] = (int*) malloc(n_max_res * sizeof(int));
            stix_data->nres[i] = n_max_res; // All the waves should look for the same number of resonances.
            for(int j = 0; j < n_max_res; j++) {
                stix_data->R_resonances[i][j] = 0.0; // Dummy value
                stix_data->res_nums[i][j] = j + 1; // Resonance number
            }
        }

        return 0;
    }

    real Rzaxis[2];
    B_field_get_axis_rz(&Rzaxis[0], bdata, 0.0); // Get the axis R and Z coordinates
    real zaxis = Rzaxis[1]; 
    
    real* Rgrid = (real*) malloc(1024 * sizeof(real));
    real Btmp[3];
    real* Bnorm = (real*) malloc(1024 * sizeof(real));
    if(Rgrid == NULL || Bnorm == NULL) {
        print_err("RF2D_gc_stix_compute_cold_resonances: Memory allocation failed for Rgrid or Bnorm.");
        return 1;
    }
    for(int i = 0; i < 1024; i++) {
        Rgrid[i] = stix_data->rmin + i * (stix_data->rmax - stix_data->rmin) / 1023.0;

        // Evaluating the local magnetic field strength at the point.
        B_field_eval_B(&Btmp[0], Rgrid[i], 0.0, zaxis, 0.0, bdata);

        Bnorm[i] = sqrt(Btmp[0] * Btmp[0] + Btmp[1] * Btmp[1] + Btmp[2] * Btmp[2]);
    }

    real* Rres_tmp = (real*) malloc(n_max_res * sizeof(real));
    int*  res_indx = (int*)  malloc(n_max_res * sizeof(int));
    int maxreson = 0;
    memset(Rres_tmp, 0, n_max_res * sizeof(real));
    memset(res_indx, 0, n_max_res * sizeof(int));

    for (int iwave=0; iwave < stix_data->nwaves; iwave++){
        int nres = 0;
        for(int i = 0; i < n_max_res; i++){
            int l_res = i + 1; // Resonance number
            real prev = stix_data->omega[iwave] - l_res * qm * Bnorm[0];
            real next = stix_data->omega[iwave] - l_res * qm * Bnorm[0];

            // We loop through the radial coordinate to find the point where there is a
            // change of sign in the resonance condition -> Resonance found. Then
            // we perform linear interpolation to find the exact resonance location.
            for(int j = 1; j < 1024; j++) {
                prev = next;
                next = stix_data->omega[iwave] - l_res * qm * Bnorm[j];
                if(prev*next < 0.0){
                    // We have a root between j-1 and j.
                    real Rres = Rgrid[j-1] - fabs(prev) * (Rgrid[j] - Rgrid[j-1])  / (fabs(next) - fabs(prev));
                    Rres_tmp[i] = Rres;
                    res_indx[i] = l_res; // Store the index of the root
                    nres++;
                    break; // We found a resonance, break the loop.
                }
            }
        }
        // We store the total resonance number.
        stix_data->nres[iwave] = nres;
        maxreson = fmax(maxreson, nres);

        stix_data->R_resonances[iwave] = (real*) malloc(nres * sizeof(real));
        stix_data->res_nums[iwave] = (int*) malloc(nres * sizeof(int));

        // We fill the resonance locations and numbers.
        int idx_res = 0;
        for(int i = 0; i < n_max_res; i++) {
            if(res_indx[i] > 0) {
                stix_data->R_resonances[iwave][idx_res] = Rres_tmp[i];
                stix_data->res_nums[iwave][idx_res] = res_indx[i];
                idx_res++;
            }
        }

        // Reset the temporary arrays for the next wave.
        memset(Rres_tmp, 0, n_max_res * sizeof(real));
        memset(res_indx, 0, n_max_res * sizeof(int));
    }
    stix_data->n_max_res = 0;

    // Debugging the resonances found.
    print_out(VERBOSE_IO, "RF2D_gc_stix_compute_cold_resonances: Cold resonances found:\n");
    for (int iwave=0; iwave < stix_data->nwaves; iwave++){
        print_out(VERBOSE_IO, " > Wave %d: found %d resonances:\n", iwave, stix_data->nres[iwave]);
        for(int i = 0; i < stix_data->nres[iwave]; i++) {
            print_out(VERBOSE_IO, "   >> Resonance l=%d: R = %f\n", stix_data->res_nums[iwave][i], stix_data->R_resonances[iwave][i]);
            stix_data->n_max_res = fmax(stix_data->n_max_res, stix_data->res_nums[iwave][i]);
        }
    }

    // Free the temporary arrays.
    free(Rgrid);
    free(Bnorm);
    free(Rres_tmp);
    free(res_indx);
    return 0;
}

/**
 * @brief Computes the interaction time for a particle with the RF wave during resonance.
 * 
 * This function calculates the squared interaction time for a particle with the RF wave
 * during resonance, based on the particle's history and the Stix diffusion operator.
 * It uses the particle's gyrofrequency and parallel momentum to compute the resonance
 * condition derivatives and the interaction time. 
 * 
 * The interaction time is based on the fixed-point technique to integrate a fast oscillating
 * integral. See more details in "M. Choi, PPCF (2010)".
 * 
 * @param stix_data Pointer to the RF2D_gc_stix structure containing the Stix data.
 * @param hist Pointer to the RF_particle_history structure containing the particle's history.
 * @param iwave Index of the wave for which the interaction time is computed.
 * @return The squared interaction time for the particle with the RF wave, or -1.
 *         If the interaction time cannot be computed, it returns -1.
 */
real RF2D_gc_stix_get_interaction_time(RF2D_gc_stix* stix_data, 
                                       RF_particle_history* hist, 
                                       int iwave, int l){
    real time = -1.0; // Default value indicating no interaction time computed.

    const int curr = 0; // Current index
    const int prev = 1; // Previous index
    const int prev_prev = 2; // Previous-previous index

    // Time step calculations
    real ddt = hist->dt[curr] + hist->dt[prev];
    real ddt2 = hist->dt[curr]*hist->dt[curr] + hist->dt[prev]*hist->dt[prev];
    real diff_dt = hist->dt[curr] - hist->dt[prev];
    real kpara = stix_data->ntor[iwave] / hist->R[curr]; // Parallel wave vector

    real nu_curr = kpara * hist->rhopara[curr] * hist->bnorm[curr] + hist->qm * l * hist->bnorm[curr];
    real nu_prev = kpara * hist->rhopara[prev] * hist->bnorm[prev] + hist->qm * l * hist->bnorm[prev];
    real nu_prev_prev = kpara * hist->rhopara[prev_prev] * hist->bnorm[prev_prev] + hist->qm * l * hist->bnorm[prev_prev];

    // Resonance condition derivatives: one sided derivatives.
    real nudot = ( nu_curr - nu_prev_prev ) / ddt;
    real nudot2 = 2.0 / ddt2 * (
        nu_curr - 2.0 * nu_prev + nu_prev_prev
        - nudot * diff_dt
    );

    // Squared interaction time
    if (fabs(nudot) > (0.5 * fabs(nudot2) * hist->dt[curr])) {
        time = M_PI / fabs(nudot);
    } else {
        // Airy function argument
        real zeta = -pow(nudot, 2) / pow(2.0 * nudot2, 2.0/3.0);
        // You need to implement or link an Airy function Ai(zeta)
        // For now, assume a function real Airy_Ai(real x) is available.
        real airy_val = gsl_sf_airy_Ai(zeta, GSL_PREC_DOUBLE);
        time = 2.0 * M_PI * airy_val / pow(fabs(nudot2 / 2.0), 1.0/3.0);
        time = 0.5 * time * time;
    }

    return time;
}

/**
 * @brief Evaluates the RF fields at given coordinates and time.
 * 
 * Uses the interpolation objects to compute the electric fields, taking into
 * account the internal flags.
 * 
 * @param r Radial coordinate.
 * @param phi Toroidal angle coordinate (not used in this function).
 * @param z Vertical coordinate.
 * @param t Time coordinate (not used in this function).
 * @param stix_data Pointer to the RF2D_gc_stix structure containing the Stix data.
 * @param Eplus_2 Variable to store the squared Eplus field.
 * @param Eminus_2 Variable to store the squared Eminus field (if included).
 * @param E2cross Variable to store the cross term of the electric fields.
 * @param kperp Variable to store the perpendicular wave vector.
 * @return a5err Error code indicating success or failure of the evaluation.
 */
a5err RF2D_gc_stix_eval_fields(real r, real phi, real z, real t,
                              RF2D_gc_stix* stix_data,
                              real* Eplus_2, real* Eminus_2,
                              real* E2cross, real* kperp){
    // Setting the values to 0.
    *Eplus_2 = 0.0;
    *Eminus_2 = 0.0;
    *E2cross = 0.0;
    *kperp = 0.0;

    a5err err = 0;

    // Evaluate the fields at the given coordinates and time.
    if (!stix_data->enabled) return err;

    // Evaluating the fields using the interpolation objects.
    err += interp2Dcomp_eval_f(Eplus_2, stix_data->Eplus_2, r, z);
    if(stix_data->include_Eminus) err += interp2Dcomp_eval_f(Eminus_2, stix_data->Eminus_2, r, z);
    err += interp2Dcomp_eval_f(E2cross, stix_data->E2cross, r, z);
    err += interp2Dcomp_eval_f(kperp, stix_data->kperp, r, z);

    return err;
}


/**
 * @brief Applies the Stix scattering operator to a particle and updates its internal state.
 * 
 * @param stix Array with the Stix data for the RF field.
 * @param hist Pointer to the RF_particle_history structure containing the particle's history.
 * @param p Pointer to the particle_simd_gc structure containing the current particle state.
 */
void RF2D_gc_stix_scatter(RF2D_gc_stix* stix, RF_particle_history* hist,
                          particle_simd_gc* p, real* h, real* rnd, uint8* used){
    // Apply the Stix scattering operator to the particle and update its internal state.
    if (stix == NULL || hist == NULL || p == NULL) {
        print_err("RF2D_gc_stix_scatter: One or more input pointers are NULL.\n");
        return;
    }

    if (!stix->enabled) return; // Stix operator is not enabled.

    // Updating the internal status of the particles.
    #pragma omp simd
    for(int imrk=0; imrk < NSIMD; imrk++) {
        if(p->running[imrk] == 0) continue; // Particle is not running, skip.
        int imrk_rnd_idx = 2 * stix->n_max_res * stix->nwaves * imrk;
        RF_particle_history_update(&hist[imrk], p, imrk, h[imrk]);

        // Checking for all the waves whether there's been some resonance.
        for(int iwave=0; iwave < stix->nwaves; iwave++) {
            // Evaluating the number of kicks.
            int nkicks, lres;
            RF_particle_eval_nkicks(&hist[imrk], p, imrk, iwave, &nkicks, &lres);

            if(nkicks <= 0 || lres == -1) continue; // No kicks, go to the next wave.

            // We have a resonance crossing. Let's evaluate the fields.
            real Eplus_2 = 0.0;
            real Eminus_2 = 0.0;
            real E2cross = 0.0;
            real kperp = 0.0;

            // Getting the resonance time.
            real t_inter = RF2D_gc_stix_get_interaction_time(stix, &hist[imrk], iwave, lres);
            if(t_inter < 0.0 || !isfinite(t_inter)) continue; // No interaction time, go to the next wave.

            RF2D_gc_stix_eval_fields(p->r[imrk], p->phi[imrk], p->z[imrk], p->time[imrk],
                                     stix, &Eplus_2, &Eminus_2, &E2cross, &kperp);
            if(kperp == 0.0) continue; // Outside domain.

            // Evaluating so prefactors.
            for(int ikick = 0; ikick < nkicks; ikick++){
                real rhoL = 1/ p->charge[imrk] * sqrt(2.0 * p->mu[imrk] * p->mass[imrk] / hist[imrk].bnorm[0]);
                real arg = fabs(kperp * rhoL);

                // Evaluating the Bessel functions
                gsl_sf_result rJm1, rJl, rJp1;
                real Jm1, Jl, Jp1; // Bessel functions of order lres-1, lres and lres+1
                int err1, err2, err3;
                err1 = gsl_sf_bessel_Jn_e(lres - 1, arg, &rJm1);
                err2 = gsl_sf_bessel_Jn_e(lres, arg, &rJl);
                err3 = gsl_sf_bessel_Jn_e(lres + 1, arg, &rJp1);
                if(err1 || err2 || err3){
                    print_err("RF2D_gc_stix_scatter: Error computing Bessel function J_%d(%f): kperp = %.3e, mu = %.3e, Bnorm = %.2f, mass=%.3e, charge=%.2e.\n", \
                        lres, arg, kperp, p->mu[imrk], hist[imrk].bnorm[0], p->mass[imrk], p->charge[imrk]);
                    p->err[imrk] = error_raise( ERR_BESSEL_EVALUATION, __LINE__, EF_RF_GC2D );
                    continue;
                }
                Jm1 = rJm1.val;
                Jl = rJl.val;
                Jp1 = rJp1.val;

                // We evaluate the value of the derivatives using the recurrence relations.
                real dbessel_l_m1 = 0.5 * (lres - 1) * Jm1 - 0.5*arg*Jl;
                real dbessel_l_p1 = 0.5 * arg * Jl - 0.5*(lres + 1) * Jp1;

                real term1 = Eplus_2 * Jm1 * Jm1 + \
                             Eminus_2 * Jp1 * Jp1 + \
                             E2cross * Jm1 * Jp1;
                real term2 = 2 * Eplus_2 * Jm1 * dbessel_l_m1 + \
                             2 * Eminus_2 * Jp1 * dbessel_l_p1 + \
                             E2cross * (Jm1 * dbessel_l_p1 + Jp1 * dbessel_l_m1);
                
                real omega_cycl = hist[imrk].qm * hist[imrk].bnorm[0];
                real scaling_para = lres * omega_cycl / stix->omega[iwave] * p->charge[imrk];
                real scaling2 = scaling_para * scaling_para;

                // Kick in the magnetic moment.
                real dmu = scaling2 * (term1 + term2) * t_inter / (4.0 * hist[imrk].bnorm[0] * p->mass[imrk]);
                real kpara = stix->ntor[iwave] / hist[imrk].R[0];
                real drhopara = kpara / (lres * omega_cycl) * dmu;

                // Computing the stochastic contribution.
                real irnd = cos(lres * rnd[imrk_rnd_idx] * 2.0 * M_PI);
                used[imrk_rnd_idx] = 0; // Mark as used.

                // The following part is prone to underflow since the mu values are always
                // scaled by the mass of the particle. We temporarily remove the charge.
                dmu /= p->charge[imrk];
                real tmp_mu = p->mu[imrk] / p->charge[imrk];
                real dmu_stoch = sqrt(fabs(2 * tmp_mu * term1 / (term1 + term2) * dmu)) * irnd;
                imrk_rnd_idx++; // For next fetch.
                real drhopara_stoch = kpara / (lres * omega_cycl) * dmu_stoch;

                // We recover the original scaling.
                dmu *= p->charge[imrk];
                dmu_stoch *= p->charge[imrk];
                drhopara_stoch *= p->charge[imrk];

                // Updating the values of the parallel momentum and the magnetic moment.
                p->mu[imrk] +=  dmu_stoch + dmu;
                p->ppar[imrk] += (drhopara + drhopara_stoch) * hist[imrk].bnorm[0];
            }
        }
    }
}

/**
 * @brief From the particle queue, it gets the value of the charge
 * over mass ratio, that it will be used to compute the cold resonances.
 * 
 * This function will check whether all particles in the queue have the same
 * charge-over-mass ratio only for valid particles.
 * 
 * @param pq Particle queue to check.
 * @return real Value of the charge-over-mass ratio, or -1 if there's an error.
 */
real guess_qm(particle_queue* pq){
    if (pq == NULL) {
        print_err("guess_qm: pq pointer is NULL.\n");
        return -1.0;
    }
    real qm = -1.0;
    for(int i = 0; i < pq->n; i++) {
        if(pq->p[i]->mass <= 0.0) continue; // Invalid particle, skip.
        if(pq->p[i]->charge == 0.0) continue; // Neutral particle, skip.
        real qm_i = fabs(pq->p[i]->charge / pq->p[i]->mass);

        if(qm < 0.0) {
            qm = qm_i; // First valid particle, set qm.
            continue;
        }
        if(fabs(qm - qm_i) > 1e-3 * qm){
            print_err("guess_qm: Particles in the queue have different charge-over-mass ratios: %e vs %e.\n", qm, qm_i);
            return -1.0;
        }
        // They are the same, continue.
        qm = qm_i;
        
    }
    return qm;
}


/**
 * @brief Sets the usage flag for random numbers for a given particle.
 * 
 * This function will set the usage flag for all random numbers associated
 * with a specific particle in the Stix data structure.
 * 
 * @param stix_data Structure with the Stix data.
 * @param used Pointer to the array of usage flags.
 * @param imrk Index of the particle in the queue.
 * @param value Value to set the usage flag.
 */
void set_rndusage(RF2D_gc_stix* stix_data, uint8* used, int imrk, uint8 value){
    if(stix_data == NULL || used == NULL) {
        print_err("set_rndusage: stix_data or used pointer is NULL.\n");
        return;
    }
    int nflags = 2 * stix_data->n_max_res * stix_data->nwaves;
    int start_idx = nflags * imrk;
    for(int i = 0; i < nflags; i++) {
        used[start_idx + i] = value;
    }
}