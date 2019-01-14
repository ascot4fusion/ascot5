#include <stdlib.h>
#include "../ascot5.h"
#include "../diag/dist_5D.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_rho5D.h"
#include "../diag/dist_rho6D.h"
#include "../print.h"
#include "../math.h"
#include "hdf5_helpers.h"
#include "hdf5_histogram.h"

void hdf5_dist_write_5D(dist_5D_offload_data* dist, real* hist, char* filename,
                        char* qid) {

    #if VERBOSE > 0
    printf("\nWriting distributions to HDF5 file...\n");
    #endif

    int abscissa_dim = 7;
    int ordinate_length = 1;

    int abscissa_n_slots[7];
    abscissa_n_slots[0] = dist->n_r;
    abscissa_n_slots[1] = dist->n_phi;
    abscissa_n_slots[2] = dist->n_z;
    abscissa_n_slots[3] = dist->n_vpara;
    abscissa_n_slots[4] = dist->n_vperp;
    abscissa_n_slots[5] = dist->n_time;
    abscissa_n_slots[6] = dist->n_q;

    double abscissa_min[7];
    abscissa_min[0] = dist->min_r;
    abscissa_min[1] = math_rad2deg(dist->min_phi);
    abscissa_min[2] = dist->min_z;
    abscissa_min[3] = dist->min_vpara;
    abscissa_min[4] = dist->min_vperp;
    abscissa_min[5] = dist->min_time;
    abscissa_min[6] = dist->min_q;

    double abscissa_max[7];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = math_rad2deg(dist->max_phi);
    abscissa_max[2] = dist->max_z;
    abscissa_max[3] = dist->max_vpara;
    abscissa_max[4] = dist->max_vperp;
    abscissa_max[5] = dist->max_time;
    abscissa_max[6] = dist->max_q;

    char* abscissa_names[] = { "R", "phi", "z", "vpa", "vpe", "time", "charge" };
    char* abscissa_units[] = { "m", "deg", "m", "m/s", "m/s", "s", "e" };
    char* ordinate_names[] = { "density" };
    char* ordinate_units[] = { "s/m^5*deg*e" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run-XXXXXXXXXX/", qid, path);

    int retval;
    retval =  hdf5_histogram_write_uniform_double(
              filename,
              path,
              "R_phi_z_vpa_vpe_t_q",
              abscissa_dim,
              ordinate_length,
              abscissa_n_slots,
              abscissa_min,
              abscissa_max,
              abscissa_units,
              abscissa_names,
              ordinate_units,
              ordinate_names,
              hist);

    if(retval) {
        print_err("Error: Could not write distributions.");
        return;
    }

    #if VERBOSE > 0
    printf("\nDone writing distributions to HDF5 file.\n");
    #endif
}

void hdf5_dist_write_6D(dist_6D_offload_data* dist, real* hist, char* filename,
                        char* qid) {

    #if VERBOSE > 0
    printf("\nWriting distributions to HDF5 file...\n");
    #endif

    int abscissa_dim = 8;
    int ordinate_length = 1;

    int abscissa_n_slots[8];
    abscissa_n_slots[0] = dist->n_r;
    abscissa_n_slots[1] = dist->n_phi;
    abscissa_n_slots[2] = dist->n_z;
    abscissa_n_slots[3] = dist->n_vr;
    abscissa_n_slots[4] = dist->n_vphi;
    abscissa_n_slots[5] = dist->n_vz;
    abscissa_n_slots[6] = dist->n_time;
    abscissa_n_slots[7] = dist->n_q;

    double abscissa_min[8];
    abscissa_min[0] = dist->min_r;
    abscissa_min[1] = math_rad2deg(dist->min_phi);
    abscissa_min[2] = dist->min_z;
    abscissa_min[3] = dist->min_vr;
    abscissa_min[4] = dist->min_vphi;
    abscissa_min[5] = dist->min_vz;
    abscissa_min[6] = dist->min_time;
    abscissa_min[7] = dist->min_q;

    double abscissa_max[8];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = math_rad2deg(dist->max_phi);
    abscissa_max[2] = dist->max_z;
    abscissa_max[3] = dist->max_vr;
    abscissa_max[4] = dist->max_vphi;
    abscissa_max[5] = dist->max_vz;
    abscissa_max[6] = dist->max_time;
    abscissa_max[7] = dist->max_q;

    char* abscissa_names[] = { "R", "phi", "z", "vr", "vphi", "vz", "time",
                               "charge" };
    char* abscissa_units[] = { "m", "deg", "m", "m/s", "m/s", "m/s", "s", "e" };
    char* ordinate_names[] = { "density" };
    char* ordinate_units[] = { "s/m^6*deg*e" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run-XXXXXXXXXX/", qid, path);

    int retval;
    retval =  hdf5_histogram_write_uniform_double(
              filename,
              path,
              "R_phi_z_vr_vphi_vz_t_q",
              abscissa_dim,
              ordinate_length,
              abscissa_n_slots,
              abscissa_min,
              abscissa_max,
              abscissa_units,
              abscissa_names,
              ordinate_units,
              ordinate_names,
              hist);

    if(retval) {
        print_err("Error: Could not write distributions.");
        return;
    }

    #if VERBOSE > 0
    printf("\nDone writing distributions to HDF5 file.\n");
    #endif
}

void hdf5_dist_write_rho5D(dist_rho5D_offload_data* dist, real* hist, char* filename,
                        char* qid) {

    #if VERBOSE > 0
    printf("\nWriting distributions to HDF5 file...\n");
    #endif

    int abscissa_dim = 7;
    int ordinate_length = 1;

    int abscissa_n_slots[7];
    abscissa_n_slots[0] = dist->n_rho;
    abscissa_n_slots[1] = dist->n_pol;
    abscissa_n_slots[2] = dist->n_phi;
    abscissa_n_slots[3] = dist->n_vpara;
    abscissa_n_slots[4] = dist->n_vperp;
    abscissa_n_slots[5] = dist->n_time;
    abscissa_n_slots[6] = dist->n_q;

    double abscissa_min[7];
    abscissa_min[0] = dist->min_rho;
    abscissa_min[1] = dist->min_pol;
    abscissa_min[2] = dist->min_phi;
    abscissa_min[3] = dist->min_vpara;
    abscissa_min[4] = dist->min_vperp;
    abscissa_min[5] = dist->min_time;
    abscissa_min[6] = dist->min_q;

    double abscissa_max[7];
    abscissa_max[0] = dist->max_rho;
    abscissa_max[1] = dist->max_pol;
    abscissa_max[2] = dist->max_phi;
    abscissa_max[3] = dist->max_vpara;
    abscissa_max[4] = dist->max_vperp;
    abscissa_max[5] = dist->max_time;
    abscissa_max[6] = dist->max_q;

    char* abscissa_names[] = { "rho", "pol", "phi", "vpa", "vpe", "time", "charge" };
    char* abscissa_units[] = { " ", "deg", "deg", "m/s", "m/s", "s", "e" };
    char* ordinate_names[] = { "density" };
    char* ordinate_units[] = { " " };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run-XXXXXXXXXX/", qid, path);

    int retval;
    retval =  hdf5_histogram_write_uniform_double(
              filename,
              path,
              "rho_pol_phi_vpa_vpe_t_q",
              abscissa_dim,
              ordinate_length,
              abscissa_n_slots,
              abscissa_min,
              abscissa_max,
              abscissa_units,
              abscissa_names,
              ordinate_units,
              ordinate_names,
              hist);

    if(retval) {
        print_err("Error: Could not write distributions.");
        return;
    }

    #if VERBOSE > 0
    printf("\nDone writing distributions to HDF5 file.\n");
    #endif
}

void hdf5_dist_write_rho6D(dist_rho6D_offload_data* dist, real* hist, char* filename,
                        char* qid) {

    #if VERBOSE > 0
    printf("\nWriting distributions to HDF5 file...\n");
    #endif

    int abscissa_dim = 8;
    int ordinate_length = 1;

    int abscissa_n_slots[8];
    abscissa_n_slots[0] = dist->n_rho;
    abscissa_n_slots[1] = dist->n_pol;
    abscissa_n_slots[2] = dist->n_phi;
    abscissa_n_slots[3] = dist->n_vr;
    abscissa_n_slots[4] = dist->n_vphi;
    abscissa_n_slots[5] = dist->n_vz;
    abscissa_n_slots[6] = dist->n_time;
    abscissa_n_slots[7] = dist->n_q;

    double abscissa_min[8];
    abscissa_min[0] = dist->min_rho;
    abscissa_min[1] = dist->min_pol;
    abscissa_min[2] = dist->min_phi;
    abscissa_min[3] = dist->min_vr;
    abscissa_min[4] = dist->min_vphi;
    abscissa_min[5] = dist->min_vz;
    abscissa_min[6] = dist->min_time;
    abscissa_min[7] = dist->min_q;

    double abscissa_max[8];
    abscissa_max[0] = dist->max_rho;
    abscissa_max[1] = dist->max_pol;
    abscissa_max[2] = dist->max_phi;
    abscissa_max[3] = dist->max_vr;
    abscissa_max[4] = dist->max_vphi;
    abscissa_max[5] = dist->max_vz;
    abscissa_max[6] = dist->max_time;
    abscissa_max[7] = dist->max_q;

    char* abscissa_names[] = { "rho", "pol", "phi", "vr", "vphi", "vz", "time",
                               "charge" };
    char* abscissa_units[] = { " ", "deg", "deg", "m/s", "m/s", "m/s", "s", "e" };
    char* ordinate_names[] = { "density" };
    char* ordinate_units[] = { " " };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run-XXXXXXXXXX/", qid, path);

    int retval;
    retval =  hdf5_histogram_write_uniform_double(
              filename,
              path,
              "rho_pol_phi_vr_vphi_vz_t_q",
              abscissa_dim,
              ordinate_length,
              abscissa_n_slots,
              abscissa_min,
              abscissa_max,
              abscissa_units,
              abscissa_names,
              ordinate_units,
              ordinate_names,
              hist);

    if(retval) {
        print_err("Error: Could not write distributions.");
        return;
    }

    #if VERBOSE > 0
    printf("\nDone writing distributions to HDF5 file.\n");
    #endif
}
