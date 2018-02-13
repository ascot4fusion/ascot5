#include <stdlib.h>
#include "../ascot5.h"
#include "hdf5_histogram.h"
#include "hdf5_helpers.h"
#include "../distributions.h"

void hdf5_dist_write_rzvv(dist_rzvv_offload_data* dist, real* hist,
			  char* filename, char* qid) {

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
    abscissa_min[1] = dist->min_phi;
    abscissa_min[2] = dist->min_z;
    abscissa_min[3] = dist->min_vpara;
    abscissa_min[4] = dist->min_vperp;
    abscissa_min[5] = dist->min_time;
    abscissa_min[6] = dist->min_q;

    double abscissa_max[7];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = dist->max_phi;
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
    
    #if VERBOSE > 0
    printf("\nDone writing distributions to HDF5 file.\n");
    #endif
}
