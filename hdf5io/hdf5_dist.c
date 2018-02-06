#include <stdlib.h>
#include "../ascot5.h"
#include "hdf5_histogram.h"
#include "hdf5_helpers.h"
#include "../distributions.h"

void hdf5_dist_write_rzvv(dist_rzvv_offload_data* dist, real* hist,
			  char* filename, char* qid) {
    int abscissa_dim = 6;
    int ordinate_length = 1;

    /* transpose the histogram data for ascot4 rzVDist ordinate */
    double* ordinate = (double*) malloc(dist->n_r * dist->n_z * dist->n_vpara
                                        * dist->n_vperp * sizeof(double));
    int i, j, k, l;
    for(i = 0; i < dist->n_r; i++) {
        for(j = 0; j < dist->n_z; j++) {
            for(k = 0; k < dist->n_vpara; k++) {
                for(l = 0; l < dist->n_vperp; l++) {
                    ordinate[  l * (dist->n_vpara * dist->n_z * dist->n_r)
                             + k * (dist->n_z * dist->n_r)
                             + j * (dist->n_r)
                             + i] =
                    hist[  i * (dist->n_z * dist->n_vpara * dist->n_vperp)
                         + j * (dist->n_vpara * dist->n_vperp)
                         + k * (dist->n_vperp)
                         + l];
                }
            }
        }
    }

    int abscissa_n_slots[6];
    abscissa_n_slots[0] = dist->n_r;
    abscissa_n_slots[1] = dist->n_z;
    abscissa_n_slots[2] = dist->n_vpara;
    abscissa_n_slots[3] = dist->n_vperp;
    abscissa_n_slots[5] = 1;
    abscissa_n_slots[4] = 1;

    double abscissa_min[6];
    abscissa_min[0] = dist->min_r;
    abscissa_min[1] = dist->min_z;
    abscissa_min[2] = dist->min_vpara;
    abscissa_min[3] = dist->min_vperp;
    abscissa_min[4] = 0;
    abscissa_min[5] = 0.5;

    double abscissa_max[6];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = dist->max_z;
    abscissa_max[2] = dist->max_vpara;
    abscissa_max[3] = dist->max_vperp;
    abscissa_max[4] = 100;
    abscissa_max[5] = 1.5;

    char* abscissa_names[] = { "R", "z", "vpa", "vpe", "time", "species" };
    char* abscissa_units[] = { "m", "m", "m/s", "m/s", "s", "" };
    char* ordinate_names[] = { "density" };
    char* ordinate_units[] = { "s^2/m^5" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run-XXXXXXXXXX/", qid, path);
    
    int retval;
    retval =  hdf5_histogram_write_uniform_double(
		      filename,
		      path,
		      "r-phi-z-vpa-vpe-t-q-dist",
		      abscissa_dim,
		      ordinate_length,
		      abscissa_n_slots,
		      abscissa_min,
		      abscissa_max,
		      abscissa_units,
		      abscissa_names,
		      ordinate_units,
		      ordinate_names,
		      ordinate);
}
