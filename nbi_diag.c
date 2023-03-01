#include <stdlib.h>
#include <hdf5.h>
#include "ascot5.h"
#include "nbi_diag.h"
#include "hdf5io/hdf5_histogram.h"
#include "hdf5io/hdf5_helpers.h"

void nbi_diag_init(nbi_diag_data* data, int n_x, real min_x, real max_x,
                   int n_y, real min_y, real max_y, int n_z, real min_z,
                   real max_z, real full_energy) {
    data->n_x = n_x;
    data->min_x = min_x;
    data->max_x = max_x;
    data->n_y = n_y;
    data->min_y = min_y;
    data->max_y = max_y;
    data->n_z = n_z;
    data->min_z = min_z;
    data->max_z = max_z;
    data->full_energy = full_energy;

    data->hist_neutral = (real*) malloc(3 * n_x * n_y * n_z * sizeof(real));
}

void nbi_diag_free(nbi_diag_data* data) {
    free(data->hist_neutral);
    data->hist_neutral = NULL;
}

void nbi_diag_update_dist(nbi_diag_data* data, real x, real y, real z,
                          real energy, real weight) {
    int i_energy;
    if(energy > 0.75*data->full_energy) {
        i_energy = 0;
    }
    else if(energy > 0.4*data->full_energy) {
        i_energy = 1;
    }
    else {
        i_energy = 2;
    }

    int i_x = floor((x - data->min_x)
                    / ((data->max_x - data-> min_x)/data->n_x));
    int i_y = floor((y - data->min_y)
                    / ((data->max_y - data-> min_y)/data->n_y));
    int i_z = floor((z - data->min_z)
                    / ((data->max_z - data-> min_z)/data->n_z));

    if(   i_x >= 0 && i_x < data->n_x && i_y >= 0 && i_y < data->n_y
       && i_z >= 0 && i_z < data->n_z) {
        #pragma omp atomic
        data->hist_neutral[i_energy*data->n_x*data->n_y*data->n_z
                     + i_x*data->n_y*data->n_z + i_y*data->n_z + i_z] += weight;
    }
}

void nbi_diag_normalize_dist(nbi_diag_data* data, real weight) {
    for(int i = 0; i < 3*data->n_x*data->n_y*data->n_z; i++) {
        data->hist_neutral[i] *= weight;
    }
}

void nbi_diag_write(char* fn, nbi_diag_data* data) {
    hid_t f = hdf5_create(fn);
    hdf5_close(f);

    f = hdf5_open(fn);

    if(hdf5_find_group(f, "/bbnbi_results/")) {
        hdf5_create_group(f, "/bbnbi_results/");
    }

    hid_t grp = H5Gcreate2(f, "/bbnbi_results/run_1234567890", H5P_DEFAULT,
                           H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(grp);

    int abscissa_dim = 4;
    int ordinate_dim = 1;
    int abscissa_nslots[] = {3, data->n_x, data->n_y, data->n_z};
    double abscissa_min[] = {0.5, data->min_x, data->min_y, data->min_z};
    double abscissa_max[] = {3.5, data->max_x, data->max_y, data->max_z};
    char* abscissa_units[] = {"", "m", "m", "m"};
    char* abscissa_names[] = {"fraction", "x", "y", "z"};
    char* ordinate_units[] = {"1/m^3"};
    char* ordinate_names[] = {"neutral density"};

    hdf5_histogram_write_uniform_double(f, "/bbnbi_results/run_1234567890/neutraldist", abscissa_dim, ordinate_dim, abscissa_nslots, abscissa_min, abscissa_max, abscissa_units, abscissa_names, ordinate_units, ordinate_names, data->hist_neutral);
}
