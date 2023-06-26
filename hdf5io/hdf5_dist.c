/**
 * @file hdf5_dist.c
 * @brief Distribution HDF5 IO
 */
#include <stdlib.h>
#include <hdf5.h>
#include "../ascot5.h"
#include "../diag/dist_5D.h"
#include "../diag/dist_6D.h"
#include "../diag/dist_rho5D.h"
#include "../diag/dist_rho6D.h"
#include "../diag/dist_com.h"
#include "../math.h"
#include "hdf5_histogram.h"
#include "hdf5_helpers.h"
#include "hdf5_dist.h"

/**
 * @brief Write 5D distribution to an existing result group
 *
 * @param f HDF5 file id
 * @param qid run QID where distribution is written
 * @param dist pointer to distribution data struct
 * @param hist pointer to distribution data
 */
int hdf5_dist_write_5D(hid_t f, char* qid, dist_5D_offload_data* dist,
                       real* hist) {

    int abscissa_dim = 7;
    int ordinate_dim = 1;

    int abscissa_n_slots[7];
    abscissa_n_slots[0] = dist->n_r;
    abscissa_n_slots[1] = dist->n_phi;
    abscissa_n_slots[2] = dist->n_z;
    abscissa_n_slots[3] = dist->n_ppara;
    abscissa_n_slots[4] = dist->n_pperp;
    abscissa_n_slots[5] = dist->n_time;
    abscissa_n_slots[6] = dist->n_q;

    double abscissa_min[7];
    abscissa_min[0] = dist->min_r;
    abscissa_min[1] = math_rad2deg(dist->min_phi);
    abscissa_min[2] = dist->min_z;
    abscissa_min[3] = dist->min_ppara;
    abscissa_min[4] = dist->min_pperp;
    abscissa_min[5] = dist->min_time;
    abscissa_min[6] = dist->min_q;

    double abscissa_max[7];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = math_rad2deg(dist->max_phi);
    abscissa_max[2] = dist->max_z;
    abscissa_max[3] = dist->max_ppara;
    abscissa_max[4] = dist->max_pperp;
    abscissa_max[5] = dist->max_time;
    abscissa_max[6] = dist->max_q;

    char* abscissa_names[] = { "r", "phi", "z", "ppar", "pperp", "time",
                               "charge" };
    char* abscissa_units[] = { "m", "deg", "m", "kg*m/s", "kg*m/s", "s", "e" };
    char* ordinate_names[] = { "distribution" };
    char* ordinate_units[] = { "s/m^5*kg^2*deg*e" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run_XXXXXXXXXX/dist5d", qid, path);

    int retval = hdf5_histogram_write_uniform_double(f, path,
                                                     abscissa_dim, ordinate_dim,
                                                     abscissa_n_slots,
                                                     abscissa_min,
                                                     abscissa_max,
                                                     abscissa_units,
                                                     abscissa_names,
                                                     ordinate_units,
                                                     ordinate_names,
                                                     hist);

    return retval;
}

/**
 * @brief Write 6D distribution to an existing result group
 *
 * @param f HDF5 file id
 * @param qid run QID where distribution is written
 * @param dist pointer to distribution data struct
 * @param hist pointer to distribution data
 */
int hdf5_dist_write_6D(hid_t f, char* qid, dist_6D_offload_data* dist,
                       real* hist) {

    int abscissa_dim = 8;
    int ordinate_dim = 1;

    int abscissa_n_slots[8];
    abscissa_n_slots[0] = dist->n_r;
    abscissa_n_slots[1] = dist->n_phi;
    abscissa_n_slots[2] = dist->n_z;
    abscissa_n_slots[3] = dist->n_pr;
    abscissa_n_slots[4] = dist->n_pphi;
    abscissa_n_slots[5] = dist->n_pz;
    abscissa_n_slots[6] = dist->n_time;
    abscissa_n_slots[7] = dist->n_q;

    double abscissa_min[8];
    abscissa_min[0] = dist->min_r;
    abscissa_min[1] = math_rad2deg(dist->min_phi);
    abscissa_min[2] = dist->min_z;
    abscissa_min[3] = dist->min_pr;
    abscissa_min[4] = dist->min_pphi;
    abscissa_min[5] = dist->min_pz;
    abscissa_min[6] = dist->min_time;
    abscissa_min[7] = dist->min_q;

    double abscissa_max[8];
    abscissa_max[0] = dist->max_r;
    abscissa_max[1] = math_rad2deg(dist->max_phi);
    abscissa_max[2] = dist->max_z;
    abscissa_max[3] = dist->max_pr;
    abscissa_max[4] = dist->max_pphi;
    abscissa_max[5] = dist->max_pz;
    abscissa_max[6] = dist->max_time;
    abscissa_max[7] = dist->max_q;

    char* abscissa_names[] = { "r", "phi", "z", "pr", "pphi", "pz", "time",
                               "charge" };
    char* abscissa_units[] = { "m", "deg", "m", "kg m/s", "kg m/s", "kg m/s",
                               "s", "e" };
    char* ordinate_names[] = { "distribution" };
    char* ordinate_units[] = { "s/m^6*kg^3*deg*e" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run_XXXXXXXXXX/dist6d", qid, path);

    int retval = hdf5_histogram_write_uniform_double(f, path,
                                                     abscissa_dim, ordinate_dim,
                                                     abscissa_n_slots,
                                                     abscissa_min,
                                                     abscissa_max,
                                                     abscissa_units,
                                                     abscissa_names,
                                                     ordinate_units,
                                                     ordinate_names,
                                                     hist);

    return retval;
}

/**
 * @brief Write rho 5D distribution to an existing result group
 *
 * @param f HDF5 file id
 * @param qid run QID where distribution is written
 * @param dist pointer to distribution data struct
 * @param hist pointer to distribution data
 */
int hdf5_dist_write_rho5D(hid_t f, char* qid, dist_rho5D_offload_data* dist,
                          real* hist) {

    int abscissa_dim = 7;
    int ordinate_dim = 1;

    int abscissa_n_slots[7];
    abscissa_n_slots[0] = dist->n_rho;
    abscissa_n_slots[1] = dist->n_theta;
    abscissa_n_slots[2] = dist->n_phi;
    abscissa_n_slots[3] = dist->n_ppara;
    abscissa_n_slots[4] = dist->n_pperp;
    abscissa_n_slots[5] = dist->n_time;
    abscissa_n_slots[6] = dist->n_q;

    double abscissa_min[7];
    abscissa_min[0] = dist->min_rho;
    abscissa_min[1] = math_rad2deg(dist->min_theta);
    abscissa_min[2] = math_rad2deg(dist->min_phi);
    abscissa_min[3] = dist->min_ppara;
    abscissa_min[4] = dist->min_pperp;
    abscissa_min[5] = dist->min_time;
    abscissa_min[6] = dist->min_q;

    double abscissa_max[7];
    abscissa_max[0] = dist->max_rho;
    abscissa_max[1] = math_rad2deg(dist->max_theta);
    abscissa_max[2] = math_rad2deg(dist->max_phi);
    abscissa_max[3] = dist->max_ppara;
    abscissa_max[4] = dist->max_pperp;
    abscissa_max[5] = dist->max_time;
    abscissa_max[6] = dist->max_q;

    char* abscissa_names[] = { "rho", "theta", "phi", "ppar", "pperp", "time",
                               "charge" };
    char* abscissa_units[] = { "1", "deg", "deg", "kg m/s", "kg m/s",
                               "s", "e" };
    char* ordinate_names[] = { "distribution" };
    char* ordinate_units[] = { "s/m^2*kg^2*deg^2e" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run_XXXXXXXXXX/distrho5d", qid, path);

    int retval = hdf5_histogram_write_uniform_double(f, path,
                                                     abscissa_dim, ordinate_dim,
                                                     abscissa_n_slots,
                                                     abscissa_min,
                                                     abscissa_max,
                                                     abscissa_units,
                                                     abscissa_names,
                                                     ordinate_units,
                                                     ordinate_names,
                                                     hist);

    return retval;
}

/**
 * @brief Write rho 6D distribution to an existing result group
 *
 * @param f HDF5 file id
 * @param qid run QID where distribution is written
 * @param dist pointer to distribution data struct
 * @param hist pointer to distribution data
 */
int hdf5_dist_write_rho6D(hid_t f, char* qid, dist_rho6D_offload_data* dist,
                          real* hist) {

    int abscissa_dim = 8;
    int ordinate_dim = 1;

    int abscissa_n_slots[8];
    abscissa_n_slots[0] = dist->n_rho;
    abscissa_n_slots[1] = dist->n_theta;
    abscissa_n_slots[2] = dist->n_phi;
    abscissa_n_slots[3] = dist->n_pr;
    abscissa_n_slots[4] = dist->n_pphi;
    abscissa_n_slots[5] = dist->n_pz;
    abscissa_n_slots[6] = dist->n_time;
    abscissa_n_slots[7] = dist->n_q;

    double abscissa_min[8];
    abscissa_min[0] = dist->min_rho;
    abscissa_min[1] = math_rad2deg(dist->min_theta);
    abscissa_min[2] = math_rad2deg(dist->min_phi);
    abscissa_min[3] = dist->min_pr;
    abscissa_min[4] = dist->min_pphi;
    abscissa_min[5] = dist->min_pz;
    abscissa_min[6] = dist->min_time;
    abscissa_min[7] = dist->min_q;

    double abscissa_max[8];
    abscissa_max[0] = dist->max_rho;
    abscissa_max[1] = math_rad2deg(dist->max_theta);
    abscissa_max[2] = math_rad2deg(dist->max_phi);
    abscissa_max[3] = dist->max_pr;
    abscissa_max[4] = dist->max_pphi;
    abscissa_max[5] = dist->max_pz;
    abscissa_max[6] = dist->max_time;
    abscissa_max[7] = dist->max_q;

    char* abscissa_names[] = { "rho", "theta", "phi", "pr", "pphi", "pz",
                               "time", "charge" };
    char* abscissa_units[] = { "1", "deg", "deg", "kg m/s", "kg m/s", "kg m/s",
                               "s", "e" };
    char* ordinate_names[] = { "distribution" };
    char* ordinate_units[] = { "s^2/m^3*kg^3*deg^2*e" };

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run_XXXXXXXXXX/distrho6d", qid, path);

    int retval = hdf5_histogram_write_uniform_double(f, path,
                                                     abscissa_dim, ordinate_dim,
                                                     abscissa_n_slots,
                                                     abscissa_min,
                                                     abscissa_max,
                                                     abscissa_units,
                                                     abscissa_names,
                                                     ordinate_units,
                                                     ordinate_names,
                                                     hist);

    return retval;
}

/**
 * @brief Write constants-of-motion distribution to an existing result group
 *
 * @param f HDF5 file id
 * @param qid run QID where distribution is written
 * @param dist pointer to distribution data struct
 * @param hist pointer to distribution data
 */
int hdf5_dist_write_COM(hid_t f, char* qid, dist_COM_offload_data* dist,
                        real* hist) {

    int abscissa_dim = 3;
    int ordinate_dim = 1;

    int abscissa_n_slots[3];
    abscissa_n_slots[0] = dist->n_mu;
    abscissa_n_slots[1] = dist->n_Ekin;
    abscissa_n_slots[2] = dist->n_Ptor;

    double abscissa_min[3];
    abscissa_min[0] = dist->min_mu;
    abscissa_min[1] = dist->min_Ekin;
    abscissa_min[2] = dist->min_Ptor;

    double abscissa_max[3];
    abscissa_max[0] = dist->max_mu;
    abscissa_max[1] = dist->max_Ekin;
    abscissa_max[2] = dist->max_Ptor;

    char* abscissa_names[] = { "mu", "ekin", "ptor"};
    char* abscissa_units[] = { "J/T", "J", "J s"};
    char* ordinate_names[] = { "distribution" };
    char* ordinate_units[] = { "T / J**3 s"};

    /* Create a group for this distribution and write the data in it */
    char path[256];
    hdf5_generate_qid_path("/results/run_XXXXXXXXXX/distcom", qid, path);

    int retval = hdf5_histogram_write_uniform_double(f, path,
                                                     abscissa_dim, ordinate_dim,
                                                     abscissa_n_slots,
                                                     abscissa_min,
                                                     abscissa_max,
                                                     abscissa_units,
                                                     abscissa_names,
                                                     ordinate_units,
                                                     ordinate_names,
                                                     hist);

    return retval;
}
