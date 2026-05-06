/**
 * @file neural_network.h
 * @brief Neural network interpolation library
 *
 */
#include "../error.h"

// TODO: add the gpu magic words

/**
 * @brief Neural network parameters for 2D data
 *
 * Assumes a 2d input (x, y) that gets mapped to 1d output (f1). In practice used to map from (r, z) to (psi).
 */
typedef struct {
    int n_layers;           /**< number of layers */
    int* matrix_dimensions; /**< dimensions of each layer (n_layers-1) */
    real* weights;          /**< weights */
    real x_mean;            /**< mean of x */
    real y_mean;            /**< mean of y */
    real f1_mean;           /**< mean of f1 */
    real x_std;             /**< standard deviation of x */
    real y_std;             /**< standard deviation of y */
    real f1_std;            /**< standard deviation of f1 */
} neural2D_data;

/**
 * @brief Neural network parameters for 3D data
 *
 * Assumes a 3d input (x, y, z) that gets mapped to 3d output (f1, f2, f3). In practice used to map from (r, phi, z) to (Br, Bphi, Bz).
 */
typedef struct {
    int n_layers;           /**< number of layers */
    int* matrix_dimensions; /**< dimensions of each layer (n_layers-1) */
    real* weights;          /**< weights */
    real x_mean;            /**< mean of x */
    real y_mean;            /**< mean of y */
    real z_mean;            /**< mean of z */
    real x_std;             /**< standard deviation of x */
    real y_std;             /**< standard deviation of y */
    real z_std;             /**< standard deviation of z */
    real f1_mean;           /**< mean of f1 */
    real f2_mean;           /**< mean of f2 */
    real f3_mean;           /**< mean of f3 */
    real f1_std;            /**< standard deviation of f1 */
    real f2_std;            /**< standard deviation of f2 */
    real f3_std;            /**< standard deviation of f3 */
} neural3D_data;

/* Setuon functions */
int neural2Dsetup(neural2D_data* neural2D_data, real* weights, int n_layers,
    int* matrix_dimensions, real x_mean, real y_mean, real f1_mean, real x_std,
    real y_std, real f1_std);

int neural3Dsetup(neural3D_data* neural3D_data, real* weights, int n_layers,
    int* matrix_dimensions, real x_mean, real y_mean, real z_mean, real f1_mean,
    real f2_mean, real f3_mean, real x_std, real y_std, real z_std,
    real f1_std, real f2_std, real f3_std);

/* Evaluation functions */
a5err neural_network2Deval_f(real* f, neural2D_data* str, real x, real y);

a5err neural_network2Deval_df(real* f, real* partial_f, neural2D_data* str,
    real x, real y);

a5err neural_network3Deval_f(real* f, neural3D_data* str,
                         real x, real y, real z);

a5err neural_network3Deval_df(real* f, real* partial_f, neural3D_data* str,
                           real x, real y, real z);