/**
 * @file neural_network.c
 * @brief Neural network interpolation library
 *
 */

#include "neural_network.h"

/**
 * @brief Setup network parameters for 2D data
 *
 * @param neural2D_data data struct containing network weights for 2D
 * @param weights weights for the neural network
 * @param n_weights number of weights
 */
int neural2Dsetup(neural2D_data* neural2D_data, real* weights,
    int n_layers, int* matrix_dimensions, real x_mean, real y_mean, real f1_mean,
    real x_std, real y_std, real f1_std) {

    int n_weights = 1*1;

    for (int i = 0; i < n_layers; i++) {
        n_weights *= matrix_dimensions[i]*matrix_dimensions[i];
    }

    neural2D_data->n_layers = n_layers;

    neural2D_data->matrix_dimensions = (int*)malloc((n_layers-1) * sizeof(int));
    for (int i = 0; i < n_layers-1; i++) {
        neural2D_data->matrix_dimensions[i] = matrix_dimensions[i];
    }

    neural2D_data->weights = (real*)malloc(n_weights * sizeof(real));

    if(neural2D_data->weights == NULL) {
        return -1;
    }

    for (int i = 0; i < n_weights; i++) {
        neural2D_data->weights[i] = weights[i];
    }

    neural2D_data->x_mean = x_mean;
    neural2D_data->y_mean = y_mean;
    neural2D_data->f1_mean = f1_mean;
    neural2D_data->x_std = x_std;
    neural2D_data->y_std = y_std;
    neural2D_data->f1_std = f1_std;

    return 0;
}

/**
 * @brief Setup network parameters for 3D data
 *
 * @param neural3D_data data struct containing network weights for 3D
 * @param weights weights for the neural network
 * @param n_weights number of weights
 */
int neural3Dsetup(neural3D_data* neural3D_data, real* weights,
    int n_layers, int* matrix_dimensions, real x_mean, real y_mean, real z_mean,
    real f1_mean, real f2_mean, real f3_mean, real x_std, real y_std, real z_std,
    real f1_std, real f2_std, real f3_std) {

    int n_weights = 3*3;

    for (int i = 0; i < n_layers; i++) {
        n_weights *= matrix_dimensions[i]*matrix_dimensions[i];
    }

    neural3D_data->n_layers = n_layers;

    neural3D_data->matrix_dimensions = (int*)malloc((n_layers-1) * sizeof(int));
    for (int i = 0; i < n_layers-1; i++) {
        neural3D_data->matrix_dimensions[i] = matrix_dimensions[i];
    }

    neural3D_data->weights = (real*)malloc(n_weights * sizeof(real));

    if(neural3D_data->weights == NULL) {
        return -1;
    }

    for (int i = 0; i < n_weights; i++) {
        neural3D_data->weights[i] = weights[i];
    }

    neural3D_data->x_mean = x_mean;
    neural3D_data->y_mean = y_mean;
    neural3D_data->z_mean = z_mean;
    neural3D_data->x_std = x_std;
    neural3D_data->y_std = y_std;
    neural3D_data->z_std = z_std;
    neural3D_data->f1_mean = f1_mean;
    neural3D_data->f2_mean = f2_mean;
    neural3D_data->f3_mean = f3_mean;
    neural3D_data->f1_std = f1_std;
    neural3D_data->f2_std = f2_std;
    neural3D_data->f3_std = f3_std;

    return 0;
}

/**
 * @brief 2D neural network evaluation
 *
 * @param f variable in which to place the evaluated value
 * @param neural2D_data data struct containing network weights for 2D
 * @param x x-coordinate
 * @param y y-coordinate
 */
a5err neural_network2Deval_f(real* f, neural2D_data* neural2D_data, real x, real y) {
    // TODO: implement
}

/**
 * @brief 2D neural network evaluation
 *
 * @param f variable in which to place the evaluated value
 * @param neural2D_data data struct containing network weights for 2D
 * @param x x-coordinate
 * @param y y-coordinate
 */
a5err neural_network2Deval_df(real* f_df, neural2D_data* neural2D_data, real x, real y) {
    // TODO: implement
}

/**
 * @brief 3D neural network evaluation
 *
 * @param f variable in which to place the evaluated value
 * @param neural3D_data data struct containing network weights for 3D
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 */
a5err neural_network3Deval_f(real* f, neural3D_data* neural3D_data, real x, real y, real z) {
    // TODO: implement
}

/**
 * @brief 3D neural network evaluation
 *
 * @param f variable in which to place the evaluated value
 * @param neural3D_data data struct containing network weights for 3D
 * @param x x-coordinate
 * @param y y-coordinate
 * @param z z-coordinate
 */
a5err neural_network3Deval_df(real* f_df, neural3D_data* neural3D_data, real x, real y, real z) {
    // TODO: implement
}