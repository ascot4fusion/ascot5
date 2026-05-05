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
int neural2Dsetup(neural2D_data* neural2D_data, real* weights, int n_weights) {

    neural2D_data->weights = (real*)malloc(n_weights * sizeof(real));

    if(neural2D_data->weights == NULL) {
        return -1;
    }

    for (int i = 0; i < n_weights; i++) {
        neural2D_data->weights[i] = weights[i];
    }

    return 0;
}

/**
 * @brief Setup network parameters for 3D data
 *
 * @param neural3D_data data struct containing network weights for 3D
 * @param weights weights for the neural network
 * @param n_weights number of weights
 */
int neural3Dsetup(neural3D_data* neural3D_data, real* weights, int n_weights) {

    neural3D_data->weights = (real*)malloc(n_weights * sizeof(real));

    if(neural3D_data->weights == NULL) {
        return -1;
    }

    for (int i = 0; i < n_weights; i++) {
        neural3D_data->weights[i] = weights[i];
    }

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