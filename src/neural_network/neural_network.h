/**
 * @file neural_network.h
 * @brief Neural network interpolation library
 *
 */
#include "../error.h"

// TODO: add the gpu magic words

/**
 * @brief Neural network parameters for 2D data
 */
typedef struct {
    int n;     /**< number                        */
} neural2D_data;

/**
 * @brief Neural network parameters for 3D data
 */
typedef struct {
    int n;     /**< number                        */
} neural3D_data;


a5err neural_network2Deval_f(real* f, neural2D_data* str, real x, real y);

a5err neural_network2Deval_df(real* f_df, neural2D_data* str, real x, real y);

a5err neural_network3Deval_f(real* f, neural3D_data* str,
                         real x, real y, real z);

a5err neural_network3Deval_df(real* f_df, neural3D_data* str,
                           real x, real y, real z);