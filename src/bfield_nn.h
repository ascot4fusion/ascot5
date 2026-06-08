#ifndef BFIELD_NN_H
#define BFIELD_NN_H

void bfield_nn_init(const char* config_file, int device);

void bfield_nn_train_gpu(
    float* nn_input,
    float* nn_label,
    int n_mrk
);

void bfield_nn_finalize(void);

#endif
