#ifndef BFIELD_NN_H
#define BFIELD_NN_H

void bfield_nn_init(const char* config_B, int device);

void bfield_nn_train_gpu(
    float* nn_input,
    float* nn_label_B,
    float* nn_output_B,
    int n_mrk
);

void bfield_nn_finalize(void);

int bfield_nn_get_step(void);
int bfield_nn_validation_ready(void);

#endif
