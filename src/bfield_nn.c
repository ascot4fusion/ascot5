#include "bfield_nn.h"

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <cuda_runtime.h>

typedef enum {
    TORCHFORT_RESULT_SUCCESS = 0,
    TORCHFORT_RESULT_INVALID_USAGE = 1,
    TORCHFORT_RESULT_NOT_SUPPORTED = 2,
    TORCHFORT_RESULT_INTERNAL_ERROR = 3,
    TORCHFORT_RESULT_CUDA_ERROR = 4,
    TORCHFORT_RESULT_MPI_ERROR = 5,
    TORCHFORT_RESULT_NCCL_ERROR = 6
} torchfort_result_t;

typedef enum {
    TORCHFORT_FLOAT = -1,
    TORCHFORT_DOUBLE = -2,
    TORCHFORT_INT32 = -3,
    TORCHFORT_INT64 = -4
} torchfort_datatype_t;

extern torchfort_result_t torchfort_set_cudnn_benchmark(const bool flag);
extern torchfort_result_t torchfort_create_model(
    const char* name, const char* config_fname, int device);
extern torchfort_result_t torchfort_train(
    const char* name, void* input, size_t input_dim, int64_t* input_shape,
    void* label, size_t label_dim, int64_t* label_shape, void* loss_val,
    torchfort_datatype_t dtype, cudaStream_t stream);
extern torchfort_result_t torchfort_save_model(const char* name, const char* fname);
extern torchfort_result_t torchfort_save_checkpoint(const char* name, const char* checkpoint_dir);
extern torchfort_result_t torchfort_inference(
    const char* name, void* input, size_t input_dim, int64_t* input_shape,
    void* output, size_t output_dim, int64_t* output_shape,
    torchfort_datatype_t dtype, cudaStream_t stream);

static int initialized = 0;
static int train_step = 0;
static int validation_ready = 0;
static const char* model_name = "bfield_B_dB";

int bfield_nn_validation_ready(void) {
    return validation_ready;
}

int bfield_nn_get_step(void) {
    return train_step;
}

void bfield_nn_init(const char* config_file, int device) {
    torchfort_set_cudnn_benchmark(true);
    if (torchfort_create_model(model_name, config_file, device) != TORCHFORT_RESULT_SUCCESS) {
        fprintf(stderr, "[bfield_nn] create %s failed\n", model_name);
        exit(1);
    }
    initialized = 1;
    printf("[bfield_nn] initialized periodic SiLU B+dB model on GPU %d\n", device);
}

void bfield_nn_train_gpu(float* input, float* label, float* output, int n_mrk) {
    if (!initialized || n_mrk <= 0) return;

    int64_t input_shape[2] = {n_mrk, 3};
    int64_t label_shape[2] = {n_mrk, 13};
    int64_t output_shape[2] = {n_mrk, 13};
    float loss = 0.0f;
    validation_ready = 0;

    torchfort_result_t ret = torchfort_train(
        model_name, input, 2, input_shape, label, 2, label_shape,
        &loss, TORCHFORT_FLOAT, 0);
    if (ret != TORCHFORT_RESULT_SUCCESS) {
        fprintf(stderr, "[bfield_nn] train failed at step %d\n", train_step);
    }
    if (train_step % 10 == 0) {
        printf("[bfield_nn] step %d combined_loss %.6e\n", train_step, loss);
    }

    if ((train_step + 1) % 1000 == 0) {
        ret = torchfort_inference(
            model_name, input, 2, input_shape, output, 2, output_shape,
            TORCHFORT_FLOAT, 0);
        if (ret != TORCHFORT_RESULT_SUCCESS) {
            fprintf(stderr, "[bfield_nn] inference failed at step %d\n", train_step + 1);
        } else {
            validation_ready = 1;
        }
    }
    train_step++;
}

void bfield_nn_finalize(void) {
    if (initialized) {
        torchfort_save_model(model_name, "bfield_B_periodic_dB_trained.pt");
        torchfort_save_checkpoint(model_name, "bfield_B_periodic_dB_checkpoint");
    }
    initialized = 0;
}
