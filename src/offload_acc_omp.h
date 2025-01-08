0;9;55M0;9;55m0;9;55M0;9;55m/**
 * @brief Macros for parallelizing the code for GPUs (with OpenACC or OpenMP) or
 * for CPUs with OpenMP
 */
#ifndef OFFLOAD_ACC_OMP_H
#define OFFLOAD_ACC_OMP_H

#include "ascot5.h"

/**
 * @brief Applies parallel execution to loops
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_PARALLEL_LOOP_ALL_LEVELS \
    str_pragma(omp target teams distribute parallel for simd)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_PARALLEL_LOOP_ALL_LEVELS str_pragma(acc parallel loop)
#else
#define GPU_PARALLEL_LOOP_ALL_LEVELS str_pragma(omp simd)
#endif

/**
 * @brief Applies parallel execution to loops with reduction
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(...) \
    str_pragma(omp target teams distribute parallel for simd \
    reduction(+:__VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(...) \
    str_pragma(acc parallel loop reduction(+:__VA_ARGS__))
#else
#define GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(...)
#endif

/**
 * @brief Declares routines or functions for target execution
 */
#if defined(GPU) && defined(_OPENMP)
#define DECLARE_TARGET str_pragma(omp declare target)
#elif defined(GPU) && defined(_OPENACC)
#define DECLARE_TARGET str_pragma(acc routine seq)
#else
#define DECLARE_TARGET
#endif

/**
 * @brief Ends the target declaration for routines or functions
 */
#if defined(GPU) && defined(_OPENMP)
#define DECLARE_TARGET_END str_pragma(omp end declare target)
#elif defined(GPU) && defined(_OPENACC)
#define DECLARE_TARGET_END
#else
#define DECLARE_TARGET_END
#endif

/**
 * @brief Maps variables to the target device
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_MAP_TO_DEVICE(...) \
    str_pragma(omp target enter data map (to: __VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_MAP_TO_DEVICE(...) str_pragma(acc enter data copyin (__VA_ARGS__))
#else
#define GPU_MAP_TO_DEVICE(...)
#endif

/**
 * @brief Updates host variables from the target device
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_UPDATE_FROM_DEVICE(...) \
    str_pragma(omp target update from (__VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_UPDATE_FROM_DEVICE(...) str_pragma(acc update host (__VA_ARGS__))
#else
#define GPU_UPDATE_FROM_DEVICE(...)
#endif

/**
 * @brief Maps variables from the target device back to the host
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_MAP_FROM_DEVICE(...) \
    str_pragma(omp target exit data map (from: __VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_MAP_FROM_DEVICE(...) \
    str_pragma(acc exit data copyout (__VA_ARGS__))
#else
#define GPU_MAP_FROM_DEVICE(...)
#endif

/**
 * @brief Deletes mappings of variables from the target device
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_MAP_DELETE_DEVICE(...) \
    str_pragma(omp target exit data (delete: __VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_MAP_DELETE_DEVICE(...) \
    str_pragma(acc exit data delete (__VA_ARGS__))
#else
#define GPU_MAP_DELETE_DEVICE(...)
#endif

/**
 * @brief Wrapper for "omp atomic"
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_ATOMIC str_pragma(omp atomic)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_ATOMIC str_pragma(omp atomic)
#else
#define GPU_ATOMIC str_pragma(omp atomic)
#endif

/**
 * @brief Declares SIMD loop
 */
#if defined(GPU) && defined(_OPENMP)
#define DECLARE_TARGET_SIMD
#elif defined(GPU) && defined(_OPENACC)
#define DECLARE_TARGET_SIMD
#else
#define DECLARE_TARGET_SIMD str_pragma(omp declare simd)
#endif

/**
 * @brief Declares SIMD loop and specifies variables which are uniform
 */
#if defined(GPU) && defined(_OPENMP)
#define DECLARE_TARGET_SIMD_UNIFORM(...)
#elif defined(GPU) && defined(_OPENACC)
#define DECLARE_TARGET_SIMD_UNIFORM(...) \
    str_pragma(omp declare simd uniform (__VA_ARGS__))
#else
#define DECLARE_TARGET_SIMD_UNIFORM(...)
#endif

/**
 * @brief Declares target or SIMD routines
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_DECLARE_TARGET_SIMD str_pragma(omp declare target)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_DECLARE_TARGET_SIMD str_pragma(acc routine seq)
#else
#define GPU_DECLARE_TARGET_SIMD str_pragma(omp declare simd)
#endif

/**
 * @brief Ends declaring target or SIMD routines
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_DECLARE_TARGET_SIMD_END str_pragma(omp end declare target)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_DECLARE_TARGET_SIMD_END
#else
#define GPU_DECLARE_TARGET_SIMD_END
#endif

/**
 * @brief Declares target or SIMD routines with uniform variables
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_DECLARE_TARGET_SIMD_UNIFORM(...) str_pragma(omp declare target)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_DECLARE_TARGET_SIMD_UNIFORM(...) str_pragma(acc routine seq)
#else
#define GPU_DECLARE_TARGET_SIMD_UNIFORM(...) \
    str_pragma(omp declare simd uniform (__VA_ARGS__))
#endif

/**
 * @brief Ends declaring target or SIMD routines (with uniform variables)
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_DECLARE_TARGET_SIMD_UNIFORM_END str_pragma(omp end declare target)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_DECLARE_TARGET_SIMD_UNIFORM_END
#else
#define GPU_DECLARE_TARGET_SIMD_UNIFORM_END
#endif

/**
 * @brief Makes a parallel region (CPU only)
 */
#if defined(GPU)
#define OMP_PARALLEL_CPU_ONLY
#else
#define OMP_PARALLEL_CPU_ONLY str_pragma(omp parallel)
#endif

/**
 * @brief Ensures the following loop within a parallel region is executed by
 *        a single thread only (OpenACC only)
 */
#if defined(GPU) && defined(_OPENACC)
#define GPU_SEQUENTIAL_LOOP str_pragma(acc loop seq)
#else
#define GPU_SEQUENTIAL_LOOP
#endif

/**
 * @brief Applies omp task for GPU
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_TASK \
    str_pragma(omp task)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_TASK str_pragma(omp task)
#else
#define GPU_OMP_TASK //
#endif

/**
 * @brief Applies omp task depend(in:) for GPU
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_TASK_DEPEND_IN(...) \
  str_pragma(omp task depend(in:__VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_TASK_DEPEND_IN(...) str_pragma(omp task depend(in:__VA_ARGS__))
#else
#define GPU_OMP_TASK_DEPEND_IN(...) //
#endif

/**
 * @brief Applies omp task depend(out:) for GPU
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_TASK_DEPEND_OUT(...) \
  str_pragma(omp task depend(out:__VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_TASK_DEPEND_OUT(...) str_pragma(omp task depend(out:__VA_ARGS__))
#else
#define GPU_OMP_TASK_DEPEND_OUT(...) //
#endif

/**
 * @brief Applies omp parallel for GPU
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_PARALLEL \
    str_pragma(omp parallel)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_PARALLEL str_pragma(omp parallel)
#else
#define GPU_OMP_PARALLEL //
#endif

/**
 * @brief Applies omp parallel do for GPU runs
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_PARALLEL_DO \
    str_pragma(omp parallel do)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_PARALLEL_DO str_pragma(omp parallel do )
#else
#define GPU_OMP_PARALLEL_DO str_pragma(omp simd)//
#endif

/**
 * @brief Applies omp parallel num_threads() for GPU
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_PARALLEL_NUM_THREADS(...) \
    str_pragma(omp parallel num_threads(__VA_ARGS__))
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_PARALLEL_NUM_THREADS(...) str_pragma(omp parallel num_threads(__VA_ARGS__))
#else
#define GPU_OMP_PARALLEL_NUM_THREADS(...) //
#endif

/**
 * @brief Applies omp master for GPU
 */
#if defined(GPU) && defined(_OPENMP)
#define GPU_OMP_MASTER \
    str_pragma(omp master)
#elif defined(GPU) && defined(_OPENACC)
#define GPU_OMP_MASTER str_pragma(omp master)
#else
#define GPU_OMP_MASTER //
#endif

/**
 * @brief Hints compiler that the following data is already present in the GPU
 *        (OpenACC only)
 */
#if defined(GPU) && defined(_OPENACC)
#define GPU_DATA_IS_MAPPED(...) str_pragma(acc data present(__VA_ARGS__))
#else
#define GPU_DATA_IS_MAPPED(...)
#endif

#ifdef _OPENACC
/**
 * @brief Number of gangs (OpenACC parallel execution units).
 */
#define NGANGS 1024

/**
 * @brief Number of workers per gang (OpenACC).
 */
#define NWORKERS 8

/**
 * @brief Number of vectors per worker (OpenACC).
 */
#define NVECTORS 32
#endif // _OPENACC

#endif
