#ifndef OFFLOAD_ACC_OMP_H
#define OFFLOAD_ACC_OMP_H

#define STRINGIFY(X) #X
#define MY_PRAGMA(X) _Pragma(STRINGIFY(X))

#ifndef GPU

#define GPU_PARALLEL_LOOP_ALL_LEVELS MY_PRAGMA(omp simd)
#define GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(x ...)
#define DECLARE_TARGET
#define DECLARE_TARGET_END
#define DECLARE_TARGET_COPYIN
#define GPU_MAP_TO_DEVICE(x ...)
#define GPU_UPDATE_FROM_DEVICE(x ...)
#define GPU_MAP_FROM_DEVICE(x ...)
#define GPU_MAP_DELETE_DEVICE(x ...)
#define GPU_ATOMIC MY_PRAGMA(omp atomic)

#else

#ifdef _OPENMP
#define GPU_PARALLEL_LOOP_ALL_LEVELS MY_PRAGMA(omp target teams distribute parallel for simd )
#define GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(x ...) MY_PRAGMA(omp target teams distribute parallel for simd reduction(+:x))
//#define GPU_PARALLEL_LOOP_ALL_LEVELS 
#define DECLARE_TARGET     MY_PRAGMA(omp declare target)
#define DECLARE_TARGET_END MY_PRAGMA(omp end declare target)
#define GPU_MAP_TO_DEVICE(x ...) MY_PRAGMA(omp target enter data map (to: x))
#define GPU_UPDATE_FROM_DEVICE(x ...) MY_PRAGMA(omp target update from ( x))
#define GPU_MAP_FROM_DEVICE(x ...) MY_PRAGMA(omp target exit data map (from: x))
#define GPU_MAP_DELETE_DEVICE(x ...)  MY_PRAGMA(omp target exit data (delete: x))
#define GPU_ATOMIC MY_PRAGMA(acc atomic)
#endif
#ifdef _OPENACC
//#warning "OpenACC"
#define NGANGS 1024
#define NWORKERS 8
#define NVECTORS 32
#define GPU_PARALLEL_LOOP_ALL_LEVELS MY_PRAGMA(acc parallel loop)
#define GPU_PARALLEL_LOOP_ALL_LEVELS_REDUCTION(x ...) MY_PRAGMA(acc parallel loop reduction(+:x))
//#define OMP_L0 MY_PRAGMA()
#define DECLARE_TARGET     MY_PRAGMA(acc routine seq)
#define DECLARE_TARGET_END
#define DECLARE_TARGET_COPYIN     MY_PRAGMA(acc declare target copyin)
#define GPU_MAP_TO_DEVICE(x ...)  MY_PRAGMA(acc enter data copyin ( x))
#define GPU_UPDATE_FROM_DEVICE(x ...) MY_PRAGMA(acc update host ( x))
#define GPU_MAP_FROM_DEVICE(x ...)  MY_PRAGMA(acc exit data copyout ( x))
#define GPU_MAP_DELETE_DEVICE(x ...)  MY_PRAGMA(acc exit data delete ( x))
#define GPU_ATOMIC MY_PRAGMA(omp atomic)
#endif

#endif

#endif 
