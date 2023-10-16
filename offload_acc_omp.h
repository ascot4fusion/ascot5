#ifndef offload_acc_omp.h
#define offload_acc_omp.h

#define STRINGIFY(X) #X
#define MY_PRAGMA(X) _Pragma(STRINGIFY(X))
#define GPU_PARALLEL_LOOP_ALL_LEVELS MY_PRAGMA(omp simd) 
#define DECLARE_TARGET     
#define DECLARE_TARGET_END 
#define DECLARE_TARGET_COPYIN
#define GPU_MAP_TO_DEVICE
//#define OMP_L1 MY_PRAGMA(omp distribute parallel for) 
#ifdef _OPENMP
#define GPU_PARALLEL_LOOP_ALL_LEVELS MY_PRAGMA(omp target teams distribute parallel for)
//#define GPU_PARALLEL_LOOP_ALL_LEVELS 
#define DECLARE_TARGET     MY_PRAGMA(omp declare target)
#define DECLARE_TARGET_END MY_PRAGMA(omp end declare target)
#endif
#ifdef _OPENACC
#warning "OpenACC"
#define NGANGS 1024
#define NWORKERS 8
#define NVECTORS 32
#define GPU_PARALLEL_LOOP_ALL_LEVELS MY_PRAGMA(acc parallel loop)
//#define OMP_L0 MY_PRAGMA()
#define DECLARE_TARGET     MY_PRAGMA(acc routine seq)
#define DECLARE_TARGET_END MY_PRAGMA()
#define DECLARE_TARGET_COPYIN     MY_PRAGMA(acc declare target copyin)
#define GPU_MAP_TO_DEVICE MY_PRAGMA(acc enter data copyin)
#endif

#endif 
