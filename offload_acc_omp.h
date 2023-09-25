#define STRINGIFY(X) #X
#define MY_PRAGMA(X) _Pragma(STRINGIFY(X))
#define OMP_L0 
#define OMP_L1 
#define OMP_L2 
#define DECLARE_TARGET     
#define DECLARE_TARGET_END 
#define DECLARE_TARGET_COPYIN 
//#define OMP_L1 MY_PRAGMA(omp distribute parallel for) 
#ifdef _OPENMP
#define OMP_L1 MY_PRAGMA(omp distribute)
#define OMP_L2 MY_PRAGMA(omp parallel for simd)
#define DECLARE_TARGET     MY_PRAGMA(omp declare target)
#define DECLARE_TARGET_END MY_PRAGMA(omp end declare target)
//#define OMP_L2 MY_PRAGMA(omp parallel for)
//#define OMP_L2 MY_PRAGMA(omp simd)
//#define OMP_L2 MY_PRAGMA()
#endif
#ifdef _OPENACC
#warning "OpenACC"
#define NGANGS 1024
#define NWORKERS 8
#define NVECTORS 32
#define OMP_L1 MY_PRAGMA()
//#define OMP_L2 MY_PRAGMA(acc loop worker vector)
#define OMP_L2 MY_PRAGMA()
#define OMP_L0 MY_PRAGMA()
#define OMP_L0_DBG MY_PRAGMA(acc parallel loop)
//#define OMP_L0 MY_PRAGMA()
#define DECLARE_TARGET     MY_PRAGMA(acc routine seq)
#define DECLARE_TARGET_END MY_PRAGMA()
#define DECLARE_TARGET_COPYIN     MY_PRAGMA(acc declare target copyin)
#endif

