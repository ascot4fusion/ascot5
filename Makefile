CC=mpicc 

ascot5_main:
	$(MAKE) -C src ascot5_main
	mkdir -p build
	mv src/ascot5_main build/ascot5_main

libascot:
	$(MAKE) -C src libascot
	mkdir -p build
	mv src/libascot.so build/libascot.so

bbnbi5:
	$(MAKE) -C src bbnbi5
	mkdir -p build
	mv src/bbnbi5 build/bbnbi5

ifdef GPU
        DEFINES+=-DGPU
endif

ifdef VERBOSE
	DEFINES+=-DVERBOSE=$(VERBOSE)
else
	DEFINES+=-DVERBOSE=1
endif

ifeq ($(SINGLEPRECISION),1)
	DEFINES+=-DSINGLEPRECISION
endif

ifeq ($(MPI),1)
	DEFINES+=-DMPI
	#CC=h5pcc
endif

ifeq ($(RANDOM),MKL)
	DEFINES+=-DRANDOM_MKL
	CFLAGS+=-mkl
	ifdef RANDOM_MKL_RNG
		DEFINES+=-DRANDOM_MKL_RNG=$(RANDOM_MKL_RNG)
	endif
else ifeq ($(RANDOM),GSL)
	DEFINES+=-DRANDOM_GSL
	CFLAGS+=-lgsl -lgslcblas
else ifeq ($(RANDOM),LCG)
	DEFINES+=-DRANDOM_LCG
endif

ifneq ($(CC),h5cc)
	ifneq ($(CC),h5pcc)
		CFLAGS+=-L$(HDF5_ROOT)/lib -lhdf5 -lhdf5_hl -I$(HDF5_ROOT)/include
	endif
endif

ifdef TARGET
    DEFINES+=-DTARGET=$(TARGET)
endif

ifeq ($(OMP),1)
# uncomment this to compile for CPU only without SIMD
        #CFLAGS+=-fopenmp -foffload=disable
# uncomment this to compile for CPU only with SIMD
        #CFLAGS+=-fopenmp -DSIMD -foffload=disable
# uncomment this to compile for GPU 
       CFLAGS+=-fopenmp 
ifeq ($(GPU), 1)
	CFLAGS+=-foffload="nvptx-none=-latomic"
        CFLAGS+=-foffload="-mgomp" 
	CFLAGS+="-g"
	CFLAGS+=-foffload="-lm"
else
	CFLAGS+=-fopenmp -DSIMD -foffload=disable -DNSIMD=8
endif
endif

CFLAGS+= -lm -fPIC -std=c11 $(DEFINES) $(FLAGS) #-DFULLMCCC

# Write CFLAGS and CC to a file to be included into output
$(shell echo "#define CFLAGS " $(CFLAGS) > compiler_flags.h)
$(shell echo "#define CC " $(CC) >> compiler_flags.h)

SIMDIR = simulate/
SIMHEADERS = $(wildcard $(SIMDIR)simulate*.h $(SIMDIR)atomic.h)
SIMOBJS = $(patsubst %.c,%.o,$(wildcard $(SIMDIR)simulate*.c) $(SIMDIR)atomic.c)

STEPDIR = $(SIMDIR)step/
STEPHEADERS = $(wildcard $(STEPDIR)step*.h)
STEPOBJS = $(patsubst %.c,%.o,$(wildcard $(STEPDIR)step*.c))

MCCCDIR = $(SIMDIR)mccc/
MCCCHEADERS = $(wildcard $(MCCCDIR)mccc*.h)
MCCCOBJS = $(patsubst %.c,%.o,$(wildcard $(MCCCDIR)mccc*.c))

DIAGDIR = diag/
DIAGHEADERS = $(wildcard $(DIAGDIR)diag*.h $(DIAGDIR)dist*.h)
DIAGOBJS = $(patsubst %.c,%.o,$(wildcard $(DIAGDIR)diag*.c $(DIAGDIR)dist*.c))

BFDIR = Bfield/
BFHEADERS = $(wildcard $(BFDIR)B_*.h)
BFOBJS = $(patsubst %.c,%.o,$(wildcard $(BFDIR)B_*.c))

EFDIR = Efield/
EFHEADERS = $(wildcard $(EFDIR)E_*.h)
EFOBJS = $(patsubst %.c,%.o,$(wildcard $(EFDIR)E_*.c))

PLSDIR = plasma/
PLSHEADERS =  $(wildcard $(PLSDIR)plasma_*.h)
PLSOBJS = $(patsubst %.c,%.o,$(wildcard $(PLSDIR)plasma_*.c))

WALLDIR = wall/
WALLHEADERS = $(wildcard $(WALLDIR)wall_*.h)
WALLOBJS = $(patsubst %.c,%.o,$(wildcard $(WALLDIR)wall_*.c))

HDF5IODIR = hdf5io/
HDF5IOHEADERS = $(wildcard $(HDF5IODIR)hdf5*.h)
HDF5IOOBJS = $(patsubst %.c,%.o,$(wildcard $(HDF5IODIR)hdf5*.c))

N0DIR = neutral/
N0HEADERS =  $(wildcard $(N0DIR)N0_*.h)
N0OBJS = $(patsubst %.c,%.o,$(wildcard $(N0DIR)N0_*.c))

MHDDIR = mhd/
MHDHEADERS =  $(wildcard $(MHDDIR)mhd_*.h)
MHDOBJS = $(patsubst %.c,%.o,$(wildcard $(MHDDIR)mhd_*.c))

ASIGMADIR = asigma/
ASIGMAHEADERS = $(wildcard $(ASIGMADIR)asigma_*.h)
ASIGMAOBJS = $(patsubst %.c,%.o,$(wildcard $(ASIGMADIR)asigma_*.c))

LINTDIR = linint/
LINTHEADERS =  $(wildcard $(LINTDIR)linint*.h)
LINTOBJS = $(patsubst %.c,%.o,$(wildcard $(LINTDIR)linint*.c))

SPLINEDIR = spline/
SPLINEHEADERS  = $(wildcard $(SPLINEDIR)spline.h $(SPLINEDIR)interp.h)
SPLINEOBJS  = $(patsubst %.c,%.o,$(wildcard $(SPLINEDIR)spline*.c \
						$(SPLINEDIR)interp*.c))

UTESTDIR = unit_tests/
DOCDIR = doc/

HEADERS=ascot5.h math.h consts.h list.h octree.h physlib.h error.h \
	$(DIAGHEADERS) $(BFHEADERS) $(EFHEADERS) $(WALLHEADERS) \
	$(MCCCHEADERS) $(STEPHEADERS) $(SIMHEADERS) $(HDF5IOHEADERS) \
	$(PLSHEADERS) $(N0HEADERS) $(MHDHEADERS) $(ASIGMAHEADERS) \
	$(LINTHEADERS) $(SPLINEHEADERS) \
	neutral.h plasma.h particle.h endcond.h B_field.h gctransform.h \
	E_field.h wall.h simulate.h diag.h offload.h boozer.h mhd.h \
	random.h print.h hdf5_interface.h suzuki.h nbi.h biosaw.h \
	asigma.h boschhale.h mpi_interface.h libascot_mem.h

OBJS= math.o list.o octree.o error.c \
	$(DIAGOBJS)  $(BFOBJS) $(EFOBJS) $(WALLOBJS) \
	$(MCCCOBJS) $(STEPOBJS) $(SIMOBJS) $(HDF5IOOBJS) \
	$(PLSOBJS) $(N0OBJS) $(MHDOBJS) $(ASIGMAOBJS) $(LINTOBJS) \
	$(SPLINEOBJS) \
	neutral.o plasma.o particle.o endcond.o B_field.o gctransform.o \
	E_field.o wall.o simulate.o diag.o offload.o boozer.o mhd.o \
	random.o print.c hdf5_interface.o suzuki.o nbi.o biosaw.o \
	asigma.o mpi_interface.o boschhale.o

#BINS=test_math test_nbi test_bsearch \
#	test_wall_2d test_plasma test_random \
#	test_wall_3d test_B test_offload test_E \
#	test_interp1Dcomp test_linint3D test_N0 test_N0_1D \
#	test_spline ascot5_main bbnbi5 test_diag_orb test_asigma \
#	test_afsi

BINS=ascot5_main

ifdef NOGIT
	DUMMY_GIT_INFO := $(shell touch gitver.h)
else
	GET_GIT_INFO := $(shell bash gitver.sh)
endif

all: $(BINS)

libascot: libascot.so
	true

libascot.so: CFLAGS+=-fPIC -shared

ifeq ($(CC),h5cc)
libascot.so: CFLAGS+=-shlib
endif

ifeq ($(CC),h5pcc)
libascot.so: CFLAGS+=-shlib
endif

libascot.so: libascot.o ascot5_main.o libascot_mem.o afsi.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

ascot5_main: ascot5_main.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

bbnbi5: bbnbi5.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

doc:
	$(MAKE) -C src doc
	$(MAKE) -C doc

clean:
	$(MAKE) -C src clean
	rm -rf build
