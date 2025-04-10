CC=h5cc
LFLAGS=

ifdef TRAP_FPE
	DEFINES+=-DTRAP_FPE=$(TRAP_FPE)
	CFLAGS+= -fsignaling-nans -ftrapping-math
endif

ifdef NSIMD
	DEFINES+=-DNSIMD=$(NSIMD)
endif

ifdef TARGET
	DEFINES+=-DTARGET=$(TARGET)
endif

ifdef VERBOSE
	DEFINES+=-DVERBOSE=$(VERBOSE)
else
	DEFINES+=-DVERBOSE=1
endif

ifeq ($(SINGLEPRECISION),1)
	DEFINES+=-DSINGLEPRECISION
endif

ifeq ($(RFOF),1)
	DEFINES+=-DRFOF
	CFLAGS+=-lgfortran -lquadmath -lxml2 -lrfof
endif

ifeq ($(MPI),1)
	DEFINES+=-DMPI
	CC=h5pcc
endif

ifneq ($(CC),h5cc)
	ifneq ($(CC),h5pcc)
		CFLAGS+=-lhdf5 -lhdf5_hl
	endif
endif

ifeq ($(GPU),1)
	DEFINES+=-DGPU
	RANDOM=LCG
	CFLAGS+=-c++libs -cuda
	ifeq ($(ACC), 1)
		CFLAGS+=-acc -Mnoopenmp
	else ifeq ($(OMP), 1)
		CFLAGS+=-mp=gpu -fopenmp
		DEFINES+=-D_OPENMP
	endif
else
	CFLAGS+=-fopenmp
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

CFLAGS+=-Wall -fPIC -std=c11 -D_XOPEN_SOURCE=700 $(DEFINES) $(FLAGS)

ifeq ($(DEBUG),1)
	CFLAGS+=-g -fsanitize=address,undefined
else
	CFLAGS+=-O2
endif

ifeq ($(CC), $(filter $(CC), h5cc h5pcc))
	CFLAGS += -shlib
endif

UNAME_S := $(shell uname -s)

# macOS-specific adjustments
ifeq ($(UNAME_S),Darwin)
	CFLAGS+=-D_DARWIN_C_SOURCE
	H5_INCLUDES=$(shell dirname $(shell which h5cc))/../include
	CLANG2PY_ARGS=$(CFLAGS) -I$(H5_INCLUDES)
endif
# Linux/other OS adjustments
ifeq ($(UNAME_S),Linux)
libascot.so: CFLAGS+=-fPIC -shared
	CFLAGS+=-lm
	ifdef CONDA_PREFIX
		CLANG2PY_ARGS=-I$(CONDA_PREFIX)/include/ \
		  -I$(CONDA_PREFIX)/x86_64-conda-linux-gnu/sysroot/usr/include/
	endif
endif

# Write CFLAGS and CC to a file to be included into output
$(shell echo "#define CFLAGS " $(CFLAGS) > compiler_flags.h)
$(shell echo "#define CC " $(CC) >> compiler_flags.h)

# Also store the version from current git tag or from VERSION file
$(shell bash gitver.sh)

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
DIAGHEADERS = $(wildcard $(DIAGDIR)diag*.h $(DIAGDIR)dist*.h $(DIAGDIR)hist*.h)
DIAGOBJS = $(patsubst %.c,%.o,$(wildcard $(DIAGDIR)diag*.c $(DIAGDIR)dist*.c \
	$(DIAGDIR)hist.c))

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
N0HEADERS = $(wildcard $(N0DIR)N0_*.h)
N0OBJS = $(patsubst %.c,%.o,$(wildcard $(N0DIR)N0_*.c))

MHDDIR = mhd/
MHDHEADERS = $(wildcard $(MHDDIR)mhd_*.h)
MHDOBJS = $(patsubst %.c,%.o,$(wildcard $(MHDDIR)mhd_*.c))

ASIGMADIR = asigma/
ASIGMAHEADERS = $(wildcard $(ASIGMADIR)asigma_*.h)
ASIGMAOBJS = $(patsubst %.c,%.o,$(wildcard $(ASIGMADIR)asigma_*.c))

LINTDIR = linint/
LINTHEADERS = $(wildcard $(LINTDIR)linint*.h)
LINTOBJS = $(patsubst %.c,%.o,$(wildcard $(LINTDIR)linint*.c))

SPLINEDIR = spline/
SPLINEHEADERS = $(wildcard $(SPLINEDIR)spline.h $(SPLINEDIR)interp.h)
SPLINEOBJS = $(patsubst %.c,%.o,$(wildcard $(SPLINEDIR)spline*.c \
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
	asigma.h boschhale.h mpi_interface.h libascot_mem.h bbnbi5.h rfof.h

OBJS= math.o list.o octree.o error.o \
	$(DIAGOBJS)  $(BFOBJS) $(EFOBJS) $(WALLOBJS) \
	$(MCCCOBJS) $(STEPOBJS) $(SIMOBJS) $(HDF5IOOBJS) \
	$(PLSOBJS) $(N0OBJS) $(MHDOBJS) $(ASIGMAOBJS) $(LINTOBJS) \
	$(SPLINEOBJS) \
	neutral.o plasma.o particle.o endcond.o B_field.o gctransform.o \
	E_field.o wall.o simulate.o diag.o boozer.o mhd.o \
	random.o print.o hdf5_interface.o suzuki.o nbi.o biosaw.o \
	asigma.o mpi_interface.o boschhale.o bbnbi5.o rfof.o \
	asigma.o mpi_interface.o boschhale.o bbnbi5.o

BINS=test_math test_nbi test_bsearch \
	test_wall_2d test_plasma test_random \
	test_wall_3d test_B test_offload test_E \
	test_interp1Dcomp test_linint3D test_N0 test_N0_1D \
	test_spline ascot5_main bbnbi5 test_diag_orb test_asigma \
	test_afsi

all: $(BINS)

libascot: libascot.so
	true

libascot.so: libascot.o ascot5_main.o libascot_mem.o afsi.o bbnbi5.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LFLAGS)

ascot5_main: ascot5_main.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

bbnbi5: bbnbi5_main.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

doc:
	doxygen Doxyfile

test_B: $(UTESTDIR)test_B.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_diag_orb: $(UTESTDIR)test_diag_orb.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_wall_3d: $(UTESTDIR)test_wall_3d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_math: $(UTESTDIR)test_math.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_wall_2d: $(UTESTDIR)test_wall_2d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_E: $(UTESTDIR)test_E.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp1Dcomp: $(UTESTDIR)test_interp1Dcomp.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_plasma: $(UTESTDIR)test_plasma.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_random: $(UTESTDIR)test_random.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_linint3D: $(UTESTDIR)test_linint3D.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_N0: $(UTESTDIR)test_N0.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_bsearch: $(UTESTDIR)test_bsearch.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_afsi: $(UTESTDIR)test_afsi.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_nbi: $(UTESTDIR)test_nbi.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_spline: $(UTESTDIR)test_spline.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_asigma: $(UTESTDIR)test_asigma.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(HEADERS) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

ASCOT2PY_HEADERFILES=ascot5.h particle.h mpi_interface.h endcond.h simulate.h  \
	$(DIAGHEADERS) $(BFHEADERS) $(EFHEADERS) $(WALLHEADERS) $(PLSHEADERS)  \
	$(ASIGMAHEADERS) $(N0HEADERS) $(MHDHEADERS) hdf5_interface.h           \
	libascot_mem.h diag.h B_field.h E_field.h wall.h plasma.h neutral.h    \
	boozer.h mhd.h nbi.h asigma.h afsi.h ascot5_main.h biosaw.h boschhale.h\
	bbnbi5.h rfof.h

ascot2py.py : libascot.so
	$(eval CLANGOMP=$(shell clang --print-resource-dir)) \
	clang2py -l libascot.so -o $@ $(ASCOT2PY_HEADERFILES) \
		--clang-args="$(CLANG2PY_ARGS) -I$(CLANGOMP)/include/"

clean:
	@rm -f *.o *.so *.test *.optrpt $(BINS) $(SIMDIR)*.o $(STEPDIR)*.o \
		$(MCCCDIR)*.o $(HDF5IODIR)*.o $(PLSDIR)*.o $(DIAGDIR)*.o \
		$(BFDIR)*.o $(EFDIR)*.o $(WALLDIR)*.o $(MHDDIR)*.o \
		$(N0DIR)*.o $(ASIGMADIR)*.o $(LINTDIR)*.o $(SPLINEDIR)*.o \
	        $(UTESTDIR)*.o *.pyc
	@rm -rf $(DOCDIR)
	@rm -f gitver.h
