CC=h5cc


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

ifeq ($(MPI),1)
	DEFINES+=-DMPI
	CC=h5pcc
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
		CFLAGS+=-lhdf5 -lhdf5_hl
	endif
endif

CFLAGS+=-O2 -lm -Wall -fopenmp -fPIC -std=c11 $(DEFINES) $(FLAGS)

# Write CFLAGS and CC to a file to be included into output
$(shell echo "#define CFLAGS " $(CFLAGS) > compiler_flags.h)
$(shell echo "#define CC " $(CC) >> compiler_flags.h)

SIMDIR = simulate/
SIMHEADERS = $(wildcard $(SIMDIR)simulate*.h)
SIMOBJS = $(patsubst %.c,%.o,$(wildcard $(SIMDIR)simulate*.c))

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
	$(PLSHEADERS) $(N0HEADERS) $(MHDHEADERS) $(LINTHEADERS) $(SPLINEHEADERS) \
	neutral.h plasma.h particle.h endcond.h B_field.h gctransform.h \
	E_field.h wall.h simulate.h diag.h offload.h boozer.h mhd.h \
	random.h print.h hdf5_interface.h suzuki.h nbi.h biosaw.h \
	mpi_interface.h libascot_mem.h

OBJS= math.o list.o octree.o error.c \
	$(DIAGOBJS)  $(BFOBJS) $(EFOBJS) $(WALLOBJS) \
	$(MCCCOBJS) $(STEPOBJS) $(SIMOBJS) $(HDF5IOOBJS) \
	$(PLSOBJS) $(N0OBJS) $(MHDOBJS) $(LINTOBJS) $(SPLINEOBJS) \
	neutral.o plasma.o particle.o endcond.o B_field.o gctransform.o \
	E_field.o wall.o simulate.o diag.o offload.o boozer.o mhd.o \
	random.o print.c hdf5_interface.o suzuki.o nbi.o biosaw.o \
	mpi_interface.o

BINS=test_math test_nbi test_bsearch \
	test_wall_2d test_plasma test_random \
	test_wall_3d test_B test_offload test_E \
	test_interp1Dcomp test_linint3D test_N0 \
	test_spline ascot5_main bbnbi5 test_diag_orb

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

libascot.so: libascot.o ascot5_main.o libascot_mem.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lhdf5 -lhdf5_hl

ascot5_main: ascot5_main.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

bbnbi5: bbnbi5.o $(OBJS)
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

test_offload: $(UTESTDIR)test_offload.o
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

test_nbi: $(UTESTDIR)test_nbi.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_spline: $(UTESTDIR)test_spline.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(HEADERS) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

ASCOTPY2_HEADERFILES=particle.h hdf5_interface.h ascot5.h mpi_interface.h \
	simulate.h ascot5_main.h offload.h diag.h libascot_mem.h \
	wall.h Bfield/B_STS.h B_field.h

ascotpy2.py : libascot.so
	clang2py -l libascot.so -o $@  \
		$(ASCOTPY2_HEADERFILES) \
		--clang-args="-I/usr/include/hdf5/serial"
# The above hdf5-include folder should not be hardcoded...

clean:
	@rm -f *.o *.so *.test *.optrpt $(BINS) $(SIMDIR)*.o $(STEPDIR)*.o \
		$(MCCCDIR)*.o $(HDF5IODIR)*.o $(PLSDIR)*.o $(DIAGDIR)*.o \
		$(BFDIR)*.o $(EFDIR)*.o $(WALLDIR)*.o $(MHDDIR)*.o \
		$(N0DIR)*.o $(LINTDIR)*.o $(SPLINEDIR)*.o $(UTESTDIR)*.o *.pyc
	@rm -rf $(DOCDIR)
	@rm -f gitver.h
