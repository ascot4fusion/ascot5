CC=h5cc

ifdef NSIMD
	DEFINES+=-DNSIMD=$(NSIMD)
endif

ifneq ($(TARGET),1)
	DEFINES+=-DNOTARGET
#	CFLAGS+=-qno-openmp-offload -diag-disable 3180
endif

ifdef VERBOSE
	DEFINES+=-DVERBOSE=$(VERBOSE)
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

else ifeq ($(RANDOM),GSL)
	DEFINES+=-DRANDOM_GSL
	CFLAGS+=-lgsl -lgslcblas
endif

ifneq ($(CC),h5cc)
	ifneq ($(CC),h5pcc)
		CFLAGS+=-lhdf5 -lhdf5_hl
	endif
endif

CFLAGS+=-lm -Wall -fopenmp -std=c99 $(DEFINES) $(FLAGS) 

SIMDIR = simulate/
SIMHEADERS = $(wildcard $(SIMDIR)simulate*.h)
SIMOBJS = $(patsubst %.c,%.o,$(wildcard $(SIMDIR)simulate*.c))

STEPDIR = $(SIMDIR)step/
STEPHEADERS = $(wildcard $(STEPDIR)step*.h)
STEPOBJS = $(patsubst %.c,%.o,$(wildcard $(STEPDIR)step*.c))

MCCCDIR = $(SIMDIR)mccc/
MCCCHEADERS = $(wildcard $(MCCCDIR)mccc*.h)
MCCCOBJS = $(patsubst %.c,%.o,$(wildcard $(MCCCDIR)mccc*.c))

DIAGHEADERS = diag.h diag_orb.h dist_5D.h dist_6D.h
DIAGOBJS = diag.o diag_orb.o dist_5D.o dist_6D.o

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

HEADERS=ascot5.h math.h consts.h list.h octree.h physlib.h error.h \
	$(DIAGHEADERS) $(BFHEADERS) $(EFHEADERS) $(WALLHEADERS) \
	$(MCCCHEADERS) $(STEPHEADERS) $(SIMHEADERS) $(HDF5IOHEADERS) \
	$(PLSHEADERS) plasma.h particle.h endcond.h B_field.h E_field.h \
	wall.h simulate.h diag.h diag_orb.h offload.h \
	spline/interp2D.h spline/interp3D.h spline/spline1D.h \
	spline/interp2Dexpl.h spline/interp3Dexpl.h \
	spline/interp2Detoc.h spline/interp3Detoc.h \
	spline/interp2Dcomp.h spline/interp3Dcomp.h \
	spline/spline1Dcomp.h spline/interp1Dcomp.h \
	random.h

OBJS= math.o list.o octree.o physlib.o \
	$(DIAGOBJS)  $(BFOBJS) $(EFOBJS) $(WALLOBJS) \
	$(MCCCOBJS) $(STEPOBJS) $(SIMOBJS) $(HDF5IOOBJS) \
	$(PLSOBJS) plasma.o particle.o endcond.o B_field.o E_field.o \
	wall.o simulate.o diag.o diag_orb.o offload.o \
	spline/interp2D.o spline/interp3D.o spline/spline1D.o \
	spline/interp2Dexpl.o spline/interp3Dexpl.o \
	spline/interp2Detoc.o spline/interp3Detoc.o \
	spline/interp2Dcomp.o spline/interp3Dcomp.o \
	spline/spline1Dcomp.o spline/interp1Dcomp.o \
	random.o

BINS=test_math \
	 test_wall_2d test_plasma_1d test_random \
	 test_hdf5 test_wall_3d test_particle \
	 test_B test_simulate_orbit test_offload test_E \
	 test_mccc test_interp1Dcomp ascot5_main

all: $(BINS)

ascotpy: ascotpy.so
	true

ascotpy.so: CFLAGS+=-shlib -fPIC -shared

ascotpy.so: ascotpy.o $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

ascot5_main: ascot5_main.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_B: test_B.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_particle: test_particle.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_wall_3d: test_wall_3d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_plasma_1d: test_plasma_1d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_math: test_math.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_wall_2d: test_wall_2d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_simulate_orbit: test_simulate_orbit.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_offload: test_offload.o 
	$(CC) -o $@ $^ $(CFLAGS)

test_E: test_E.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_mccc: simulate/mccc/test_mccc.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp1Dcomp: test_interp1Dcomp.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_random: test_random.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(HEADERS) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *.o *.so *.test *.optrpt $(BINS) $(SIMDIR)*.o $(STEPDIR)*.o \
		$(MCCCDIR)*.o $(HDF5IODIR)*.o $(PLSDIR)*.o \
		$(BFDIR)*.o $(EFDIR)*.o $(WALLDIR)*.o \
		spline/*.o *.pyc
