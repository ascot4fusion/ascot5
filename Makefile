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

DIAGHEADERS = diag.h diag_orb.h distributions.h
DIAGOBJS = diag.o diag_orb.o distributions.o

BFDIR = Bfield/
BFHEADERS = $(wildcard $(BFDIR)B_*.h)
BFOBJS = $(patsubst %.c,%.o,$(wildcard $(BFDIR)B_*.c))

EFDIR = Efield/
EFHEADERS = $(wildcard $(EFDIR)E_*.h)
EFOBJS = $(patsubst %.c,%.o,$(wildcard $(EFDIR)E_*.c))

WALLDIR = wall/
WALLHEADERS = $(wildcard $(WALLDIR)wall_*.h)
WALLOBJS = $(patsubst %.c,%.o,$(wildcard $(WALLDIR)wall_*.c))

HDF5IODIR = hdf5io/
HDF5IOHEADERS = $(wildcard $(HDF5IODIR)hdf5*.h)
HDF5IOOBJS = $(patsubst %.c,%.o,$(wildcard $(HDF5IODIR)hdf5*.c))

UIDIR = ui/
ASCOT4IFDIR = ascot4_interface/


HEADERS=ascot5.h math.h consts.h list.h octree.h physlib.h error.h \
	$(DIAGHEADERS) $(BFHEADERS) $(EFHEADERS) $(WALLHEADERS) \
	$(MCCCHEADERS) $(STEPHEADERS) $(SIMHEADERS) $(HDF5IOHEADERS) \
	plasma_1d.h particle.h endcond.h B_field.h E_field.h \
	wall.h simulate.h diag.h diag_orb.h offload.h \
	spline/interp2D.h spline/interp3D.h spline/spline1D.h \
	spline/interp2Dexpl.h spline/interp3Dexpl.h \
	spline/interp2Detoc.h spline/interp3Detoc.h \
	spline/interp2Dcomp.h spline/interp3Dcomp.h spline/spline1Dcomp.h \

OBJS= math.o consts.o list.o octree.o physlib.o \
     $(DIAGOBJS)  $(BFOBJS) $(EFOBJS) $(WALLOBJS) \
	$(MCCCOBJS) $(STEPOBJS) $(SIMOBJS) $(HDF5IOOBJS) \
	plasma_1d.o particle.o endcond.o B_field.o E_field.o \
	wall.o simulate.o diag.o diag_orb.o offload.o \
	spline/interp2D.o spline/interp3D.o spline/spline1D.o \
	spline/interp2Dexpl.o spline/interp3Dexpl.o \
	spline/interp2Detoc.o spline/interp3Detoc.o \
	spline/interp2Dcomp.o spline/interp3Dcomp.o spline/spline1Dcomp.o \

BINS=test_math \
	 test_wall_2d test_ascot4_interface test_plasma_1d \
	 test_interact test_hdf5 test_wall_3d test_particle filip5 \
	 test_B test_simulate_orbit test_offload test_E \
	 test_mccc ascot5_main \

all: $(BINS)

ascotpy: CFLAGS+=-fPIC
ascotpy: ascotpy.o math.o B_field.o $(HDF5IOOBJS) $(BFOBJS)
	f2py ascotpy.pyf ascotpy.c math.o B_field.o -c $(DEFINES)

ascot5_main: ascot5_main.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_B: test_B.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

filip5: filip5.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_particle: test_particle.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_wall_3d: test_wall_3d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_hdf5: test_hdf5.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_plasma_1d: test_plasma_1d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_ascot4_interface: test_ascot4_interface.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_math: test_math.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_wall_2d: test_wall_2d.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_simulate_orbit: test_simulate_orbit.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interact: test_interact.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_offload: test_offload.o 
	$(CC) -o $@ $^ $(CFLAGS)

test_diag_offload: test_diag_offload.o simulate.o $(DIAGOBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_E: test_E.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_mccc: simulate/mccc/test_mccc.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp2D: test_interp2D.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp3D: test_interp3D.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp2Dexpl: test_interp2Dexpl.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp3Dexpl: test_interp3Dexpl.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp2Detoc: test_interp2Detoc.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp3Detoc: test_interp3Detoc.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp2Dcomp: test_interp2Dcomp.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp3Dcomp: test_interp3Dcomp.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(HEADERS) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *.o *.test *.optrpt $(BINS) $(SIMDIR)*.o $(STEPDIR)*.o \
		$(MCCCDIR)*.o $(HDF5IODIR)*.o *.pyc $(ASCOT4IFDIR)*.pyc \
		$(UIDIR)*.pyc $(BFDIR)*.o $(EFDIR)*.o $(WALLDIR)*.o \
		spline/*.o
