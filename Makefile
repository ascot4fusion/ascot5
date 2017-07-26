ifndef CC
	CC=icc
endif

ifdef NSIMD
	DEFINES+=-DNSIMD=$(NSIMD)
endif

ifneq ($(TARGET),1)
	DEFINES+=-DNOTARGET
    CFLAGS+=-qno-openmp-offload -diag-disable 3180
endif

ifdef VERBOSE
	DEFINES+=-DVERBOSE=$(VERBOSE)
endif

ifdef REPORTORBIT
	DEFINES+=-DREPORTORBIT=$(REPORTORBIT)
endif

ifdef COULOMBCOLL
	DEFINES+=-DCOULOMBCOLL=$(COULOMBCOLL)
else
    DEFINES+=-DCOULOMBCOLL=0
endif

ifdef ORBITFOLLOWING
	DEFINES+=-DORBITFOLLOWING=$(ORBITFOLLOWING)
else
    DEFINES+=-DORBITFOLLOWING=1
endif

ifeq ($(SINGLEPRECISION),1)
	DEFINES+=-DSINGLEPRECISION
endif

ifeq ($(MPI),1)
	DEFINES+=-DMPI
	CC=mpiicc
endif

ifeq ($(CC),h5cc)
	CFLAGS=-Wall -fopenmp -std=c99 $(DEFINES) $(FLAGS) 
else
	CFLAGS=-lm -lhdf5 -lhdf5_hl -fopenmp -std=c99 $(DEFINES) $(FLAGS) 
endif

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


HEADERS=ascot5.h math.h consts.h \
	   	$(DIAGHEADERS)  \
		plasma_1d.h simulate.h \
		list.h octree.h $(BFHEADERS) $(EFHEADERS) \
		particle.h filip5.h endcond.h $(WALLHEADERS) \
		B_field.h E_field.h wall.h phys_orbit.h \
		$(MCCCHEADERS) $(STEPHEADERS) $(SIMHEADERS) \
		diag.h diag_orb.h \
		$(HDF5IOHEADERS) \
        offload.h

OBJS= math.o consts.o  \
     $(DIAGOBJS)  $(BFOBJS) $(EFOBJS) \
	plasma_1d.o \
	 $(WALLOBJS) list.o octree.o \
     particle.o endcond.o B_field.o E_field.o wall.o simulate.o \
	phys_orbit.o \
	$(MCCCOBJS) $(STEPOBJS) $(SIMOBJS)  \
	diag.o diag_orb.o \
	$(HDF5IOOBJS) \
    offload.c \
	splinePatrik/interp2D.o splinePatrik/interp3D.o splinePatrik/spline1D.o \
	splinePatrik/interp2Dexpl.o splinePatrik/interp3Dexpl.o

BINS=test_math \
	 test_wall_2d test_ascot4_interface test_plasma_1d \
	 test_interact test_hdf5 test_wall_3d test_particle filip5 \
	 test_B ascot5_gc test_simulate_orbit test_offload test_E \
	 test_mccc ascot5_main\

all: $(BINS)

ascotpy: CFLAGS+=-fPIC
ascotpy: ascotpy.o math.o B_field.o $(HDF5IOOBJS) $(BFOBJS)
	f2py ascotpy.pyf ascotpy.c math.o B_field.o -c $(DEFINES)

ascot5_gc: ascot5_gc.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

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

test_mccc: mccc/test_mccc.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp2D: test_interp2D.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp3D: test_interp3D.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp2Dexpl: test_interp2Dexpl.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

test_interp3Dexpl: test_interp3Dexpl.o $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

%.o: %.c $(HEADERS) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *.o *.test *.optrpt $(BINS) $(SIMDIR)*.o $(STEPDIR)*.o \
		$(MCCCDIR)*.o $(HDF5IODIR)*.o *.pyc $(ASCOT4IFDIR)*.pyc \
		$(UIDIR)*.pyc $(BFDIR)*.o $(EFDIR)*.o $(WALLDIR)*.o
