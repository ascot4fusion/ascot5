CC=icc

ifdef NSIMD
	DEFINES+=-DNSIMD=$(NSIMD)
endif

ifeq ($(NOTARGET),1)
	DEFINES+=-DNOTARGET
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

ifeq ($(SINGLEPRECISION),1)
	DEFINES+=-DSINGLEPRECISION
endif

ifeq ($(MPI),1)
	DEFINES+=-DMPI
	CC=mpiicc
endif

CFLAGS=-lm -lhdf5 -lhdf5_hl -fopenmp -std=c99 $(DEFINES) $(FLAGS) 

HEADERS=ascot5.h B_GS.h math.h consts.h \
	   	wall_2d.h ascot4_interface.h distributions.h B_2D.h \
		plasma_1d.h interact.h step_gc_rk4.h step_fo_lf.h simulate.h \
		hdf5_helpers.h hdf5_histogram.h B_3D.h simulate_fo_lf.h \
		simulate_gc_rk4.h wall_3d.h list.h octree.h hdf5_particlestate.h \
		step_fo_vpa.h B_ST.h B_TC.h \
		particle.h filip5.h endcond.h orbit_write.h \
		B_field.h E_field.h wall.h

OBJS=ascot4_interface.o B_GS.o math.o consts.o  \
     wall_2d.o distributions.o B_2D.o B_ST.o B_TC.o  \
	 plasma_1d.o interact.o step_gc_rk4.o step_fo_lf.o step_fo_vpa.o \
	 hdf5_helpers.o hdf5_histogram.o B_3D.o simulate_fo_lf.o \
	 simulate_gc_rk4.o wall_3d.o list.o octree.o hdf5_particlestate.o \
     particle.o endcond.o B_field.o E_field.o wall.o simulate.o \

BINS=test_math \
	 test_wall_2d test_ascot4_interface test_plasma_1d \
	 test_interact test_hdf5 test_wall_3d test_particle filip5 \
	 test_B ascot5_gc test_simulate_orbit

all: $(BINS)

ascotpy: CFLAGS+=-fPIC
ascotpy: ascotpy.c B_none.o B_GS.o B_2Dlin.o B_2D.o B_3D.o
	f2py ascotpy.pyf ascotpy.c B_none.o B_GS.o B_2Dlin.o B_2D.o B_3D.o -c $(DEFINES)

ascot5_gc: ascot5_gc.o $(OBJS)
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

%.o: %.c $(HEADERS) Makefile
	$(CC) -c -o $@ $< $(CFLAGS)

clean:
	rm -f *.o *.optrpt $(BINS)
