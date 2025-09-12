CC=cc

SRC_DIR :=.
BUILD_DIR := ../build
OBJ_DIR := $(BUILD_DIR)/obj

DEPFLAGS = -MMD -MP

ifdef TRAP_FPE
	DEFINES+=-DTRAP_FPE=$(TRAP_FPE)
	CFLAGS+= -fsignaling-nans -ftrapping-math
endif

ifdef NSIMD
	DEFINES+=-DNSIMD=$(NSIMD)
endif

ifeq ($(RFOF),1)
	DEFINES+=-DRFOF
	LDFLAGS+=-lgfortran -lquadmath -lxml2 -lrfof
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
	LDFLAGS+=-mkl
	ifdef RANDOM_MKL_RNG
		DEFINES+=-DRANDOM_MKL_RNG=$(RANDOM_MKL_RNG)
	endif
else ifeq ($(RANDOM),GSL)
	DEFINES+=-DRANDOM_GSL
	LDFLAGS+=-lgsl -lgslcblas
else ifeq ($(RANDOM),LCG)
	DEFINES+=-DRANDOM_LCG
endif

CFLAGS+=-Wall -Wextra -fPIC -std=c11 -D_XOPEN_SOURCE=700 $(DEFINES) $(FLAGS)

ifeq ($(DEBUG),1)
	CFLAGS+=-g
else
	CFLAGS+=-O2
endif

UNAME_S := $(shell uname -s)

# macOS-specific adjustments
ifeq ($(UNAME_S),Darwin)
	CFLAGS+=-D_DARWIN_C_SOURCE
endif
# Linux/other OS adjustments
ifeq ($(UNAME_S),Linux)
libascot.so: CFLAGS+=-fPIC
	LDFLAGS+=-lm -shared
endif

# Write CFLAGS and CC to a file to be included into output
$(file >compiler_flags.h,#define CFLAGS "$(CFLAGS)")

# Also store the version from current git tag or from VERSION file
$(shell bash gitver.sh)

SRCS := $(shell find $(SRC_DIR) -name '*.c')
SRCS := $(filter-out $(SRC_DIR)/test_ascot5.c, $(SRCS))

OBJS := $(patsubst $(SRC_DIR)/%, $(OBJ_DIR)/%, $(SRCS:.c=.o))
DEPS := $(OBJS:.o=.d)

DOCDIR = doc/

all: libascot doc

libascot: $(BUILD_DIR)/libascot.so
	true

$(BUILD_DIR)/libascot.so: $(OBJS)
	@mkdir -p $(dir $@)
	$(CC) -shared -o $@ $^ $(LDFLAGS) -fopenmp

$(BUILD_DIR)/test_ascot5: $(OBJS) $(OBJ_DIR)/test_ascot5.o
	$(CC) -o $@ $^ $(LDFLAGS)

doc:
	doxygen Doxyfile

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c Makefile
	@mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(CFLAGS) $(DEPFLAGS)

clean:
	@rm -f $(BUILD_DIR)/libascot.so
	@rm -rf $(OBJ_DIR)
	@rm -rf $(DOCDIR)
	@rm -f gitver.h

.PHONY: all clean $(BUILD_DIR)/libascot.so

-include $(DEPS)
