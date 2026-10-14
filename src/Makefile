CC=cc

SRC_DIR :=.
BUILD_DIR := ../build
OBJ_DIR := $(BUILD_DIR)/obj

DEPFLAGS = -MMD -MP

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
	ifeq ($(ACC), 1)
		CFLAGS+=-c++libs -cuda -acc -Mnoopenmp
	else ifeq ($(OMP), 1)
		CFLAGS+=-fopenmp
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

ifeq ($(DEBUG),1)
	CFLAGS+=-g -fsanitize=address,undefined
	LDFLAGS+=-fsanitize=address,undefined
else
	CFLAGS+=-O2
endif

CFLAGS+=-Wall -Wextra -fPIC -std=c11 -D_XOPEN_SOURCE=700 -Iinclude -Iascot \
	-Irandom123 $(DEFINES) $(FLAGS)

SRCS := $(shell find $(SRC_DIR) -name '*.c')
OBJS := $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%,$(SRCS:.c=.o))
DEPS := $(OBJS:.o=.d)

all: ascot test_ascot5

ascot: $(BUILD_DIR)/libascot.so

test_ascot5: $(BUILD_DIR)/test_ascot5

$(BUILD_DIR)/libascot.so: $(OBJS)
	@mkdir -p $(dir $@)
	$(CC) -shared -o $@ $^ $(LDFLAGS) -lm -fopenmp

$(BUILD_DIR)/test_ascot5: $(OBJS)
	@mkdir -p $(dir $@)
	$(CC) -o $@ $^ $(LDFLAGS) -lm -fopenmp

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c Makefile
	@mkdir -p $(dir $@)
	$(CC) -c -o $@ $< $(CFLAGS) $(DEPFLAGS)

clean:
	@rm -f $(BUILD_DIR)/libascot.so
	@rm -f $(BUILD_DIR)/test_ascot5
	@rm -rf $(OBJ_DIR)

.PHONY: all clean test_ascot5 libascot.so

-include $(DEPS)
