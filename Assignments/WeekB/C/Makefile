#HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O3 -Wall -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 48#     smallest size to be timed.
NLAST_SMALL:= 1000#    largest size to be timed for slow implementations.
NLAST      := 2000#   largest size to be timed for fast implementations.
NINC       := 48#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# ---------------------

TEST_OBJS_IJP := driver.o Gemm_JI_4x4Kernel.o

driver_JI_4x4Kernel.x: $(TEST_OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_4x4Kernel.x $(LDFLAGS) 

test_JI_4x4Kernel: driver_JI_4x4Kernel.x
	echo "your_version = 'JI\_4x4Kernel';" > output_JI_4x4Kernel.m
	echo " " >> output_JI_4x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_4x4Kernel.x >> output_JI_4x4Kernel.m

# ---------------------

TEST_OBJS_IJP := driver.o Gemm_JI_8x4Kernel.o

driver_JI_8x4Kernel.x: $(TEST_OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_8x4Kernel.x $(LDFLAGS) 

test_JI_8x4Kernel: driver_JI_8x4Kernel.x
	echo "your_version = 'JI\_8x4Kernel';" > output_JI_8x4Kernel.m
	echo " " >> output_JI_8x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_8x4Kernel.x >> output_JI_8x4Kernel.m

# ---------------------

TEST_OBJS_IJP := driver.o Gemm_JI_12x4Kernel.o

driver_JI_12x4Kernel.x: $(TEST_OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_12x4Kernel.x $(LDFLAGS) 

test_JI_12x4Kernel: driver_JI_12x4Kernel.x
	echo "your_version = 'JI\_12x4Kernel';" > output_JI_12x4Kernel.m
	echo " " >> output_JI_12x4Kernel.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_12x4Kernel.x >> output_JI_12x4Kernel.m

# ---------------------

clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
