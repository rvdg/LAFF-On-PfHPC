#HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O2 -Wall -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 24#     smallest size to be timed.
NLAST      := 1200#     largest size to be timed.
NINC       := 24#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $n@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

# ---------------------

TEST_OBJS_GEMMJI_KERNEL_MRNR   := driver.o GemmJI_Kernel_MRxNR.o

driver_GemmJI_Kernel_MRxNR.x: $(TEST_OBJS_GEMMJI_KERNEL_MRNR) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMJI_KERNEL_MRNR) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmJI_Kernel_MRxNR.x $(LDFLAGS) 

test_GemmJI_Kernel_MRxNR: driver_GemmJI_Kernel_MRxNR.x
	echo "your_version = 'GemmJI\_Kernel\_4x4';" > output_data_GemmJI_Kernel_MRxNR.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmJI_Kernel_MRxNR.x >> output_data_GemmJI_Kernel_MRxNR.m

# ---------------------

TEST_OBJS_GEMMIJP_JI_KERNEL_MRNR   := driver.o GemmIJP_JI_Kernel_MRxNR.o 

driver_GemmIJP_JI_Kernel_MRxNR.x: $(TEST_OBJS_GEMMIJP_JI_KERNEL_MRNR) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMIJP_JI_KERNEL_MRNR) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmIJP_JI_Kernel_MRxNR.x $(LDFLAGS) 

test_GemmIJP_JI_Kernel_MRxNR: driver_GemmIJP_JI_Kernel_MRxNR.x
	echo "your_version = 'GemmIJP\_JI\_Kernel\_4x4';" > output_data_GemmIJP_JI_Kernel_MRxNR.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmIJP_JI_Kernel_MRxNR.x >> output_data_GemmIJP_JI_Kernel_MRxNR.m

# ---------------------

TEST_OBJS_GEMMFIVELOOPS_4x4   := driver.o GemmFiveLoops_4x4.o Pack.o 

driver_GemmFiveLoops_4x4.x: $(TEST_OBJS_GEMMFIVELOOPS_4x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMFIVELOOPS_4x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmFiveLoops_4x4.x $(LDFLAGS) 

test_GemmFiveLoops_4x4: driver_GemmFiveLoops_4x4.x
	echo "your_version = 'GemmFiveLoops4x4';" > output_data_GemmFiveLoops_4x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmFiveLoops_4x4.x >> output_data_GemmFiveLoops_4x4.m

# ---------------------

TEST_OBJS_GEMMFIVELOOPS_12x4   := driver.o GemmFiveLoops_12x4.o 

driver_GemmFiveLoops_12x4.x: $(TEST_OBJS_GEMMFIVELOOPS_12x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMFIVELOOPS_12x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmFiveLoops_12x4.x $(LDFLAGS) 

test_GemmFiveLoops_12x4: driver_GemmFiveLoops_12x4.x
	echo "your_version = 'GemmFiveLoops12x4';" > output_data_GemmFiveLoops_12x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmFiveLoops_12x4.x >> output_data_GemmFiveLoops_12x4.m

# ---------------------

clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
