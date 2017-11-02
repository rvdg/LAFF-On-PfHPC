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
NFIRST     := 128#     smallest size to be timed.
NLAST      := 1500#     largest size to be timed.
NINC       := 128#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

TEST_OBJS_GEMMJI_KERNEL_MRNR   := driver_GemmJI_Kernel_MRxNR.o GemmJI_Kernel_MRxNR.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $n@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@


driver_GemmJI_Kernel_MRxNR.x: $(TEST_OBJS_GEMMJI_KERNEL_MRNR) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMJI_KERNEL_MRNR) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmJI_Kernel_MRxNR.x $(LDFLAGS) 

test_GemmJI_Kernel_MRxNR: driver_GemmJI_Kernel_MRxNR.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmJI_Kernel_MRxNR.x > output_data_GemmJI_Kernel_MRxNR.m


clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
