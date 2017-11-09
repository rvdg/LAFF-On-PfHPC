#HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O2 -Wall -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L -fopenmp
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 120#     smallest size to be timed.
NLAST      := 3600#     largest size to be timed.
NINC       := 120#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm -fopenmp

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $n@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

# ---------------------

TEST_OBJS_PARALLEL_LOOP1_12x4   := driver.o Gemm_Parallel_Loop1_12x4.o 

driver_Gemm_Parallel_Loop1_12x4.x: $(TEST_OBJS_PARALLEL_LOOP1_12x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PARALLEL_LOOP1_12x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Gemm_Parallel_Loop1_12x4.x $(LDFLAGS) 

test_Gemm_Parallel_Loop1_12x4: driver_Gemm_Parallel_Loop1_12x4.x
	echo "your_version = 'Parallel Loop 1';" > output_data_Gemm_Parallel_Loop1_12x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Gemm_Parallel_Loop1_12x4.x >> output_data_Gemm_Parallel_Loop1_12x4.m

# ---------------------

TEST_OBJS_PARALLEL_LOOP2_12x4   := driver.o Gemm_Parallel_Loop2_12x4.o 

driver_Gemm_Parallel_Loop2_12x4.x: $(TEST_OBJS_PARALLEL_LOOP2_12x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PARALLEL_LOOP2_12x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Gemm_Parallel_Loop2_12x4.x $(LDFLAGS) 

test_Gemm_Parallel_Loop2_12x4: driver_Gemm_Parallel_Loop2_12x4.x
	echo "your_version = 'Parallel Loop 2';" > output_data_Gemm_Parallel_Loop2_12x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Gemm_Parallel_Loop2_12x4.x >> output_data_Gemm_Parallel_Loop2_12x4.m

# ---------------------

TEST_OBJS_PARALLEL_LOOP3_12x4   := driver.o Gemm_Parallel_Loop3_12x4.o 

driver_Gemm_Parallel_Loop3_12x4.x: $(TEST_OBJS_PARALLEL_LOOP3_12x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PARALLEL_LOOP3_12x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Gemm_Parallel_Loop3_12x4.x $(LDFLAGS) 

test_Gemm_Parallel_Loop3_12x4: driver_Gemm_Parallel_Loop3_12x4.x
	echo "your_version = 'Parallel Loop 3';" > output_data_Gemm_Parallel_Loop3_12x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Gemm_Parallel_Loop3_12x4.x >> output_data_Gemm_Parallel_Loop3_12x4.m

# ---------------------

TEST_OBJS_PARALLEL_LOOP4_12x4   := driver.o Gemm_Parallel_Loop4_12x4.o 

driver_Gemm_Parallel_Loop4_12x4.x: $(TEST_OBJS_PARALLEL_LOOP4_12x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PARALLEL_LOOP4_12x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Gemm_Parallel_Loop4_12x4.x $(LDFLAGS) 

test_Gemm_Parallel_Loop4_12x4: driver_Gemm_Parallel_Loop4_12x4.x
	echo "your_version = 'Parallel Loop 4';" > output_data_Gemm_Parallel_Loop4_12x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Gemm_Parallel_Loop4_12x4.x >> output_data_Gemm_Parallel_Loop4_12x4.m

# ---------------------

TEST_OBJS_PARALLEL_LOOP5_12x4   := driver.o Gemm_Parallel_Loop5_12x4.o 

driver_Gemm_Parallel_Loop5_12x4.x: $(TEST_OBJS_PARALLEL_LOOP5_12x4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PARALLEL_LOOP5_12x4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Gemm_Parallel_Loop5_12x4.x $(LDFLAGS) 

test_Gemm_Parallel_Loop5_12x4: driver_Gemm_Parallel_Loop5_12x4.x
	echo "your_version = 'Parallel Loop 5';" > output_data_Gemm_Parallel_Loop5_12x4.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Gemm_Parallel_Loop5_12x4.x >> output_data_Gemm_Parallel_Loop5_12x4.m

# ---------------------

clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
