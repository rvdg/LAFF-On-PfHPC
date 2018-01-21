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

TEST_OBJS_IJP := driver.o Gemm_IJP.o

driver_IJP.x: $(TEST_OBJS_IJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJP.x $(LDFLAGS) 

test_IJP: driver_IJP.x
	echo "your_version = 'IJP';" > output_IJP.m
	echo " " >> output_IJP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJP.x >> output_IJP.m

# ---------------------

TEST_OBJS_IPJ := driver.o Gemm_IPJ.o

driver_IPJ.x: $(TEST_OBJS_IPJ) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IPJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IPJ.x $(LDFLAGS) 

test_IPJ: driver_IPJ.x
	echo "your_version = 'IPJ';" > output_IPJ.m
	echo " " >> output_IPJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IPJ.x >> output_IPJ.m

# ---------------------

TEST_OBJS_JIP := driver.o Gemm_JIP.o

driver_JIP.x: $(TEST_OBJS_JIP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JIP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JIP.x $(LDFLAGS) 

test_JIP: driver_JIP.x
	echo "your_version = 'JIP';" > output_JIP.m
	echo " " >> output_JIP.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JIP.x >> output_JIP.m

# ---------------------

TEST_OBJS_JPI := driver.o Gemm_JPI.o

driver_JPI.x: $(TEST_OBJS_JPI) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JPI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JPI.x $(LDFLAGS) 

test_JPI: driver_JPI.x
	echo "your_version = 'JPI';" > output_JPI.m
	echo " " >> output_JPI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JPI.x >> output_JPI.m

# ---------------------

TEST_OBJS_PIJ := driver.o Gemm_PIJ.o

driver_PIJ.x: $(TEST_OBJS_PIJ) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PIJ) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PIJ.x $(LDFLAGS) 

test_PIJ: driver_PIJ.x
	echo "your_version = 'PIJ';" > output_PIJ.m
	echo " " >> output_PIJ.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PIJ.x >> output_PIJ.m

# ---------------------

TEST_OBJS_PJI := driver.o Gemm_PJI.o

driver_PJI.x: $(TEST_OBJS_PJI) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_PJI) $(UTIL_OBJS) $(BLAS_LIB) -o driver_PJI.x $(LDFLAGS) 

test_PJI: driver_PJI.x
	echo "your_version = 'PJI';" > output_PJI.m
	echo " " >> output_PJI.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_PJI.x >> output_PJI.m

# ---------------------

TEST_OBJS_IJ_DOTS := driver.o Gemm_IJ_Dots.o

driver_IJ_Dots.x: $(TEST_OBJS_IJ_DOTS) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IJ_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_IJ_Dots.x $(LDFLAGS) 

test_IJ_Dots: driver_IJ_Dots.x
	echo "your_version = 'IJ_Dots';" > output_IJ_Dots.m
	echo " " >> output_IJ_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_IJ_Dots.x >> output_IJ_Dots.m

# ---------------------

TEST_OBJS_JI_DOTS := driver.o Gemm_JI_Dots.o

driver_JI_Dots.x: $(TEST_OBJS_JI_DOTS) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_JI_DOTS) $(UTIL_OBJS) $(BLAS_LIB) -o driver_JI_Dots.x $(LDFLAGS) 

test_JI_Dots: driver_JI_Dots.x
	echo "your_version = 'JI_Dots';" > output_JI_Dots.m
	echo " " >> output_JI_Dots.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_JI_Dots.x >> output_JI_Dots.m

# ---------------------

TEST_OBJS_P_GER := driver.o Gemm_P_Ger.o

driver_P_Ger.x: $(TEST_OBJS_P_GER) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_P_GER) $(UTIL_OBJS) $(BLAS_LIB) -o driver_P_Ger.x $(LDFLAGS) 

test_P_Ger: driver_P_Ger.x
	echo "your_version = 'P_Ger';" > output_P_Ger.m
	echo " " >> output_P_Ger.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_P_Ger.x >> output_P_Ger.m

# ---------------------

TEST_OBJS_J_GEMV := driver.o Gemm_J_Gemv.o

driver_J_Gemv.x: $(TEST_OBJS_J_GEMV) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_J_GEMV) $(UTIL_OBJS) $(BLAS_LIB) -o driver_J_Gemv.x $(LDFLAGS) 

test_J_Gemv: driver_J_Gemv.x
	echo "your_version = 'J_Gemv';" > output_J_Gemv.m
	echo " " >> output_J_Gemv.m
	echo "$(NREPEATS) $(NFIRST) $(NLAST_SMALL) $(NINC)" | ./driver_J_Gemv.x >> output_J_Gemv.m

# ---------------------

clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
