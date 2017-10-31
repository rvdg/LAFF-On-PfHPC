#HOME      := /Users/rvdg
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O2 -Wall -I$(BLAS_INC) -m64 -msse3 -std=c99 -D_POSIX_C_SOURCE=200112L
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 32#     smallest size to be timed.
NLAST      := 500#     largest size to be timed.
NINC       := 32#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

TEST_OBJS_GEMMIJP      := driver_GemmIJP.o GemmIJP.o 
TEST_OBJS_GEMMALL      := driver_GemmAll.o GemmIJP.o GemmJIP.o GemmIPJ.o GemmPIJ.o GemmJPI.o GemmPJI.o
TEST_OBJS_AXPY         := driver_Axpy.o Axpy.o
TEST_OBJS_IPJ_AXPY     := driver_GemmIPJ_Axpy.o GemmIPJ_Axpy.o Axpy.o
TEST_OBJS_GERJI_AXPY   := driver_GerJI_Axpy.o GerJI_Axpy.o Axpy.o
TEST_OBJS_GEMMJPI_GEMV := driver_GemmJPI_Gemv.o GemmJPI_Gemv.o GemvJI_Axpy.o Axpy.o
TEST_OBJS_GEMMPJI_GER  := driver_GemmPJI_Ger.o GemmPJI_Ger.o GerJI_Axpy.o Axpy.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $n@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver_GemmIJP.x: $(TEST_OBJS_GEMMIJP) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMIJP) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmIJP.x $(LDFLAGS) 

test_GemmIJP:   driver_GemmIJP.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmIJP.x > output_data_GemmIJP.m

driver_GemmAll.x: $(TEST_OBJS_GEMMALL) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMALL) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmAll.x $(LDFLAGS) 

test_GemmAll: driver_GemmAll.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmAll.x > output_data_GemmAll.m

driver_Axpy.x: $(TEST_OBJS_AXPY) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_AXPY) $(UTIL_OBJS) $(BLAS_LIB) -o driver_Axpy.x $(LDFLAGS) 

test_Axpy: driver_Axpy.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_Axpy.x > output_data_Axpy.m

driver_GemmIPJ_Axpy.x: $(TEST_OBJS_IPJ_AXPY) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_IPJ_AXPY) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmIPJ_Axpy.x $(LDFLAGS) 

test_GemmIPJ_Axpy: driver_GemmIPJ_Axpy.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmIPJ_Axpy.x > output_data_GemmIPJ_Axpy.m

driver_GerJI_Axpy.x: $(TEST_OBJS_GERJI_AXPY) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GERJI_AXPY) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GerJI_Axpy.x $(LDFLAGS) 

test_GerJI_Axpy: driver_GerJI_Axpy.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GerJI_Axpy.x > output_data_GerJI_Axpy.m

driver_GemmJPI_Gemv.x: $(TEST_OBJS_GEMMJPI_GEMV) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMJPI_GEMV) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmJPI_Gemv.x $(LDFLAGS) 

test_GemmJPI_Gemv: driver_GemmJPI_Gemv.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmJPI_Gemv.x > output_data_GemmJPI_Gemv.m

driver_GemmPJI_Ger.x: $(TEST_OBJS_GEMMPJI_GER) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_GEMMPJI_GER) $(UTIL_OBJS) $(BLAS_LIB) -o driver_GemmPJI_Ger.x $(LDFLAGS) 

test_GemmPJI_Ger: driver_GemmPJI_Ger.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_GemmPJI_Ger.x > output_data_GemmPJI_Ger.m


clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
