HOME      := /Users/rvdg
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
CC         := gcc
LINKER     := $(CC)
CFLAGS     := -O2 -Wall -I$(BLAS_INC) -m64 -msse3 -std=c99
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 5#       number of times each experiment is repeated.  The best time is reported.
NFIRST     := 32#     smallest size to be timed.
NLAST      := 1000#     largest size to be timed.
NINC       := 32#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := FLA_Clock.o MaxAbsDiff.o RandomMatrix.o

TEST_OBJS_STEP_1  := driver_step_1.o GemmIJP.o 

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $n@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver_step_1.x: $(TEST_OBJS_STEP_1) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_1) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_1.x $(LDFLAGS) 

test_step_1:   driver_step_1.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_1.x > data_step_1.m

TEST_OBJS_STEP_2  := driver_step_2.o GemmIJSapdot.o 

driver_step_2.x: $(TEST_OBJS_STEP_2) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_2) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_2.x $(LDFLAGS) 

test_step_2:   driver_step_2.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_2.x > data_step_2.m

TEST_OBJS_STEP_3  := driver_step_3.o GemmIJSapdot4x1.o 

driver_step_3.x: $(TEST_OBJS_STEP_3) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_3) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_3.x $(LDFLAGS) 

test_step_3:   driver_step_3.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_3.x > data_step_3.m

TEST_OBJS_STEP_4  := driver_step_4.o GemmIJSapdot4x1Opt.o 

driver_step_4.x: $(TEST_OBJS_STEP_4) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_4) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_4.x $(LDFLAGS) 

test_step_4:   driver_step_4.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_4.x > data_step_4.m

TEST_OBJS_STEP_5  := driver_step_5.o GemmIJSapdot4x1Reg.o 

driver_step_5.x: $(TEST_OBJS_STEP_5) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_5) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_5.x $(LDFLAGS) 

test_step_5:   driver_step_5.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_5.x > data_step_5.m

TEST_OBJS_STEP_6  := driver_step_6.o GemmIJSapdot4x4Reg.o 

driver_step_6.x: $(TEST_OBJS_STEP_6) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_6) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_6.x $(LDFLAGS) 

test_step_6:   driver_step_6.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_6.x > data_step_6.m

TEST_OBJS_STEP_7  := driver_step_7.o GemmIJSapdot4x4Vec.o 

driver_step_7.x: $(TEST_OBJS_STEP_7) $(UTIL_OBJS)
	$(LINKER) $(TEST_OBJS_STEP_7) $(UTIL_OBJS) $(BLAS_LIB) -o driver_step_7.x $(LDFLAGS) 

test_step_7:   driver_step_7.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver_step_7.x > data_step_7.m




clean:
	rm -f *.o *~ core *.x data_step_*.m *.pdf
