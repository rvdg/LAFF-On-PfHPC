HOME       := ../../Resources
BLAS_LIB   := $(HOME)/lib/libopenblas_sandybridgep-r0.2.20.a
FLAME_LIB  := $(HOME)/lib/libflame-x86_64-5.1.0-36.a
FLAME_INC  := $(HOME)/include

# indicate where the object files are to be created
CC         := gcc
FC         := gfortran
LINKER     := $(CC)
CFLAGS     := -O3 -Wall -I$(FLAME_INC) -I$(CBLAS_INC) -m64 -msse3 -std=c99 
FFLAGS     := $(CFLAGS) 

# set the range of experiments to be performed
NREPEATS   := 3#       number of times each experiment is repeated.
NFIRST     := 50#     smallest size to be timed.
NLAST      := 4000#    largest size to be timed.
NINC       := 50#     increment between sizes.

LDFLAGS    := -lpthread -m64 -lm 

# indicate where the FLAME include files reside

# add files in which your implementations are coded
TEST_OBJS  := driver.o Symv_unb_var2.o Symv_unb_var4.o

# $%.o: %.c
#	$(CC) $(CFLAGS) -c $< -o $@
# $%.o: %.f
#	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)  $(LAPACK_OBJS)
	$(LINKER) $(TEST_OBJS) $(LAPACK_OBJS) $(FLAME_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(LDFLAGS) -o driver.x

test:   driver.x
	echo "$(NREPEATS) $(NFIRST) $(NLAST) $(NINC)" | ./driver.x > output.m

clean:
	rm -f *.o *~ core *.x output.m
