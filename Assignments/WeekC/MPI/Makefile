#HOME      := /Users/rvdg
# Make sure you have BLIS installed in your home directory.
BLAS_LIB  := $(HOME)/blis/lib/libblis.a
BLAS_INC  := $(HOME)/blis/include/blis

# indicate how the object files are to be created
MPICC	   := mpicc
LINKER     := $(MPICC)
CFLAGS     := -O2 -Wall -I$(BLAS_INC) -m64 -msse3 -std=c99 -march=native -D_POSIX_C_SOURCE=200112L
FFLAGS     := $(CFLAGS) 

LDFLAGS    := -lpthread -m64 -lm 

UTIL_OBJS  := 

# ---------------------

HelloWorld:	HelloWorld.c
	mpicc HelloWorld.c -o HelloWorld.x
	./HelloWorld.x

# ---------------------

clean:
	rm -f *.o *~ core *.x output_data_*.m *.pdf
