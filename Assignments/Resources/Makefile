INCLUDE=-I include/
LIBS=lib/libmpi.a lib/libflame-x86_64-5.1.0-36.a lib/libopenblas_sandybridgep-r0.2.20.a -lpthread

test : test.c
	gcc -g $(INCLUDE) $^ $(LIBS) -o $@

clean :
	rm -f test

run : test
	bin/mpirun -n 2 ./test
