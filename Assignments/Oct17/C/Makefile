INCLUDE=-I include/
LIBS=lib/libmpi.a lib/libflame-x86_64-5.1.0-36.a lib/libopenblas_sandybridgep-r0.2.20.a -lpthread

test : driver_symv_unb_var4.c
	gcc -g $(INCLUDE) $^ $(LIBS) -o $@.x

clean :
	rm -f test.x

run : test
	./test.x

