COMPILER_C = icx
BASE_C = -O3 -ipo -xHost -fast -qmkl -qopenmp

all: mod_rand_generator.o PCG_func.o

mod_rand_generator.o: mod_rand_generator.f90 PCG_func.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

PCG_func.o: PCG_func.c
	$(COMPILER_C) $(BASE_C) -c -o $@ $<	

clean:
	$(RM) *.o *.mod *.dat *.a





