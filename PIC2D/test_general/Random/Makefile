
all: mod_rand_generator.o

mod_rand_generator.o: mod_rand_generator.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

clean:
	$(RM) *.o *.mod *.dat *.a





