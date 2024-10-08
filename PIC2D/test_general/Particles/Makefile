
all: mod_particle.o

mod_particle.o: mod_particle.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

clean:
	$(RM) *.o *.mod *.dat *.a





