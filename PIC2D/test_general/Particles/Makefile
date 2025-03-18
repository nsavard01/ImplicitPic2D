
all: mod_particle_per_cell.o mod_particle_contiguous.o mod_particle.o

mod_particle_per_cell.o: mod_particle_per_cell.f90 mod_particle.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_particle_contiguous.o: mod_particle_contiguous.f90 mod_particle.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_particle.o: mod_particle.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

clean:
	$(RM) *.o *.mod *.dat *.a





