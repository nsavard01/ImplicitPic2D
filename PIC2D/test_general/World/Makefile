
all: mod_domain_uniform.o mod_domain_curv.o

mod_domain_curv.o: mod_domain_curv.f90 mod_domain_base.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_domain_uniform.o: mod_domain_uniform.f90 mod_domain_base.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_domain_base.o: mod_domain_base.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

clean:
	$(RM) *.o *.mod *.dat *.a





