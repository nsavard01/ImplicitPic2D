
all: constants.o

constants.o: constants.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

clean:
	$(RM) *.o *.mod *.dat *.a





