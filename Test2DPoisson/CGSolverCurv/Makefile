FC = ifx
CC=mpiifort
BASE= -O3 -ipo -qopenmp -qmkl -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= main.o mod_GSSolver.o
RM = rm -f

all: clean main

debug: clean setDebug main

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all)
	
main: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

main.o: main.f90 mod_GSSolver.o
	$(FC) $(BASE) -c  -o $@ $<

# mod_PreCondCGSolver.o: mod_PreCondCGSolver.f90 mod_MGSolver.o mod_GSSolver.o mod_pardisoSolver.o mod_CSRMAtrix.o
# 	$(FC) $(BASE) -c  -o $@ $<

# mod_MGSolver.o: mod_MGSolver.f90 mod_GSSolver.o mod_pardisoSolver.o mod_CSRMAtrix.o
# 	$(FC) $(BASE) -c  -o $@ $<

# mod_CSRMAtrix.o: mod_CSRMAtrix.f90
# 	$(FC) $(BASE) -c  -o $@ $<

mod_GSSolver.o: mod_GSSolver.f90
	$(FC) $(BASE) -c  -o $@ $<

# mod_pardisoSolver.o: mod_pardisoSolver.f90
# 	$(FC) $(BASE) -c  -o $@ $<


clean:
	$(RM) main *.o *.mod *.dat





