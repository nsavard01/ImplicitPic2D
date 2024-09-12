FC = ifx
CC=mpiifort
BASE= -Ofast -ipo -qopenmp -qmkl -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
OBJ= main.o lib.a
RM = rm -f

all: clean main

debug: clean setDebug main

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all)
	
main: $(OBJ)
	$(FC) $(BASE) -o $@ $(OBJ)

lib.a: mod_PreCondBiCGSTABSolver.o mod_PreCondCGSolver.o mod_MGSolver.o mod_MG_Stage.o mod_RedBlackSolverCurv.o mod_ZebraSolverCurv.o mod_RedBlackSolverEven.o mod_ZebraSolverEven.o \
mod_GS_Base_Curv.o mod_GS_Base_Even.o mod_GS_Base.o mod_CSRMAtrix.o mod_pardisoSolver.o
	ar rs $@ $?

main.o: main.f90 mod_PreCondCGSolver.o mod_PreCondBiCGSTABSolver.o
	$(FC) $(BASE) -c  -o $@ $<

mod_PreCondBiCGSTABSolver.o: mod_PreCondBiCGSTABSolver.f90 mod_MGSolver.o
	$(FC) $(BASE) -c  -o $@ $<

mod_PreCondCGSolver.o: mod_PreCondCGSolver.f90 mod_MGSolver.o
	$(FC) $(BASE) -c  -o $@ $<

mod_MGSolver.o: mod_MGSolver.f90 mod_MG_Stage.o mod_pardisoSolver.o mod_CSRMAtrix.o
	$(FC) $(BASE) -c  -o $@ $<


mod_MG_Stage.o: mod_MG_Stage.f90 mod_RedBlackSolverCurv.o mod_ZebraSolverCurv.o mod_RedBlackSolverEven.o mod_ZebraSolverEven.o mod_GS_Base_Curv.o mod_GS_Base_Even.o mod_GS_Base.o
	$(FC) $(BASE) -c  -o $@ $<

mod_RedBlackSolverCurv.o: mod_RedBlackSolverCurv.f90 mod_GS_Base_Curv.o
	$(FC) $(BASE) -c  -o $@ $<

mod_RedBlackSolverEven.o: mod_RedBlackSolverEven.f90 mod_GS_Base_Even.o
	$(FC) $(BASE) -c  -o $@ $<

mod_ZebraSolverCurv.o: mod_ZebraSolverCurv.f90 mod_GS_Base_Curv.o
	$(FC) $(BASE) -c  -o $@ $<

mod_ZebraSolverEven.o: mod_ZebraSolverEven.f90 mod_GS_Base_Even.o
	$(FC) $(BASE) -c  -o $@ $<

mod_GS_Base_Even.o: mod_GS_Base_Even.f90 mod_GS_Base.o
	$(FC) $(BASE) -c  -o $@ $<

mod_GS_Base_Curv.o: mod_GS_Base_Curv.f90 mod_GS_Base.o
	$(FC) $(BASE) -c  -o $@ $<

mod_GS_Base.o: mod_GS_Base.f90
	$(FC) $(BASE) -c  -o $@ $<

mod_CSRMAtrix.o: mod_CSRMAtrix.f90
	$(FC) $(BASE) -c  -o $@ $<

mod_pardisoSolver.o: mod_pardisoSolver.f90
	$(FC) $(BASE) -c  -o $@ $<


clean:
	$(RM) main *.o *.mod *.dat





