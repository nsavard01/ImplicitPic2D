FC = ifx
BASE= -Ofast -ipo -qopenmp -qmkl -xHost #-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)
include_dir = $(addprefix -I, $(solver_dir))
solver_dir = PoissonSolver/
OBJ_solver = mod_PreCondBiCGSTABSolver.o mod_PreCondCGSolver.o mod_MGSolver.o mod_MG_Stage.o mod_RedBlackSolverCurv.o mod_ZebraSolverCurv.o mod_RedBlackSolverEven.o mod_ZebraSolverEven.o \
mod_GS_Base_Curv.o mod_GS_Base_Even.o mod_GS_Base.o mod_CSRMAtrix.o mod_pardisoSolver.o
solver_name_list = $(addprefix $(solver_dir), $(OBJ_solver))
RM = rm -f

all: main

debug: clean setDebug main

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all)
	
main: main.o $(solver_name_list)
	$(FC) $(BASE) -o $@ main.o $(solver_name_list)

main.o: main.f90 solver
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

solver: 
	cd $(solver_dir); $(MAKE) "FC=$(FC)" "BASE=$(BASE)"


clean:
	$(RM) main *.o *.mod *.dat *.a
	cd $(solver_dir); $(MAKE) clean





