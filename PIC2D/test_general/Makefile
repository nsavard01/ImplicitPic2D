FC = ifx
BASE= -Ofast -ipo -qopenmp -qmkl -xHost -heap-arrays 8192 #-xCORE-AVX512 -mt_mpi 
DEBUG = -qmkl -qopenmp -g -traceback -check all -heap-arrays 8192 -fp-stack-check -warn all #-traceback -CB -fp-stack-check -warn all -check uninit -check pointers -check output_conversion -check format -check stack -fpe0 -check arg_temp_created -g
OPT= $(BASE)

constants_dir = Constants/
OBJ_constants = constants.o
constants_obj_list = $(addprefix $(constants_dir), $(OBJ_constants))

rand_dir = Random/
OBJ_rand = mod_rand_generator.o
rand_obj_list = $(addprefix $(rand_dir), $(OBJ_rand))

world_dir = World/
OBJ_world = mod_domain_curv.o mod_domain_uniform.o mod_domain_base.o
world_obj_list = $(addprefix $(world_dir), $(OBJ_world))

solver_dir = PoissonSolver/
OBJ_solver = mod_PreCondBiCGSTABSolver.o mod_PreCondCGSolver.o mod_MGSolver.o mod_MG_Stage.o mod_ZebraSolverCurv.o mod_RedBlackSolverCurv.o mod_ZebraSolverEven.o mod_RedBlackSolverEven.o mod_GS_Base_Curv.o \
mod_GS_Base_Even.o mod_GS_Base.o mod_CSRMAtrix.o mod_pardisoSolver.o
solver_obj_list = $(addprefix $(solver_dir), $(OBJ_solver))

particle_dir = Particles/
OBJ_particle = mod_particle.o
particle_obj_list = $(addprefix $(particle_dir), $(OBJ_particle))

total_obj_list = $(solver_obj_list) $(particle_obj_list) $(world_obj_list) $(rand_obj_list) $(constants_obj_list)
total_dir_list = $(solver_dir) $(particle_dir) $(world_dir) $(rand_dir) $(constants_dir)

include_dir = $(addprefix -I,  $(total_dir_list))

RM = rm -f

all: main

debug: clean setDebug main

setDebug: 
	@echo "Setting up Debug flags"
	$(eval BASE = -qmkl -qopenmp -g -traceback -check all -fp-stack-check -warn all)
	
main: main.o $(total_obj_list)
	$(FC) $(BASE) -o $@ main.o $(total_obj_list)

main.o: main.f90 solver particle world random constants
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

solver: particle world constants
	cd $(solver_dir); $(MAKE) "FC=$(FC)" "BASE=$(BASE)" "include_dir=$(addprefix -I../, $(total_dir_list))"

particle: world constants random
	cd $(particle_dir); $(MAKE) "FC=$(FC)" "BASE=$(BASE)" "include_dir=$(addprefix -I../, $(total_dir_list))"

world: constants random
	cd $(world_dir); $(MAKE) "FC=$(FC)" "BASE=$(BASE)" "include_dir=$(addprefix -I../, $(total_dir_list))"

random: constants
	cd $(rand_dir); $(MAKE) "FC=$(FC)" "BASE=$(BASE)" "include_dir=$(addprefix -I../, $(total_dir_list))"

constants:
	cd $(constants_dir); $(MAKE) "FC=$(FC)" "BASE=$(BASE)" "include_dir=$(addprefix -I../, $(total_dir_list))"

clean:
	$(RM) main *.o *.mod *.dat *.a
	cd $(constants_dir); $(MAKE) clean
	cd $(solver_dir); $(MAKE) clean
	cd $(world_dir); $(MAKE) clean
	cd $(particle_dir); $(MAKE) clean
	cd $(rand_dir); $(MAKE) clean





