
all: mod_pardisoSolver.o mod_CSRMAtrix.o mod_GS_Base.o mod_RedBlackSolverEven.o

# mod_PreCondBiCGSTABSolver.o: mod_PreCondBiCGSTABSolver.f90 mod_MGSolver.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_PreCondCGSolver.o: mod_PreCondCGSolver.f90 mod_MGSolver.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_MGSolver.o: mod_MGSolver.f90 mod_MG_Stage.o mod_pardisoSolver.o mod_CSRMAtrix.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_MG_Stage.o: mod_MG_Stage.f90 mod_RedBlackSolverCurv.o mod_ZebraSolverCurv.o mod_RedBlackSolverEven.o mod_ZebraSolverEven.o mod_GS_Base_Curv.o mod_GS_Base_Even.o mod_GS_Base.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_RedBlackSolverCurv.o: mod_RedBlackSolverCurv.f90 mod_GS_Base_Curv.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_RedBlackSolverEven.o: mod_RedBlackSolverEven.f90 mod_GS_Base_Even.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_ZebraSolverCurv.o: mod_ZebraSolverCurv.f90 mod_GS_Base_Curv.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_ZebraSolverEven.o: mod_ZebraSolverEven.f90 mod_GS_Base_Even.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_GS_Base_Even.o: mod_GS_Base_Even.f90 mod_GS_Base.o
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

# mod_GS_Base_Curv.o: mod_GS_Base_Curv.f90 mod_GS_Base.o
# 	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_GS_Base.o: mod_GS_Base.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_CSRMAtrix.o: mod_CSRMAtrix.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<

mod_pardisoSolver.o: mod_pardisoSolver.f90
	$(FC) $(BASE) $(include_dir) -c  -o $@ $<


clean:
	$(RM) *.o *.mod *.dat *.a





