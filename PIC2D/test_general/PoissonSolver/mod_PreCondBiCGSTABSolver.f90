module mod_PreCondBiCGSTABSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_pardisoSolver
    use mod_MGSolver
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: BiCGSTAB_Solver

    type, extends(MGSolver) :: BiCGSTAB_Solver
        real(real64), allocatable :: t(:,:), v(:,:), z(:,:), p_vector(:,:), initial_res(:,:), solution(:,:), residual(:,:)
    contains
        procedure, public, pass(self) :: solve_BiCGSTAB
        procedure, public, pass(self) :: solve => solve_PreCond_BiCGSTAB
    end type

    interface BiCGSTAB_Solver
        module procedure :: BiCGSTAB_Solver_constructor
    end interface BiCGSTAB_Solver

contains

    type(BiCGSTAB_Solver) function BiCGSTAB_Solver_constructor(N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper) result(self)
        ! Construct object, set initial variables, same as MG_Solver since it extends it
        integer(int32), intent(in) :: N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper
        integer :: N_y_coarse, N_x_coarse
        self%numberStages = numberStages
        self%smoothNumber = numberStages-1 ! number smoothing GS stages (excluse coarsest grid)
        self%maxIter = maxIter
        self%numberPreSmoothOper = numberPreSmoothOper
        self%numberPostSmoothOper = numberPostSmoothOper
        self%numIter = 0
        self%N_x = N_x
        self%N_y = N_y
        allocate(self%MG_smoothers(self%smoothNumber))
        N_x_coarse = (self%N_x + (2**(self%numberStages-1) - 1))/(2**(self%numberStages-1))
        N_y_coarse = (self%N_y + (2**(self%numberStages-1) - 1))/(2**(self%numberStages-1))
        self%directSolver = pardisoSolver(N_x_coarse * N_y_coarse)
        allocate(self%t(self%N_x, self%N_y), self%p_vector(self%N_x, self%N_y), self%v(self%N_x, self%N_y), &
        self%initial_res(self%N_x, self%N_y), self%solution(self%N_x, self%N_y), self%residual(self%N_x, self%N_y), &
        self%z(self%N_x, self%N_y))
        self%t = 0
        self%p_vector = 0
        self%v = 0
        self%z = 0
        self%initial_res = 0
        self%solution = 0.0d0
        self%residual = 0.0d0

    end function BiCGSTAB_Solver_constructor

    subroutine solve_PreCond_BiCGSTAB(self, stepTol, relTol, intSelect)
        ! Solve MG preconditioned CG with F-cycles (seems to work best for different sized cells in x,y)
        class(BiCGSTAB_Solver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        integer, intent(in) :: intSelect
        real(real64) :: R2_init, R2_new, param_old, param_new, alpha, beta, w
        integer(int32) :: i
        print *, 'entering BiCGSTAB'

        associate(stageOne => self%MG_smoothers(1)%GS_Smoother)
        
        ! Since an extended type of MG_solver, should store solution always in stageOne solution at the end
        ! Calculate initial residual
        call stageOne%calcResidual()
        !$OMP parallel
        !$OMP workshare
        !set intial residual vector
        self%initial_res = stageOne%residual
        !$OMP end workshare
        !$OMP workshare
        !get p_0
        param_old = SUM(stageOne%residual**2)
        !$OMP end workshare
        !$OMP workshare
        !get p_vector
        self%p_vector = stageOne%residual
        !$OMP end workshare
        !$OMP workshare
        !store initial solution
        self%solution = stageOne%solution
        !$OMP end workshare
        !$OMP workshare
        ! set initial residual
        self%residual = stageOne%residual
        !$OMP end workshare
        !$OMP end parallel
        
        R2_init = param_old
        do i = 1, self%maxIter

            ! First preconditioning on p_vector
            !$OMP parallel
            !$OMP workshare
            ! Set residual to source term
            stageOne%sourceTerm = self%p_vector
            !$OMP end workshare
            !$OMP workshare
            ! Set solution to 0
            stageOne%solution = 0.0d0
            !$OMP end workshare
            !$OMP end parallel
            call stageOne%smoothIterations(self%numberPreSmoothOper) 
            call stageOne%calcResidual()
            if (intSelect == 0) then
                call self%V_Cycle()
            else
                call self%F_Cycle()
            end if
            call stageOne%smoothIterations(self%numberPostSmoothOper)
        
            ! Calculate v from preconditioned p_vector
            call stageOne%AX_Mult(stageOne%solution, self%v)
            !$OMP parallel
            !$OMP workshare
            ! calculate alpha denominator
            alpha = SUM(self%initial_res * self%v)
            !$OMP end workshare
            !$OMP end parallel
            alpha = param_old / alpha

            !$OMP parallel
            !$OMP workshare
            ! update solution
            self%solution = self%solution + alpha * stageOne%solution
            !$OMP end workshare
            !$OMP workshare
            ! get current step residual
            self%stepResidual = SUM(stageOne%solution**2)
            !$OMP end workshare
            !$OMP end parallel
            self%stepResidual = alpha * SQRT(self%stepResidual/stageOne%numberBoundNodes)

            if (self%stepResidual < stepTol) then
                print *, 'step tolerance low enough first time'
                !$OMP parallel
                !$OMP workshare
                ! store solution into stageOne
                stageOne%solution = self%solution
                !$OMP end workshare
                !$OMP end parallel
                exit
            end if
            
            !$OMP parallel
            !$OMP workshare
            ! update residual
            self%residual = self%residual - alpha * self%v
            !$OMP end workshare
            !$OMP workshare
            ! get new residual norm
            R2_new = SUM(self%residual**2)
            !$OMP end workshare
            !$OMP end parallel
           
            if (R2_new/R2_init < relTol) then
                print *, 'relative tolerance low enough first time'
                !$OMP parallel
                !$OMP workshare
                ! store solution into stageOne
                stageOne%solution = self%solution
                !$OMP end workshare
                !$OMP end parallel
                exit
            end if

            ! ------------------------ Second round of operations ----------------------
            ! Second preconditioning on new residual
            !$OMP parallel
            !$OMP workshare
            ! Set residual to source term
            stageOne%sourceTerm = self%residual
            !$OMP end workshare
            !$OMP workshare
            ! Set solution to 0
            stageOne%solution = 0.0d0
            !$OMP end workshare
            !$OMP end parallel
            call stageOne%smoothIterations(self%numberPreSmoothOper) 
            call stageOne%calcResidual()
            if (intSelect == 0) then
                call self%V_Cycle()
            else
                call self%F_Cycle()
            end if
            call stageOne%smoothIterations(self%numberPostSmoothOper)
            !$OMP parallel workshare
            ! store preconditioned residual into z vector
            self%z = stageOne%solution
            !$OMP end parallel workshare
            ! Get t vector from A*Z
            call stageOne%AX_Mult(self%z, self%t)

            ! preconditioning on t vector
            !$OMP parallel
            !$OMP workshare
            ! Set t vector to source term
            stageOne%sourceTerm = self%t
            !$OMP end workshare
            !$OMP workshare
            ! Set solution to 0
            stageOne%solution = 0.0d0
            !$OMP end workshare
            !$OMP end parallel
            call stageOne%smoothIterations(self%numberPreSmoothOper) 
            call stageOne%calcResidual()
            if (intSelect == 0) then
                call self%V_Cycle()
            else
                call self%F_Cycle()
            end if
            call stageOne%smoothIterations(self%numberPostSmoothOper)

            ! get reals to computer w
            !$OMP parallel
            !$OMP workshare
            ! numerator of w is preconditioned t * z vector, use param_new as temp variable
            param_new = SUM(stageOne%solution * self%z)
            !$OMP end workshare
            !$OMP workshare
            ! denominator is preconditioned t squared
            w = SUM(stageOne%solution**2)
            !$OMP end workshare
            !$OMP end parallel
            w = param_new/w

            !$OMP parallel
            !$OMP workshare
            ! update solution
            self%solution = self%solution + w * self%z
            !$OMP end workshare
            !$OMP workshare
            ! step residual calculation
            self%stepResidual = SUM(self%z**2)
            !$OMP end workshare
            !$OMP end parallel
            self%stepResidual = w * SQRT(self%stepResidual/stageOne%numberBoundNodes)

            if (self%stepResidual < stepTol) then
                print *, 'step tolerance low enough second time'
                !$OMP parallel
                !$OMP workshare
                ! store solution into stageOne
                stageOne%solution = self%solution
                !$OMP end workshare
                !$OMP end parallel
                exit
            end if

            !$OMP parallel
            !$OMP workshare
            ! update residual
            self%residual = self%residual - w * self%t
            !$OMP end workshare
            !$OMP workshare
            ! get new residual norm
            R2_new = SUM(self%residual**2)
            !$OMP end workshare
            !$OMP end parallel
        
            if (R2_new/R2_init < relTol) then
                print *, 'relative tolerance low enough second time'
                !$OMP parallel
                !$OMP workshare
                ! store solution into stageOne
                stageOne%solution = self%solution
                !$OMP end workshare
                !$OMP end parallel
                exit
            end if
            
            ! ------------------------ Final operations before next iteration ---------------------
            !$OMP parallel
            !$OMP workshare
            param_new = SUM(self%initial_res*self%residual)
            !$OMP end workshare
            !$OMP end parallel

            ! Solve for beta
            beta = (param_new/param_old) * (alpha/w)

            !$OMP parallel workshare
            ! Update p_vector
            self%p_vector = self%residual + beta * (self%p_vector - w * self%v)
            !$OMP end parallel workshare

            ! set new param to old
            param_old = param_new

        end do
        end associate
        self%numIter = i-1
    end subroutine solve_PreCond_BiCGSTAB

    subroutine solve_BiCGSTAB(self, stepTol, relTol)
        ! Solve MG preconditioned CG with F-cycles (seems to work best for different sized cells in x,y)
        class(BiCGSTAB_Solver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        real(real64) :: R2_init, R2_new, param_old, param_new, alpha, beta, w
        integer(int32) :: i
        print *, 'entering BiCGSTAB'

        associate(stageOne => self%MG_smoothers(1)%GS_Smoother)
        
        ! Since an extended type of MG_solver, should store solution always in stageOne solution
        ! Calculate initial residual (CG and first smoother solution should already be the same!), store into CG residual
        call stageOne%calcResidual()
        !$OMP parallel
        !$OMP workshare
        !set intial residual vector
        self%initial_res = stageOne%residual
        !$OMP end workshare
        !$OMP workshare
        !get p_0
        param_old = SUM(stageOne%residual**2)
        !$OMP end workshare
        !$OMP workshare
        !get p_vector
        self%p_vector = stageOne%residual
        !$OMP end workshare
        !$OMP end parallel
        
        R2_init = param_old
        do i = 1, self%maxIter
            ! Calculate v
            call stageOne%AX_Mult(self%p_vector, self%v)
            !$OMP parallel
            !$OMP workshare
            ! calculate alpha
            alpha = SUM(self%initial_res * self%v)
            !$OMP end workshare
            !$OMP end parallel
            alpha = param_old / alpha

            !$OMP parallel
            !$OMP workshare
            ! update solution
            stageOne%solution = stageOne%solution + alpha * self%p_vector
            !$OMP end workshare
            !$OMP workshare
            ! update residual
            stageOne%residual = stageOne%residual - alpha * self%v
            !$OMP end workshare
            !$OMP workshare
            ! get new residual norm
            R2_new = SUM(stageOne%residual**2)
            !$OMP end workshare
            !$OMP end parallel
           
            if (R2_new/R2_init < relTol) then
                print *, 'relative tolerance low enough first time'
                exit
            end if

            ! store As into t
            call stageOne%AX_Mult(stageOne%residual, self%t)
            ! get reals to computer w
            !$OMP parallel
            !$OMP workshare
            param_new = SUM(self%t * stageOne%residual)
            !$OMP end workshare
            !$OMP workshare
            w = SUM(self%t**2)
            !$OMP end workshare
            !$OMP end parallel
            w = param_new/w

            !$OMP parallel
            !$OMP workshare
            ! update solution
            stageOne%solution = stageOne%solution + w * stageOne%residual
            !$OMP end workshare
            !$OMP workshare
            ! update residual
            stageOne%residual = stageOne%residual - w * self%t
            !$OMP end workshare
            !$OMP workshare
            ! get new residual norm
            R2_new = SUM(stageOne%residual**2)
            !$OMP end workshare
            !$OMP end parallel
        
            if (R2_new/R2_init < relTol) then
                print *, 'relative tolerance low enough second time'
                exit
            end if
            
            !$OMP parallel
            !$OMP workshare
            param_new = SUM(self%initial_res*stageOne%residual)
            !$OMP end workshare
            !$OMP end parallel

            ! Solve for beta
            beta = (param_new/param_old) * (alpha/w)

            !$OMP parallel workshare
            ! Update p_vector
            self%p_vector = stageOne%residual + beta * (self%p_vector - w * self%v)
            !$OMP end parallel workshare

            ! set new param to old
            param_old = param_new

        end do
        end associate
        self%numIter = i-1
    end subroutine solve_BiCGSTAB

    




end module mod_PreCondBiCGSTABSolver