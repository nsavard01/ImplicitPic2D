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
        real(real64), allocatable :: t(:,:), v(:,:), p_vector(:,:), initial_res(:,:)
    contains
        procedure, public, pass(self) :: solve_BiCGSTAB
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
        allocate(self%t(self%N_x, self%N_y), self%p_vector(self%N_x, self%N_y), self%v(self%N_x, self%N_y), self%initial_res(self%N_x, self%N_y))
        self%t = 0
        self%p_vector = 0
        self%v = 0
        self%initial_res = 0
        
    end function BiCGSTAB_Solver_constructor


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