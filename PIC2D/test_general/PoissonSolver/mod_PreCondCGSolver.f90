module mod_PreCondCGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_pardisoSolver
    use mod_MGSolver
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: PreCondCGSolver

    type, extends(MGSolver) :: PreCondCGSolver
        real(real64), allocatable :: D_Vector(:,:), solution(:,:), residual(:,:)
    contains
        procedure, public, pass(self) :: solve => solve_PCG
        procedure, public, pass(self) :: solve_CG
    end type

    interface PreCondCGSolver
        module procedure :: PreCondCGSolver_constructor
    end interface PreCondCGSolver

contains

    type(PreCondCGSolver) function PreCondCGSolver_constructor(N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper) result(self)
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
        allocate(self%D_Vector(self%N_x, self%N_y), self%solution(self%N_x, self%N_y), self%residual(self%N_x, self%N_y))

        !$OMP parallel
        !$OMP workshare
        self%D_Vector = 0.0d0
        !$OMP end workshare
        !$OMP workshare
        self%solution = 0.0d0
        !$OMP end workshare
        !$OMP workshare
        self%residual = 0.0d0
        !$OMP end workshare
        !$OMP end parallel
    end function PreCondCGSolver_constructor

    subroutine solve_PCG(self, stepTol, relTol, intSelect)
        ! Solve MG preconditioned CG with F-cycles (seems to work best for different sized cells in x,y)
        class(PreCondCGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        integer, intent(in) :: intSelect
        real(real64) :: R2_init, resProduct_old, resProduct_new, alpha, beta
        integer(int32) :: i
        
        associate(stageOne => self%MG_smoothers(1)%GS_Smoother)
        
        ! Since an extended type of MG_solver, should store solution always in stageOne solution
        !$OMP parallel workshare
        self%solution = stageOne%solution
        !$OMP end parallel workshare
        ! Calculate initial residual (CG and first smoother solution should already be the same!), store into CG residual
        call stageOne%calcResidual()
        !$OMP parallel workshare
        self%residual = stageOne%residual
        !$OMP end parallel workshare
        
        ! Undergo F_cycle MG with pre, post smoothing
        call stageOne%smoothIterations(self%numberPreSmoothOper) 
        call stageOne%calcResidual()
        if (intSelect == 0) then
            call self%V_Cycle()
        else
            call self%F_Cycle()
        end if
        call stageOne%smoothIterations(self%numberPostSmoothOper) 
        !$OMP parallel
        !$OMP workshare
        ! Store preconditioned residual into smoother solution
        stageOne%solution = stageOne%solution - self%solution
        !$OMP end workshare
        !$OMP workshare
        ! D_vector initialization
        self%D_Vector = stageOne%solution
        !$OMP end workshare
        !$OMP end parallel
        


        !$OMP parallel
        !$OMP workshare
        ! Find initial residual
        R2_init = SUM(self%residual**2)
        !$OMP end workshare
        !$OMP workshare
        ! initial residual x preconditioned residual
        resProduct_old = SUM(self%residual * stageOne%solution)
        !$OMP end workshare
        !$OMP end parallel
        R2_init = sqrt(R2_init)

        do i = 1, self%maxIter
            ! Calculate denominator of alpha = D^T * A  * D
            alpha = stageOne%YAX_Mult(self%D_Vector, self%D_Vector)
            alpha = resProduct_old/alpha

            !$OMP parallel
            !$OMP workshare
            ! Update solution
            self%solution = self%solution + alpha * self%D_Vector
            !$OMP end workshare
            !$OMP workshare
            ! Transfer to smoother and calculate residuals
            stageOne%solution = self%solution
            !$OMP end workshare
            !$OMP end parallel

            call stageOne%calcResidual()
            !$OMP parallel
            !$OMP workshare
            self%R2_current = SUM(stageOne%residual**2)
            !$OMP end workshare
            !$OMP end parallel
            self%R2_current = sqrt(self%R2_current)
            if (self%R2_current/R2_init < relTol) then
                ! exit for reaching relative tolerance
                exit
            end if

            ! Save non-preconditioned residual 
            !$OMP parallel workshare
            self%residual = stageOne%residual
            !$OMP end parallel workshare

            ! MG preconditioning
            call stageOne%smoothIterations(self%numberPreSmoothOper) 
            call stageOne%calcResidual()
            if (intSelect == 0) then
                call self%V_Cycle()
            else
                call self%F_Cycle()
            end if
            call stageOne%smoothIterations(self%numberPostSmoothOper-1)
            self%stepResidual = stageOne%smoothWithRes()

            if (self%stepResidual < stepTol) then
                exit
            end if

            
            !$OMP parallel 
            !$OMP workshare
            ! Store preconditioned residual into smoother solution
            stageOne%solution = stageOne%solution - self%solution
            !$OMP end workshare
            !$OMP workshare
            ! Solve for new residual product
            resProduct_new = SUM(self%residual * stageOne%solution)
            !$OMP end workshare
            !$OMP end parallel

            ! Solve for beta
            beta = resProduct_new/resProduct_old

            !$OMP parallel workshare
            ! Update D_vector
            self%D_Vector = stageOne%solution + beta * self%D_Vector
            !$OMP end parallel workshare

            ! set old resProduct to new
            resProduct_old = resProduct_new

        end do
        end associate
        self%numIter = i-1
    end subroutine solve_PCG

    subroutine solve_CG(self, stepTol, relTol)
        ! Solve MG preconditioned CG with F-cycles (seems to work best for different sized cells in x,y)
        class(PreCondCGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        real(real64) :: R2_init, resProduct_old, resProduct_new, alpha, beta
        integer(int32) :: i
        print *, 'entering CG'
        associate(stageOne => self%MG_smoothers(1)%GS_Smoother)
        
        ! Since an extended type of MG_solver, should store solution always in stageOne solution
        ! Calculate initial residual (CG and first smoother solution should already be the same!), store into CG residual
        call stageOne%calcResidual()
        !$OMP parallel workshare
        self%residual = stageOne%residual
        !$OMP end parallel workshare
        
        !$OMP parallel
        !$OMP workshare
        ! D_vector initialization
        self%D_Vector = stageOne%residual
        !$OMP end workshare
        !$OMP end parallel
        
        !$OMP parallel
        !$OMP workshare
        ! Find initial residual
        R2_init = SUM(stageOne%residual**2)
        !$OMP end workshare
        !$OMP end parallel
        self%R2_current = R2_init
        do i = 1, self%maxIter
            ! Calculate denominator of alpha = D^T * A  * D
            alpha = stageOne%YAX_Mult(self%D_Vector, self%D_Vector)
            alpha = self%R2_current/alpha

            !$OMP parallel
            !$OMP workshare
            ! Update solution
            stageOne%solution = stageOne%solution + alpha * self%D_Vector
            !$OMP end workshare
            !$OMP end parallel
            self%stepResidual = abs(alpha) * sqrt(SUM(self%D_Vector**2)/self%N_x/self%N_y)
            if (self%stepResidual < stepTol) then
                print *, 'reached step tol'
                exit
            end if

            call stageOne%calcResidual()
            !$OMP parallel
            !$OMP workshare
            resProduct_new = SUM(stageOne%residual**2)
            !$OMP end workshare
            !$OMP end parallel

            if (sqrt(resProduct_new/R2_init) < relTol) then
                ! exit for reaching relative tolerance
                print *, 'reached rel tolerance'
                exit
            end if

            ! Solve for beta
            beta = resProduct_new/self%R2_current

            !$OMP parallel workshare
            ! Update D_vector
            self%D_Vector = stageOne%residual + beta * self%D_Vector
            !$OMP end parallel workshare

            ! set old resProduct to new
            self%R2_current = resProduct_new

        end do
        end associate
        self%numIter = i-1
    end subroutine solve_CG

    




end module mod_PreCondCGSolver