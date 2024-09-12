module mod_PreCondCGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_PardisoSolver
    use mod_GSSolver
    use mod_MGSolver
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: PreCondCGSolver

    type, extends(MGSolver) :: PreCondCGSolver
        real(real64), allocatable :: D_Vector(:), solution(:), residual(:)
    contains
        procedure, public, pass(self) :: solve
    end type

    interface PreCondCGSolver
        module procedure :: PreCondCGSolver_constructor
    end interface PreCondCGSolver

contains

    type(PreCondCGSolver) function PreCondCGSolver_constructor(N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper) result(self)
        ! Construct object, set initial variables, same as MG_Solver since it extends it
        integer(int32), intent(in) :: N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper
        integer :: i, matDimension
        self%numberStages = numberStages
        self%smoothNumber = numberStages-1 ! number smoothing GS stages (excluse coarsest grid)
        self%maxIter = maxIter
        self%numberPreSmoothOper = numberPreSmoothOper
        self%numberPostSmoothOper = numberPostSmoothOper
        self%numIter = 0
        ! allocate smoothers and dimension
        allocate(self%GS_smoothers(self%smoothNumber), self%N_x(numberStages), self%N_y(numberStages))
        do i = 1, self%numberStages
            self%N_x(i) = (N_x + (2**(i-1) - 1))/(2**(i-1))
            self%N_y(i) = (N_y + (2**(i-1) - 1))/(2**(i-1))
            print *, 'Stage ', i, ':', self%N_x(i), self%N_y(i)
        end do
        matDimension = N_x * N_y
        ! allocate specific vectors to CG method
        allocate(self%D_Vector(matDimension), self%solution(matDimension), self%residual(matDimension))
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
        ! initialize direct Solver
        self%directSolver = pardisoSolver(self%N_x(self%numberStages) * self%N_y(self%numberStages))
    end function PreCondCGSolver_constructor

    subroutine solve(self, stepTol, relTol)
        ! Solve MG preconditioned CG with F-cycles (seems to work best for different sized cells in x,y)
        class(PreCondCGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        real(real64) :: R2_init, resProduct_old, resProduct_new, alpha, beta
        integer(int32) :: i

        ! Calculate initial residual (CG and first smoother solution should already be the same!), store into CG residual
        call self%GS_smoothers(1)%calcResidual()
        !$OMP parallel workshare
        self%residual = self%GS_smoothers(1)%residual
        !$OMP end parallel workshare
        
        ! Undergo F_cycle MG with pre, post smoothing
        call self%GS_smoothers(1)%smoothIterations(self%numberPreSmoothOper, .false.) 
        call self%GS_smoothers(1)%calcResidual()
        call self%F_Cycle()
        self%stepResidual = self%GS_smoothers(1)%smoothIterationsWithRes(self%numberPostSmoothOper) 
        !$OMP parallel
        !$OMP workshare
        ! Store preconditioned residual into smoother solution
        self%GS_smoothers(1)%solution = self%GS_smoothers(1)%solution - self%solution
        !$OMP end workshare
        !$OMP workshare
        ! D_vector initialization
        self%D_Vector = self%GS_smoothers(1)%solution
        !$OMP end workshare
        !$OMP end parallel
        


        !$OMP parallel
        !$OMP workshare
        ! Find initial residual
        R2_init = SUM(self%residual**2)
        !$OMP end workshare
        !$OMP workshare
        ! initial residual x preconditioned residual
        resProduct_old = SUM(self%residual * self%GS_smoothers(1)%solution)
        !$OMP end workshare
        !$OMP end parallel

        do i = 1, self%maxIter
            ! Calculate denominator of alpha = D^T * A  * D
            alpha = self%GS_smoothers(1)%XAX_Mult(self%D_Vector)
            alpha = resProduct_old/alpha

            !$OMP parallel
            !$OMP workshare
            ! Update solution
            self%solution = self%solution + alpha * self%D_Vector
            !$OMP end workshare
            !$OMP workshare
            ! Transfer to smoother and calculate residuals
            self%GS_smoothers(1)%solution = self%solution
            !$OMP end workshare
            !$OMP end parallel

            call self%GS_smoothers(1)%calcResidual
            !$OMP parallel
            !$OMP workshare
            self%R2_current = SUM(self%GS_smoothers(1)%residual**2)
            !$OMP end workshare
            !$OMP end parallel

            
            if (self%R2_current/R2_init < relTol) then
                ! exit for reaching relative tolerance
                exit
            end if

            ! Save non-preconditioned residual 
            !$OMP parallel workshare
            self%residual = self%GS_smoothers(1)%residual
            !$OMP end parallel workshare

            ! MG preconditioning
            call self%GS_smoothers(1)%smoothIterations(self%numberPreSmoothOper, .false.) 
            call self%GS_smoothers(1)%calcResidual()
            call self%F_Cycle()
            self%stepResidual = self%GS_smoothers(1)%smoothIterationsWithRes(self%numberPostSmoothOper) 

            if (self%stepResidual < stepTol) then
                ! If step size in smoother is sufficiently small, store into solution and exit
                !$OMP parallel workshare
                self%solution = self%GS_smoothers(1)%solution
                !$OMP end parallel workshare
                exit
            end if

            
            !$OMP parallel 
            !$OMP workshare
            ! Store preconditioned residual into smoother solution
            self%GS_smoothers(1)%solution = self%GS_smoothers(1)%solution - self%solution
            !$OMP end workshare
            !$OMP workshare
            ! Solve for new residual product
            resProduct_new = SUM(self%residual * self%GS_smoothers(1)%solution)
            !$OMP end workshare
            !$OMP end parallel

            ! Solve for beta
            beta = resProduct_new/resProduct_old

            !$OMP parallel workshare
            ! Update D_vector
            self%D_Vector = self%GS_smoothers(1)%solution + beta * self%D_Vector
            !$OMP end parallel workshare

            ! set old resProduct to new
            resProduct_old = resProduct_new

        end do
        self%numIter = i-1
    end subroutine solve

    




end module mod_PreCondCGSolver