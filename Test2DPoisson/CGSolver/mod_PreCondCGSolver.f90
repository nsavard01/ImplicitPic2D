module mod_PreCondCGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
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
        allocate(self%D_Vector(matDimension), self%solution(matDimension), self%residual(matDimension))
        ! initialize direct Solver
        self%directSolver = pardisoSolver(self%N_x(self%numberStages) * self%N_y(self%numberStages))
    end function PreCondCGSolver_constructor

    




end module mod_PreCondCGSolver