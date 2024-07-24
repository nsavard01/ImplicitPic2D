module mod_MGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_GSSolver
    use mod_PardisoSolver
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: MGSolver

    type :: MGSolver
        type(GSSolver), allocatable :: GS_smoothers(:)
        integer :: numberStages, smoothNumber ! number of stages is total, so includes finest and coarsest grid
        integer, allocatable :: matDimensions(:), N_x(:), N_y(:)
    contains
        procedure, public, pass(self) :: makeSmootherStages
    end type

    interface MGSolver
        module procedure :: MGSolver_constructor
    end interface MGSolver

contains

    type(MGSolver) function MGSolver_constructor(N_x, N_y, numberStages) result(self)
        ! Construct object, set initial variables
        integer(int32), intent(in) :: N_x, N_y, numberStages
        integer :: i
        self%numberStages = numberStages
        self%smoothNumber = numberStages
        allocate(self%GS_smoothers(self%smoothNumber), self%N_x(numberStages), self%N_y(numberStages))
        do i = 1, self%numberStages
            self%N_x(i) = (N_x + (2**(i-1) - 1))/(2**(i-1))
            self%N_y(i) = (N_y + (2**(i-1) - 1))/(2**(i-1))
            print *, 'Stage 1:', self%N_x(i), self%N_y(i)
        end do
    end function MGSolver_constructor

    subroutine makeSmootherStages(self, delX, delY, NESW_wallBoundaries, NESW_phiValues, boundaryConditions, omega)
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x(1)*self%N_y(1))
        real(real64), intent(in) :: NESW_phiValues(4), delX, delY, omega
        real(real64) :: delY_temp, delX_temp
        integer :: stageIter, i, j, k_orig, k_temp, j_temp, i_temp, matDimension
        integer, allocatable :: boundaryConditionsTemp(:)

        ! Construct first stage GS smoother (finest grid)
        self%GS_smoothers(1) = GSSolver(omega)
        call self%GS_smoothers(1)%constructPoissonEven(self%N_x(1), self%N_y(1), delX, delY, NESW_wallBoundaries, NESW_phiValues, boundaryConditions)
        delY_temp = delY
        delX_temp = delX
    
        allocate(boundaryConditionsTemp(self%N_x(2) * self%N_y(2)))
        print *, 'assembling new boundaries and generating new GS_smoothers'
        do stageIter = 2, self%numberStages
            delY_temp = delY_temp * 2.0d0
            delX_temp = delX_temp * 2.0d0
            do j = 1, self%N_y(1), 2**(stageIter-1)
                do i = 1, self%N_x(1), 2**(stageIter-1)
                    k_orig = (j-1) * self%N_x(1) + i
                    i_temp = (i-1)/(2**(stageIter-1)) + 1
                    j_temp = (j-1)/(2**(stageIter-1)) + 1
                    k_temp = (j_temp - 1) * self%N_x(stageIter) + i_temp
                    boundaryConditionsTemp(k_temp) = boundaryConditions(k_orig)
                end do
            end do
            self%GS_smoothers(stageIter) = GSSolver(omega)
            call self%GS_smoothers(stageIter)%constructPoissonEven(self%N_x(stageIter), self%N_y(stageIter), &
                delX_temp, delY_temp, NESW_wallBoundaries, NESW_phiValues, boundaryConditionsTemp(1:(self%N_x(stageIter) * self%N_y(stageIter))))
        end do


    end subroutine makeSmootherStages




end module mod_MGSolver