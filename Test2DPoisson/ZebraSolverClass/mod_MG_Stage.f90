module mod_MG_Stage
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base
    use mod_ZebraSolverEven
    use mod_ZebraSolverCurv
    use mod_RedBlackSolverEven
    use mod_RedBlackSolverCurv
    implicit none

    ! Stage for each multigrid, contains gauss-seidel smoother

    private
    public :: MG_Stage

    type :: MG_Stage
        class(GS_Base), allocatable :: GS_Smoother
    end type

    interface MG_Stage
        module procedure :: MG_Stage_constructor
    end interface MG_Stage

contains 
    type(MG_Stage) function MG_Stage_constructor(omega, N_x, N_y, diffX, diffY, NESW_wallBoundaries, boundaryConditions, evenGridBool, redBlackBool) result(self)
    ! Construct object, set initial variables
    real(real64), intent(in) :: omega, diffX(N_x-1), diffY(N_y-1)
    integer(int32), intent(in) :: N_x, N_y, NESW_wallBoundaries(4), boundaryConditions(N_x, N_y)
    logical, intent(in) :: evenGridBool, redBlackBool
    if (evenGridBool) then
        if (redBlackBool) then
            self%GS_Smoother = RedBlackSolverEven(omega, N_x, N_y)
        else
            self%GS_Smoother = ZebraSolverEven(omega, N_x, N_y)
        end if
    else
        if (RedBlackBool) then
            self%GS_Smoother = RedBlackSolverCurv(omega, N_x, N_y)
        else
            self%GS_Smoother = ZebraSolverCurv(omega, N_x, N_y)
        end if
    end if
    call self%GS_Smoother%constructPoissonOrthogonal(diffX, diffY, NESW_wallBoundaries, boundaryConditions)
    end function MG_Stage_constructor

end module mod_MG_Stage