module mod_MG_Stage
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_ZebraSolverEven
    use mod_RedBlackSolverEven
    use mod_ZebraSolverCurv
    use mod_RedBlackSolverCurv
    implicit none

    ! Stage for each multigrid, contains gauss-seidel smoother for even grid

    private
    public :: MG_Stage

    type :: MG_Stage
        class(GS_Base), allocatable :: GS_Smoother
    end type

    interface MG_Stage
        module procedure :: MG_Stage_constructor
    end interface MG_Stage

contains 
    type(MG_Stage) function MG_Stage_constructor(omega, world, N_x, N_y, redBlackBool) result(self)
    ! Construct object, set initial variables
    real(real64), intent(in) :: omega
    integer(int32), intent(in) :: N_x, N_y
    class(domain_base), intent(in), target :: world
    logical, intent(in) :: redBlackBool
    select type (world)
    type is (domain_uniform)
        if (redBlackBool) then
            self%GS_smoother = RedBlackSolverEven(omega, world, N_x, N_y)
        else
            self%GS_smoother = ZebraSolverEven(omega, world, N_x, N_y)
        end if
    type is (domain_curv)
        if (redBlackBool) then
            self%GS_smoother = RedBlackSolverCurv(omega, world, N_x, N_y)
        else
            self%GS_smoother = ZebraSolverCurv(omega, world, N_x, N_y)
        end if
    end select
    end function MG_Stage_constructor


end module mod_MG_Stage