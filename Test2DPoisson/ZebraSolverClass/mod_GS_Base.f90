module mod_GS_Base
    use iso_fortran_env, only: int32, int64, real64
    implicit none


    type :: GS_Base
        ! store grid quantities
        real(real64), allocatable :: sourceTerm(:,:), solution(:,:), residual(:,:)
        integer(int32), allocatable :: horzIndx(:,:), vertIndx(:,:)
        integer(int32) :: numberRows, numberColumns, numberBoundNodes, startRow, endRow, startCol, endCol, iterNumber, N_x, N_y
        integer(int32) :: startRowCoarse, endRowCoarse, startColCoarse, endColCoarse
        real(real64) :: omega
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: solveGS
        procedure, public, pass(self) :: smoothIterations
        procedure, public, pass(self) :: smoothWithRes
        procedure, public, pass(self) :: calcResidual
        procedure, public, pass(self) :: restriction
        procedure, public, pass(self) :: prolongation
        procedure, public, pass(self) :: AX_Mult
        procedure, public, pass(self) :: YAX_Mult
    end type

contains

    subroutine constructPoissonOrthogonal(self, diffX, diffY, NESW_wallBoundaries, boundaryConditions)
        ! Construct orthogonal grid solver
        class(GS_Base), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
        real(real64), intent(in) :: diffX(self%N_x-1), diffY(self%N_y-1)
    end subroutine constructPoissonOrthogonal
    
    subroutine solveGS(self, tol)
        ! Solve GS down to some tolerance
        class(GS_Base), intent(in out) :: self
        real(real64), intent(in) :: tol
        real(real64) :: Res
        Res = 1.0
        self%iterNumber = 0
        do while (Res > tol)
            ! single iterations slightly faster
            call self%smoothIterations(100)
            Res = self%smoothWithRes()
            self%iterNumber = self%iterNumber + 101
        end do

    end subroutine solveGS

    subroutine smoothIterations(self, iterNum)
        ! Solve GS down to some tolerance
        class(GS_Base), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
    end subroutine smoothIterations

    function smoothWithRes(self) result(Res)
        ! Solve GS down to some tolerance
        class(GS_Base), intent(in out) :: self
        real(real64) :: Res
        Res = 0.0d0
    end function smoothWithRes


    subroutine calcResidual(self)
        ! Solve GS down to some tolerance
        class(GS_Base), intent(in out) :: self
    end subroutine calcResidual

    subroutine restriction(self, fineGrid, coarseGrid)
        ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
        class(GS_Base), intent(in) :: self
        real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
        real(real64), intent(in out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    end subroutine restriction

    subroutine prolongation(self, fineGrid, coarseGrid)
        ! Prolongate operator from coarse to fine grid using gauss seidel data
        class(GS_Base), intent(in) :: self
        real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
        real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    end subroutine prolongation

    function YAX_Mult(self, x, y) result(res)
        ! Use gauss-seidel to calculate y^T * A * x
        class(GS_Base), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y), y(self%N_x, self%N_y)
        real(real64) :: res
        res = 0.0d0

    end function YAX_Mult

    subroutine AX_Mult(self, x, y)
        ! Calculate y = A * x
        class(GS_Base), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y)
        real(real64), intent(in out) :: y(self%N_x, self%N_y)

    end subroutine AX_Mult


end module mod_GS_Base