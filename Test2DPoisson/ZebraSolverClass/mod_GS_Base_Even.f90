module mod_GS_Base_Even
    use iso_fortran_env, only: int32, int64, real64
    implicit none


    type, extends(GS_Base) :: GS_Base_Even
        ! store grid quantities
        real(real64) :: coeffX, coeffY, centerCoeff
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: solveGS
        procedure, public, pass(self) :: smoothIterations
        procedure, public, pass(self) :: smoothWithRes
        procedure, public, pass(self) :: calcResidual
        procedure, public, pass(self) :: restriction
        procedure, public, pass(self) :: prolongation
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type

contains


    subroutine restriction(self, fineGrid, coarseGrid)
        ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
        class(GS_Base_Even), intent(in) :: self
        real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
        real(real64), intent(in out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    end subroutine restriction

    subroutine prolongation(self, fineGrid, coarseGrid)
        ! Prolongate operator from coarse to fine grid using gauss seidel data
        class(GS_Base_Even), intent(in) :: self
        real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
        real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    end subroutine prolongation


end module mod_GS_Base_Even