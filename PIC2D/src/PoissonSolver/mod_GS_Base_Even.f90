module mod_GS_Base_Even
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base
    implicit none


    type, extends(GS_Base) :: GS_Base_Even
        ! store grid quantities
        real(real64) :: coeffX, coeffY, centerCoeff
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: restriction => restriction_even
        procedure, public, pass(self) :: prolongation => prolongation_even
        procedure, public, pass(self) :: calcResidual => calcResidual_even
        procedure, public, pass(self) :: AX_Mult => AX_Mult_even
        procedure, public, pass(self) :: YAX_Mult => YAX_Mult_even
    end type

contains

    subroutine constructPoissonOrthogonal(self, diffX, diffY, NESW_wallBoundaries, boundaryConditions)
        ! Construct orthogonal grid solver
        class(GS_Base_Even), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
        real(real64), intent(in) :: diffX(self%N_x-1), diffY(self%N_y-1)
    end subroutine constructPoissonOrthogonal

    subroutine AX_Mult_even(self, x, y)
        ! Use gauss-seidel to calculate x^T * A * x
        class(GS_Base_Even), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y)
        real(real64), intent(in out) :: y(self%N_x, self%N_y)
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j, p
        real(real64) :: inv_centerCoeff
        inv_centerCoeff = 1.0d0/self%centerCoeff
        !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx)
        ! loop through inner nodes
        !$OMP do collapse(2)
        do k = 1, self%numberRows
            do p = 1, self%numberColumns
                i = self%startCol + p - 1
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1,k)
                E_indx = self%horzIndx(1,p)
                S_indx = self%vertIndx(2,k)
                W_indx = self%horzIndx(2,p)
                y(i,j) = self%coeffY * (x(i,N_indx) + x(i,S_indx)) + &
                    self%coeffX * (x(E_indx, j) + x(W_indx, j)) + x(i,j)*inv_centerCoeff
            end do
        end do
        !$OMP end do
        !$OMP end parallel 
    end subroutine AX_Mult_even

    function YAX_Mult_even(self, x, y) result(res)
        ! Use gauss-seidel to calculate x^T * A * x
        class(GS_Base_Even), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y), y(self%N_x, self%N_y)
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j, p
        real(real64) :: res, inv_centerCoeff
        inv_centerCoeff = 1.0d0/self%centerCoeff
        res = 0.0d0
        !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx) reduction(+:res)
        ! loop through inner nodes
        !$OMP do collapse(2)
        do k = 1, self%numberRows
            do p = 1, self%numberColumns
                i = self%startCol + p - 1
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1,k)
                E_indx = self%horzIndx(1,p)
                S_indx = self%vertIndx(2,k)
                W_indx = self%horzIndx(2,p)
                res = res + y(i,j) * (self%coeffY * (x(i,N_indx) + x(i,S_indx)) + &
                    self%coeffX * (x(E_indx, j) + x(W_indx, j)) + x(i,j)*inv_centerCoeff)
            end do
        end do
        !$OMP end do
        !$OMP end parallel 
    end function YAX_Mult_even

subroutine calcResidual_even(self)
    ! Solve GS down to some tolerance
    class(GS_Base_Even), intent(in out) :: self
    integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p

  
    !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
    !$OMP&  W_indx)
    ! loop through inner nodes
    !$OMP do collapse(2)
    do k = 1, self%numberRows
        do p = 1, self%numberColumns
            i = self%startCol + p - 1
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1,k)
            E_indx = self%horzIndx(1,p)
            S_indx = self%vertIndx(2,k)
            W_indx = self%horzIndx(2,p)
            self%residual(i, j) = self%sourceTerm(i,j) - (self%solution(i,N_indx) + self%solution(i,S_indx)) * self%coeffY &
                - (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX - self%solution(i,j)/self%centerCoeff
        end do
    end do
    !$OMP end do
    !$OMP end parallel
end subroutine calcResidual_even

subroutine restriction_even(self, fineGrid, coarseGrid)
    ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
    class(GS_Base_Even), intent(in) :: self
    real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
    real(real64), intent(in out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k, p

    !$OMP parallel private(k, p, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx)
    !$OMP do collapse(2)
    ! loop over fine grid nodes which overlap with coarse grid nodes
    do k = self%startRowCoarse, self%endRowCoarse, 2
        do p = self%startColCoarse, self%endColCoarse, 2
            ! calculate fine indices around overlapping index
            j_fine = self%startRow + k - 1
            j_coarse = (j_fine + 1)/2
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            i_fine = self%startCol + p - 1
            i_coarse = (i_fine + 1)/2
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)

            ! Bilinear interpolation
            coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
            0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
            fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
            0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
            fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
        end do
    end do
    !$OMP end do
    !$OMP end parallel
end subroutine restriction_even

subroutine prolongation_even(self, fineGrid, coarseGrid)
    ! Prolongate operator from coarse to fine grid using gauss seidel data
    class(GS_Base_Even), intent(in) :: self
    real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
    real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse
    N_x_coarse = (self%N_x+1)/2
    N_y_coarse = (self%N_y+1)/2
    ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
    !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine)
    ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
    !$OMP do collapse(2)
    ! set similar nodes
    do j_coarse = 1, N_y_coarse
        do i_coarse = 1, N_x_coarse
            i_fine = i_coarse*2 - 1
            j_fine = j_coarse*2-1
            ! Add overlapping coarse nodes
            fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
        end do
    end do
    !$OMP end do nowait
    !$OMP do collapse(2)
    ! go horizontal each row,
    do j_coarse = 1, N_y_coarse
        do i_coarse = 1, N_x_coarse-1
            i_fine = i_coarse*2 - 1
            j_fine = j_coarse*2 - 1

            ! simple interpolation
            fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
        end do
    end do
    !$OMP end do nowait
    !$OMP do collapse(2)
    ! go vertical
    do j_coarse = 1, N_y_coarse-1
        do i_coarse = 1, N_x_coarse
            i_fine = i_coarse*2 - 1
            j_fine = j_coarse*2 - 1
            ! Simple interpolation
            fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
        end do
    end do
    !$OMP end do nowait
    !$OMP do collapse(2)
    ! add offset fine nodes which require 4 point interpolation
    ! first loop through SW corner coarse nodes
    do j_coarse = 1, N_y_coarse-1
        do i_coarse = 1, N_x_coarse-1
            i_fine = i_coarse*2 - 1
            j_fine = j_coarse*2 - 1
            ! simple interpolation
            fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse) * 0.25d0 &
                + coarseGrid(i_coarse, j_coarse+1) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse+1) * 0.25d0
        end do
    end do
    !$OMP end do
    !$OMP end parallel
end subroutine prolongation_even


end module mod_GS_Base_Even