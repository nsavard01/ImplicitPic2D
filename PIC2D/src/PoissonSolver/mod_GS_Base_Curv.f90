module mod_GS_Base_Curv
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base
    implicit none


    type, extends(GS_Base) :: GS_Base_Curv
        ! store grid quantities
        real(real64), allocatable :: horzCoeffs(:,:), vertCoeffs(:,:), centerCoeffs(:,:) !, rowBoundCoeffs(:,:,:), colBoundCoeffs(:,:, :) ! lower, then upper diagonals
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: restriction => restriction_curv
        procedure, public, pass(self) :: prolongation => prolongation_curv
        procedure, public, pass(self) :: calcResidual => calcResidual_curv
        procedure, public, pass(self) :: YAX_Mult => YAX_Mult_curv
        procedure, public, pass(self) :: AX_Mult => AX_Mult_curv
    end type

contains

subroutine constructPoissonOrthogonal(self, diffX, diffY, NESW_wallBoundaries, boundaryConditions)
    ! Construct orthogonal grid solver
    class(GS_Base_Curv), intent(in out) :: self
    integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
    real(real64), intent(in) :: diffX(self%N_x-1), diffY(self%N_y-1)
end subroutine constructPoissonOrthogonal

subroutine AX_Mult_curv(self, x, y)
    ! Solve GS down to some tolerance
    class(GS_Base_Curv), intent(in out) :: self
    real(real64), intent(in) :: x(self%N_x, self%N_y)
    real(real64), intent(in out) :: y(self%N_x, self%N_y)
    real(real64) :: C_N, C_E, C_O, C_W, C_S
    integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
    !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
    !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W)
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
            C_O = self%centerCoeffs(p,k)
            C_N = self%vertCoeffs(1,k)
            C_E = self%horzCoeffs(1,p)
            C_S = self%vertCoeffs(2,k)
            C_W = self%horzCoeffs(2,p)
            y(i,j) = C_N * x(i,N_indx) + C_S * x(i,S_indx) + &
                    C_E * x(E_indx, j) + C_W * x(W_indx, j) + x(i,j)/C_O
        end do
    end do
    !$OMP end do
    !$OMP end parallel
end subroutine AX_Mult_curv

function YAX_Mult_curv(self, x, y) result(res)
    ! Solve GS down to some tolerance
    class(GS_Base_Curv), intent(in out) :: self
    real(real64), intent(in) :: x(self%N_x, self%N_y), y(self%N_x, self%N_y)
    real(real64) :: res, C_N, C_E, C_O, C_W, C_S
    integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
    res = 0.0d0
    !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
    !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W) reduction(+:res)
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
            C_O = self%centerCoeffs(p,k)
            C_N = self%vertCoeffs(1,k)
            C_E = self%horzCoeffs(1,p)
            C_S = self%vertCoeffs(2,k)
            C_W = self%horzCoeffs(2,p)
            res = res + y(i,j) * (C_N * x(i,N_indx) + C_S * x(i,S_indx) + &
                    C_E * x(E_indx, j) + C_W * x(W_indx, j) + x(i,j)/C_O)
        end do
    end do
    !$OMP end do
    !$OMP end parallel
end function YAX_Mult_curv

subroutine calcResidual_curv(self)
    ! Solve GS down to some tolerance
    class(GS_Base_Curv), intent(in out) :: self
    real(real64) :: C_N, C_E, C_O, C_W, C_S
    integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p

    !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
    !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W)
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
            C_O = self%centerCoeffs(p,k)
            C_N = self%vertCoeffs(1,k)
            C_E = self%horzCoeffs(1,p)
            C_S = self%vertCoeffs(2,k)
            C_W = self%horzCoeffs(2,p)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
    end do
    !$OMP end do nowait
    !$OMP end parallel
end subroutine calcResidual_curv

subroutine restriction_curv(self, fineGrid, coarseGrid)
    ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
    class(GS_Base_Curv), intent(in) :: self
    real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
    real(real64), intent(in out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    real(real64) :: a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S
    integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k, p

    !$OMP parallel private(k, p, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, &
    !$OMP& a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S)
    !$OMP do collapse(2)
    ! loop over fine grid nodes which overlap with coarse grid nodes
    do k = self%startRowCoarse, self%endRowCoarse, 2
        do p = self%startColCoarse, self%endColCoarse, 2
            ! calculate fine indices around overlapping index
            i_fine = self%startCol + p - 1
            j_fine = self%startRow + k - 1
            N_indx = self%vertIndx(1,k)
            E_indx = self%horzIndx(1,p)
            S_indx = self%vertIndx(2,k)
            W_indx = self%horzIndx(2,p)
            i_coarse = (i_fine + 1)/2
            j_coarse = (j_fine + 1)/2

            ! Transpose of prolongation, seems to work best for red black
            if (E_indx /= W_indx) then
                ! interpolation from east point
                C_E = self%horzCoeffs(1, E_indx + 1 - self%startCol)
                C_W = self%horzCoeffs(2, E_indx + 1 - self%startCol)
                a_E = C_W/(C_E + C_W)

                ! interpolation from west point
                C_E = self%horzCoeffs(1, W_indx + 1 - self%startCol)
                C_W = self%horzCoeffs(2, W_indx + 1 - self%startCol)
                a_W = C_E/(C_E + C_W)
            else if (p == self%startColCoarse) then
                ! left neumann
                C_E = self%horzCoeffs(1, E_indx + 1 - self%startCol)
                C_W = self%horzCoeffs(2, E_indx + 1 - self%startCol)
                a_E = C_W/(C_E + C_W)
                a_W = a_E
            else
                ! right neumann
                C_E = self%horzCoeffs(1, W_indx + 1 - self%startCol)
                C_W = self%horzCoeffs(2, W_indx + 1 - self%startCol)
                a_W = C_E/(C_E + C_W)
                a_E = a_W
            end if

            if (N_indx /= S_indx) then
                ! interpolation from north point
                C_N = self%vertCoeffs(1, N_indx + 1 - self%startRow)
                C_S = self%vertCoeffs(2, N_indx + 1 - self%startRow)
                a_N = C_S/(C_N + C_S)

                ! interpolation from west point
                C_N = self%vertCoeffs(1, S_indx + 1 - self%startRow)
                C_S = self%vertCoeffs(2, S_indx + 1 - self%startRow)
                a_S = C_N/(C_N + C_S)
            else if (k == self%startRowCoarse) then
                ! South neumann
                C_N = self%vertCoeffs(1, N_indx + 1 - self%startRow)
                C_S = self%vertCoeffs(2, N_indx + 1 - self%startRow)
                a_N = C_S/(C_N + C_S)
                a_S = a_N
            else
                ! top neumann
                C_N = self%vertCoeffs(1, S_indx + 1 - self%startRow)
                C_S = self%vertCoeffs(2, S_indx + 1 - self%startRow)
                a_S = C_N/(C_N + C_S)
                a_N = a_S
            end if


            !node based restriction, works best for zebra
            ! C_N = self%vertCoeffs(1,k)
            ! C_E = self%horzCoeffs(1,p)
            ! C_S = self%vertCoeffs(2,k)
            ! C_W = self%horzCoeffs(2,p)
            ! a_N = C_N/(C_N + C_S)
            ! a_S = C_S/(C_N + C_S)
            ! a_E = C_E/(C_E + C_W)
            ! a_W = C_W/(C_E + C_W)


            ! Interpolate residual to coarse grid
            coarseGrid(i_coarse, j_coarse) = 0.25d0 * (fineGrid(i_fine, j_fine) + &
            fineGrid(E_indx, j_fine) * a_E + fineGrid(W_indx, j_fine) * a_W +  &
            fineGrid(i_fine, N_indx) * a_N + fineGrid(i_fine, S_indx) * a_S + &
            fineGrid(E_indx, N_indx) * a_E * a_N + fineGrid(E_indx, S_indx) * a_E * a_S + &
            fineGrid(W_indx, N_indx) * a_W * a_N + fineGrid(W_indx, S_indx) * a_W * a_S)
        end do
    end do
    !$OMP end do
    !$OMP end parallel

end subroutine restriction_curv

subroutine prolongation_curv(self, fineGrid, coarseGrid)
    ! Prolongate operator from coarse to fine grid using gauss seidel data
    class(GS_Base_Curv), intent(in) :: self
    real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
    real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    real(real64) :: w_1, w_2, w_3, w_4, w_tot
    integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse, k, p
    N_x_coarse = (self%N_x+1)/2
    N_y_coarse = (self%N_y+1)/2

    !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine, w_1, w_2, w_3, w_4, w_tot)
    ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
    !$OMP do collapse(2)
    ! set similar nodes
    do k = self%startRowCoarse, self%endRowCoarse, 2
        do p = self%startColCoarse, self%endColCoarse, 2
            i_fine = p + self%startCol -1
            j_fine = k + self%startRow - 1
            j_coarse = (j_fine+1)/2
            i_coarse = (i_fine+1)/2
            ! Add overlapping coarse nodes
            fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
        end do
    end do
    !$OMP end do nowait
    !$OMP do collapse(2)
    ! go horizontal each row
    do k = self%startRowCoarse, self%endRowCoarse, 2
        do i_coarse = 1, N_x_coarse-1
            i_fine = 2*i_coarse ! i_fine right of overlapping coarse node
            j_fine = k + self%startRow-1
            j_coarse = (j_fine+1)/2
            p = i_fine + 1 - self%startCol ! p at i_fine
            w_1 = self%horzCoeffs(1, p) !E
            w_2 = self%horzCoeffs(2, p) !W
            w_tot = w_1 + w_2
            w_1 = w_1 / w_tot
            w_2 = w_2 / w_tot

            ! Average fine nodes between coarse nodes horizontally
            fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse+1, j_coarse) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
        end do
    end do
    !$OMP end do nowait
    !$OMP do collapse(2)
    ! go vertical
    do j_coarse = 1, N_y_coarse-1
        do p = self%startColCoarse, self%endColCoarse, 2
            i_fine = p + self%startCol-1
            i_coarse = (i_fine + 1)/2
            j_fine = j_coarse*2 ! j_fine north of overlapping coarse node
            k = j_fine + 1 - self%startRow
            w_1 = self%vertCoeffs(1, k) !N
            w_2 = self%vertCoeffs(2, k) !S
            w_tot = w_1 + w_2
            w_1 = w_1 / w_tot
            w_2 = w_2 / w_tot
            ! Average fine nodes between coarse nodes vertical direction
            fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse+1) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
        end do
    end do
    !$OMP end do nowait
    !$OMP do collapse(2)
    ! add offset fine nodes which require 4 point interpolation
    ! first loop through SW corner coarse nodes
    do j_coarse = 1, N_y_coarse-1
        do i_coarse = 1, N_x_coarse-1
            i_fine = i_coarse*2 ! i_fine east i_coarse
            j_fine = j_coarse*2 ! j_fine north j_fine
            p = i_fine + 1 - self%startCol
            k = j_fine + 1 - self%startRow
            w_1 = self%vertCoeffs(1, k) !N
            w_2 = self%vertCoeffs(2, k) !S
            w_tot = w_1 + w_2
            w_1 = w_1/w_tot
            w_2 = w_2/w_tot

            w_3 = self%horzCoeffs(1,p) !E
            w_4 = self%horzCoeffs(2,p) !W
            w_tot = w_3 + w_4
            w_3 = w_3/w_tot
            w_4 = w_4/w_tot

            fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse) * w_2*w_4 + coarseGrid(i_coarse+1, j_coarse) * w_3 * w_2 &
                + coarseGrid(i_coarse, j_coarse+1) * w_1 * w_4 + coarseGrid(i_coarse+1, j_coarse+1) * w_1 * w_3
        end do
    end do
    !$OMP end do
    !$OMP end parallel

end subroutine prolongation_curv


end module mod_GS_Base_Curv