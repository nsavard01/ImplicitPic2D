module mod_GS_Base_Curv
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_GS_Base
    use mod_domain_curv
    implicit none


    type, extends(GS_Base) :: GS_Base_Curv
        ! store grid quantities
        real(real64), allocatable :: centerCoeff(:,:) ! just store all central coefficients on grid, assume most will be used, out of simplicity
        real(real64), allocatable :: inner_node_coeff_north(:), inner_node_coeff_south(:), &
            inner_node_coeff_east(:), inner_node_coeff_west(:) 
        real(real64), allocatable :: bottom_row_coeff_south(:), bottom_row_coeff_north(:)!, bottom_row_coeff_center(:)
        real(real64), allocatable :: top_row_coeff_south(:), top_row_coeff_north(:)!, top_row_coeff_center(:)
        real(real64), allocatable :: left_column_coeff_east(:), left_column_coeff_west(:)!, left_column_coeff_center(:)
        real(real64), allocatable :: right_column_coeff_east(:), right_column_coeff_west(:)!, right_column_coeff_center(:)
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: initialize_GS_Curv
        ! procedure, public, pass(self) :: restriction => restriction_curv
        ! procedure, public, pass(self) :: prolongation => prolongation_curv
        procedure, public, pass(self) :: calcResidual => calcResidual_curv
        ! procedure, public, pass(self) :: YAX_Mult => YAX_Mult_curv
        ! procedure, public, pass(self) :: AX_Mult => AX_Mult_curv
    end type

contains

subroutine constructPoissonOrthogonal(self)
    ! Construct orthogonal grid solver
    class(GS_Base_Curv), intent(in out) :: self
    ! integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
    ! real(real64), intent(in) :: del_x, del_y
end subroutine constructPoissonOrthogonal

subroutine initialize_GS_Curv(self, world)
    ! Construct derived object for everything needed from GS_Base_Even class
    class(GS_Base_Curv), intent(in out) :: self
    type(domain_curv), intent(in), target :: world
    real(real64) :: del_x_west, del_x_east, del_y_north, del_y_south
    integer :: i, j, p, i_fine, j_fine

    self%world => world

    ! Allocate inner node coefficients
    allocate(self%centerCoeff(self%N_x, self%N_y), self%inner_node_coeff_north(self%N_y-2), self%inner_node_coeff_south(self%N_y-2), &
    self%inner_node_coeff_east(self%N_x-2), self%inner_node_coeff_west(self%N_x-2))
    do i = 1, self%N_x-2
        i_fine = i * self%x_indx_step - self%x_indx_step + 1
        del_x_west = SUM(world%del_x(i_fine:i_fine + self%x_indx_step-1))
        del_x_east = SUM(world%del_x(i_fine+self%x_indx_step:i_fine + 2*self%x_indx_step-1))
        self%inner_node_coeff_west(i) = 1.0d0 / (del_x_west * 0.5d0 * (del_x_west + del_x_east))
        self%inner_node_coeff_east(i) = 1.0d0 / (del_x_east * 0.5d0 * (del_x_west + del_x_east))
    end do
    do j = 1, self%N_y-2
        j_fine = j * self%y_indx_step - self%y_indx_step + 1
        del_y_south = SUM(world%del_y(j_fine:j_fine+self%y_indx_step-1))
        del_y_north = SUM(world%del_y(j_fine + self%y_indx_step:j_fine + 2 * self%y_indx_step -1))
        self%inner_node_coeff_north(j) = 1.0d0 / (del_y_north * 0.5d0 * (del_y_north + del_y_south))
        self%inner_node_coeff_south(j) = 1.0d0 / (del_y_south * 0.5d0 * (del_y_north + del_y_south))
    end do
    

    !$OMP parallel private(i,j)
    !$OMP do collapse(2)
    do j = 1, self%N_y-2
        do i = 1, self%N_x-2
            self%centerCoeff(i+1,j+1) = -1.0d0 / (self%inner_node_coeff_north(j) + self%inner_node_coeff_south(j) + &
            self%inner_node_coeff_east(i) + self%inner_node_coeff_west(i))
        end do
    end do
    !$OMP end do
    !$OMP end parallel

    del_x_west = SUM(world%del_x(1:self%x_indx_step))
    del_x_east = SUM(world%del_x(self%N_x-self%x_indx_step:self%N_x-1))
    del_y_south = SUM(world%del_y(1:self%y_indx_step))
    del_y_north = SUM(world%del_y(self%N_y-self%y_indx_step:self%N_y-1))
    self%centerCoeff(1,1) = -0.5d0 / (1.0d0 / (del_x_west**2) + 1.0d0 / (del_y_south**2))
    self%centerCoeff(self%N_x, 1) = -0.5d0 / (1.0d0 / (del_x_east**2) + 1.0d0 / (del_y_south**2))
    self%centerCoeff(1, self%N_y) = -0.5d0 / (1.0d0 / (del_x_west**2) + 1.0d0 / (del_y_north**2))
    self%centerCoeff(self%N_x, self%N_y) = -0.5d0 / (1.0d0 / (del_x_east**2) + 1.0d0 / (del_y_north**2))

    ! allocate boundary values
    ! lower boundary
    allocate(self%bottom_row_coeff_north(self%number_bottom_row_sections), self%bottom_row_coeff_south(self%number_bottom_row_sections))
    do i = 1, self%number_bottom_row_sections
        if (self%bottom_row_boundary_type(i) == 2) then
            !Neumann
            self%bottom_row_coeff_north(i) = 1.0d0 / (del_y_south**2)
            self%bottom_row_coeff_south(i) = 1.0d0 / (del_y_south**2)
        else
            !Periodic
            self%bottom_row_coeff_north(i) = 1.0d0 / (del_y_south * 0.5d0 * (del_y_south + del_y_north))
            self%bottom_row_coeff_south(i) = 1.0d0 / (del_y_north * 0.5d0 * (del_y_south + del_y_north))
        end if
        do p = self%start_bottom_row_indx(i), self%end_bottom_row_indx(i)
            self%centerCoeff(p, 1) = -1.0d0 / (self%inner_node_coeff_east(p-1) + &
            self%inner_node_coeff_west(p-1) + self%bottom_row_coeff_south(i) + self%bottom_row_coeff_north(i))
        end do
    end do

    ! upper boundary
    allocate(self%top_row_coeff_north(self%number_top_row_sections), self%top_row_coeff_south(self%number_top_row_sections))
    do i = 1, self%number_top_row_sections
        if (self%top_row_boundary_type(i) == 2) then
            !Neumann
            self%top_row_coeff_north(i) = 1.0d0 / (del_y_north**2)
            self%top_row_coeff_south(i) = 1.0d0 / (del_y_north**2)
        else
            !Periodic
            self%top_row_coeff_north(i) = 1.0d0 / (del_y_south * 0.5d0 * (del_y_south + del_y_north))
            self%top_row_coeff_south(i) = 1.0d0 / (del_y_north * 0.5d0 * (del_y_south + del_y_north))
        end if
        do p = self%start_top_row_indx(i), self%end_top_row_indx(i)
            self%centerCoeff(p, self%N_y) = -1.0d0 / (self%inner_node_coeff_east(p-1) + &
            self%inner_node_coeff_west(p-1) + self%top_row_coeff_south(i) + self%top_row_coeff_north(i))
        end do
    end do

    ! left boundary
    allocate(self%left_column_coeff_east(self%number_left_column_sections), self%left_column_coeff_west(self%number_left_column_sections))
    do i = 1, self%number_left_column_sections
        if (self%left_column_boundary_type(i) == 2) then
            !Neumann
            self%left_column_coeff_east(i) = 1.0d0 / (del_x_west**2)
            self%left_column_coeff_west(i) = 1.0d0 / (del_x_west**2)
        else
            !Periodic
            self%left_column_coeff_east(i) = 1.0d0 / (del_x_west * 0.5d0 * (del_x_west + del_x_east))
            self%left_column_coeff_west(i) = 1.0d0 / (del_x_east * 0.5d0 * (del_x_west + del_x_east))
        end if
        do p = self%start_left_column_indx(i), self%end_left_column_indx(i)
            self%centerCoeff(1, p) = -1.0d0 / (self%inner_node_coeff_south(p-1) + &
            self%inner_node_coeff_north(p-1) + self%left_column_coeff_east(i) + self%left_column_coeff_west(i))
        end do
    end do

    ! right boundary
    allocate(self%right_column_coeff_east(self%number_right_column_sections), self%right_column_coeff_west(self%number_right_column_sections))
    do i = 1, self%number_right_column_sections
        if (self%right_column_boundary_type(i) == 2) then
            !Neumann
            self%right_column_coeff_east(i) = 1.0d0 / (del_x_east**2)
            self%right_column_coeff_west(i) = 1.0d0 / (del_x_east**2)
        else
            !Periodic
            self%right_column_coeff_east(i) = 1.0d0 / (del_x_west * 0.5d0 * (del_x_west + del_x_east))
            self%right_column_coeff_west(i) = 1.0d0 / (del_x_east * 0.5d0 * (del_x_west + del_x_east))
        end if
        do p = self%start_right_column_indx(i), self%end_right_column_indx(i)
            self%centerCoeff(self%N_x, p) = -1.0d0 / (self%inner_node_coeff_south(p-1) + &
            self%inner_node_coeff_north(p-1) + self%right_column_coeff_east(i) + self%right_column_coeff_west(i))
        end do
    end do


end subroutine initialize_GS_Curv

! subroutine AX_Mult_curv(self, x, y)
!     ! Solve GS down to some tolerance
!     class(GS_Base_Curv), intent(in out) :: self
!     real(real64), intent(in) :: x(self%N_x, self%N_y)
!     real(real64), intent(in out) :: y(self%N_x, self%N_y)
!     real(real64) :: C_N, C_E, C_O, C_W, C_S
!     integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
!     !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
!     !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W)
!     ! loop through inner nodes
!     !$OMP do collapse(2)
!     do k = 1, self%numberRows
!         do p = 1, self%numberColumns
!             i = self%startCol + p - 1
!             j = self%startRow + k - 1
!             N_indx = self%vertIndx(1,k)
!             E_indx = self%horzIndx(1,p)
!             S_indx = self%vertIndx(2,k)
!             W_indx = self%horzIndx(2,p)
!             C_O = self%centerCoeffs(p,k)
!             C_N = self%vertCoeffs(1,k)
!             C_E = self%horzCoeffs(1,p)
!             C_S = self%vertCoeffs(2,k)
!             C_W = self%horzCoeffs(2,p)
!             y(i,j) = C_N * x(i,N_indx) + C_S * x(i,S_indx) + &
!                     C_E * x(E_indx, j) + C_W * x(W_indx, j) + x(i,j)/C_O
!         end do
!     end do
!     !$OMP end do
!     !$OMP end parallel
! end subroutine AX_Mult_curv

! function YAX_Mult_curv(self, x, y) result(res)
!     ! Solve GS down to some tolerance
!     class(GS_Base_Curv), intent(in out) :: self
!     real(real64), intent(in) :: x(self%N_x, self%N_y), y(self%N_x, self%N_y)
!     real(real64) :: res, C_N, C_E, C_O, C_W, C_S
!     integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
!     res = 0.0d0
!     !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
!     !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W) reduction(+:res)
!     ! loop through inner nodes
!     !$OMP do collapse(2)
!     do k = 1, self%numberRows
!         do p = 1, self%numberColumns
!             i = self%startCol + p - 1
!             j = self%startRow + k - 1
!             N_indx = self%vertIndx(1,k)
!             E_indx = self%horzIndx(1,p)
!             S_indx = self%vertIndx(2,k)
!             W_indx = self%horzIndx(2,p)
!             C_O = self%centerCoeffs(p,k)
!             C_N = self%vertCoeffs(1,k)
!             C_E = self%horzCoeffs(1,p)
!             C_S = self%vertCoeffs(2,k)
!             C_W = self%horzCoeffs(2,p)
!             res = res + y(i,j) * (C_N * x(i,N_indx) + C_S * x(i,S_indx) + &
!                     C_E * x(E_indx, j) + C_W * x(W_indx, j) + x(i,j)/C_O)
!         end do
!     end do
!     !$OMP end do
!     !$OMP end parallel
! end function YAX_Mult_curv

    subroutine calcResidual_curv(self)
        ! Solve GS down to some tolerance
        class(GS_Base_Curv), intent(in out) :: self
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
        real(real64) :: C_N, C_E, C_S, C_W, C_O


        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, &
        !$OMP& C_N, C_E, C_S, C_W, C_O)

        ! Inner Nodes
        !$OMP do
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            C_N = self%inner_node_coeff_north(j-1)
            C_S = self%inner_node_coeff_south(j-1)
            do p = 1, self%number_row_sections(k)
                do i = self%start_inner_indx_x(p, k), self%end_inner_indx_x(p,k)  
                    E_indx = i+1
                    W_indx = i-1
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = 1.0d0/self%centerCoeff(i, j)
                    self%residual(i,j) = self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx)*C_S - &
                        self%solution(E_indx, j)*C_E - self%solution(W_indx, j) * C_W - self%solution(i,j) * C_O
                end do
            end do
        end do
        !$OMP end do nowait

        !first do corners
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            C_N = self%bottom_row_coeff_north(1)
            C_E = self%left_column_coeff_east(1)
            C_O = 1.0d0/self%centerCoeff(1,1)
            self%residual(1,1) = self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * C_N - &
                2.0d0 * self%solution(2, 1) * C_E  - self%solution(1,1)*C_O
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            C_N = self%bottom_row_coeff_north(self%number_bottom_row_sections)
            C_W = self%right_column_coeff_west(1)
            C_O = 1.0d0/self%centerCoeff(self%N_x,1)
            self%residual(self%N_x, 1) = self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * C_N - &
                2.0d0 * self%solution(self%N_x-1, 1) * C_W -  self%solution(self%N_x, 1)*C_O
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            C_S = self%top_row_coeff_south(1)
            C_E = self%left_column_coeff_east(self%number_left_column_sections)
            C_O = 1.0d0/self%centerCoeff(1,self%N_y)
            self%residual(1, self%N_y) = self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * C_S - &
                2.0d0 * self%solution(2, self%N_y) * C_E - self%solution(1, self%N_y)*C_O
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            C_S = self%top_row_coeff_south(self%number_top_row_sections)
            C_W = self%right_column_coeff_west(self%number_right_column_sections)
            C_O = 1.0d0/self%centerCoeff(self%N_x, self%N_y)
            self%residual(self%N_x, self%N_y) = self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * C_S - &
                2.0d0 * self%solution(self%N_x-1, self%N_y) * C_W - self%solution(self%N_x, self%N_y)*C_O
        end if
        !$OMP end sections nowait

        !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            if (self%bottom_row_boundary_type(p) == 2) then
                S_indx = 2
            else
                S_indx = self%N_y-1
            end if
            C_S = self%bottom_row_coeff_south(p)
            C_N = self%bottom_row_coeff_north(p)
            do i = self%start_bottom_row_indx(p), self%end_bottom_row_indx(p)
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                C_O = 1.0d0/self%centerCoeff(i,1)
                self%residual(i,1) = self%sourceTerm(i,1) - self%solution(i, 2) * C_N - self%solution(i, S_indx) * C_S - &
                    self%solution(i-1, 1) * C_W - self%solution(i+1, 1) * C_E -  self%solution(i,1)*C_O
            end do
        end do
        !$OMP end do nowait

        ! ! upper boundary
        !$OMP do
        do p = 1, self%number_top_row_sections
            if (self%top_row_boundary_type(p) == 2) then
                N_indx = self%N_y-1
            else
                N_indx = 2
            end if
            C_S = self%top_row_coeff_south(p)
            C_N = self%top_row_coeff_north(p)
            do i = self%start_top_row_indx(p), self%end_top_row_indx(p)
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                C_O = 1.0d0/self%centerCoeff(i, self%N_y)
                self%residual(i,self%N_y) = self%sourceTerm(i,self%N_y) - self%solution(i, N_indx) * C_N - self%solution(i, self%N_y-1) * C_S - &
                    self%solution(i-1, self%N_y) * C_W - self%solution(i+1, self%N_y) * C_E -  self%solution(i, self%N_y)*C_O
            end do
        end do
        !$OMP end do nowait

        ! left boundary
        !$OMP do
        do p = 1, self%number_left_column_sections
            if (self%left_column_boundary_type(p) == 2) then
                W_indx = 2
            else
                W_indx = self%N_x-1
            end if
            C_E = self%left_column_coeff_east(p)
            C_W = self%left_column_coeff_west(p)
            do j = self%start_left_column_indx(p), self%end_left_column_indx(p)
                C_S = self%inner_node_coeff_south(j-1)
                C_N = self%inner_node_coeff_north(j-1)
                C_O = 1.0d0 / self%centerCoeff(1,j)
                self%residual(1,j) = self%sourceTerm(1,j) - self%solution(1, j+1) * C_N - self%solution(1, j-1) * C_S - &
                    self%solution(2, j) * C_E - self%solution(W_indx, j) * C_W - self%solution(1,j)*C_O
            end do
        end do
        !$OMP end do nowait

        ! right boundary
        !$OMP do
        do p = 1, self%number_right_column_sections
            if (self%right_column_boundary_type(p) == 2) then
                E_indx = self%N_x-1
            else
                E_indx = 2
            end if
            C_E = self%right_column_coeff_east(p)
            C_W = self%right_column_coeff_west(p)
            do j = self%start_right_column_indx(p), self%end_right_column_indx(p)
                C_S = self%inner_node_coeff_south(j-1)
                C_N = self%inner_node_coeff_north(j-1)
                C_O = 1.0d0/ self%centerCoeff(self%N_x,j)
                self%residual(self%N_x,j) = self%sourceTerm(self%N_x,j) - self%solution(self%N_x, j+1) * C_N - self%solution(self%N_x, j-1) * C_S - &
                    self%solution(self%N_x-1, j) * C_W - self%solution(E_indx, j) * C_E - self%solution(self%N_x,j)*C_O
            end do
        end do
        !$OMP end do
        !$OMP end parallel
    end subroutine calcResidual_curv

! subroutine restriction_curv(self, fineGrid, coarseGrid)
!     ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
!     class(GS_Base_Curv), intent(in) :: self
!     real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
!     real(real64), intent(in out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
!     real(real64) :: a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S
!     integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k, p

!     !$OMP parallel private(k, p, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, &
!     !$OMP& a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S)
!     !$OMP do collapse(2)
!     ! loop over fine grid nodes which overlap with coarse grid nodes
!     do k = self%startRowCoarse, self%endRowCoarse, 2
!         do p = self%startColCoarse, self%endColCoarse, 2
!             ! calculate fine indices around overlapping index
!             i_fine = self%startCol + p - 1
!             j_fine = self%startRow + k - 1
!             N_indx = self%vertIndx(1,k)
!             E_indx = self%horzIndx(1,p)
!             S_indx = self%vertIndx(2,k)
!             W_indx = self%horzIndx(2,p)
!             i_coarse = (i_fine + 1)/2
!             j_coarse = (j_fine + 1)/2

!             ! Transpose of prolongation, seems to work best for red black
!             if (E_indx /= W_indx) then
!                 ! interpolation from east point
!                 C_E = self%horzCoeffs(1, E_indx + 1 - self%startCol)
!                 C_W = self%horzCoeffs(2, E_indx + 1 - self%startCol)
!                 a_E = C_W/(C_E + C_W)

!                 ! interpolation from west point
!                 C_E = self%horzCoeffs(1, W_indx + 1 - self%startCol)
!                 C_W = self%horzCoeffs(2, W_indx + 1 - self%startCol)
!                 a_W = C_E/(C_E + C_W)
!             else if (p == self%startColCoarse) then
!                 ! left neumann
!                 C_E = self%horzCoeffs(1, E_indx + 1 - self%startCol)
!                 C_W = self%horzCoeffs(2, E_indx + 1 - self%startCol)
!                 a_E = C_W/(C_E + C_W)
!                 a_W = a_E
!             else
!                 ! right neumann
!                 C_E = self%horzCoeffs(1, W_indx + 1 - self%startCol)
!                 C_W = self%horzCoeffs(2, W_indx + 1 - self%startCol)
!                 a_W = C_E/(C_E + C_W)
!                 a_E = a_W
!             end if

!             if (N_indx /= S_indx) then
!                 ! interpolation from north point
!                 C_N = self%vertCoeffs(1, N_indx + 1 - self%startRow)
!                 C_S = self%vertCoeffs(2, N_indx + 1 - self%startRow)
!                 a_N = C_S/(C_N + C_S)

!                 ! interpolation from west point
!                 C_N = self%vertCoeffs(1, S_indx + 1 - self%startRow)
!                 C_S = self%vertCoeffs(2, S_indx + 1 - self%startRow)
!                 a_S = C_N/(C_N + C_S)
!             else if (k == self%startRowCoarse) then
!                 ! South neumann
!                 C_N = self%vertCoeffs(1, N_indx + 1 - self%startRow)
!                 C_S = self%vertCoeffs(2, N_indx + 1 - self%startRow)
!                 a_N = C_S/(C_N + C_S)
!                 a_S = a_N
!             else
!                 ! top neumann
!                 C_N = self%vertCoeffs(1, S_indx + 1 - self%startRow)
!                 C_S = self%vertCoeffs(2, S_indx + 1 - self%startRow)
!                 a_S = C_N/(C_N + C_S)
!                 a_N = a_S
!             end if


!             !node based restriction, works best for zebra
!             ! C_N = self%vertCoeffs(1,k)
!             ! C_E = self%horzCoeffs(1,p)
!             ! C_S = self%vertCoeffs(2,k)
!             ! C_W = self%horzCoeffs(2,p)
!             ! a_N = C_N/(C_N + C_S)
!             ! a_S = C_S/(C_N + C_S)
!             ! a_E = C_E/(C_E + C_W)
!             ! a_W = C_W/(C_E + C_W)


!             ! Interpolate residual to coarse grid
!             coarseGrid(i_coarse, j_coarse) = 0.25d0 * (fineGrid(i_fine, j_fine) + &
!             fineGrid(E_indx, j_fine) * a_E + fineGrid(W_indx, j_fine) * a_W +  &
!             fineGrid(i_fine, N_indx) * a_N + fineGrid(i_fine, S_indx) * a_S + &
!             fineGrid(E_indx, N_indx) * a_E * a_N + fineGrid(E_indx, S_indx) * a_E * a_S + &
!             fineGrid(W_indx, N_indx) * a_W * a_N + fineGrid(W_indx, S_indx) * a_W * a_S)
!         end do
!     end do
!     !$OMP end do
!     !$OMP end parallel

! end subroutine restriction_curv

! subroutine prolongation_curv(self, fineGrid, coarseGrid)
!     ! Prolongate operator from coarse to fine grid using gauss seidel data
!     class(GS_Base_Curv), intent(in) :: self
!     real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
!     real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
!     real(real64) :: w_1, w_2, w_3, w_4, w_tot
!     integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse, k, p
!     N_x_coarse = (self%N_x+1)/2
!     N_y_coarse = (self%N_y+1)/2

!     !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine, w_1, w_2, w_3, w_4, w_tot)
!     ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
!     !$OMP do collapse(2)
!     ! set similar nodes
!     do k = self%startRowCoarse, self%endRowCoarse, 2
!         do p = self%startColCoarse, self%endColCoarse, 2
!             i_fine = p + self%startCol -1
!             j_fine = k + self%startRow - 1
!             j_coarse = (j_fine+1)/2
!             i_coarse = (i_fine+1)/2
!             ! Add overlapping coarse nodes
!             fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
!         end do
!     end do
!     !$OMP end do nowait
!     !$OMP do collapse(2)
!     ! go horizontal each row
!     do k = self%startRowCoarse, self%endRowCoarse, 2
!         do i_coarse = 1, N_x_coarse-1
!             i_fine = 2*i_coarse ! i_fine right of overlapping coarse node
!             j_fine = k + self%startRow-1
!             j_coarse = (j_fine+1)/2
!             p = i_fine + 1 - self%startCol ! p at i_fine
!             w_1 = self%horzCoeffs(1, p) !E
!             w_2 = self%horzCoeffs(2, p) !W
!             w_tot = w_1 + w_2
!             w_1 = w_1 / w_tot
!             w_2 = w_2 / w_tot

!             ! Average fine nodes between coarse nodes horizontally
!             fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse+1, j_coarse) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
!         end do
!     end do
!     !$OMP end do nowait
!     !$OMP do collapse(2)
!     ! go vertical
!     do j_coarse = 1, N_y_coarse-1
!         do p = self%startColCoarse, self%endColCoarse, 2
!             i_fine = p + self%startCol-1
!             i_coarse = (i_fine + 1)/2
!             j_fine = j_coarse*2 ! j_fine north of overlapping coarse node
!             k = j_fine + 1 - self%startRow
!             w_1 = self%vertCoeffs(1, k) !N
!             w_2 = self%vertCoeffs(2, k) !S
!             w_tot = w_1 + w_2
!             w_1 = w_1 / w_tot
!             w_2 = w_2 / w_tot
!             ! Average fine nodes between coarse nodes vertical direction
!             fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse+1) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
!         end do
!     end do
!     !$OMP end do nowait
!     !$OMP do collapse(2)
!     ! add offset fine nodes which require 4 point interpolation
!     ! first loop through SW corner coarse nodes
!     do j_coarse = 1, N_y_coarse-1
!         do i_coarse = 1, N_x_coarse-1
!             i_fine = i_coarse*2 ! i_fine east i_coarse
!             j_fine = j_coarse*2 ! j_fine north j_fine
!             p = i_fine + 1 - self%startCol
!             k = j_fine + 1 - self%startRow
!             w_1 = self%vertCoeffs(1, k) !N
!             w_2 = self%vertCoeffs(2, k) !S
!             w_tot = w_1 + w_2
!             w_1 = w_1/w_tot
!             w_2 = w_2/w_tot

!             w_3 = self%horzCoeffs(1,p) !E
!             w_4 = self%horzCoeffs(2,p) !W
!             w_tot = w_3 + w_4
!             w_3 = w_3/w_tot
!             w_4 = w_4/w_tot

!             fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse) * w_2*w_4 + coarseGrid(i_coarse+1, j_coarse) * w_3 * w_2 &
!                 + coarseGrid(i_coarse, j_coarse+1) * w_1 * w_4 + coarseGrid(i_coarse+1, j_coarse+1) * w_1 * w_3
!         end do
!     end do
!     !$OMP end do
!     !$OMP end parallel

! end subroutine prolongation_curv


end module mod_GS_Base_Curv