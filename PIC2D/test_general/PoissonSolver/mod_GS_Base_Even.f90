module mod_GS_Base_Even
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_GS_Base
    use mod_domain_uniform
    implicit none


    type, extends(GS_Base) :: GS_Base_Even
        ! store grid quantities
        real(real64) :: coeffX, coeffY, centerCoeff, inv_centerCoeff
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: initialize_GS_Even
        procedure, public, pass(self) :: restriction => restriction_even
        procedure, public, pass(self) :: prolongation => prolongation_even
        procedure, public, pass(self) :: calcResidual => calcResidual_even
        procedure, public, pass(self) :: AX_Mult => AX_Mult_even
        procedure, public, pass(self) :: YAX_Mult => YAX_Mult_even
    end type

contains

    subroutine initialize_GS_Even(self, world)
        ! Construct derived object for everything needed from GS_Base_Even class
        class(GS_Base_Even), intent(in out) :: self
        type(domain_uniform), intent(in), target :: world

        self%world => world
        self%coeffX = 1.0d0 / (self%x_indx_step * world%del_x)**2
        self%coeffY = 1.0d0 / (self%y_indx_step * world%del_y)**2
        self%centerCoeff = -1.0d0 / (2.0d0 * self%coeffX + 2.0d0 * self%coeffY)
        self%inv_centerCoeff = 1.0d0 / self%centerCoeff

    end subroutine initialize_GS_Even

    subroutine constructPoissonOrthogonal(self)
        ! Construct orthogonal grid solver
        class(GS_Base_Even), intent(in out) :: self
        ! integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
        ! real(real64), intent(in) :: del_x, del_y
    end subroutine constructPoissonOrthogonal
  

    subroutine AX_Mult_even(self, x, y)
        ! Use gauss-seidel to calculate y = x^T * A * x
        class(GS_Base_Even), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y)
        real(real64), intent(in out) :: y(self%N_x, self%N_y)
        integer :: N_indx, E_indx, S_indx, W_indx, k, i, j, p
    
        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx)

        !first do corners
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            y(1,1) = 2.0d0 * x(1, 2) * self%coeffY + &
                2.0d0 * x(2, 1) * self%coeffX + x(1,1)*self%inv_centerCoeff
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            y(self%N_x, 1) = 2.0d0 * x(self%N_x, 2) * self%coeffY + &
                2.0d0 * x(self%N_x-1, 1) * self%coeffX + x(self%N_x, 1)*self%inv_centerCoeff
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            y(1, self%N_y) = 2.0d0 * x(1, self%N_y-1) * self%coeffY + &
                2.0d0 * x(2, self%N_y) * self%coeffX + x(1, self%N_y)*self%inv_centerCoeff
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            y(self%N_x, self%N_y) = 2.0d0 * x(self%N_x, self%N_y-1) * self%coeffY + &
                2.0d0 * x(self%N_x-1, self%N_y) * self%coeffX + x(self%N_x, self%N_y)*self%inv_centerCoeff
        end if
        !$OMP end sections nowait

        ! Inner Nodes
        !$OMP do
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            do p = 1, self%number_row_sections(k)
                do i = self%start_inner_indx_x(p, k), self%end_inner_indx_x(p,k)  
                    E_indx = i+1
                    W_indx = i-1
                    y(i,j) = (x(i, N_indx) + x(i, S_indx)) * self%coeffY + &
                        (x(E_indx, j) + x(W_indx, j)) * self%coeffX + x(i,j) * self%inv_centerCoeff
                end do
            end do
        end do
        !$OMP end do nowait

        !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            if (self%bottom_row_boundary_type(p) == 2) then
                S_indx = 2
            else
                S_indx = self%N_y-1
            end if
            do i = self%start_bottom_row_indx(p), self%end_bottom_row_indx(p)
                y(i,1) = (x(i, 2) + x(i, S_indx)) * self%coeffY + &
                    (x(i-1, 1) + x(i+1, 1)) * self%coeffX + x(i,1)*self%inv_centerCoeff
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
            do i = self%start_top_row_indx(p), self%end_top_row_indx(p)
                y(i,self%N_y) = (x(i, self%N_y-1) + x(i, N_indx)) * self%coeffY + &
                    (x(i-1, self%N_y) + x(i+1, self%N_y)) * self%coeffX + x(i, self%N_y)*self%inv_centerCoeff
            end do
        end do
        !$OMP end do nowait

        ! !left boundary
        !$OMP do
        do p = 1, self%number_left_column_sections
            if (self%left_column_boundary_type(p) == 2) then
                W_indx = 2
            else
                W_indx = self%N_x-1
            end if
            do j = self%start_left_column_indx(p), self%end_left_column_indx(p)
                y(1,j) = (x(1, j+1) + x(1, j-1)) * self%coeffY + &
                    (x(2, j) + x(W_indx, j)) * self%coeffX + x(1,j)*self%inv_centerCoeff
            end do
        end do
        !$OMP end do nowait

        ! !right boundary
        !$OMP do
        do p = 1, self%number_right_column_sections
            if (self%right_column_boundary_type(p) == 2) then
                E_indx = self%N_x-1
            else
                E_indx = 2
            end if
            do j = self%start_right_column_indx(p), self%end_right_column_indx(p)
                y(self%N_x,j) = (x(self%N_x, j-1) + x(self%N_x, j+1)) * self%coeffY + &
                (x(self%N_x-1, j) + x(E_indx, j)) * self%coeffX + x(self%N_x,j)*self%inv_centerCoeff
            end do
        end do
        !$OMP end do
        !$OMP end parallel
    end subroutine AX_Mult_even


    function YAX_Mult_even(self, x, y) result(res)
        ! Use gauss-seidel to calculate x^T * A * x
        class(GS_Base_Even), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y), y(self%N_x, self%N_y)
        integer :: N_indx, E_indx, S_indx, W_indx, k, i, j, p
        real(real64) :: res
        
        res = 0.0d0
        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx) reduction(+:res)

        !first do corners
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            res = res + y(1,1) * (2.0d0 * x(1, 2) * self%coeffY + &
                2.0d0 * x(2, 1) * self%coeffX + x(1,1)*self%inv_centerCoeff)
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            res = res + y(self%N_x, 1) * (2.0d0 * x(self%N_x, 2) * self%coeffY + &
                2.0d0 * x(self%N_x-1, 1) * self%coeffX + x(self%N_x, 1)*self%inv_centerCoeff)
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            res = res + y(1, self%N_y) * (2.0d0 * x(1, self%N_y-1) * self%coeffY + &
                2.0d0 * x(2, self%N_y) * self%coeffX + x(1, self%N_y)*self%inv_centerCoeff)
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            res = res + y(self%N_x, self%N_y) * (2.0d0 * x(self%N_x, self%N_y-1) * self%coeffY + &
                2.0d0 * x(self%N_x-1, self%N_y) * self%coeffX + x(self%N_x, self%N_y)*self%inv_centerCoeff)
        end if
        !$OMP end sections nowait

        ! Inner Nodes
        !$OMP do
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            do p = 1, self%number_row_sections(k)
                do i = self%start_inner_indx_x(p, k), self%end_inner_indx_x(p,k)  
                    E_indx = i+1
                    W_indx = i-1
                    res = res + y(i,j) * ((x(i, N_indx) + x(i, S_indx)) * self%coeffY + &
                        (x(E_indx, j) + x(W_indx, j)) * self%coeffX + x(i,j) * self%inv_centerCoeff)
                end do
            end do
        end do
        !$OMP end do nowait

        !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            if (self%bottom_row_boundary_type(p) == 2) then
                S_indx = 2
            else
                S_indx = self%N_y-1
            end if
            do i = self%start_bottom_row_indx(p), self%end_bottom_row_indx(p)
                res = res + y(i,1) * ( (x(i, 2) + x(i, S_indx)) * self%coeffY + &
                    (x(i-1, 1) + x(i+1, 1)) * self%coeffX + x(i,1)*self%inv_centerCoeff)
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
            do i = self%start_top_row_indx(p), self%end_top_row_indx(p)
                res = res + y(i,self%N_y) * ((x(i, self%N_y-1) + x(i, N_indx)) * self%coeffY + &
                    (x(i-1, self%N_y) + x(i+1, self%N_y)) * self%coeffX + x(i, self%N_y)*self%inv_centerCoeff)
            end do
        end do
        !$OMP end do nowait

        ! !left boundary
        !$OMP do
        do p = 1, self%number_left_column_sections
            if (self%left_column_boundary_type(p) == 2) then
                W_indx = 2
            else
                W_indx = self%N_x-1
            end if
            do j = self%start_left_column_indx(p), self%end_left_column_indx(p)
                res = res + y(1,j) * ((x(1, j+1) + x(1, j-1)) * self%coeffY + &
                    (x(2, j) + x(W_indx, j)) * self%coeffX + x(1,j)*self%inv_centerCoeff)
            end do
        end do
        !$OMP end do nowait

        ! !right boundary
        !$OMP do
        do p = 1, self%number_right_column_sections
            if (self%right_column_boundary_type(p) == 2) then
                E_indx = self%N_x-1
            else
                E_indx = 2
            end if
            do j = self%start_right_column_indx(p), self%end_right_column_indx(p)
                res = res + y(self%N_x,j) * ((x(self%N_x, j-1) + x(self%N_x, j+1)) * self%coeffY + &
                (x(self%N_x-1, j) + x(E_indx, j)) * self%coeffX + x(self%N_x,j)*self%inv_centerCoeff)
            end do
        end do
        !$OMP end do
        !$OMP end parallel
    end function YAX_Mult_even

    subroutine calcResidual_even(self)
        ! Solve GS down to some tolerance
        class(GS_Base_Even), intent(in out) :: self
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, i_finest, j_finest


        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx)

        !first do corners
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            self%residual(1,1) = self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * self%coeffY - &
                2.0d0 * self%solution(2, 1) * self%coeffX  - self%solution(1,1)*self%inv_centerCoeff
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            self%residual(self%N_x, 1) = self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, 1) * self%coeffX -  self%solution(self%N_x, 1)*self%inv_centerCoeff
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            self%residual(1, self%N_y) = self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(2, self%N_y) * self%coeffX - self%solution(1, self%N_y)*self%inv_centerCoeff
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            self%residual(self%N_x, self%N_y) = self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, self%N_y) * self%coeffX - self%solution(self%N_x, self%N_y)*self%inv_centerCoeff
        end if
        !$OMP end sections nowait

        ! Inner Nodes
        !$OMP do
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            do p = 1, self%number_row_sections(k)
                do i = self%start_inner_indx_x(p, k), self%end_inner_indx_x(p,k)  
                    E_indx = i+1
                    W_indx = i-1
                    self%residual(i,j) = self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX - self%solution(i,j) / self%centerCoeff
                end do
            end do
        end do
        !$OMP end do nowait

        !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            if (self%bottom_row_boundary_type(p) == 2) then
                S_indx = 2
            else
                S_indx = self%N_y-1
            end if
            do i = self%start_bottom_row_indx(p), self%end_bottom_row_indx(p)
                self%residual(i,1) = self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX -  self%solution(i,1)*self%inv_centerCoeff
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
            do i = self%start_top_row_indx(p), self%end_top_row_indx(p)
                self%residual(i,self%N_y) = self%sourceTerm(i,self%N_y) - (self%solution(i, self%N_y-1) + self%solution(i, N_indx)) * self%coeffY - &
                    (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX -  self%solution(i, self%N_y)*self%inv_centerCoeff
            end do
        end do
        !$OMP end do nowait

        ! !left boundary
        !$OMP do
        do p = 1, self%number_left_column_sections
            if (self%left_column_boundary_type(p) == 2) then
                W_indx = 2
            else
                W_indx = self%N_x-1
            end if
            do j = self%start_left_column_indx(p), self%end_left_column_indx(p)
                self%residual(1,j) = self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                    (self%solution(2, j) + self%solution(W_indx, j)) * self%coeffX - self%solution(1,j)*self%inv_centerCoeff
            end do
        end do
        !$OMP end do nowait

        ! !right boundary
        !$OMP do
        do p = 1, self%number_right_column_sections
            if (self%right_column_boundary_type(p) == 2) then
                E_indx = self%N_x-1
            else
                E_indx = 2
            end if
            do j = self%start_right_column_indx(p), self%end_right_column_indx(p)
                self%residual(self%N_x,j) = self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j-1) + self%solution(self%N_x, j+1)) * self%coeffY - &
                (self%solution(self%N_x-1, j) + self%solution(E_indx, j)) * self%coeffX - self%solution(self%N_x,j)*self%inv_centerCoeff
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
        integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k, p, start_indx


        !$OMP parallel private(k, p, i_fine, j_fine, N_indx, W_indx, E_indx, S_indx, j_coarse, i_coarse, start_indx)

        ! Inner Nodes
        !$OMP do
        do k = self%start_coarse_row_indx, self%number_inner_rows, 2
            j_fine = self%start_row_indx + k - 1
            j_coarse = (j_fine + 1)/2
            N_indx = j_fine + 1
            S_indx = j_fine - 1
            do p = 1, self%number_row_sections(k)
                start_indx = self%start_inner_indx_x(p, k)
                if (MOD(start_indx,2) == 0) start_indx = start_indx + 1
                do i_fine = start_indx, self%end_inner_indx_x(p,k), 2  
                    i_coarse = (i_fine+1)/2
                    E_indx = i_fine+1
                    W_indx = i_fine-1
                    ! Bilinear interpolation
                    coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
                    0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
                    fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
                    0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
                    fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
                end do
            end do
        end do
        !$OMP end do nowait

        !first do corners
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            coarseGrid(1,1) = 0.25d0 * fineGrid(1,1) + &
            0.25d0 * (fineGrid(2, 1) +  fineGrid(1, 2)) + &
            0.25d0 * (fineGrid(2,2))
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            coarseGrid((self%N_x+1)/2,1) = 0.25d0 * fineGrid(self%N_x,1) + &
            0.25d0 * (fineGrid(self%N_x-1, 1) +  fineGrid(self%N_x, 2)) + &
            0.25d0 * (fineGrid(self%N_x-1,2))
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            coarseGrid(1,(self%N_y+1)/2) = 0.25d0 * fineGrid(1,self%N_y) + &
            0.25d0 * (fineGrid(2, self%N_y) +  fineGrid(1, self%N_y-1)) + &
            0.25d0 * (fineGrid(2, self%N_y-1))
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            coarseGrid((self%N_x+1)/2,(self%N_y+1)/2) = 0.25d0 * fineGrid(self%N_x,self%N_y) + &
            0.25d0 * (fineGrid(self%N_x-1, self%N_y) +  fineGrid(self%N_x, self%N_y-1)) + &
            0.25d0 * (fineGrid(self%N_x-1, self%N_y-1))
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
            N_indx = 2
            j_fine = 1
            j_coarse = 1
            start_indx = self%start_bottom_row_indx(p)
            if (MOD(start_indx,2) == 0) start_indx = start_indx + 1
            do i_fine = start_indx, self%end_bottom_row_indx(p), 2
                i_coarse = (i_fine+1)/2
                E_indx = i_fine+1
                W_indx = i_fine-1
                coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
                    0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
                    fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
                    0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
                    fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
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
            S_indx = self%N_y-1
            j_fine = self%N_y
            j_coarse = (self%N_y+1)/2
            start_indx = self%start_top_row_indx(p)
            if (MOD(start_indx,2) == 0) start_indx = start_indx + 1
            do i_fine = start_indx, self%end_top_row_indx(p), 2
                i_coarse = (i_fine+1)/2
                E_indx = i_fine+1
                W_indx = i_fine-1
                coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
                    0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
                    fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
                    0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
                    fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
            end do
        end do
        !$OMP end do nowait

        ! !left boundary
        !$OMP do
        do p = 1, self%number_left_column_sections
            if (self%left_column_boundary_type(p) == 2) then
                W_indx = 2
            else
                W_indx = self%N_x-1
            end if
            E_indx = 2
            i_fine = 1
            i_coarse = 1
            start_indx = self%start_left_column_indx(p)
            if (MOD(start_indx,2) == 0) start_indx = start_indx + 1
            do j_fine = start_indx, self%end_left_column_indx(p), 2
                j_coarse = (j_fine+1)/2
                N_indx = j_fine + 1
                S_indx = j_fine - 1
                coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
                    0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
                    fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
                    0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
                    fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
            end do
        end do
        !$OMP end do nowait

        ! !right boundary
        !$OMP do
        do p = 1, self%number_right_column_sections
            if (self%right_column_boundary_type(p) == 2) then
                E_indx = self%N_x-1
            else
                E_indx = 2
            end if
            W_indx = self%N_x-1
            i_fine = self%N_x
            i_coarse = (self%N_x+1)/2
            start_indx = self%start_right_column_indx(p)
            if (MOD(start_indx,2) == 0) start_indx = start_indx + 1
            do j_fine = start_indx, self%end_right_column_indx(p), 2
                j_coarse = (j_fine+1)/2
                N_indx = j_fine + 1
                S_indx = j_fine - 1
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
        integer :: i_fine, j_fine, i_coarse, j_coarse, k, p, start_indx, N_y_coarse, N_x_coarse
        N_y_coarse = (self%N_y+1)/2
        N_x_coarse = (self%N_x+1)/2

        ! Need to avoid interpolating onto dirichlet fine nodes, so can't generalize 
        !$OMP parallel private(i_fine, j_fine, i_coarse, j_coarse, k, p, start_indx)
        ! do each do loop individually because can't think of fewer loops which don't do unnecessary calculations

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

        !$OMP do
        ! inner fine nodes horizontally between coarse nodes
        do k = self%start_coarse_row_indx, self%number_inner_rows, 2
            j_fine = self%start_row_indx + k - 1
            j_coarse = (j_fine + 1)/2
            do p = 1, self%number_row_sections(k)
                start_indx = self%start_inner_indx_x(p, k)
                if (MOD(start_indx,2) == 1) start_indx = start_indx + 1
                do i_fine = start_indx, self%end_inner_indx_x(p,k), 2  
                    i_coarse = (i_fine)/2 ! i_coarse to the left of fine node
                    fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
                end do
            end do
        end do
        !$OMP end do nowait

        !$OMP do
        ! inner fine nodes vertically between coarse nodes
        do k = self%start_fine_row_indx, self%number_inner_rows, 2
            j_fine = self%start_row_indx + k - 1
            j_coarse = (j_fine)/2 ! this is coarse lower than j_fine
            do p = 1, self%number_row_sections(k)
                start_indx = self%start_inner_indx_x(p, k)
                if (MOD(start_indx,2) == 0) start_indx = start_indx + 1
                do i_fine = start_indx, self%end_inner_indx_x(p,k), 2   ! i_fine odd
                    i_coarse = (i_fine + 1)/2
                    fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse+1) * 0.5d0
                end do
            end do
        end do
        !$OMP end do nowait

        !$OMP do
        ! inner fine nodes with bilinear interpolation
        do k = self%start_fine_row_indx, self%number_inner_rows, 2
            j_fine = self%start_row_indx + k - 1
            j_coarse = (j_fine)/2 ! this is coarse lower than j_fine
            do p = 1, self%number_row_sections(k)
                start_indx = self%start_inner_indx_x(p, k)
                if (MOD(start_indx,2) == 1) start_indx = start_indx + 1
                do i_fine = start_indx, self%end_inner_indx_x(p,k), 2  ! i_fine even
                    i_coarse = (i_fine)/2 ! coarse node lower i_fine
                    fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse) * 0.25d0 &
                    + coarseGrid(i_coarse, j_coarse+1) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse+1) * 0.25d0
                end do
            end do
        end do
        !$OMP end do nowait

        !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            j_fine = 1
            j_coarse = 1
            start_indx = self%start_bottom_row_indx(p)
            if (MOD(start_indx,2) == 1) start_indx = start_indx + 1
            do i_fine = start_indx, self%end_bottom_row_indx(p), 2
                i_coarse = (i_fine)/2 ! i_coarse to the left of fine node
                fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
            end do
        end do
        !$OMP end do nowait

        ! ! upper boundary
        !$OMP do
        do p = 1, self%number_top_row_sections
            j_fine = self%N_y
            j_coarse = (self%N_y+1)/2
            start_indx = self%start_top_row_indx(p)
            if (MOD(start_indx,2) == 1) start_indx = start_indx + 1
            do i_fine = start_indx, self%end_top_row_indx(p), 2
                i_coarse = (i_fine)/2 ! i_coarse to the left of fine node
                fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
            end do
        end do
        !$OMP end do nowait

        ! !left boundary
        !$OMP do
        do p = 1, self%number_left_column_sections
            i_fine = 1
            i_coarse = 1
            start_indx = self%start_left_column_indx(p)
            if (MOD(start_indx,2) == 1) start_indx = start_indx + 1
            do j_fine = start_indx, self%end_left_column_indx(p), 2
                j_coarse = (j_fine)/2
                fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse+1) * 0.5d0
            end do
        end do
        !$OMP end do nowait

        ! !right boundary
        !$OMP do
        do p = 1, self%number_right_column_sections
            i_fine = self%N_x
            i_coarse = (self%N_x+1)/2
            start_indx = self%start_right_column_indx(p)
            if (MOD(start_indx,2) == 1) start_indx = start_indx + 1
            do j_fine = start_indx, self%end_right_column_indx(p), 2
                j_coarse = (j_fine)/2
                fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse+1) * 0.5d0
            end do
        end do
        !$OMP end do
        !$OMP end parallel
    end subroutine prolongation_even


end module mod_GS_Base_Even