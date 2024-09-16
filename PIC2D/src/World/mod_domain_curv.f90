module mod_domain_curv
    use iso_fortran_env, only: int32, int64, real64
    use mod_domain_base
    use constants, only: pi_const
    implicit none

    ! Base world type which will contain information about the world
    ! Computational grid is rectangular, so set N_x and N_y node numbers
    ! Will extend to curvilinear grid or curv grid
    ! Assume all inner nodes nodes are boundaries
    ! For orthogonal grid so will only need grid information in x and y


    type, extends(domain_base) :: domain_curv
        ! store grid quantities
        real(real64), allocatable :: del_x(:), del_y(:)
    contains
        procedure, public, pass(self) :: get_xi_from_X => get_xi_from_X_curv
        procedure, public, pass(self) :: get_eta_from_Y => get_eta_from_Y_curv
        procedure, public, pass(self) :: get_cell_del_x => get_cell_del_x_curv
        procedure, public, pass(self) :: get_cell_del_y => get_cell_del_y_curv
    end type

    interface domain_curv
        module procedure :: constructor_domain_curv
    end interface domain_curv

contains
    type(domain_curv) function constructor_domain_curv(N_x, N_y, NESW_wall_boundaries, length_x, length_y, &
    curv_grid_type_x, curv_grid_type_y, del_x_small, del_y_small) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: N_x, N_y, NESW_wall_boundaries(4), curv_grid_type_x, curv_grid_type_y
        real(real64), intent(in) :: length_x, length_y, del_x_small, del_y_small
        integer(int32) :: left_bound, right_bound, top_bound, bottom_bound, i
        self%N_x = N_x
        self%N_y = N_y
        self%N_x_cells = N_x-1
        self%N_y_cells = N_y-1
        self%start_X = 0.0d0
        self%start_Y = 0.0d0
        self%end_X = length_x
        self%end_Y = length_y
        
        allocate(self%grid_X(self%N_x), self%grid_Y(self%N_y), &
        self%del_x(self%N_x-1), self%del_y(self%N_y-1), self%boundary_conditions(self%N_x, self%N_y))

        self%NESW_wall_boundaries = NESW_wall_boundaries ! N, E, S, W indices
        top_bound = NESW_wall_boundaries(1)
        right_bound = NESW_wall_boundaries(2)
        bottom_bound = NESW_wall_boundaries(3)
        left_bound = NESW_wall_boundaries(4)

        !$OMP parallel
        !$OMP workshare
        self%boundary_conditions = 0
        !$OMP end workshare
        !$OMP end parallel
        self%boundary_conditions(2:self%N_x-1, self%N_y) = top_bound
        self%boundary_conditions(self%N_x, 2:self%N_y-1) = right_bound
        self%boundary_conditions(2:self%N_x-1, 1) = bottom_bound
        self%boundary_conditions(1, 2:self%N_y-1) = left_bound


        ! Upper right corner
        self%boundary_conditions(self%N_x, self%N_y) = MIN(top_bound, right_bound)

        ! lower right corner
        self%boundary_conditions(self%N_x, 1) = MIN(bottom_bound, right_bound)

        ! lower left corner
        self%boundary_conditions(1, 1) = MIN(bottom_bound, left_bound)

        ! upper left corner
        self%boundary_conditions(1, self%N_y) = MIN(top_bound, left_bound)

        call create_sin_grid(self%N_x, length_x, self%grid_X, self%del_x, del_x_small)
        call create_sin_grid(self%N_y, length_y, self%grid_y, self%del_y, del_y_small)


    end function constructor_domain_curv
    
    function get_xi_from_X_curv(self, x) result(xi)
        ! get xi in computational space from x
        class(domain_curv), intent(in) :: self
        real(real64), intent(in) :: x
        real(real64) :: xi
        integer(int32) :: idx_lower
        call binary_search_left_idx(self%grid_X, self%N_x, x, idx_lower)
        xi = idx_lower + (x - self%grid_X(idx_lower))/self%del_x(idx_lower)
    end function get_xi_from_X_curv

    function get_eta_from_Y_curv(self, y) result(eta)
        ! get eta in computational space from y
        class(domain_curv), intent(in) :: self
        real(real64), intent(in) :: y
        real(real64) :: eta
        integer(int32) :: idx_lower
        call binary_search_left_idx(self%grid_Y, self%N_y, y, idx_lower)
        eta = idx_lower + (y - self%grid_Y(idx_lower))/self%del_y(idx_lower)
    end function get_eta_from_Y_curv

    function get_cell_del_x_curv(self, xi_cell) result(del_x)
        ! from xi cell value get del_x
        class(domain_curv), intent(in) :: self
        integer(int32), intent(in) :: xi_cell
        real(real64) :: del_x
        del_x = self%del_x(xi_cell)
    end function get_cell_del_x_curv

    function get_cell_del_y_curv(self, eta_cell) result(del_y)
        ! from eta cell value get del_y
        class(domain_curv), intent(in) :: self
        integer(int32), intent(in) :: eta_cell
        real(real64) :: del_y
        del_y = self%del_y(eta_cell)
    end function get_cell_del_y_curv

    subroutine binary_search_left_idx(array, n, x, idx_lower)
        real(real64), intent(in) :: array(n), x
        integer(int32), intent(in) :: n
        integer(int32), intent(in out) :: idx_lower
        integer(int32) :: idx_higher, idx_middle
        idx_lower = 1
        idx_higher = n
        ! if ((x<self%grid_X(1)) .or. (x > self%grid_X(self%N_x))) then
        !     print *, 'x value outside of grid range in getLFromX!'
        !     stop
        ! end if
        do while (idx_lower /= idx_higher-1)
            idx_middle = (idx_lower + idx_higher)/2
            if (array(idx_middle) <= x) then
                idx_lower = idx_middle
            else
                idx_higher = idx_middle
            end if
        end do
    end subroutine binary_search_left_idx

    subroutine create_sin_grid(n, length, grid, del, smallest_del)
        ! generate sinusoidal symmetric grid with smallest cell at boundary
        integer, intent(in) :: n
        real(real64), intent(in) :: length, smallest_del ! smallest cell size
        real(real64), intent(in out) :: del(n-1), grid(n) ! cell sizes to place in
        integer(int32) :: i
        grid(1) = 0.0d0
        grid(n) = length
        do i = 2,n-1
            grid(i) = length * ((real(i)-1.0d0)/(real(n) - 1.0d0) - (1.0d0/(real(n) - 1.0d0) - smallest_del/length) &
            * SIN(2.0d0 * pi_const * (i-1) / real(n - 1)) / SIN(2.0d0 * pi_const / real(n - 1)) )
        end do
        do i = 1, n-1
            del(i) = grid(i+1) - grid(i)
        end do

    end subroutine create_sin_grid

end module mod_domain_curv