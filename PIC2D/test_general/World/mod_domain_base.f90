module mod_domain_base
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    implicit none

    ! Base world type which will contain information about the world
    ! Computational grid is rectangular, so set N_x and N_y node numbers
    ! Will extend to curvilinear grid or uniform grid
    ! Assume all inner nodes nodes are boundaries
    ! For orthogonal grid so will only need grid information in x and y


    type :: domain_base
        ! store grid quantities
        real(real64), allocatable :: grid_X(:), grid_Y(:) ! spatial location in grid
        integer(int32), allocatable :: boundary_conditions(:,:) ! type of node in each direction 
        real(real64) :: start_X, end_X, start_Y, end_Y! start and end locations of grid as reference
        integer(int32) :: N_x, N_y, N_x_cells, N_y_cells  ! amount of nodes in x and y direction
    contains
        procedure, public, pass(self) :: get_xi_from_X
        procedure, public, pass(self) :: get_eta_from_Y
        procedure, public, pass(self) :: get_cell_del_x
        procedure, public, pass(self) :: get_cell_del_y
    end type

    interface domain_base
        module procedure :: constructor
    end interface domain_base

contains
    type(domain_base) function constructor(N_x, N_y, length_x, length_y) result(self)
        ! Construct domain object, initialize grid, dx_dl, and nodeVol.
        integer(int32), intent(in) :: N_x, N_y
        real(real64), intent(in) :: length_x, length_y
        integer(int32) :: left_bound, right_bound, top_bound, bottom_bound
        self%N_x = N_x
        self%N_y = N_y
        self%start_X = 0.0d0
        self%start_Y = 0.0d0
        self%end_X = length_x
        self%end_Y = length_y
        allocate(self%grid_X(self%N_x), self%grid_Y(self%N_y),self%boundary_conditions(self%N_x, self%N_y))

        ! self%NESW_wall_boundaries = NESW_wall_boundaries ! N, E, S, W indices
        ! top_bound = NESW_wall_boundaries(1)
        ! right_bound = NESW_wall_boundaries(2)
        ! bottom_bound = NESW_wall_boundaries(3)
        ! left_bound = NESW_wall_boundaries(4)

        !$OMP parallel
        !$OMP workshare
        self%boundary_conditions = 0
        !$OMP end workshare
        !$OMP end parallel
        ! self%boundary_conditions(2:self%N_x-1, self%N_y) = top_bound
        ! self%boundary_conditions(self%N_x, 2:self%N_y-1) = right_bound
        ! self%boundary_conditions(2:self%N_x-1, 1) = bottom_bound
        ! self%boundary_conditions(1, 2:self%N_y-1) = left_bound


        ! ! Upper right corner
        ! self%boundary_conditions(self%N_x, self%N_y) = MIN(top_bound, right_bound)

        ! ! lower right corner
        ! self%boundary_conditions(self%N_x, 1) = MIN(bottom_bound, right_bound)

        ! ! lower left corner
        ! self%boundary_conditions(1, 1) = MIN(bottom_bound, left_bound)

        ! ! upper left corner
        ! self%boundary_conditions(1, self%N_y) = MIN(top_bound, left_bound)
        
    end function constructor
    
    function get_xi_from_X(self, x) result(xi)
        ! get xi in computational space from x
        class(domain_base), intent(in) :: self
        real(real64), intent(in) :: x
        real(real64) :: xi
        xi = 0.0d0
    end function get_xi_from_X

    function get_eta_from_Y(self, y) result(eta)
        ! get eta in computational space from y
        class(domain_base), intent(in) :: self
        real(real64), intent(in) :: y
        real(real64) :: eta
        eta = 0.0d0
    end function get_eta_from_Y

    function get_cell_del_x(self, xi_cell) result(del_x)
        ! from xi cell value get del_x
        class(domain_base), intent(in) :: self
        integer(int32), intent(in) :: xi_cell
        real(real64) :: del_x
        del_x = 0.0d0
    end function get_cell_del_x

    function get_cell_del_y(self, eta_cell) result(del_y)
        ! from eta cell value get del_y
        class(domain_base), intent(in) :: self
        integer(int32), intent(in) :: eta_cell
        real(real64) :: del_y
        del_y = 0.0d0
    end function get_cell_del_y


end module mod_domain_base