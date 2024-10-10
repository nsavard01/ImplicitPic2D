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
        integer(int32) :: N_x, N_y, N_x_cells, N_y_cells ! amount of nodes in x and y direction
        integer(int64) :: number_total_cells 
    contains
        procedure, public, pass(self) :: get_xi_from_X
        procedure, public, pass(self) :: get_eta_from_Y
        procedure, public, pass(self) :: get_cell_del_x
        procedure, public, pass(self) :: get_cell_del_y
        procedure, public, pass(self) :: get_number_cells
        procedure, public, pass(self) :: generate_node_volume
        procedure, public, pass(self) :: form_boundary_conditions
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

        !$OMP parallel
        !$OMP workshare
        self%boundary_conditions = 0
        !$OMP end workshare
        !$OMP end parallel
        
    end function constructor

    subroutine get_number_cells(self)
        class(domain_base), intent(in out) :: self
        integer :: j, i
        integer(int64) :: num
        num = 0
        !$OMP parallel private(j,i) reduction(+:num)
        !$OMP do
        do j = 1, self%N_y-1
            do i = 1, self%N_x-1
                if (self%boundary_conditions(i,j) == 0 .or. self%boundary_conditions(i+1,j) == 0 &
                .or. self%boundary_conditions(i,j+1) == 0 .or. self%boundary_conditions(i+1,j+1) == 0) then
                    num = num + 1
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        self%number_total_cells = num
    end subroutine get_number_cells

    subroutine form_boundary_conditions(self, upperBound, rightBound, lowerBound, leftBound, inner_box_first_x, inner_box_last_x, inner_box_first_y, inner_box_last_y, center_box_bool)
        class(domain_base), intent(in out) :: self
        integer, intent(in) :: upperBound, rightBound, lowerBound, leftBound
        integer, intent(in out) :: inner_box_first_x, inner_box_last_x, inner_box_first_y, inner_box_last_y
        logical, intent(in) :: center_box_bool

        self%boundary_conditions(2:self%N_x-1, self%N_y) = upperBound
        self%boundary_conditions(self%N_x, 2:self%N_y-1) = rightBound
        self%boundary_conditions(2:self%N_x-1, 1) = lowerBound
        self%boundary_conditions(1, 2:self%N_y-1) = leftBound


        ! Upper right corner
        self%boundary_conditions(self%N_x, self%N_y) = MIN(upperBound, rightBound)

        ! lower right corner
        self%boundary_conditions(self%N_x, 1) = MIN(lowerBound, rightBound)

        ! lower left corner
        self%boundary_conditions(1, 1) = MIN(lowerBound, leftBound)

        ! upper left corner
        self%boundary_conditions(1, self%N_y) = MIN(upperBound, leftBound)

        if (center_box_bool) then
            inner_box_first_x = 3* self%N_x/8
            inner_box_last_x = 5 * self%N_x / 8
            inner_box_first_y = 3 * self%N_y / 8
            inner_box_last_y = 5 * self%N_y / 8
            !$OMP parallel
            !$OMP workshare
            self%boundary_conditions(inner_box_first_x:inner_box_last_x, inner_box_first_y:inner_box_last_y) = 1
            !$OMP end workshare
            !$OMP end parallel
        end if
    end subroutine form_boundary_conditions

    subroutine generate_node_volume(self)
        class(domain_base), intent(in out) :: self

    end subroutine generate_node_volume
    
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