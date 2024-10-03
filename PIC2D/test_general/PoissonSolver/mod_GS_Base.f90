module mod_GS_Base
    use iso_fortran_env, only: int32, int64, real64
    use mod_domain_base
    use omp_lib
    implicit none

    ! Have tested, having indexed references is faster than using if statements for entire domain using boundary conditions (factor of 2)
    type :: GS_Base
        ! store grid quantities
        real(real64), allocatable :: sourceTerm(:,:), solution(:,:), residual(:,:)
        integer(int32), allocatable :: number_row_sections(:), start_inner_indx_x(:,:), end_inner_indx_x(:,:)
        integer(int32) :: number_inner_rows, number_solve_nodes, max_number_row_sections, iterNumber, N_x, N_y, start_row_indx, end_row_indx
        integer(int32) :: start_coarse_row_indx, start_fine_row_indx
        ! Get non-dirichlet indices on boundaries
        integer :: number_bottom_row_sections, number_top_row_sections, number_left_column_sections, number_right_column_sections
        integer, allocatable :: start_bottom_row_indx(:), end_bottom_row_indx(:), bottom_row_boundary_type(:)
        integer, allocatable :: start_top_row_indx(:), end_top_row_indx(:), top_row_boundary_type(:)
        integer, allocatable :: start_left_column_indx(:), end_left_column_indx(:), left_column_boundary_type(:)
        integer, allocatable :: start_right_column_indx(:), end_right_column_indx(:), right_column_boundary_type(:)
        integer :: x_indx_step, y_indx_step
        real(real64) :: omega, inv_omega
        class(domain_base), pointer :: world
    contains
        procedure, public, pass(self) :: solveGS
        procedure, public, pass(self) :: smoothIterations
        procedure, public, pass(self) :: smoothWithRes
        procedure, public, pass(self) :: calcResidual
        procedure, public, pass(self) :: restriction
        procedure, public, pass(self) :: prolongation
        procedure, public, pass(self) :: AX_Mult
        procedure, public, pass(self) :: YAX_Mult
        procedure, public, pass(self) :: initialize_GS_orthogonal_base
    end type


contains

    subroutine initialize_GS_orthogonal_base(self, omega, world, N_x, N_y)
        ! Construct derived object for everything needed from GS_Base_Even class
        class(GS_Base), intent(in out) :: self
        integer, intent(in) :: N_x, N_y
        real(real64), intent(in) :: omega
        class(domain_base), intent(in), target :: world
        integer :: i_coarse, i_fine, j_coarse, j_fine, k, numberBound, min_indx
        logical :: found_first_indx
        self%omega = omega
        self%inv_omega = 1.0d0 - omega
        self%iterNumber = 0
        self%N_x = N_x
        self%N_y = N_y
        self%world => world
        allocate(self%sourceTerm(self%N_x, self%N_y), self%solution(self%N_x, self%N_y), self%residual(self%N_x, self%N_y))   
        self%x_indx_step = (self%world%N_x-1)/(self%N_x-1) ! equivalent index step sizes in finest grid
        self%y_indx_step = (self%world%N_y-1)/(self%N_y-1)
        
        found_first_indx = .false.
        ! Set source terms and solutions to 0
        
        !$OMP parallel
        !$OMP workshare
        self%residual = 0.0d0
        !$OMP end workshare
        !$OMP workshare
        self%solution = 0.0d0
        !$OMP end workshare
        !$OMP workshare
        self%sourceTerm = 0.0d0
        !$OMP end workshare
        !$OMP end parallel

        numberBound = 0
        ! Count number red and black nodes
        ! Also fill in any dirichlet boundaries solution vector
        !$OMP parallel reduction(+:numberBound)
        !$OMP do collapse(2)
        do j_fine = 1, self%world%N_y, self%y_indx_step
            do i_fine = 1, self%world%N_x, self%x_indx_step
                if (self%world%boundary_conditions(i_fine,j_fine) /= 1) then
                    numberBound = numberBound + 1
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel

        self%number_solve_nodes = numberBound

        ! Find number of rows
        k = 0
        min_indx = self%N_y
        !$OMP parallel private(i_fine, j_fine) reduction(+:k) reduction(min:min_indx)
        !$OMP do
        do j_fine = 1, self%world%N_y, self%y_indx_step
            do i_fine = 1, self%world%N_x, self%x_indx_step
                if (self%world%boundary_conditions(i_fine,j_fine) == 0) then
                    k = k + 1
                    min_indx = min(min_indx, j_fine)
                    exit
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        self%number_inner_rows = k
        self%start_row_indx = (min_indx + (self%y_indx_step-1))/self%y_indx_step! needs to be converted to coarse grid
        self%end_row_indx = self%start_row_indx + self%number_inner_rows-1

        allocate(self%number_row_sections(self%number_inner_rows))

        ! Find amount of different inner node sections for each row
        k = 0
        !$OMP parallel private(i_coarse, j_coarse, j_fine, i_fine, k, found_first_indx)
        !$OMP do
        do j_coarse = self%start_row_indx, self%end_row_indx
            j_fine = j_coarse * self%y_indx_step - (self%y_indx_step-1)
            k = 0
            found_first_indx = .false.
            do i_fine = 1, self%world%N_x, self%x_indx_step
                if (.not. found_first_indx) then
                    if (self%world%boundary_conditions(i_fine, j_fine) == 0) then
                        found_first_indx = .true.
                    end if
                else
                    if (self%world%boundary_conditions(i_fine, j_fine) /= 0) then
                        found_first_indx = .false.
                        k = k + 1
                    end if
                end if
            end do
            self%number_row_sections(j_coarse - self%start_row_indx + 1) = k
        end do
        !$OMP end do
        !$OMP end parallel
        self%max_number_row_sections = MAXVAL(self%number_row_sections) ! Use to allocate array of maximum, probably not that wasteful in memory, could also create object of allocatable 1D arrays
        allocate(self%start_inner_indx_x(self%max_number_row_sections, self%number_inner_rows), &
        self%end_inner_indx_x(self%max_number_row_sections, self%number_inner_rows))

        ! get start and end indx for each inner node section
        k = 0
        !$OMP parallel private(i_coarse, j_coarse, j_fine, i_fine, k, found_first_indx)
        !$OMP do
        do j_coarse = self%start_row_indx, self%end_row_indx
            j_fine = j_coarse * self%y_indx_step - (self%y_indx_step-1)
            k = 0
            found_first_indx = .false.
            do i_fine = 1, self%world%N_x, self%x_indx_step
                i_coarse = (i_fine + (self%x_indx_step-1))/self%x_indx_step
                if (.not. found_first_indx) then
                    if (self%world%boundary_conditions(i_fine, j_fine) == 0) then
                        found_first_indx = .true.
                        k = k + 1
                        self%start_inner_indx_x(k, j_coarse - self%start_row_indx + 1) = i_coarse
                    end if
                else
                    if (self%world%boundary_conditions(i_fine, j_fine) /= 0) then
                        found_first_indx = .false.
                        self%end_inner_indx_x(k, j_coarse - self%start_row_indx + 1) = i_coarse - 1
                    end if
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel

        ! ------------------------ get numbers for each boundary side -------------------- 

        ! lower boundary
        k = 0
        found_first_indx = .false.
        do i_coarse = 2, self%N_x-1
            i_fine = i_coarse * self%x_indx_step - self%x_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(i_fine,1) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                end if
            else
                if (self%world%boundary_conditions(i_fine,1) == 1) then
                    found_first_indx = .false.
                end if
            end if
        end do
        self%number_bottom_row_sections = k
        allocate(self%start_bottom_row_indx(self%number_bottom_row_sections), self%end_bottom_row_indx(self%number_bottom_row_sections), &
        self%bottom_row_boundary_type(self%number_bottom_row_sections))

        k = 0
        found_first_indx = .false.
        do i_coarse = 2, self%N_x-1
            i_fine = i_coarse * self%x_indx_step - self%x_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(i_fine,1) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                    self%start_bottom_row_indx(k) = i_coarse
                    self%bottom_row_boundary_type(k) = self%world%boundary_conditions(i_fine,1)
                end if
            else
                if (self%world%boundary_conditions(i_fine,1) == 1) then
                    found_first_indx = .false.
                    self%end_bottom_row_indx(k) = i_coarse-1
                end if
            end if
        end do
        if (found_first_indx) self%end_bottom_row_indx(k) = self%N_x-1
        

        ! upper boundary
        k = 0
        found_first_indx = .false.
        do i_coarse = 2, self%N_x-1
            i_fine = i_coarse * self%x_indx_step - self%x_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(i_fine,self%world%N_y) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                end if
            else
                if (self%world%boundary_conditions(i_fine,self%world%N_y) == 1) then
                    found_first_indx = .false.
                end if
            end if
        end do
        self%number_top_row_sections = k
        allocate(self%start_top_row_indx(self%number_top_row_sections), self%end_top_row_indx(self%number_top_row_sections), &
        self%top_row_boundary_type(self%number_top_row_sections))

        k = 0
        found_first_indx = .false.
        do i_coarse = 2, self%N_x-1
            i_fine = i_coarse * self%x_indx_step - self%x_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(i_fine,self%world%N_y) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                    self%start_top_row_indx(k) = i_coarse
                    self%top_row_boundary_type(k) = self%world%boundary_conditions(i_fine,self%world%N_y)
                end if
            else
                if (self%world%boundary_conditions(i_fine,self%world%N_y) == 1) then
                    found_first_indx = .false.
                    self%end_top_row_indx(k) = i_coarse-1
                end if
            end if
        end do
        if (found_first_indx) self%end_top_row_indx(k) = self%N_x-1


        ! left boundary
        k = 0
        found_first_indx = .false.
        do j_coarse = 2, self%N_y-1
            j_fine = j_coarse * self%y_indx_step - self%y_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(1,j_fine) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                end if
            else
                if (self%world%boundary_conditions(1, j_fine) == 1) then
                    found_first_indx = .false.
                end if
            end if
        end do
        self%number_left_column_sections = k
        allocate(self%start_left_column_indx(self%number_left_column_sections), self%end_left_column_indx(self%number_left_column_sections), &
        self%left_column_boundary_type(self%number_left_column_sections))

        k = 0
        found_first_indx = .false.
        do j_coarse = 2, self%N_y-1
            j_fine = j_coarse * self%y_indx_step - self%y_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(1,j_fine) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                    self%start_left_column_indx(k) = j_coarse
                    self%left_column_boundary_type(k) = self%world%boundary_conditions(1,j_fine)
                end if
            else
                if (self%world%boundary_conditions(1, j_fine) == 1) then
                    found_first_indx = .false.
                    self%end_left_column_indx(k) = j_coarse - 1
                end if
            end if
        end do
        if (found_first_indx) self%end_left_column_indx(k) = self%N_y-1

        ! right boundary
        k = 0
        found_first_indx = .false.
        do j_coarse = 2, self%N_y-1
            j_fine = j_coarse * self%y_indx_step - self%y_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(self%world%N_x,j_fine) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                end if
            else
                if (self%world%boundary_conditions(self%world%N_x, j_fine) == 1) then
                    found_first_indx = .false.
                end if
            end if
        end do
        self%number_right_column_sections = k
        allocate(self%start_right_column_indx(self%number_right_column_sections), self%end_right_column_indx(self%number_right_column_sections), &
        self%right_column_boundary_type(self%number_right_column_sections))

        k = 0
        found_first_indx = .false.
        do j_coarse = 2, self%N_y-1
            j_fine = j_coarse * self%y_indx_step - self%y_indx_step + 1
            if (.not. found_first_indx) then
                if (self%world%boundary_conditions(self%world%N_x,j_fine) /= 1) then
                    found_first_indx = .true.
                    k = k + 1
                    self%start_right_column_indx(k) = j_coarse
                    self%right_column_boundary_type(k) = self%world%boundary_conditions(self%world%N_x,j_fine)
                end if
            else
                if (self%world%boundary_conditions(self%world%N_x, j_fine) == 1) then
                    found_first_indx = .false.
                    self%end_right_column_indx(k) = j_coarse - 1
                end if
            end if
        end do
        if (found_first_indx) self%end_right_column_indx(k) = self%N_y-1

        ! ----------- for restriction/prolongation ---------
        ! get starting row index for coarse nodes
        if (MOD(self%start_row_indx,2) == 1) then
            self%start_coarse_row_indx = 1
            self%start_fine_row_indx = 2
        else
            self%start_coarse_row_indx = 2
            self%start_fine_row_indx = 1
        end if


    end subroutine initialize_GS_orthogonal_base
    
    subroutine solveGS(self, stepTol, relTol)
        ! Solve GS down to some tolerance
        class(GS_Base), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        real(real64) :: stepRes, relRes, R2_init, R2
        stepRes = 1.0d0
        relRes = 1.0d0
        self%iterNumber = 0
        call self%calcResidual()
        !$OMP parallel workshare
        R2_init = SUM(self%residual**2)
        !$OMP end parallel workshare
        R2_init = SQRT(R2_init)
        do while (stepRes > stepTol .and. relRes > relTol)
            ! single iterations slightly faster
            call self%smoothIterations(100)
            stepRes = self%smoothWithRes()
            call self%calcResidual()
            !$OMP parallel workshare
            R2 = SUM(self%residual**2)
            !$OMP end parallel workshare
            relRes = sqrt(R2)/R2_init
            self%iterNumber = self%iterNumber + 101
        end do
        if (stepRes > stepTol) then
            print *, 'Converged due to step size'
        else
            print *, 'Converged due to R2 decrease'
        end if

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