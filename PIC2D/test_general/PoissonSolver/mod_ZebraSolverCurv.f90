module mod_ZebraSolverCurv
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base_Curv
    use mod_domain_curv
    use omp_lib
    implicit none

    ! Class for a Black Red gauss siedel smoother
    ! Since will be for stage of multigrid, and to work with Black red points, N_x, N_y must be odd
    ! Example order for black (b) and red (r) points for 5x5 matrix is:
    ! 
    ! b - r - b - r - b
    ! r - b - r - b - r
    ! b - r - b - r - b
    ! r - b - r - b - r
    ! b - r - b - r - b
    ! Must be orthogonal grid


    private
    public :: ZebraSolverCurv

    type, extends(GS_Base_Curv) :: ZebraSolverCurv
        ! store grid quantities
        integer, allocatable :: number_column_sections(:), start_inner_indx_y(:,:), end_inner_indx_y(:,:)
        integer :: number_inner_columns, start_column_indx, end_column_indx, max_number_column_sections
        integer :: inner_row_idx_white_start, inner_row_idx_black_start
        integer :: inner_column_idx_white_start, inner_column_idx_black_start
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal => constructPoissonOrthogonal_ZebraCurv
        procedure, public, pass(self) :: smoothIterations => smoothIterations_ZebraCurv
        procedure, public, pass(self) :: solve_row
        procedure, public, pass(self) :: solve_row_res
        procedure, public, pass(self) :: solve_column
        procedure, public, pass(self) :: smoothWithRes => smoothWithRes_ZebraCurv
        ! procedure, public, pass(self) :: restriction
        ! procedure, public, pass(self) :: prolongation
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type

    interface ZebraSolverCurv
        module procedure :: ZebraSolverCurv_constructor
    end interface ZebraSolverCurv

contains

    type(ZebraSolverCurv) function ZebraSolverCurv_constructor(omega, world, N_x, N_y) result(self)
        ! Construct object, set initial variables
        real(real64), intent(in) :: omega
        integer, intent(in) :: N_x, N_y
        type(domain_curv), intent(in), target :: world
        call self%initialize_GS_orthogonal_base(omega, world, N_x, N_y)
        call self%initialize_GS_Curv(world)
        call self%constructPoissonOrthogonal()
    end function ZebraSolverCurv_constructor

    subroutine constructPoissonOrthogonal_ZebraCurv(self)
        ! Construct orthogonal grid solver
        class(ZebraSolverCurv), intent(in out) :: self
        integer :: min_indx, k, j_fine, i_fine, i_coarse, j_coarse
        logical :: found_first_indx

        ! Find number of columns and minimum indx
        k = 0
        min_indx = self%N_x
        !$OMP parallel private(i_fine, j_fine) reduction(+:k) reduction(min:min_indx)
        !$OMP do
        do i_fine = 1, self%world%N_x, self%x_indx_step
            do j_fine = 1, self%world%N_y, self%y_indx_step    
                if (self%world%boundary_conditions(i_fine,j_fine) == 0) then
                    k = k + 1
                    min_indx = min(min_indx, i_fine)
                    exit
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        self%number_inner_columns = k
        self%start_column_indx = (min_indx + (self%x_indx_step-1))/self%x_indx_step! needs to be converted to coarse grid
        self%end_column_indx = self%start_column_indx + self%number_inner_columns-1

        allocate(self%number_column_sections(self%number_inner_columns))

        ! Find amount of different inner node sections for each row
        k = 0
        !$OMP parallel private(i_coarse, j_coarse, j_fine, i_fine, k, found_first_indx)
        !$OMP do
        do i_coarse = self%start_column_indx, self%end_column_indx
            i_fine = i_coarse * self%x_indx_step - self%x_indx_step + 1
            k = 0
            found_first_indx = .false.
            do j_fine = 1, self%world%N_y, self%y_indx_step
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
            self%number_column_sections(i_coarse - self%start_column_indx + 1) = k
        end do
        !$OMP end do
        !$OMP end parallel
        self%max_number_column_sections = MAXVAL(self%number_column_sections) ! Use to allocate array of maximum, probably not that wasteful in memory, could also create object of allocatable 1D arrays
        allocate(self%start_inner_indx_y(self%max_number_column_sections, self%number_inner_columns), &
        self%end_inner_indx_y(self%max_number_column_sections, self%number_inner_columns))

        ! get start and end indx for each inner node section
        k = 0
        !$OMP parallel private(i_coarse, j_coarse, j_fine, i_fine, k, found_first_indx)
        !$OMP do
        do i_coarse = self%start_column_indx, self%end_column_indx
            i_fine = i_coarse * self%x_indx_step - self%x_indx_step+1
            k = 0
            found_first_indx = .false.
            do j_fine = 1, self%world%N_y, self%y_indx_step
                j_coarse = (j_fine + (self%y_indx_step-1))/self%y_indx_step
                if (.not. found_first_indx) then
                    if (self%world%boundary_conditions(i_fine, j_fine) == 0) then
                        found_first_indx = .true.
                        k = k + 1
                        self%start_inner_indx_y(k, i_coarse - self%start_column_indx + 1) = j_coarse
                    end if
                else
                    if (self%world%boundary_conditions(i_fine, j_fine) /= 0) then
                        found_first_indx = .false.
                        self%end_inner_indx_y(k, i_coarse - self%start_column_indx + 1) = j_coarse - 1
                    end if
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel

        ! ------------ find start/end indices for white/black rows/columns ------------------

        ! white starts on first row
        if (MOD(self%start_row_indx,2) == 1) then
            ! if odd, white row is starting point
            self%inner_row_idx_white_start = 1
            self%inner_row_idx_black_start = 2
        else
            self%inner_row_idx_black_start = 1
            self%inner_row_idx_white_start = 2
        end if

        ! white starts on first column
        if (MOD(self%start_column_indx,2) == 1) then
            ! if odd, white row is starting point
            self%inner_column_idx_white_start = 1
            self%inner_column_idx_black_start = 2
        else
            self%inner_column_idx_black_start = 1
            self%inner_column_idx_white_start = 2
        end if
        
    end subroutine constructPoissonOrthogonal_ZebraCurv


    ! subroutine matMult(self, x, y)
    !     ! Use gauss-seidel to calculate A*x, store solution into y
    !     class(ZebraSolver), intent(in out) :: self
    !     real(real64), intent(in) :: x(self%N_x, self%N_y)
    !     real(real64), intent(in out) :: y(self%N_x, self%N_y)
    !     integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j
    !     real(real64) :: C_N, C_E, C_O, C_W, C_S

    !     !$OMP parallel private(i, j, N_indx, E_indx, S_indx, W_indx, C_N, C_E, C_O, C_W, C_S)
    !     !$OMP do
    !     do k = 1, self%numberBlackInnerNodes
    !         i = self%black_InnerIndx(1, k)
    !         j = self%black_InnerIndx(2, k)
    !         N_indx = j + 1
    !         E_indx = i+1
    !         S_indx = j-1
    !         W_indx = i-1
    !         C_O = self%matCoeffsInnerBlack(1, k)
    !         C_N = self%matCoeffsInnerBlack(2, k)
    !         C_E = self%matCoeffsInnerBlack(3, k)
    !         C_S = self%matCoeffsInnerBlack(4, k)
    !         C_W = self%matCoeffsInnerBlack(5, k)
    !         y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
    !     end do
    !     !$OMP end do nowait
    !     !$OMP do
    !     do k = 1, self%numberBlackBoundNodes
    !         i = self%black_NESW_BoundIndx(1, k)
    !         j = self%black_NESW_BoundIndx(2, k)
    !         E_indx = self%black_NESW_BoundIndx(3, k)
    !         W_indx = self%black_NESW_BoundIndx(4, k)
    !         N_indx = self%black_NESW_BoundIndx(5, k)
    !         S_indx = self%black_NESW_BoundIndx(6, k)
    !         C_O = self%matCoeffsBoundBlack(1, k)
    !         C_N = self%matCoeffsBoundBlack(2, k)
    !         C_E = self%matCoeffsBoundBlack(3, k)
    !         C_S = self%matCoeffsBoundBlack(4, k)
    !         C_W = self%matCoeffsBoundBlack(5, k)
    !         y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
    !     end do
    !     !$OMP end do nowait
    !     !$OMP do
    !     do k = 1, self%numberRedInnerNodes
    !         i = self%red_InnerIndx(1, k)
    !         j = self%red_InnerIndx(2, k)
    !         N_indx = j + 1
    !         E_indx = i+1
    !         S_indx = j-1
    !         W_indx = i-1
    !         C_O = self%matCoeffsInnerRed(1, k)
    !         C_N = self%matCoeffsInnerRed(2, k)
    !         C_E = self%matCoeffsInnerRed(3, k)
    !         C_S = self%matCoeffsInnerRed(4, k)
    !         C_W = self%matCoeffsInnerRed(5, k)
    !         y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
    !     end do
    !     !$OMP end do nowait
    !     !$OMP do
    !     do k = 1, self%numberRedBoundNodes
    !         i = self%red_NESW_BoundIndx(1, k)
    !         j = self%red_NESW_BoundIndx(2, k)
    !         E_indx = self%red_NESW_BoundIndx(3, k)
    !         W_indx = self%red_NESW_BoundIndx(4, k)
    !         N_indx = self%red_NESW_BoundIndx(5, k)
    !         S_indx = self%red_NESW_BoundIndx(6, k)
    !         C_O = self%matCoeffsBoundRed(1, k)
    !         C_N = self%matCoeffsBoundRed(2, k)
    !         C_E = self%matCoeffsBoundRed(3, k)
    !         C_S = self%matCoeffsBoundRed(4, k)
    !         C_W = self%matCoeffsBoundRed(5, k)
    !         y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
    !     end do
    !     !$OMP end do
    !     !$OMP end parallel
         
    ! end subroutine matMult

    ! function XAX_Mult(self, x) result(res)
    !     ! Use gauss-seidel to calculate x^T * A * x
    !     class(ZebraSolver), intent(in out) :: self
    !     real(real64), intent(in) :: x(self%N_x, self%N_y)
    !     integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j
    !     real(real64) :: res, C_N, C_E, C_O, C_W, C_S
    !     res = 0.0d0
    !     !$OMP parallel private(i, j, N_indx, E_indx, S_indx, W_indx, C_N, C_E, C_O, C_W, C_S)
    !     !$OMP do
    !     do k = 1, self%numberBlackInnerNodes
    !         i = self%black_InnerIndx(1, k)
    !         j = self%black_InnerIndx(2, k)
    !         N_indx = j + 1
    !         E_indx = i+1
    !         S_indx = j-1
    !         W_indx = i-1
    !         C_O = self%matCoeffsInnerBlack(1, k)
    !         C_N = self%matCoeffsInnerBlack(2, k)
    !         C_E = self%matCoeffsInnerBlack(3, k)
    !         C_S = self%matCoeffsInnerBlack(4, k)
    !         C_W = self%matCoeffsInnerBlack(5, k)
    !         res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
    !     end do
    !     !$OMP end do nowait
    !     !$OMP do
    !     do k = 1, self%numberBlackBoundNodes
    !         i = self%black_NESW_BoundIndx(1, k)
    !         j = self%black_NESW_BoundIndx(2, k)
    !         E_indx = self%black_NESW_BoundIndx(3, k)
    !         W_indx = self%black_NESW_BoundIndx(4, k)
    !         N_indx = self%black_NESW_BoundIndx(5, k)
    !         S_indx = self%black_NESW_BoundIndx(6, k)
    !         C_O = self%matCoeffsBoundBlack(1, k)
    !         C_N = self%matCoeffsBoundBlack(2, k)
    !         C_E = self%matCoeffsBoundBlack(3, k)
    !         C_S = self%matCoeffsBoundBlack(4, k)
    !         C_W = self%matCoeffsBoundBlack(5, k)
    !         res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
    !     end do
    !     !$OMP end do nowait
    !     !$OMP do
    !     do k = 1, self%numberRedInnerNodes
    !         i = self%red_InnerIndx(1, k)
    !         j = self%red_InnerIndx(2, k)
    !         N_indx = j + 1
    !         E_indx = i+1
    !         S_indx = j-1
    !         W_indx = i-1
    !         C_O = self%matCoeffsInnerRed(1, k)
    !         C_N = self%matCoeffsInnerRed(2, k)
    !         C_E = self%matCoeffsInnerRed(3, k)
    !         C_S = self%matCoeffsInnerRed(4, k)
    !         C_W = self%matCoeffsInnerRed(5, k)
    !         res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
    !     end do
    !     !$OMP end do nowait
    !     !$OMP do
    !     do k = 1, self%numberRedBoundNodes
    !         i = self%red_NESW_BoundIndx(1, k)
    !         j = self%red_NESW_BoundIndx(2, k)
    !         E_indx = self%red_NESW_BoundIndx(3, k)
    !         W_indx = self%red_NESW_BoundIndx(4, k)
    !         N_indx = self%red_NESW_BoundIndx(5, k)
    !         S_indx = self%red_NESW_BoundIndx(6, k)
    !         C_O = self%matCoeffsBoundRed(1, k)
    !         C_N = self%matCoeffsBoundRed(2, k)
    !         C_E = self%matCoeffsBoundRed(3, k)
    !         C_S = self%matCoeffsBoundRed(4, k)
    !         C_W = self%matCoeffsBoundRed(5, k)
    !         res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
    !             + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
    !     end do
    !     !$OMP end do
    !     !$OMP end parallel
         
    ! end function XAX_Mult


    subroutine smoothIterations_ZebraCurv(self, iterNum)
        ! Solve GS down to some tolerance
        class(ZebraSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, iter, i_fine, j_fine, idx_start, idx_end, boundary_start, boundary_end, vector_size
        real(real64) :: C_N, C_E, C_S, C_W
        ! p indexs in horizontal maps to i, k in vertical maps to j

        do iter = 1, iterNum 
            
            !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, &
            !$OMP& i_fine, j_fine, idx_start, idx_end, vector_size, boundary_start, boundary_end, &
            !$OMP& C_N, C_E, C_S, C_W)

            ! white rows 

            !$OMP do
            do k = self%inner_row_idx_white_start, self%number_inner_rows, 2
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                j_fine = j * self%y_indx_step - self%y_indx_step + 1
                C_N = self%inner_node_coeff_north(j-1)
                C_S = self%inner_node_coeff_south(j-1)
                do p = 1, self%number_row_sections(k)

                    ! Determine boundary properties
                    idx_start = self%start_inner_indx_x(p,k)
                    idx_end = self%end_inner_indx_x(p,k)
                    vector_size = idx_end - idx_start + 3
                    
                    boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, j_fine)
                    boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, j_fine)


                    call self%solve_row(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                end do
            end do
            !$OMP end do nowait

            !$OMP sections
            !$OMP section
            ! bottom boundary
            if (self%number_bottom_row_sections /= 0) then
                N_indx = 2
                j = 1
                do p = 1, self%number_bottom_row_sections
                    if (self%bottom_row_boundary_type(p) == 2) then
                        S_indx = N_indx
                    else
                        S_indx = self%N_y-1
                    end if
                    C_S = self%bottom_row_coeff_south(p)
                    C_N = self%bottom_row_coeff_north(p)
                     ! Determine boundary properties
                    idx_start = self%start_bottom_row_indx(p)
                    idx_end = self%end_bottom_row_indx(p)
                    vector_size = idx_end - idx_start + 3
                    
                    boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, 1)
                    boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, 1)


                    call self%solve_row(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                end do
            end if
            !$OMP section
            ! top boundary
            if (self%number_top_row_sections /= 0) then
                S_indx = self%N_y-1
                j = self%N_y
                do p = 1, self%number_top_row_sections
                    if (self%top_row_boundary_type(p) == 2) then
                        N_indx = S_indx
                    else
                        N_indx = 2
                    end if
                    C_S = self%top_row_coeff_south(p)
                    C_N = self%top_row_coeff_north(p)
                     ! Determine boundary properties
                    idx_start = self%start_top_row_indx(p)
                    idx_end = self%end_top_row_indx(p)
                    vector_size = idx_end - idx_start + 3
                    
                    boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, self%world%N_y)
                    boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, self%world%N_y)

                    call self%solve_row(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                end do
            end if
            !$OMP end sections

            !$OMP barrier
            
            ! black columns

            !$OMP do
            do k = self%inner_column_idx_black_start, self%number_inner_columns, 2
                i = self%start_column_indx + k - 1
                E_indx = i + 1
                W_indx = i - 1
                i_fine = i * self%x_indx_step - self%x_indx_step + 1
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                do p = 1, self%number_column_sections(k)
                    ! Determine boundary properties
                    idx_start = self%start_inner_indx_y(p,k)
                    idx_end = self%end_inner_indx_y(p,k)
                    vector_size = idx_end - idx_start + 3
                    boundary_start = self%world%boundary_conditions(i_fine, (idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                    boundary_end = self%world%boundary_conditions(i_fine, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                    call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
                end do
            end do
            !$OMP end do
            
            !$OMP barrier

            ! black row
            !$OMP do
            do k = self%inner_row_idx_black_start, self%number_inner_rows, 2
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                j_fine = j * self%y_indx_step - self%y_indx_step + 1
                C_N = self%inner_node_coeff_north(j-1)
                C_S = self%inner_node_coeff_south(j-1)
                do p = 1, self%number_row_sections(k)

                    ! Determine boundary properties
                    idx_start = self%start_inner_indx_x(p,k)
                    idx_end = self%end_inner_indx_x(p,k)
                    vector_size = idx_end - idx_start + 3
                    boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, j_fine)
                    boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, j_fine)

                    call self%solve_row(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                end do
            end do
            !$OMP end do

            !$OMP barrier

            ! white columns
            !$OMP do
            do k = self%inner_column_idx_white_start, self%number_inner_columns, 2
                i = self%start_column_indx + k - 1
                E_indx = i + 1
                W_indx = i - 1
                i_fine = i * self%x_indx_step - self%x_indx_step + 1
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                do p = 1, self%number_column_sections(k)
                    ! Determine boundary properties
                    idx_start = self%start_inner_indx_y(p,k)
                    idx_end = self%end_inner_indx_y(p,k)
                    vector_size = idx_end - idx_start + 3
                    boundary_start = self%world%boundary_conditions(i_fine, (idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                    boundary_end = self%world%boundary_conditions(i_fine, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                    call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
                end do
            end do
            !$OMP end do nowait

            !$OMP sections
            !$OMP section
            ! left boundary
            if (self%number_left_column_sections /= 0) then
                E_indx = 2
                i = 1
                do p = 1, self%number_left_column_sections
                    if (self%left_column_boundary_type(p) == 2) then
                        W_indx = E_indx
                    else
                        W_indx = self%N_x-1
                    end if
                    C_E = self%left_column_coeff_east(p)
                    C_W = self%left_column_coeff_west(p)
                     ! Determine boundary properties
                    idx_start = self%start_left_column_indx(p)
                    idx_end = self%end_left_column_indx(p)
                    vector_size = idx_end - idx_start + 3
                    
                    boundary_start = self%world%boundary_conditions(1,(idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                    boundary_end = self%world%boundary_conditions(1, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                    call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
                end do
            end if
            !$OMP section
            ! right boundary
            if (self%number_right_column_sections /= 0) then
                W_indx = self%N_x-1
                i = self%N_x
                do p = 1, self%number_right_column_sections
                    if (self%right_column_boundary_type(p) == 2) then
                        E_indx = W_indx
                    else
                        E_indx = 2
                    end if
                    C_E = self%right_column_coeff_east(p)
                    C_W = self%right_column_coeff_west(p)
                     ! Determine boundary properties
                    idx_start = self%start_right_column_indx(p)
                    idx_end = self%end_right_column_indx(p)
                    vector_size = idx_end - idx_start + 3
                    
                    boundary_start = self%world%boundary_conditions(self%world%N_x,(idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                    boundary_end = self%world%boundary_conditions(self%world%N_x, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                    call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
                end do
            end if
            !$OMP end sections
            
            !$OMP end parallel
        end do
    end subroutine smoothIterations_ZebraCurv

    function smoothWithRes_ZebraCurv(self) result(Res)
        ! Solve GS down to some tolerance
        class(ZebraSolverCurv), intent(in out) :: self
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, iter, i_fine, j_fine, idx_start, idx_end, boundary_start, boundary_end, vector_size
        real(real64) :: C_N, C_E, C_S, C_W, Res, res_rows
        ! p indexs in horizontal maps to i, k in vertical maps to j

        Res = 0.0d0
        !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, &
        !$OMP& i_fine, j_fine, idx_start, idx_end, vector_size, boundary_start, boundary_end, &
        !$OMP& C_N, C_E, C_S, C_W, res_rows) reduction(+:Res)

        ! white rows 

        !$OMP do
        do k = self%inner_row_idx_white_start, self%number_inner_rows, 2
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            j_fine = j * self%y_indx_step - self%y_indx_step + 1
            C_N = self%inner_node_coeff_north(j-1)
            C_S = self%inner_node_coeff_south(j-1)
            do p = 1, self%number_row_sections(k)

                ! Determine boundary properties
                idx_start = self%start_inner_indx_x(p,k)
                idx_end = self%end_inner_indx_x(p,k)
                vector_size = idx_end - idx_start + 3
                
                boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, j_fine)
                boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, j_fine)


                res_rows = self%solve_row_res(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                Res = Res + res_rows
            end do
        end do
        !$OMP end do nowait

        !$OMP sections
        !$OMP section
        ! bottom boundary
        if (self%number_bottom_row_sections /= 0) then
            N_indx = 2
            j = 1
            do p = 1, self%number_bottom_row_sections
                if (self%bottom_row_boundary_type(p) == 2) then
                    S_indx = N_indx
                else
                    S_indx = self%N_y-1
                end if
                C_S = self%bottom_row_coeff_south(p)
                C_N = self%bottom_row_coeff_north(p)
                    ! Determine boundary properties
                idx_start = self%start_bottom_row_indx(p)
                idx_end = self%end_bottom_row_indx(p)
                vector_size = idx_end - idx_start + 3
                
                boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, 1)
                boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, 1)


                res_rows = self%solve_row_res(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                Res = Res + res_rows
            end do
        end if
        !$OMP section
        ! top boundary
        if (self%number_top_row_sections /= 0) then
            S_indx = self%N_y-1
            j = self%N_y
            do p = 1, self%number_top_row_sections
                if (self%top_row_boundary_type(p) == 2) then
                    N_indx = S_indx
                else
                    N_indx = 2
                end if
                C_S = self%top_row_coeff_south(p)
                C_N = self%top_row_coeff_north(p)
                    ! Determine boundary properties
                idx_start = self%start_top_row_indx(p)
                idx_end = self%end_top_row_indx(p)
                vector_size = idx_end - idx_start + 3
                
                boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, self%world%N_y)
                boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, self%world%N_y)

                res_rows = self%solve_row_res(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                Res = Res + res_rows
            end do
        end if
        !$OMP end sections

        !$OMP barrier
        
        ! black columns

        !$OMP do
        do k = self%inner_column_idx_black_start, self%number_inner_columns, 2
            i = self%start_column_indx + k - 1
            E_indx = i + 1
            W_indx = i - 1
            i_fine = i * self%x_indx_step - self%x_indx_step + 1
            C_E = self%inner_node_coeff_east(i-1)
            C_W = self%inner_node_coeff_west(i-1)
            do p = 1, self%number_column_sections(k)
                ! Determine boundary properties
                idx_start = self%start_inner_indx_y(p,k)
                idx_end = self%end_inner_indx_y(p,k)
                vector_size = idx_end - idx_start + 3
                boundary_start = self%world%boundary_conditions(i_fine, (idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                boundary_end = self%world%boundary_conditions(i_fine, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
            end do
        end do
        !$OMP end do
        
        !$OMP barrier

        ! black row
        !$OMP do
        do k = self%inner_row_idx_black_start, self%number_inner_rows, 2
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            j_fine = j * self%y_indx_step - self%y_indx_step + 1
            C_N = self%inner_node_coeff_north(j-1)
            C_S = self%inner_node_coeff_south(j-1)
            do p = 1, self%number_row_sections(k)

                ! Determine boundary properties
                idx_start = self%start_inner_indx_x(p,k)
                idx_end = self%end_inner_indx_x(p,k)
                vector_size = idx_end - idx_start + 3
                boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, j_fine)
                boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, j_fine)

                res_rows = self%solve_row_res(idx_start, idx_end, boundary_start, boundary_end, vector_size, j, N_indx, S_indx, C_N, C_S)
                Res = Res + res_rows
            end do
        end do
        !$OMP end do

        !$OMP barrier

        ! white columns
        !$OMP do
        do k = self%inner_column_idx_white_start, self%number_inner_columns, 2
            i = self%start_column_indx + k - 1
            E_indx = i + 1
            W_indx = i - 1
            i_fine = i * self%x_indx_step - self%x_indx_step + 1
            C_E = self%inner_node_coeff_east(i-1)
            C_W = self%inner_node_coeff_west(i-1)
            do p = 1, self%number_column_sections(k)
                ! Determine boundary properties
                idx_start = self%start_inner_indx_y(p,k)
                idx_end = self%end_inner_indx_y(p,k)
                vector_size = idx_end - idx_start + 3
                boundary_start = self%world%boundary_conditions(i_fine, (idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                boundary_end = self%world%boundary_conditions(i_fine, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
            end do
        end do
        !$OMP end do nowait

        !$OMP sections
        !$OMP section
        ! left boundary
        if (self%number_left_column_sections /= 0) then
            E_indx = 2
            i = 1
            do p = 1, self%number_left_column_sections
                if (self%left_column_boundary_type(p) == 2) then
                    W_indx = E_indx
                else
                    W_indx = self%N_x-1
                end if
                C_E = self%left_column_coeff_east(p)
                C_W = self%left_column_coeff_west(p)
                    ! Determine boundary properties
                idx_start = self%start_left_column_indx(p)
                idx_end = self%end_left_column_indx(p)
                vector_size = idx_end - idx_start + 3
                
                boundary_start = self%world%boundary_conditions(1,(idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                boundary_end = self%world%boundary_conditions(1, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
            end do
        end if
        !$OMP section
        ! right boundary
        if (self%number_right_column_sections /= 0) then
            W_indx = self%N_x-1
            i = self%N_x
            do p = 1, self%number_right_column_sections
                if (self%right_column_boundary_type(p) == 2) then
                    E_indx = W_indx
                else
                    E_indx = 2
                end if
                C_E = self%right_column_coeff_east(p)
                C_W = self%right_column_coeff_west(p)
                    ! Determine boundary properties
                idx_start = self%start_right_column_indx(p)
                idx_end = self%end_right_column_indx(p)
                vector_size = idx_end - idx_start + 3
                
                boundary_start = self%world%boundary_conditions(self%world%N_x,(idx_start-1) * self%y_indx_step - self%y_indx_step + 1)
                boundary_end = self%world%boundary_conditions(self%world%N_x, (idx_end+1) * self%y_indx_step - self%y_indx_step + 1)

                call self%solve_column(idx_start, idx_end, boundary_start, boundary_end, vector_size, i, E_indx, W_indx, C_E, C_W)
            end do
        end if
        !$OMP end sections
        
        !$OMP end parallel
        Res = SQRT(Res/ (self%number_solve_nodes))

    end function smoothWithRes_ZebraCurv


    subroutine solve_row(self, idx_start, idx_end, boundary_start, boundary_end, vector_size, j_grid, N_indx, S_indx, C_N, C_S)
        class(ZebraSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: idx_start, idx_end, boundary_start, boundary_end, vector_size, j_grid, N_indx, S_indx
        real(real64), intent(in) :: C_N, C_S
        integer(int32) :: i_grid, i_local
        real(real64) :: source_thread(vector_size), cp_thread(vector_size-1), dp_thread(vector_size-1), m_thread, C_E, C_W, C_O
        
        !left boundary
        i_grid = idx_start - 1
        select case (boundary_start)
        case (1)
            cp_thread(1) = 0.0d0
            source_thread(1) = self%solution(i_grid,j_grid)
            dp_thread(1) = self%solution(i_grid,j_grid)
        case (2)
            ! Get Coefficients 
            do i_local = 1, self%number_left_column_sections
                if (j_grid >= self%start_left_column_indx(i_local)-1 .and. j_grid <= self%end_left_column_indx(i_local)+1) then
                    C_E = self%left_column_coeff_east(i_local)
                    exit
                end if
            end do
            C_O = self%centerCoeff(i_grid, j_grid)
            if (boundary_end /= 2) then
                ! dirichlet boundary on other side
                cp_thread(1) = 2.0d0 * C_E * C_O
                source_thread(1) = self%sourceTerm(i_grid, j_grid) - C_N * self%solution(i_grid, N_indx)  - C_S * self%solution(i_grid, S_indx)
                dp_thread(1) = source_thread(1) * C_O
            else
                ! another neumann on other side, need to do point iteration to get guess for phi at boundary
                cp_thread(1) = 0.0d0
                source_thread(1) = (self%sourceTerm(i_grid,j_grid) - C_N * self%solution(i_grid, N_indx) - C_S * self%solution(i_grid, S_indx)- &
                    2.0d0 * C_E * self%solution(2, j_grid)) * C_O * self%omega + self%inv_omega * self%solution(i_grid,j_grid)
                dp_thread(1) = source_thread(1)
            end if
        case (3)
            ! Point iteration
            do i_local = 1, self%number_left_column_sections
                if (j_grid >= self%start_left_column_indx(i_local) .and. j_grid <= self%end_left_column_indx(i_local)) then
                    C_E = self%left_column_coeff_east(i_local)
                    C_W = self%left_column_coeff_west(i_local)
                    exit
                end if
            end do
            C_O = self%centerCoeff(i_grid,j_grid)
            cp_thread(1) = 0.0d0
            source_thread(1) = (self%sourceTerm(i_grid,j_grid) - C_N * self%solution(i_grid, N_indx) - C_S * self%solution(i_grid, S_indx) - &
                C_E * self%solution(2, j_grid) - C_W * self%solution(self%N_x-1, j_grid)) * C_O * self%omega + self%inv_omega * self%solution(i_grid,j_grid)
            dp_thread(1) = source_thread(1)
        end select
    
        ! forward substitution
        do i_grid = idx_start,idx_end
            i_local = i_grid - idx_start + 2
            C_E = self%inner_node_coeff_east(i_grid-1)
            C_W = self%inner_node_coeff_west(i_grid-1)
            C_O = 1.0d0/self%centerCoeff(i_grid, j_grid)
            m_thread = C_O - cp_thread(i_local-1) * C_W
            cp_thread(i_local) = C_E / m_thread
            source_thread(i_local) = self%sourceTerm(i_grid,j_grid) - C_N * self%solution(i_grid, N_indx) - C_S * self%solution(i_grid, S_indx)
            dp_thread(i_local) = (source_thread(i_local) - dp_thread(i_local-1) * C_W) / m_thread
        end do

        ! rightmost boundary properties
        i_grid = idx_end + 1
        select case(boundary_end)
        case (1)
            source_thread(vector_size) = self%solution(i_grid, j_grid)
        case (2)
            do i_local = 1, self%number_right_column_sections
                if (j_grid >= self%start_right_column_indx(i_local)-1 .and. j_grid <= self%end_right_column_indx(i_local)+1) then
                    C_W = self%right_column_coeff_west(i_local)
                    exit
                end if
            end do
            C_O = 1.0d0/self%centerCoeff(i_grid,j_grid)
            source_thread(vector_size) = self%sourceTerm(i_grid,j_grid) - C_N * self%solution(i_grid, N_indx) - C_S * self%solution(i_grid, S_indx)
            m_thread = 2.0d0 * C_W
            source_thread(vector_size) = (source_thread(vector_size) - dp_thread(vector_size-1) * m_thread)/&
                (C_O - cp_thread(vector_size-1)  * m_thread)
        case (3)
            source_thread(vector_size) = self%solution(1,j_grid)
        end select

        ! backwards substitution
        do i_local = vector_size-1, 1, -1
            source_thread(i_local) = dp_thread(i_local)-cp_thread(i_local)*source_thread(i_local+1)
        end do

        self%solution(idx_start-1:idx_end+1,j_grid) = source_thread

    end subroutine solve_row

    function solve_row_res(self, idx_start, idx_end, boundary_start, boundary_end, vector_size, j_grid, N_indx, S_indx, C_N, C_S) result(res)
        class(ZebraSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: idx_start, idx_end, boundary_start, boundary_end, vector_size, j_grid, N_indx, S_indx
        real(real64), intent(in) :: C_N, C_S 
        real(real64) :: oldSol(vector_size), res
        oldSol = self%solution(idx_start-1:idx_end+1,j_grid)
        call self%solve_row(idx_start, idx_end, boundary_start, boundary_end, vector_size, j_grid, N_indx, S_indx, C_N, C_S)
        res = SUM((self%solution(idx_start-1:idx_end+1,j_grid) - oldSol)**2)

    end function solve_row_res


    subroutine solve_column(self, idx_start, idx_end, boundary_start, boundary_end, vector_size, i_grid, E_indx, W_indx, C_E, C_W)
        class(ZebraSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: idx_start, idx_end, boundary_start, boundary_end, vector_size, i_grid, E_indx, W_indx
        real(real64), intent(in) :: C_E, C_W
        integer(int32) :: j_grid, j_local
        real(real64) :: source_thread(vector_size), cp_thread(vector_size-1), dp_thread(vector_size-1), m_thread, C_S, C_N, C_O
        

        !left boundary
        j_grid = idx_start - 1
        select case (boundary_start)
        case (1)
            cp_thread(1) = 0.0d0
            source_thread(1) = self%solution(i_grid,j_grid)
            dp_thread(1) = self%solution(i_grid,j_grid)
        case (2)
            ! Get Coefficients 
            do j_local = 1, self%number_bottom_row_sections
                if (i_grid >= self%start_bottom_row_indx(j_local)-1 .and. i_grid <= self%end_bottom_row_indx(j_local)+1) then
                    C_N = self%bottom_row_coeff_north(j_local)
                    exit
                end if
            end do
            C_O = self%centerCoeff(i_grid, j_grid)
            if (boundary_end /= 2) then
                ! dirichlet boundary on other side
                cp_thread(1) = 2.0d0 * C_N * C_O
                source_thread(1) = self%sourceTerm(i_grid, j_grid) - C_E * self%solution(E_indx, j_grid) - C_W * self%solution(W_indx, j_grid)
                dp_thread(1) = source_thread(1) * C_O
            else
                ! another neumann on other side, need to do point iteration to get guess for phi at boundary
                cp_thread(1) = 0.0d0
                source_thread(1) = (self%sourceTerm(i_grid,j_grid) - C_E * self%solution(E_indx, j_grid) - C_W * self%solution(W_indx, j_grid) - &
                    2.0d0 * C_N * self%solution(i_grid, 2)) * C_O * self%omega + self%inv_omega * self%solution(i_grid,j_grid)
                dp_thread(1) = source_thread(1)
            end if
        case (3)
            ! Point iteration
            do j_local = 1, self%number_bottom_row_sections
                if (i_grid >= self%start_bottom_row_indx(j_local)-1 .and. i_grid <= self%end_bottom_row_indx(j_local)+1) then
                    C_N = self%bottom_row_coeff_north(j_local)
                    C_S = self%bottom_row_coeff_south(j_local)
                    exit
                end if
            end do
            C_O = self%centerCoeff(i_grid,j_grid)
            cp_thread(1) = 0.0d0
            source_thread(1) = (self%sourceTerm(i_grid,j_grid) - C_E * self%solution(E_indx, j_grid) - C_W * self%solution(W_indx, j_grid) - &
                C_N * self%solution(i_grid, 2) - C_S * self%solution(i_grid, self%N_y-1)) * C_O * self%omega + self%inv_omega * self%solution(i_grid,j_grid)
            dp_thread(1) = source_thread(1)
        end select

        ! forward substitution
        do j_grid = idx_start,idx_end
            j_local = j_grid - idx_start + 2
            C_N = self%inner_node_coeff_north(j_grid-1)
            C_S = self%inner_node_coeff_south(j_grid-1)
            C_O = 1.0d0/self%centerCoeff(i_grid, j_grid)
            m_thread = C_O - cp_thread(j_local-1) * C_S
            cp_thread(j_local) = C_N / m_thread
            source_thread(j_local) = self%sourceTerm(i_grid,j_grid) - C_E * self%solution(E_indx, j_grid) - C_W * self%solution(W_indx, j_grid)
            dp_thread(j_local) = (source_thread(j_local) - dp_thread(j_local-1) * C_S) / m_thread
        end do

        ! rightmost boundary properties
        j_grid = idx_end + 1
        select case(boundary_end)
        case (1)
            source_thread(vector_size) = self%solution(i_grid, j_grid)
        case (2)
            do j_local = 1, self%number_top_row_sections
                if (i_grid >= self%start_top_row_indx(j_local)-1 .and. i_grid <= self%end_top_row_indx(j_local)+1) then
                    C_S = self%top_row_coeff_south(j_local)
                    exit
                end if
            end do
            C_O = 1.0d0/self%centerCoeff(i_grid, j_grid)
            source_thread(vector_size) = self%sourceTerm(i_grid,j_grid) - C_E * self%solution(E_indx, j_grid) - C_W * self%solution(W_indx, j_grid)
            m_thread = 2.0d0 * C_S
            source_thread(vector_size) = (source_thread(vector_size) - dp_thread(vector_size-1) * m_thread)/&
                (C_O - cp_thread(vector_size-1)  * m_thread)
        case (3)
            source_thread(vector_size) = self%solution(i_grid,1)
        end select

        ! backwards substitution
        do j_local = vector_size-1, 1, -1
            source_thread(j_local) = dp_thread(j_local)-cp_thread(j_local)*source_thread(j_local+1)
        end do

        self%solution(i_grid,idx_start-1:idx_end+1) = source_thread

    end subroutine solve_column


end module mod_ZebraSolverCurv