module mod_ZebraSolverEven
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base_Even
    use mod_domain_uniform
    use omp_lib
    implicit none

    ! Class for a Zebra line block gauss siedel method
    ! Implicitly solve for each row using tridiagonal direct solver
    ! separate into white and black lines which are separated so they don't interfere (allows parallelization)
    ! do alternating white rows, black columns, black rows, white columns: find it's the best for general anisotropic grids
    ! 
    ! w - w - w - w - w
    ! b - b - b - b - b
    ! w - w - w - w - w
    ! b - b - b - b - b
    ! w - w - w - w - w
    ! Must be orthogonal grid


    private
    public :: ZebraSolverEven

    type, extends(GS_Base_Even) :: ZebraSolverEven
        integer, allocatable :: number_column_sections(:), start_inner_indx_y(:,:), end_inner_indx_y(:,:)
        integer :: number_inner_columns, start_column_indx, end_column_indx, max_number_column_sections
        integer :: inner_row_idx_white_start, inner_row_idx_black_start
        integer :: inner_column_idx_white_start, inner_column_idx_black_start
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal => constructPoissonOrthogonal_ZebraEven
        procedure, public, pass(self) :: smoothIterations => smoothIterations_ZebraEven
        ! procedure, public, pass(self) :: smoothWithRes => smoothWithRes_ZebraEven
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type
    interface ZebraSolverEven
        module procedure :: ZebraSolverEven_constructor
    end interface ZebraSolverEven

contains
    type(ZebraSolverEven) function ZebraSolverEven_constructor(omega, world, N_x, N_y) result(self)
        ! Construct object, set initial variables
        real(real64), intent(in) :: omega
        integer, intent(in) :: N_x, N_y
        type(domain_uniform), intent(in), target :: world
        call self%initialize_GS_orthogonal_base(omega, world, N_x, N_y)
        call self%initialize_GS_Even(world)
        call self%constructPoissonOrthogonal()
    end function ZebraSolverEven_constructor

    subroutine constructPoissonOrthogonal_ZebraEven(self)
        ! Construct orthogonal grid solver
        class(ZebraSolverEven), intent(in out) :: self
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

    end subroutine constructPoissonOrthogonal_ZebraEven


    ! subroutine matMult(self, x, y)
    !     ! Use gauss-seidel to calculate A*x, store solution into y
    !     class(ZebraSolverEven), intent(in out) :: self
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

    subroutine smoothIterations_ZebraEven(self, iterNum)
        ! Solve GS down to some tolerance
        class(ZebraSolverEven), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        real(real64) :: omega_inv, inv_centerCoeff
        real(real64) :: cp_thread(MAX(self%N_x, self%N_y)-1), dp_thread(MAX(self%N_x, self%N_y)-1), source_thread(MAX(self%N_x, self%N_y)), m_thread
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, iter, i_fine, j_fine, idx_start, idx_end, boundary_start, boundary_end
        omega_inv = 1.0d0 - self%omega
        inv_centerCoeff = 1.0d0/self%centerCoeff
        ! p indexs in horizontal maps to i, k in vertical maps to j

        do iter = 1, iterNum 
            ! !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, cp_y, cp_x, dp_y, dp_x, sourceX, sourceY, m)

            ! ! horizontal sweeps left to right
            ! !$OMP do
            ! ! Sweep odd numbered rows
            ! do k = 1, , 2
            !     j = self%startRow + k - 1
            !     N_indx = self%vertIndx(1, k)
            !     S_indx = self%vertIndx(2, k)
            !     select case (self%leftBound)
            !     case (1)
            !         cp_x(1) = 0.0d0
            !         sourceX(1) = self%solution(1,j)
            !         dp_x(1) = sourceX(1)
            !     case (2)
            !         cp_x(1) = 2.0d0 * self%coeffX * self%centerCoeff
            !         sourceX(1) = self%sourceTerm(1, j) - self%coeffY * (self%solution(1, N_indx)  + self%solution(1, S_indx))
            !         dp_x(1) = (sourceX(1)) * self%centerCoeff
            !     case (3,4)
            !         cp_x(1) = 0.0d0
            !         E_indx = self%horzIndx(1,1)
            !         W_indx = self%horzIndx(2,1)
            !         sourceX(1) = (self%sourceTerm(1,j) - self%coeffY * (self%solution(1, N_indx) + self%solution(1, S_indx))- &
            !         self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))) * self%centerCoeff* self%omega + omega_inv * self%solution(1,j)
            !         dp_x(1) = sourceX(1)
            !     end select
            !     ! solve for vectors c-prime and d-prime
            !     do i = 2,self%N_x-1
            !         m = inv_centerCoeff - cp_x(i-1) * self%coeffX
            !         cp_x(i) = self%coeffX / m
            !         sourceX(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
            !         dp_x(i) = (sourceX(i) - dp_x(i-1) * self%coeffX) / m
            !     end do
            !     select case(self%rightBound)
            !     case (1)
            !         sourceX(self%N_x) = self%solution(self%N_x, j)
            !     case (2)
            !         sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - self%coeffY * ( self%solution(self%N_x, N_indx) + self%solution(self%N_x, S_indx))
            !         m = 2.0d0 * self%coeffX
            !         sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * m)/&
            !             (inv_centerCoeff - cp_x(self%N_x-1)  * m)
            !     case (3)
            !         sourceX(self%N_x) = sourceX(1)
            !     end select
            !     do i = self%N_x-1, 1, -1
            !         sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
            !     end do
            !     self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
            ! end do
            ! !$OMP end do
            
            ! !$OMP do
            ! ! sweep even numbered columns
            ! do p = 2, self%numberColumns, 2
            !     i = self%startCol + p - 1
            !     E_indx = self%horzIndx(1, p)
            !     W_indx = self%horzIndx(2, p)
            !     select case (self%lowerBound)
            !     case (1)
            !         cp_y(1) = 0.0d0
            !         sourceY(1) = self%solution(i,1)
            !         dp_y(1) = sourceY(1)
            !     case (2)
            !         cp_y(1) = 2.0d0 * self%coeffY * self%centerCoeff
            !         sourceY(1) = self%sourceTerm(i, 1) - self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))
            !         dp_y(1) = (sourceY(1)) * self%centerCoeff
            !     case (3,4)
            !         cp_y(1) = 0.0d0
            !         N_indx = self%vertIndx(1,1)
            !         S_indx = self%vertIndx(2,1)
            !         sourceY(1) = (self%sourceTerm(i,1) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx)) - &
            !         self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
            !         dp_y(1) = sourceY(1)
            !     end select
            !     ! solve for vectors c-prime and d-prime
            !     do j = 2,self%N_y-1
            !         m = inv_centerCoeff - cp_y(j-1) * self%coeffY
            !         cp_y(j) = self%coeffY / m
            !         sourceY(j) = self%sourceTerm(i, j) - self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))
            !         dp_y(j) = (sourceY(j) - dp_y(j-1) * self%coeffY) / m
            !     end do
            !     select case(self%upperBound)
            !     case (1)
            !         sourceY(self%N_y) = self%solution(i, self%N_y)
            !     case (2)
            !         sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - self%coeffX * ( self%solution(E_indx, self%N_y) + self%solution(W_indx, self%N_y))
            !         m = 2.0d0 * self%coeffY
            !         sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * m)/&
            !             (inv_centerCoeff - cp_y(self%N_y-1)  *m)
            !     case (3)
            !         sourceY(self%N_y) = sourceY(1)
            !     end select
            !     do j = self%N_y-1, 1, -1
            !         sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
            !     end do
            !     self%solution(i,:) = sourceY *self%omega + omega_inv * self%solution(i,:)
            ! end do
            ! !$OMP end do

            ! ! Sweep horizontal even rows
            ! !$OMP do
            ! do k = 2, self%numberRows, 2
            !     j = self%startRow + k - 1
            !     N_indx = self%vertIndx(1, k)
            !     S_indx = self%vertIndx(2, k)
            !     select case (self%leftBound)
            !     case (1)
            !         cp_x(1) = 0.0d0
            !         sourceX(1) = self%solution(1,j)
            !         dp_x(1) = sourceX(1)
            !     case (2)
            !         cp_x(1) = 2.0d0 * self%coeffX * self%centerCoeff
            !         sourceX(1) = self%sourceTerm(1, j) - self%coeffY * (self%solution(1, N_indx)  + self%solution(1, S_indx))
            !         dp_x(1) = (sourceX(1)) * self%centerCoeff
            !     case (3,4)
            !         cp_x(1) = 0.0d0
            !         E_indx = self%horzIndx(1,1)
            !         W_indx = self%horzIndx(2,1)
            !         sourceX(1) = (self%sourceTerm(1,j) - self%coeffY * (self%solution(1, N_indx) + self%solution(1, S_indx))- &
            !         self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))) * self%centerCoeff* self%omega + omega_inv * self%solution(1,j)
            !         dp_x(1) = sourceX(1)
            !     end select
            !     ! solve for vectors c-prime and d-prime
            !     do i = 2,self%N_x-1
            !         m = inv_centerCoeff - cp_x(i-1) * self%coeffX
            !         cp_x(i) = self%coeffX / m
            !         sourceX(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
            !         dp_x(i) = (sourceX(i) - dp_x(i-1) * self%coeffX) / m
            !     end do
            !     select case(self%rightBound)
            !     case (1)
            !         sourceX(self%N_x) = self%solution(self%N_x, j)
            !     case (2)
            !         sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - self%coeffY * ( self%solution(self%N_x, N_indx) + self%solution(self%N_x, S_indx))
            !         m = 2.0d0 * self%coeffX
            !         sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * m)/&
            !             (inv_centerCoeff - cp_x(self%N_x-1)  * m)
            !     case (3)
            !         sourceX(self%N_x) = sourceX(1)
            !     end select
            !     do i = self%N_x-1, 1, -1
            !         sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
            !     end do
            !     self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
            ! end do
            ! !$OMP end do

            ! ! vertical sweeps bottom to top
            ! !$OMP do
            ! ! Sweep odd numbers upwards
            ! do p = 1, self%numberColumns, 2
            !     i = self%startCol + p - 1
            !     E_indx = self%horzIndx(1, p)
            !     W_indx = self%horzIndx(2, p)
            !     select case (self%lowerBound)
            !     case (1)
            !         cp_y(1) = 0.0d0
            !         sourceY(1) = self%solution(i,1)
            !         dp_y(1) = sourceY(1)
            !     case (2)
            !         cp_y(1) = 2.0d0 * self%coeffY * self%centerCoeff
            !         sourceY(1) = self%sourceTerm(i, 1) - self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))
            !         dp_y(1) = (sourceY(1)) * self%centerCoeff
            !     case (3,4)
            !         cp_y(1) = 0.0d0
            !         N_indx = self%vertIndx(1,1)
            !         S_indx = self%vertIndx(2,1)
            !         sourceY(1) = (self%sourceTerm(i,1) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx)) - &
            !         self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
            !         dp_y(1) = sourceY(1)
            !     end select
            !     ! solve for vectors c-prime and d-prime
            !     do j = 2,self%N_y-1
            !         m = inv_centerCoeff - cp_y(j-1) * self%coeffY
            !         cp_y(j) = self%coeffY / m
            !         sourceY(j) = self%sourceTerm(i, j) - self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))
            !         dp_y(j) = (sourceY(j) - dp_y(j-1) * self%coeffY) / m
            !     end do
            !     select case(self%upperBound)
            !     case (1)
            !         sourceY(self%N_y) = self%solution(i, self%N_y)
            !     case (2)
            !         sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - self%coeffX * ( self%solution(E_indx, self%N_y) + self%solution(W_indx, self%N_y))
            !         m = 2.0d0 * self%coeffY
            !         sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * m)/&
            !             (inv_centerCoeff - cp_y(self%N_y-1)  *m)
            !     case (3)
            !         sourceY(self%N_y) = sourceY(1)
            !     end select
            !     do j = self%N_y-1, 1, -1
            !         sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
            !     end do
            !     self%solution(i,:) = sourceY *self%omega + omega_inv * self%solution(i,:)
            ! end do
            ! !$OMP end do
            ! !$OMP end parallel

            !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, cp_thread, dp_thread, &
            !$OMP& source_thread, m_thread, i_fine, j_fine, idx_start, idx_end)

            ! white rows 

            !$OMP do
            do k = self%inner_row_idx_white_start, self%number_inner_rows, 2
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                j_fine = j * self%y_indx_step - self%y_indx_step + 1
                do p = 1, self%number_row_sections(k)

                    ! Determine boundary properties
                    idx_start = self%start_inner_indx_x(p,k)
                    idx_end = self%end_inner_indx_x(p,k)
                    boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, j_fine)
                    boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, j_fine)


                    !left boundary
                    i = idx_start - 1
                    select case (boundary_start)
                    case (1)
                        cp_thread(i) = 0.0d0
                        source_thread(i) = self%solution(i,j)
                        dp_thread(i) = self%solution(i,j)
                    case (2)
                        if (boundary_end /= 2) then
                            ! dirichlet boundary on other side
                            cp_thread(i) = 2.0d0 * self%coeffX * self%centerCoeff
                            source_thread(i) = self%sourceTerm(i, j) - 2.0d0 * self%coeffY * (self%solution(i, N_indx)  + self%solution(i, S_indx))
                            dp_thread(i) = source_thread(i) * self%centerCoeff
                        else
                            ! another neumann on other side, need to do point iteration to get guess for phi at boundary
                            cp_thread(i) = 0.0d0
                            source_thread(i) = (self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))- &
                                2.0d0 * self%coeffX * self%solution(2, j)) * self%centerCoeff* self%omega + omega_inv * self%solution(i,j)
                            dp_thread(i) = source_thread(i)
                        end if
                    case (3)
                        ! Point iteration
                        cp_thread(i) = 0.0d0
                        source_thread(i) = (self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))- &
                        self%coeffX * (self%solution(2, j) + self%solution(self%N_x-1, j))) * self%centerCoeff* self%omega + omega_inv * self%solution(i,j)
                        dp_thread(i) = source_thread(i)
                    end select

                    ! forward substitution
                    do i = idx_start,idx_end
                        m_thread = inv_centerCoeff - cp_thread(i-1) * self%coeffX
                        cp_thread(i) = self%coeffX / m_thread
                        source_thread(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
                        dp_thread(i) = (source_thread(i) - dp_thread(i-1) * self%coeffX) / m_thread
                    end do

                    ! rightmost boundary properties
                    i = idx_end + 1
                    select case(boundary_end)
                    case (1)
                        source_thread(i) = self%solution(i, j)
                    case (2)
                        source_thread(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
                        m_thread = 2.0d0 * self%coeffX
                        source_thread(i) = (source_thread(i) - dp_thread(i-1) * m_thread)/&
                            (inv_centerCoeff - cp_thread(i-1)  * m_thread)
                    case (3)
                        source_thread(i) = source_thread(idx_start-1)
                    end select

                    ! backwards substitution
                    do i = idx_end, idx_start-1, -1
                        source_thread(i) = dp_thread(i)-cp_thread(i)*source_thread(i+1)
                    end do

                    self%solution(idx_start-1:idx_end+1,j) = source_thread(idx_start-1:idx_end+1)
                end do
            end do
            !$OMP end do nowait

            ! !$OMP sections
            ! !$OMP section
            ! ! bottom boundary
            ! if (self%number_bottom_row_sections /= 0) then
            !     do p = 1, self%number_bottom_row_sections
            !         if (self%bottom_row_boundary_type(p) == 2) then
            !             S_indx = 2
            !         else
            !             S_indx = self%N_y-1
            !         end if
            !         do i = self%start_bottom_row_indx(p), self%end_bottom_row_indx(p)
            !             self%solution(i,1) = 500.0d0
            !         end do
            !     end do
            ! end if
            ! !$OMP section
            ! ! top boundary
            ! if (self%number_top_row_sections /= 0) then
            !     do p = 1, self%number_top_row_sections
            !         if (self%top_row_boundary_type(p) == 2) then
            !             S_indx = 2
            !         else
            !             S_indx = self%N_y-1
            !         end if
            !         do i = self%start_top_row_indx(p), self%end_top_row_indx(p)
            !             self%solution(i,self%N_y) = 500.0d0
            !         end do
            !     end do
            ! end if
            ! !$OMP end sections

            !$OMP barrier
            
            ! black columns

            ! !$OMP do
            ! do k = self%inner_column_idx_black_start, self%number_inner_columns, 2
            !     i = self%start_column_indx + k - 1
            !     E_indx = i + 1
            !     W_indx = j - 1
            !     do p = 1, self%number_column_sections(k)
            !         do j = self%start_inner_indx_y(p,k), self%end_inner_indx_y(p,k)
            !             self%solution(i,j) = 500.0d0
            !         end do
            !     end do
            ! end do
            ! !$OMP end do
            
            !$OMP barrier

            ! black row
            !$OMP do
            do k = self%inner_row_idx_black_start, self%number_inner_rows, 2
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                j_fine = j * self%y_indx_step - self%y_indx_step + 1
                do p = 1, self%number_row_sections(k)

                    ! Determine boundary properties
                    idx_start = self%start_inner_indx_x(p,k)
                    idx_end = self%end_inner_indx_x(p,k)
                    boundary_start = self%world%boundary_conditions((idx_start-1) * self%x_indx_step - self%x_indx_step + 1, j_fine)
                    boundary_end = self%world%boundary_conditions((idx_end+1) * self%x_indx_step - self%x_indx_step + 1, j_fine)


                    !left boundary
                    i = idx_start - 1
                    select case (boundary_start)
                    case (1)
                        cp_thread(i) = 0.0d0
                        source_thread(i) = self%solution(i,j)
                        dp_thread(i) = self%solution(i,j)
                    case (2)
                        if (boundary_end /= 2) then
                            ! dirichlet boundary on other side
                            cp_thread(i) = 2.0d0 * self%coeffX * self%centerCoeff
                            source_thread(i) = self%sourceTerm(i, j) - 2.0d0 * self%coeffY * (self%solution(i, N_indx)  + self%solution(i, S_indx))
                            dp_thread(i) = source_thread(i) * self%centerCoeff
                        else
                            ! another neumann on other side, need to do point iteration to get guess for phi at boundary
                            cp_thread(i) = 0.0d0
                            source_thread(i) = (self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))- &
                                2.0d0 * self%coeffX * self%solution(2, j)) * self%centerCoeff* self%omega + omega_inv * self%solution(i,j)
                            dp_thread(i) = source_thread(i)
                        end if
                    case (3)
                        ! Point iteration
                        cp_thread(i) = 0.0d0
                        source_thread(i) = (self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))- &
                        self%coeffX * (self%solution(2, j) + self%solution(self%N_x-1, j))) * self%centerCoeff* self%omega + omega_inv * self%solution(i,j)
                        dp_thread(i) = source_thread(i)
                    end select

                    ! forward substitution
                    do i = idx_start,idx_end
                        m_thread = inv_centerCoeff - cp_thread(i-1) * self%coeffX
                        cp_thread(i) = self%coeffX / m_thread
                        source_thread(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
                        dp_thread(i) = (source_thread(i) - dp_thread(i-1) * self%coeffX) / m_thread
                    end do

                    ! rightmost boundary properties
                    i = idx_end + 1
                    select case(boundary_end)
                    case (1)
                        source_thread(i) = self%solution(i, j)
                    case (2)
                        source_thread(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
                        m_thread = 2.0d0 * self%coeffX
                        source_thread(i) = (source_thread(i) - dp_thread(i-1) * m_thread)/&
                            (inv_centerCoeff - cp_thread(i-1)  * m_thread)
                    case (3)
                        source_thread(i) = source_thread(idx_start-1)
                    end select

                    ! backwards substitution
                    do i = idx_end, idx_start-1, -1
                        source_thread(i) = dp_thread(i)-cp_thread(i)*source_thread(i+1)
                    end do

                    self%solution(idx_start-1:idx_end+1,j) = source_thread(idx_start-1:idx_end+1)
                end do
            end do
            !$OMP end do

            !$OMP barrier

            ! white columns
            ! !$OMP do
            ! do k = self%inner_column_idx_white_start, self%number_inner_columns, 2
            !     i = self%start_column_indx + k - 1
            !     E_indx = i + 1
            !     W_indx = j - 1
            !     do p = 1, self%number_column_sections(k)
            !         do j = self%start_inner_indx_y(p,k), self%end_inner_indx_y(p,k)
            !             self%solution(i,j) = 500.0d0
            !         end do
            !     end do
            ! end do
            ! !$OMP end do nowait

            ! !$OMP sections
            ! !$OMP section
            ! ! left boundary
            ! if (self%number_left_column_sections /= 0) then
            !     do p = 1, self%number_left_column_sections
            !         if (self%left_column_boundary_type(p) == 2) then
            !             W_indx = 2
            !         else
            !             W_indx = self%N_x-1
            !         end if
            !         do j = self%start_left_column_indx(p), self%end_left_column_indx(p)
            !             self%solution(1,j) = 500.0d0
            !         end do
            !     end do
            ! end if
            ! !$OMP section
            ! ! right boundary
            ! if (self%number_right_column_sections /= 0) then
            !     do p = 1, self%number_right_column_sections
            !         if (self%right_column_boundary_type(p) == 2) then
            !             E_indx = self%N_x-1
            !         else
            !             E_indx = 2
            !         end if
            !         do j = self%start_right_column_indx(p), self%end_right_column_indx(p)
            !             self%solution(self%N_x,j) = 500.0d0
            !         end do
            !     end do
            ! end if
            ! !$OMP end sections
            
            !$OMP end parallel
        end do

    
    end subroutine smoothIterations_ZebraEven

    ! function smoothWithRes_ZebraEven(self) result(Res)
    !     ! Solve GS down to some tolerance
    !     class(ZebraSolverEven), intent(in out) :: self
    !     real(real64) :: Res, omega_inv, inv_centerCoeff
    !     real(real64) :: cp_x(self%N_x-1), dp_x(self%N_x-1), cp_y(self%N_y-1), dp_y(self%N_y-1), sourceX(self%N_x), sourceY(self%N_y), m
    !     integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
    !     omega_inv = 1.0d0 - self%omega
    !     inv_centerCoeff = 1.0d0/self%centerCoeff
    !     Res = 0.0d0
    !     !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, cp_y, cp_x, dp_y, dp_x, sourceX, sourceY, m) reduction(+:Res)

    !     ! horizontal sweeps left to right
    !     !$OMP do
    !     ! Sweep odd numbers forward
    !     do k = 1, self%numberRows, 2
    !         j = self%startRow + k - 1
    !         N_indx = self%vertIndx(1, k)
    !         S_indx = self%vertIndx(2, k)
    !         select case (self%leftBound)
    !         case (1)
    !             cp_x(1) = 0.0d0
    !             sourceX(1) = self%solution(1,j)
    !             dp_x(1) = sourceX(1)
    !         case (2)
    !             cp_x(1) = 2.0d0 * self%coeffX * self%centerCoeff
    !             sourceX(1) = self%sourceTerm(1, j) - self%coeffY * (self%solution(1, N_indx)  + self%solution(1, S_indx))
    !             dp_x(1) = (sourceX(1)) * self%centerCoeff
    !         case (3,4)
    !             cp_x(1) = 0.0d0
    !             E_indx = self%horzIndx(1,1)
    !             W_indx = self%horzIndx(2,1)
    !             sourceX(1) = (self%sourceTerm(1,j) - self%coeffY * (self%solution(1, N_indx) + self%solution(1, S_indx))- &
    !             self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))) * self%centerCoeff* self%omega + omega_inv * self%solution(1,j)
    !             dp_x(1) = sourceX(1)
    !         end select
    !         ! solve for vectors c-prime and d-prime
    !         do i = 2,self%N_x-1
    !             m = inv_centerCoeff - cp_x(i-1) * self%coeffX
    !             cp_x(i) = self%coeffX / m
    !             sourceX(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
    !             dp_x(i) = (sourceX(i) - dp_x(i-1) * self%coeffX) / m
    !         end do
    !         select case(self%rightBound)
    !         case (1)
    !             sourceX(self%N_x) = self%solution(self%N_x, j)
    !         case (2)
    !             sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - self%coeffY * ( self%solution(self%N_x, N_indx) + self%solution(self%N_x, S_indx))
    !             m = 2.0d0 * self%coeffX
    !             sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * m)/&
    !                 (inv_centerCoeff - cp_x(self%N_x-1)  * m)
    !         case (3)
    !             sourceX(self%N_x) = sourceX(1)
    !         end select
    !         do i = self%N_x-1, 1, -1
    !             sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
    !         end do
    !         self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
    !     end do
    !     !$OMP end do
        
    !     !$OMP do
    !     ! sweep even numbered columns
    !     do p = 2, self%numberColumns, 2
    !         i = self%startCol + p - 1
    !         E_indx = self%horzIndx(1, p)
    !         W_indx = self%horzIndx(2, p)
    !         select case (self%lowerBound)
    !         case (1)
    !             cp_y(1) = 0.0d0
    !             sourceY(1) = self%solution(i,1)
    !             dp_y(1) = sourceY(1)
    !         case (2)
    !             cp_y(1) = 2.0d0 * self%coeffY * self%centerCoeff
    !             sourceY(1) = self%sourceTerm(i, 1) - self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))
    !             dp_y(1) = (sourceY(1)) * self%centerCoeff
    !         case (3,4)
    !             cp_y(1) = 0.0d0
    !             N_indx = self%vertIndx(1,1)
    !             S_indx = self%vertIndx(2,1)
    !             sourceY(1) = (self%sourceTerm(i,1) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx)) - &
    !             self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
    !             dp_y(1) = sourceY(1)
    !         end select
    !         ! solve for vectors c-prime and d-prime
    !         do j = 2,self%N_y-1
    !             m = inv_centerCoeff - cp_y(j-1) * self%coeffY
    !             cp_y(j) = self%coeffY / m
    !             sourceY(j) = self%sourceTerm(i, j) - self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))
    !             dp_y(j) = (sourceY(j) - dp_y(j-1) * self%coeffY) / m
    !         end do
    !         select case(self%upperBound)
    !         case (1)
    !             sourceY(self%N_y) = self%solution(i, self%N_y)
    !         case (2)
    !             sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - self%coeffX * ( self%solution(E_indx, self%N_y) + self%solution(W_indx, self%N_y))
    !             m = 2.0d0 * self%coeffY
    !             sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * m)/&
    !                 (inv_centerCoeff - cp_y(self%N_y-1)  *m)
    !         case (3)
    !             sourceY(self%N_y) = sourceY(1)
    !         end select
    !         do j = self%N_y-1, 1, -1
    !             sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
    !         end do
    !         sourceY = sourceY*self%omega + omega_inv * self%solution(i,:)
    !         Res = Res + SUM((sourceY(self%startRow:self%endRow) - self%solution(i, self%startRow:self%endRow))**2)
    !         self%solution(i,:) = sourceY
    !     end do
    !     !$OMP end do

    !     ! Sweep horizontal even rows
    !     !$OMP do
    !     do k = 2, self%numberRows, 2
    !         j = self%startRow + k - 1
    !         N_indx = self%vertIndx(1, k)
    !         S_indx = self%vertIndx(2, k)
    !         select case (self%leftBound)
    !         case (1)
    !             cp_x(1) = 0.0d0
    !             sourceX(1) = self%solution(1,j)
    !             dp_x(1) = sourceX(1)
    !         case (2)
    !             cp_x(1) = 2.0d0 * self%coeffX * self%centerCoeff
    !             sourceX(1) = self%sourceTerm(1, j) - self%coeffY * (self%solution(1, N_indx)  + self%solution(1, S_indx))
    !             dp_x(1) = (sourceX(1)) * self%centerCoeff
    !         case (3,4)
    !             cp_x(1) = 0.0d0
    !             E_indx = self%horzIndx(1,1)
    !             W_indx = self%horzIndx(2,1)
    !             sourceX(1) = (self%sourceTerm(1,j) - self%coeffY * (self%solution(1, N_indx) + self%solution(1, S_indx))- &
    !             self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))) * self%centerCoeff* self%omega + omega_inv * self%solution(1,j)
    !             dp_x(1) = sourceX(1)
    !         end select
    !         ! solve for vectors c-prime and d-prime
    !         do i = 2,self%N_x-1
    !             m = inv_centerCoeff - cp_x(i-1) * self%coeffX
    !             cp_x(i) = self%coeffX / m
    !             sourceX(i) = self%sourceTerm(i,j) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx))
    !             dp_x(i) = (sourceX(i) - dp_x(i-1) * self%coeffX) / m
    !         end do
    !         select case(self%rightBound)
    !         case (1)
    !             sourceX(self%N_x) = self%solution(self%N_x, j)
    !         case (2)
    !             sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - self%coeffY * ( self%solution(self%N_x, N_indx) + self%solution(self%N_x, S_indx))
    !             m = 2.0d0 * self%coeffX
    !             sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * m)/&
    !                 (inv_centerCoeff - cp_x(self%N_x-1)  * m)
    !         case (3)
    !             sourceX(self%N_x) = sourceX(1)
    !         end select
    !         do i = self%N_x-1, 1, -1
    !             sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
    !         end do
    !         self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
    !     end do
    !     !$OMP end do

    !     ! vertical sweeps bottom to top
    !     !$OMP do
    !     ! Sweep odd numbers upwards
    !     do p = 1, self%numberColumns, 2
    !         i = self%startCol + p - 1
    !         E_indx = self%horzIndx(1, p)
    !         W_indx = self%horzIndx(2, p)
    !         select case (self%lowerBound)
    !         case (1)
    !             cp_y(1) = 0.0d0
    !             sourceY(1) = self%solution(i,1)
    !             dp_y(1) = sourceY(1)
    !         case (2)
    !             cp_y(1) = 2.0d0 * self%coeffY * self%centerCoeff
    !             sourceY(1) = self%sourceTerm(i, 1) - self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))
    !             dp_y(1) = (sourceY(1)) * self%centerCoeff
    !         case (3,4)
    !             cp_y(1) = 0.0d0
    !             N_indx = self%vertIndx(1,1)
    !             S_indx = self%vertIndx(2,1)
    !             sourceY(1) = (self%sourceTerm(i,1) - self%coeffY * (self%solution(i, N_indx) + self%solution(i, S_indx)) - &
    !             self%coeffX * (self%solution(E_indx, 1) + self%solution(W_indx, 1))) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
    !             dp_y(1) = sourceY(1)
    !         end select
    !         ! solve for vectors c-prime and d-prime
    !         do j = 2,self%N_y-1
    !             m = inv_centerCoeff - cp_y(j-1) * self%coeffY
    !             cp_y(j) = self%coeffY / m
    !             sourceY(j) = self%sourceTerm(i, j) - self%coeffX * (self%solution(E_indx, j) + self%solution(W_indx, j))
    !             dp_y(j) = (sourceY(j) - dp_y(j-1) * self%coeffY) / m
    !         end do
    !         select case(self%upperBound)
    !         case (1)
    !             sourceY(self%N_y) = self%solution(i, self%N_y)
    !         case (2)
    !             sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - self%coeffX * ( self%solution(E_indx, self%N_y) + self%solution(W_indx, self%N_y))
    !             m = 2.0d0 * self%coeffY
    !             sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * m)/&
    !                 (inv_centerCoeff - cp_y(self%N_y-1)  *m)
    !         case (3)
    !             sourceY(self%N_y) = sourceY(1)
    !         end select
    !         do j = self%N_y-1, 1, -1
    !             sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
    !         end do
    !         sourceY = sourceY*self%omega + omega_inv * self%solution(i,:)
    !         Res = Res + SUM((sourceY(self%startRow:self%endRow) - self%solution(i, self%startRow:self%endRow))**2)
    !         self%solution(i,:) = sourceY
    !     end do
    !     !$OMP end do
    !     !$OMP end parallel
        
    !     Res = SQRT(Res/ (self%numberColumns * self%numberRows))

    ! end function smoothWithRes_ZebraEven


end module mod_ZebraSolverEven