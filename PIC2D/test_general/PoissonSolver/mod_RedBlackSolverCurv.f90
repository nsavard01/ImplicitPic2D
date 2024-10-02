module mod_RedBlackSolverCurv
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base_Curv
    use mod_domain_curv
    use omp_lib
    implicit none

    ! Class for a red black gauss siedel smoother
    ! Since will be for stage of multigrid, and to work with red black points, N_x, N_y must be odd
    ! Example order for black (b) and red (r) points for 5x5 matrix is:
    ! 
    ! r - b - r - b - r
    ! b - r - b - r - b
    ! r - b - r - b - r
    ! b - r - b - r - b
    ! r - b - r - b - r
    ! Must be orthogonal grid


    private
    public :: RedBlackSolverCurv

    type, extends(GS_Base_Curv) :: RedBlackSolverCurv
        ! store grid quantities
        integer(int32), allocatable :: start_red_indx(:,:), start_black_indx(:,:) ! Point to starting index for each section for black or red points
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal => constructPoissonOrthogonal_RedBlackCurv
        procedure, public, pass(self) :: smoothIterations => smoothIterations_RedBlackCurv
        procedure, public, pass(self) :: smoothWithRes => smoothWithRes_RedBlackCurv
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type

    interface RedBlackSolverCurv
        module procedure :: RedBlackSolverCurv_constructor
    end interface RedBlackSolverCurv

contains

    type(RedBlackSolverCurv) function RedBlackSolverCurv_constructor(omega, world, N_x, N_y) result(self)
        ! Construct object, set initial variables
        real(real64), intent(in) :: omega
        integer, intent(in) :: N_x, N_y
        type(domain_curv), intent(in), target :: world
        call self%initialize_GS_orthogonal_base(omega, world, N_x, N_y)
        call self%initialize_GS_Curv(world)
        call self%constructPoissonOrthogonal()
    end function RedBlackSolverCurv_constructor

    subroutine constructPoissonOrthogonal_RedBlackCurv(self)
        ! Construct orthogonal grid solver
        class(RedBlackSolverCurv), intent(in out) :: self
        integer :: i, j, l, k
        allocate(self%start_black_indx(self%max_number_row_sections, self%number_inner_rows), &
            self%start_red_indx(self%max_number_row_sections, self%number_inner_rows)) !
        ! Get starting point for red or black points
        !$OMP parallel private(i, j, l, k)
        !$OMP do
        do l = 1, self%number_inner_rows
            j = l - 1 + self%start_row_indx
            do k = 1, self%number_row_sections(l)
                i = self%start_inner_indx_x(k, l)
                if ((MOD(i,2) == 0) .eq. (MOD(j,2) == 0)) then
                    !Red point start
                    self%start_red_indx(k,l) = i
                    self%start_black_indx(k,l) = i + 1
                else
                    !Black point
                    self%start_red_indx(k,l) = i+1
                    self%start_black_indx(k,l) = i
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
       
        
    end subroutine constructPoissonOrthogonal_RedBlackCurv


    ! subroutine matMult(self, x, y)
    !     ! Use gauss-seidel to calculate A*x, store solution into y
    !     class(RedBlackSolver), intent(in out) :: self
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
    !     class(RedBlackSolver), intent(in out) :: self
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

    subroutine smoothIterations_RedBlackCurv(self, iterNum)
        ! Solve GS down to some tolerance
        class(RedBlackSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        ! integer(int32), intent(in) :: iterNum
        real(real64) :: oldSol, C_N, C_E, C_S, C_W, C_O
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p,iter
        ! p indexs in horizontal maps to i, k in vertical maps to j
        do iter = 1, iterNum 


            !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol, &
            !$OMP& C_N, C_E, C_S, C_W, C_O)


            !Go through all red points first
            !$OMP do
            ! Sweep rows that start with red
            do k = 1, self%number_inner_rows
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                C_N = self%inner_node_coeff_north(j-1)
                C_S = self%inner_node_coeff_south(j-1)
                do p = 1, self%number_row_sections(k)
                    do i = self%start_red_indx(p, k), self%end_inner_indx_x(p,k), 2     
                        E_indx = i+1
                        W_indx = i-1
                        C_E = self%inner_node_coeff_east(i-1)
                        C_W = self%inner_node_coeff_west(i-1)
                        C_O = self%centerCoeff(i, j)
                        oldSol = self%solution(i,j)
                        self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx)*C_S - &
                            self%solution(E_indx, j)*C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * oldSol
                    end do
                end do
            end do
            !$OMP end do nowait

            ! do corners which need to be neumann for change
            !$OMP sections
            !$OMP section
            if (self%world%boundary_conditions(1,1) == 2) then
                ! lower left corner
                C_N = self%bottom_row_coeff_north(1)
                C_E = self%left_column_coeff_east(1)
                C_O = self%centerCoeff(1,1)
                self%solution(1,1) = (self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * C_N - &
                    2.0d0 * self%solution(2, 1) * C_E) * C_O * self%omega + self%inv_omega * self%solution(1,1)
            end if
            !$OMP section
            if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
                ! lower right corner
                C_N = self%bottom_row_coeff_north(self%number_bottom_row_sections)
                C_W = self%right_column_coeff_west(1)
                C_O = self%centerCoeff(self%N_x, 1)
                self%solution(self%N_x, 1) = (self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * C_N - &
                    2.0d0 * self%solution(self%N_x-1, 1) * C_W) * C_O * self%omega + self%inv_omega * self%solution(self%N_x, 1)
            end if
            !$OMP section
            if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
                ! upper left corner
                C_S = self%top_row_coeff_south(1)
                C_E = self%left_column_coeff_east(self%number_left_column_sections)
                C_O = self%centerCoeff(1, self%N_y)
                self%solution(1, self%N_y) = (self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * C_S - &
                    2.0d0 * self%solution(2, self%N_y) * C_E) * C_O * self%omega + self%inv_omega * self%solution(1, self%N_y)
            end if
            !$OMP section
            if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
                ! upper right corner
                C_S = self%top_row_coeff_south(self%number_top_row_sections)
                C_W = self%right_column_coeff_west(self%number_right_column_sections)
                C_O = self%centerCoeff(self%N_x, self%N_y)
                self%solution(self%N_x, self%N_y) = (self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * C_S - &
                    2.0d0 * self%solution(self%N_x-1, self%N_y) * C_W) * C_O * self%omega + self%inv_omega * self%solution(self%N_x, self%N_y)
            end if
            !$OMP end sections nowait

            ! !Now go through red boundaries
            !lower boundary
            !$OMP do
            do p = 1, self%number_bottom_row_sections
                if (self%bottom_row_boundary_type(p) == 2) then
                    S_indx = 2
                else
                    S_indx = self%N_y-1
                end if
                if (MOD(self%start_bottom_row_indx(p),2) == 1) then
                    k = self%start_bottom_row_indx(p)
                else
                    k = self%start_bottom_row_indx(p) + 1
                end if
                C_S = self%bottom_row_coeff_south(p)
                C_N = self%bottom_row_coeff_north(p)
                do i = k, self%end_bottom_row_indx(p), 2
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = self%centerCoeff(i,1)
                    self%solution(i,1) = (self%sourceTerm(i,1) - self%solution(i, 2) * C_N - self%solution(i, S_indx) * C_S - &
                            self%solution(i-1, 1) * C_W - self%solution(i+1, 1) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,1)
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
                if (MOD(self%start_top_row_indx(p),2) == 1) then
                    k = self%start_top_row_indx(p)
                else
                    k = self%start_top_row_indx(p) + 1
                end if
                C_S = self%top_row_coeff_south(p)
                C_N = self%top_row_coeff_north(p)
                do i = k, self%end_top_row_indx(p), 2
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = self%centerCoeff(i,self%N_y)
                    self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - self%solution(i, N_indx) * C_N - self%solution(i, self%N_y-1) * C_S - &
                            self%solution(i-1, self%N_y) * C_W - self%solution(i+1, self%N_y) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,self%N_y)
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
                if (MOD(self%start_left_column_indx(p),2) == 1) then
                    k = self%start_left_column_indx(p)
                else
                    k = self%start_left_column_indx(p) + 1
                end if
                C_E = self%left_column_coeff_east(p)
                C_W = self%left_column_coeff_west(p)
                do j = k, self%end_left_column_indx(p), 2
                    C_S = self%inner_node_coeff_south(j-1)
                    C_N = self%inner_node_coeff_north(j-1)
                    C_O = self%centerCoeff(1,j)
                    self%solution(1,j) = (self%sourceTerm(1,j) - self%solution(1, j+1) * C_N - self%solution(1, j-1) * C_S - &
                            self%solution(2, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * self%solution(1,j)
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
                if (MOD(self%start_right_column_indx(p),2) == 1) then
                    k = self%start_right_column_indx(p)
                else
                    k = self%start_right_column_indx(p) + 1
                end if
                C_E = self%right_column_coeff_east(p)
                C_W = self%right_column_coeff_west(p)
                do j = k, self%end_right_column_indx(p), 2
                    C_S = self%inner_node_coeff_south(j-1)
                    C_N = self%inner_node_coeff_north(j-1)
                    C_O = self%centerCoeff(self%N_x,j)
                    self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - self%solution(self%N_x, j+1) * C_N - self%solution(self%N_x, j-1) * C_S - &
                            self%solution(self%N_x-1, j) * C_W - self%solution(E_indx, j) * C_E) * C_O * self%omega + self%inv_omega * self%solution(self%N_x,j)
                end do
            end do
            !$OMP end do

            !$OMP barrier
            ! ------------------------------------------------- Now black points ------------------------------------------------------------------

            ! !Go through all black points
            !$OMP do
            do k = 1, self%number_inner_rows
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                C_N = self%inner_node_coeff_north(j-1)
                C_S = self%inner_node_coeff_south(j-1)
                do p = 1, self%number_row_sections(k)
                    do i = self%start_black_indx(p, k), self%end_inner_indx_x(p,k), 2     
                        E_indx = i+1
                        W_indx = i-1
                        C_E = self%inner_node_coeff_east(i-1)
                        C_W = self%inner_node_coeff_west(i-1)
                        C_O = self%centerCoeff(i, j)
                        oldSol = self%solution(i,j)
                        self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx)*C_S - &
                            self%solution(E_indx, j)*C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * oldSol
                    end do
                end do
            end do
            !$OMP end do nowait

            ! !lower boundary
            !$OMP do
            do p = 1, self%number_bottom_row_sections
                if (self%bottom_row_boundary_type(p) == 2) then
                    S_indx = 2
                else
                    S_indx = self%N_y-1
                end if
                if (MOD(self%start_bottom_row_indx(p),2) == 0) then
                    k = self%start_bottom_row_indx(p)
                else
                    k = self%start_bottom_row_indx(p) + 1
                end if
                C_S = self%bottom_row_coeff_south(p)
                C_N = self%bottom_row_coeff_north(p)
                do i = k, self%end_bottom_row_indx(p), 2
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = self%centerCoeff(i,1)
                    self%solution(i,1) = (self%sourceTerm(i,1) - self%solution(i, 2) * C_N - self%solution(i, S_indx) * C_S - &
                            self%solution(i-1, 1) * C_W - self%solution(i+1, 1) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,1)
                end do
            end do
            !$OMP end do nowait

            ! ! ! upper boundary
            !$OMP do
            do p = 1, self%number_top_row_sections
                if (self%top_row_boundary_type(p) == 2) then
                    N_indx = self%N_y-1
                else
                    N_indx = 2
                end if
                if (MOD(self%start_top_row_indx(p),2) == 0) then
                    k = self%start_top_row_indx(p)
                else
                    k = self%start_top_row_indx(p) + 1
                end if
                C_S = self%top_row_coeff_south(p)
                C_N = self%top_row_coeff_north(p)
                do i = k, self%end_top_row_indx(p), 2
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = self%centerCoeff(i,self%N_y)
                    self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - self%solution(i, N_indx) * C_N - self%solution(i, self%N_y-1) * C_S - &
                            self%solution(i-1, self%N_y) * C_W - self%solution(i+1, self%N_y) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,self%N_y)
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
                if (MOD(self%start_left_column_indx(p),2) == 0) then
                    k = self%start_left_column_indx(p)
                else
                    k = self%start_left_column_indx(p) + 1
                end if
                C_E = self%left_column_coeff_east(p)
                C_W = self%left_column_coeff_west(p)
                do j = k, self%end_left_column_indx(p), 2
                    C_S = self%inner_node_coeff_south(j-1)
                    C_N = self%inner_node_coeff_north(j-1)
                    C_O = self%centerCoeff(1,j)
                    self%solution(1,j) = (self%sourceTerm(1,j) - self%solution(1, j+1) * C_N - self%solution(1, j-1) * C_S - &
                            self%solution(2, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * self%solution(1,j)
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
                if (MOD(self%start_right_column_indx(p),2) == 0) then
                    k = self%start_right_column_indx(p)
                else
                    k = self%start_right_column_indx(p) + 1
                end if
                C_E = self%right_column_coeff_east(p)
                C_W = self%right_column_coeff_west(p)
                do j = k, self%end_right_column_indx(p), 2
                    C_S = self%inner_node_coeff_south(j-1)
                    C_N = self%inner_node_coeff_north(j-1)
                    C_O = self%centerCoeff(self%N_x,j)
                    self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - self%solution(self%N_x, j+1) * C_N - self%solution(self%N_x, j-1) * C_S - &
                            self%solution(self%N_x-1, j) * C_W - self%solution(E_indx, j) * C_E) * C_O * self%omega + self%inv_omega * self%solution(self%N_x,j)
                end do
            end do
            !$OMP end do
            !$OMP end parallel
        end do
    end subroutine smoothIterations_RedBlackCurv

    function smoothWithRes_RedBlackCurv(self) result(Res)
        ! Solve GS down to some tolerance
        class(RedBlackSolverCurv), intent(in out) :: self
        real(real64) :: oldSol, C_N, C_E, C_S, C_W, C_O, Res
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p,iter
        Res = 0.0d0
       
        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol, &
            !$OMP& C_N, C_E, C_S, C_W, C_O) reduction(+:Res)


        !Go through all red points first
        !$OMP do
        ! Sweep rows that start with red
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            C_N = self%inner_node_coeff_north(j-1)
            C_S = self%inner_node_coeff_south(j-1)
            do p = 1, self%number_row_sections(k)
                do i = self%start_red_indx(p, k), self%end_inner_indx_x(p,k), 2     
                    E_indx = i+1
                    W_indx = i-1
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = self%centerCoeff(i, j)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx)*C_S - &
                        self%solution(E_indx, j)*C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * oldSol
                    Res = Res + (self%solution(i,j) - oldSol)**2
                end do
            end do
        end do
        !$OMP end do nowait

        ! do corners which need to be neumann for change
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            C_N = self%bottom_row_coeff_north(1)
            C_E = self%left_column_coeff_east(1)
            C_O = self%centerCoeff(1,1)
            oldSol = self%solution(1,1)
            self%solution(1,1) = (self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * C_N - &
                2.0d0 * self%solution(2, 1) * C_E) * C_O * self%omega + self%inv_omega * self%solution(1,1)
            Res = Res + (self%solution(1,1) - oldSol)**2
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            C_N = self%bottom_row_coeff_north(self%number_bottom_row_sections)
            C_W = self%right_column_coeff_west(1)
            C_O = self%centerCoeff(self%N_x,1)
            oldSol = self%solution(self%N_x, 1)
            self%solution(self%N_x, 1) = (self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * C_N - &
                2.0d0 * self%solution(self%N_x-1, 1) * C_W) * C_O * self%omega + self%inv_omega * self%solution(self%N_x, 1)
            Res = Res + (self%solution(self%N_x,1) - oldSol)**2
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            C_S = self%top_row_coeff_south(1)
            C_E = self%left_column_coeff_east(self%number_left_column_sections)
            C_O = self%centerCoeff(1,self%N_y)
            oldSol = self%solution(1, self%N_y)
            self%solution(1, self%N_y) = (self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * C_S - &
                2.0d0 * self%solution(2, self%N_y) * C_E) * C_O * self%omega + self%inv_omega * self%solution(1, self%N_y)
            Res = Res + (self%solution(1, self%N_y) - oldSol)**2
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            C_S = self%top_row_coeff_south(self%number_top_row_sections)
            C_W = self%right_column_coeff_west(self%number_right_column_sections)
            C_O = self%centerCoeff(self%N_x, self%N_y)
            oldSol = self%solution(self%N_x, self%N_y)
            self%solution(self%N_x, self%N_y) = (self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * C_S - &
                2.0d0 * self%solution(self%N_x-1, self%N_y) * C_W) * C_O * self%omega + self%inv_omega * self%solution(self%N_x, self%N_y)
            Res = Res + (self%solution(self%N_x, self%N_y) - oldSol)**2
        end if
        !$OMP end sections nowait

        ! !Now go through red boundaries
        !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            if (self%bottom_row_boundary_type(p) == 2) then
                S_indx = 2
            else
                S_indx = self%N_y-1
            end if
            if (MOD(self%start_bottom_row_indx(p),2) == 1) then
                k = self%start_bottom_row_indx(p)
            else
                k = self%start_bottom_row_indx(p) + 1
            end if
            C_S = self%bottom_row_coeff_south(p)
            C_N = self%bottom_row_coeff_north(p)
            do i = k, self%end_bottom_row_indx(p), 2
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                C_O = self%centerCoeff(i,1)
                oldSol = self%solution(i,1)
                self%solution(i,1) = (self%sourceTerm(i,1) - self%solution(i, 2) * C_N - self%solution(i, S_indx) * C_S - &
                        self%solution(i-1, 1) * C_W - self%solution(i+1, 1) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,1)
                Res = Res + (self%solution(i,1) - oldSol)**2
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
            if (MOD(self%start_top_row_indx(p),2) == 1) then
                k = self%start_top_row_indx(p)
            else
                k = self%start_top_row_indx(p) + 1
            end if
            C_S = self%top_row_coeff_south(p)
            C_N = self%top_row_coeff_north(p)
            do i = k, self%end_top_row_indx(p), 2
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                C_O = self%centerCoeff(i,self%N_y)
                oldSol = self%solution(i,self%N_y)
                self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - self%solution(i, N_indx) * C_N - self%solution(i, self%N_y-1) * C_S - &
                        self%solution(i-1, self%N_y) * C_W - self%solution(i+1, self%N_y) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,self%N_y)
                Res = Res + (self%solution(i,self%N_y) - oldSol)**2
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
            if (MOD(self%start_left_column_indx(p),2) == 1) then
                k = self%start_left_column_indx(p)
            else
                k = self%start_left_column_indx(p) + 1
            end if
            C_E = self%left_column_coeff_east(p)
            C_W = self%left_column_coeff_west(p)
            do j = k, self%end_left_column_indx(p), 2
                C_S = self%inner_node_coeff_south(j-1)
                C_N = self%inner_node_coeff_north(j-1)
                C_O = self%centerCoeff(1,j)
                oldSol = self%solution(1,j)
                self%solution(1,j) = (self%sourceTerm(1,j) - self%solution(1, j+1) * C_N - self%solution(1, j-1) * C_S - &
                        self%solution(2, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * self%solution(1,j)
                Res = Res + (self%solution(1,j) - oldSol)**2
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
            if (MOD(self%start_right_column_indx(p),2) == 1) then
                k = self%start_right_column_indx(p)
            else
                k = self%start_right_column_indx(p) + 1
            end if
            C_E = self%right_column_coeff_east(p)
            C_W = self%right_column_coeff_west(p)
            do j = k, self%end_right_column_indx(p), 2
                C_S = self%inner_node_coeff_south(j-1)
                C_N = self%inner_node_coeff_north(j-1)
                C_O = self%centerCoeff(self%N_x,j)
                oldSol = self%solution(self%N_x, j)
                self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - self%solution(self%N_x, j+1) * C_N - self%solution(self%N_x, j-1) * C_S - &
                        self%solution(self%N_x-1, j) * C_W - self%solution(E_indx, j) * C_E) * C_O * self%omega + self%inv_omega * self%solution(self%N_x,j)
                Res = Res + (self%solution(self%N_x,j) - oldSol)**2
            end do
        end do
        !$OMP end do

        !$OMP barrier
        ! ------------------------------------------------- Now black points ------------------------------------------------------------------

        ! !Go through all black points
        !$OMP do
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            C_N = self%inner_node_coeff_north(j-1)
            C_S = self%inner_node_coeff_south(j-1)
            do p = 1, self%number_row_sections(k)
                do i = self%start_black_indx(p, k), self%end_inner_indx_x(p,k), 2     
                    E_indx = i+1
                    W_indx = i-1
                    C_E = self%inner_node_coeff_east(i-1)
                    C_W = self%inner_node_coeff_west(i-1)
                    C_O = self%centerCoeff(i, j)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx)*C_S - &
                        self%solution(E_indx, j)*C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * oldSol
                    Res = Res + (self%solution(i,j) - oldSol)**2
                end do
            end do
        end do
        !$OMP end do nowait

        ! !lower boundary
        !$OMP do
        do p = 1, self%number_bottom_row_sections
            if (self%bottom_row_boundary_type(p) == 2) then
                S_indx = 2
            else
                S_indx = self%N_y-1
            end if
            if (MOD(self%start_bottom_row_indx(p),2) == 0) then
                k = self%start_bottom_row_indx(p)
            else
                k = self%start_bottom_row_indx(p) + 1
            end if
            C_S = self%bottom_row_coeff_south(p)
            C_N = self%bottom_row_coeff_north(p)
            do i = k, self%end_bottom_row_indx(p), 2
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                C_O = self%centerCoeff(i,1)
                oldSol = self%solution(i,1)
                self%solution(i,1) = (self%sourceTerm(i,1) - self%solution(i, 2) * C_N - self%solution(i, S_indx) * C_S - &
                        self%solution(i-1, 1) * C_W - self%solution(i+1, 1) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,1)
                Res = Res + (self%solution(i,1) - oldSol)**2
            end do
        end do
        !$OMP end do nowait

        ! ! ! upper boundary
        !$OMP do
        do p = 1, self%number_top_row_sections
            if (self%top_row_boundary_type(p) == 2) then
                N_indx = self%N_y-1
            else
                N_indx = 2
            end if
            if (MOD(self%start_top_row_indx(p),2) == 0) then
                k = self%start_top_row_indx(p)
            else
                k = self%start_top_row_indx(p) + 1
            end if
            C_S = self%top_row_coeff_south(p)
            C_N = self%top_row_coeff_north(p)
            do i = k, self%end_top_row_indx(p), 2
                C_E = self%inner_node_coeff_east(i-1)
                C_W = self%inner_node_coeff_west(i-1)
                C_O = self%centerCoeff(i,self%N_y)
                oldSol = self%solution(i,self%N_y)
                self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - self%solution(i, N_indx) * C_N - self%solution(i, self%N_y-1) * C_S - &
                        self%solution(i-1, self%N_y) * C_W - self%solution(i+1, self%N_y) * C_E) * C_O * self%omega + self%inv_omega * self%solution(i,self%N_y)
                Res = Res + (self%solution(i,self%N_y) - oldSol)**2
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
            if (MOD(self%start_left_column_indx(p),2) == 0) then
                k = self%start_left_column_indx(p)
            else
                k = self%start_left_column_indx(p) + 1
            end if
            C_E = self%left_column_coeff_east(p)
            C_W = self%left_column_coeff_west(p)
            do j = k, self%end_left_column_indx(p), 2
                C_S = self%inner_node_coeff_south(j-1)
                C_N = self%inner_node_coeff_north(j-1)
                C_O = self%centerCoeff(1,j)
                oldSol = self%solution(1,j)
                self%solution(1,j) = (self%sourceTerm(1,j) - self%solution(1, j+1) * C_N - self%solution(1, j-1) * C_S - &
                        self%solution(2, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + self%inv_omega * self%solution(1,j)
                Res = Res + (self%solution(1,j) - oldSol)**2
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
            if (MOD(self%start_right_column_indx(p),2) == 0) then
                k = self%start_right_column_indx(p)
            else
                k = self%start_right_column_indx(p) + 1
            end if
            C_E = self%right_column_coeff_east(p)
            C_W = self%right_column_coeff_west(p)
            do j = k, self%end_right_column_indx(p), 2
                C_S = self%inner_node_coeff_south(j-1)
                C_N = self%inner_node_coeff_north(j-1)
                C_O = self%centerCoeff(self%N_x,j)
                oldSol = self%solution(self%N_x,j)
                self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - self%solution(self%N_x, j+1) * C_N - self%solution(self%N_x, j-1) * C_S - &
                        self%solution(self%N_x-1, j) * C_W - self%solution(E_indx, j) * C_E) * C_O * self%omega + self%inv_omega * self%solution(self%N_x,j)
                Res = Res + (self%solution(self%N_x,j) - oldSol)**2
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        Res = SQRT(Res/ self%number_solve_nodes)

    end function smoothWithRes_RedBlackCurv


    


end module mod_RedBlackSolverCurv