module mod_RedBlackSolverEven
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base_Even
    use mod_domain_uniform
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
    public :: RedBlackSolverEven

    type, extends(GS_Base_Even) :: RedBlackSolverEven
        ! store grid quantities
        integer(int32), allocatable :: start_red_indx(:,:), start_black_indx(:,:) ! Point to starting index for each section for black or red points
        ! integer, allocatable :: start_bottom_row_red_indx(:), start_bottom_row_black_indx(:)
        ! integer, allocatable :: start_top_row_red_indx(:), start_top_row_black_indx(:)
        ! integer, allocatable :: start_left_column_red_indx(:), start_left_column_black_indx(:)
        ! integer, allocatable :: start_right_column_red_indx(:), start_right_column_black_indx(:)
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal => constructPoissonOrthogonal_RedBlackEven
        procedure, public, pass(self) :: smoothIterations => smoothIterations_RedBlackEven
        procedure, public, pass(self) :: smoothWithRes => smoothWithRes_RedBlackEven
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type

    interface RedBlackSolverEven
        module procedure :: RedBlackSolverEven_constructor
    end interface RedBlackSolverEven

contains
    type(RedBlackSolverEven) function RedBlackSolverEven_constructor(omega, world, N_x, N_y) result(self)
        ! Construct object, set initial variables
        real(real64), intent(in) :: omega
        integer, intent(in) :: N_x, N_y
        type(domain_uniform), intent(in), target :: world
        call self%initialize_GS_Even(omega, world, N_x, N_y)
        call self%constructPoissonOrthogonal()
    end function RedBlackSolverEven_constructor

    subroutine constructPoissonOrthogonal_RedBlackEven(self)
        ! Construct orthogonal grid solver
        class(RedBlackSolverEven), intent(in out) :: self
        integer :: i, j, l, k, numberBound, x_indx_step, y_indx_step, min_indx
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
        
        ! allocate(self%start_bottom_row_red_indx(self%number_bottom_row_sections), self%start_bottom_row_black_indx(self%number_bottom_row_sections))
        ! do k = 1, self%number_bottom_row_sections
        !     i = self%start_bottom_row_indx(k)
        !     if (i == 1) then
        !         ! don't start at corner
        !         self%start_bottom_row_red_indx(k) = 3
        !         self%start_bottom_row_black_indx(k) = 2
        !     else if (MOD(i, 2) == 1) then
        !         self%start_bottom_row_red_indx(k) = i
        !         self%start_bottom_row_black_indx(k) = i+1
        !     else
        !         self%start_bottom_row_black_indx(k) = i
        !         self%start_bottom_row_red_indx(k) = i+1
        !     end if
        ! end do
        ! allocate(self%start_top_row_red_indx(self%number_top_row_sections), self%start_top_row_black_indx(self%number_top_row_sections))
        ! do k = 1, self%number_top_row_sections
        !     i = self%start_top_row_indx(k)
        !     if (i == 1) then
        !         ! don't start at corner
        !         self%start_top_row_red_indx(k) = 3
        !         self%start_top_row_black_indx(k) = 2
        !     else if (MOD(i, 2) == 1) then
        !         self%start_top_row_red_indx(k) = i
        !         self%start_top_row_black_indx(k) = i+1
        !     else
        !         self%start_top_row_black_indx(k) = i
        !         self%start_top_row_red_indx(k) = i+1
        !     end if
        ! end do
        ! allocate(self%start_left_column_red_indx(self%number_left_column_sections), self%start_left_column_black_indx(self%number_left_column_sections))
        ! do k = 1, self%number_left_column_sections
        !     j = self%start_left_column_indx(k)
        !     if (j == 1) then
        !         ! don't start at corner
        !         self%start_left_column_red_indx(k) = 3
        !         self%start_left_column_black_indx(k) = 2
        !     else if (MOD(j, 2) == 1) then
        !         self%start_left_column_red_indx(k) = j
        !         self%start_left_column_black_indx(k) = j+1
        !     else
        !         self%start_left_column_black_indx(k) = j
        !         self%start_left_column_red_indx(k) = j+1
        !     end if
        ! end do
        ! allocate(self%start_right_column_red_indx(self%number_right_column_sections), self%start_right_column_black_indx(self%number_right_column_sections))
        ! do k = 1, self%number_right_column_sections
        !     j = self%start_right_column_indx(k)
        !     if (j == 1) then
        !         ! don't start at corner
        !         self%start_right_column_red_indx(k) = 3
        !         self%start_right_column_black_indx(k) = 2
        !     else if (MOD(j, 2) == 1) then
        !         self%start_right_column_red_indx(k) = j
        !         self%start_right_column_black_indx(k) = j+1
        !     else
        !         self%start_right_column_black_indx(k) = j
        !         self%start_right_column_red_indx(k) = j+1
        !     end if
        ! end do
    end subroutine constructPoissonOrthogonal_RedBlackEven

    subroutine smoothIterations_RedBlackEven(self, iterNum)
        ! Solve GS down to some tolerance
        class(RedBlackSolverEven), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        real(real64) :: oldSol, omega_inv
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p,iter
        omega_inv = 1.0d0 - self%omega
        ! p indexs in horizontal maps to i, k in vertical maps to j
        do iter = 1, iterNum 


            !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol)

            !first do corners which need to be neumann for change
            !$OMP sections
            !$OMP section
            if (self%world%boundary_conditions(1,1) == 2) then
                ! lower left corner
                self%solution(1,1) = (self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * self%coeffY - &
                    2.0d0 * self%solution(2, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,1)
            end if
            !$OMP section
            if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
                ! lower right corner
                self%solution(self%N_x, 1) = (self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x, 1)
            end if
            !$OMP section
            if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
                ! upper left corner
                self%solution(1, self%N_y) = (self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * self%coeffY - &
                    2.0d0 * self%solution(2, self%N_y) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1, self%N_y)
            end if
            !$OMP section
            if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
                ! upper right corner
                self%solution(self%N_x, self%N_y) = (self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, self%N_y) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x, self%N_y)
            end if
            !$OMP end sections nowait

            !Go through all red points first
            !$OMP do
            ! Sweep rows that start with red
            do k = 1, self%number_inner_rows
                j = self%start_row_indx + k - 1
                N_indx = j + 1
                S_indx = j - 1
                do p = 1, self%number_row_sections(k)
                    do i = self%start_red_indx(p, k), self%end_inner_indx_x(p,k), 2     
                        E_indx = i+1
                        W_indx = i-1
                        oldSol = self%solution(i,j)
                        self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                            (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                    end do
                end do
            end do
            !$OMP end do nowait

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
                do i = k, self%end_bottom_row_indx(p), 2
                    self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, S_indx)) * self%coeffY - &
                            (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
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
                do i = k, self%end_top_row_indx(p), 2
                    self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - (self%solution(i, N_indx) + self%solution(i, self%N_y-1)) * self%coeffY - &
                            (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,self%N_y)
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
                do j = k, self%end_left_column_indx(p), 2
                    self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                            (self%solution(2, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,j)
                end do
            end do
            !$OMP end do

            !$OMP barrier

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
                do j = k, self%end_right_column_indx(p), 2
                    self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j+1) + self%solution(self%N_x, j-1)) * self%coeffY - &
                            (self%solution(self%N_x-1, j) + self%solution(E_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x,j)
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
                do p = 1, self%number_row_sections(k)
                    do i = self%start_black_indx(p, k), self%end_inner_indx_x(p,k), 2     
                        E_indx = i+1
                        W_indx = i-1
                        oldSol = self%solution(i,j)
                        self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                            (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
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
                do i = k, self%end_bottom_row_indx(p), 2
                    self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, S_indx)) * self%coeffY - &
                            (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
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
                do i = k, self%end_top_row_indx(p), 2
                    self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - (self%solution(i, N_indx) + self%solution(i, self%N_y-1)) * self%coeffY - &
                            (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,self%N_y)
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
                do j = k, self%end_left_column_indx(p), 2
                    self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                            (self%solution(2, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,j)
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
                do j = k, self%end_right_column_indx(p), 2
                    self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j+1) + self%solution(self%N_x, j-1)) * self%coeffY - &
                            (self%solution(self%N_x-1, j) + self%solution(E_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x,j)
                end do
            end do
            !$OMP end do
            
            !$OMP end parallel
        end do

    
    end subroutine smoothIterations_RedBlackEven

    function smoothWithRes_RedBlackEven(self) result(Res)
        ! Solve GS down to some tolerance
        class(RedBlackSolverEven), intent(in out) :: self
        real(real64) :: oldSol, Res, omega_inv
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
        omega_inv = 1.0d0 - self%omega

        Res = 0.0d0
        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol) reduction(+:Res)

        !first do corners which need to be neumann for change
        !$OMP sections
        !$OMP section
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            oldSol = self%solution(1,1)
            self%solution(1,1) = (self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * self%coeffY - &
                2.0d0 * self%solution(2, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            Res = Res + (self%solution(1,1) - oldSol)**2
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            oldSol = self%solution(self%N_x, 1)
            self%solution(self%N_x, 1) = (self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            Res = Res + (self%solution(self%N_x,1) - oldSol)**2
        end if
        !$OMP section
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            oldSol = self%solution(1, self%N_y)
            self%solution(1, self%N_y) = (self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(2, self%N_y) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            Res = Res + (self%solution(1, self%N_y) - oldSol)**2
        end if
        !$OMP section
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            oldSol = self%solution(self%N_x, self%N_y)
            self%solution(self%N_x, self%N_y) = (self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, self%N_y) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            Res = Res + (self%solution(self%N_x, self%N_y) - oldSol)**2
        end if
        !$OMP end sections nowait

        !Go through all red points first
        !$OMP do
        ! Sweep rows that start with red
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            do p = 1, self%number_row_sections(k)
                do i = self%start_red_indx(p, k), self%end_inner_indx_x(p,k), 2     
                    E_indx = i+1
                    W_indx = i-1
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                    Res = Res + (self%solution(i,j) - oldSol)**2
                end do
            end do
        end do
        !$OMP end do nowait

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
            do i = k, self%end_bottom_row_indx(p), 2
                oldSol = self%solution(i,1)
                self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
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
            do i = k, self%end_top_row_indx(p), 2
                oldSol = self%solution(i,self%N_y)
                self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - (self%solution(i, N_indx) + self%solution(i, self%N_y-1)) * self%coeffY - &
                        (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,self%N_y)
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
            do j = k, self%end_left_column_indx(p), 2
                oldSol = self%solution(1,j)
                self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                        (self%solution(2, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,j)
                Res = Res + (self%solution(1,j) - oldSol)**2
            end do
        end do
        !$OMP end do

        !$OMP barrier
        
        ! Go through all black points
        !$OMP do
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            do p = 1, self%number_row_sections(k)
                do i = self%start_black_indx(p, k), self%end_inner_indx_x(p,k), 2     
                    E_indx = i+1
                    W_indx = i-1
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
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
            do i = k, self%end_bottom_row_indx(p), 2
                oldSol = self%solution(i,1)
                self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,1)
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
            do i = k, self%end_top_row_indx(p), 2
                oldSol = self%solution(i,self%N_y)
                self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - (self%solution(i, N_indx) + self%solution(i, self%N_y-1)) * self%coeffY - &
                        (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(i,self%N_y)
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
            do j = k, self%end_left_column_indx(p), 2
                oldSol = self%solution(1,j)
                self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                        (self%solution(2, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,j)
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
            do j = k, self%end_right_column_indx(p), 2
                oldSol = self%solution(self%N_x,j)
                self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j+1) + self%solution(self%N_x, j-1)) * self%coeffY - &
                        (self%solution(self%N_x-1, j) + self%solution(E_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x,j)
                Res = Res + (self%solution(self%N_x,j) - oldSol)**2
            end do
        end do
        !$OMP end do
        
        !$OMP end parallel
        
        Res = SQRT(Res/ self%number_solve_nodes)

    end function smoothWithRes_RedBlackEven


end module mod_RedBlackSolverEven