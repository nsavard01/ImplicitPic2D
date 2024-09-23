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
    end subroutine constructPoissonOrthogonal_RedBlackEven

    subroutine smoothIterations_RedBlackEven(self, iterNum)
        ! Solve GS down to some tolerance
        class(RedBlackSolverEven), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        real(real64) :: oldSol, omega_inv
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, iter, i_finest, j_finest
        omega_inv = 1.0d0 - self%omega
        ! p indexs in horizontal maps to i, k in vertical maps to j
        do iter = 1, iterNum 

            ! first do corners (due to inconvenience in vectorizing) which have to be neumann for changes
            ! these are all red points
            if (self%world%boundary_conditions(1,1) == 2) then
                ! lower left corner
                self%solution(1,1) = (self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * self%coeffY - &
                    2.0d0 * self%solution(2, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,1)
            end if
            if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
                ! lower right corner
                self%solution(self%N_x, 1) = (self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x, 1)
            end if
            if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
                ! upper left corner
                self%solution(1, self%N_y) = (self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * self%coeffY - &
                    2.0d0 * self%solution(2, self%N_y-1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1, self%N_y)
            end if
            if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
                ! upper right corner
                self%solution(self%N_x, self%N_y) = (self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, self%N_y) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x, self%N_y)
            end if

            !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol, i_finest, j_finest)

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

            ! Now go through red boundaries
            ! Red Upper/lower rows
            !$OMP do 
            do i = 3, self%N_x-2, 2
                i_finest = i * self%x_indx_step - self%x_indx_step + 1

                ! Lower row
                oldSol = self%solution(i,1)
                if (self%world%boundary_conditions(i_finest, 1) == 2) then
                    self%solution(i,1) = (self%sourceTerm(i,1) - 2.0d0 * self%solution(i, 2) * self%coeffY - &
                        (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                else if (self%world%boundary_conditions(i_finest, 1) == 3) then
                    self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, self%N_y-1)) * self%coeffY - &
                        (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                    self%solution(i, self%N_y) = self%solution(i,1)
                end if

                ! Upper row
                oldSol = self%solution(i,self%N_y)
                if (self%world%boundary_conditions(i_finest, self%world%N_y) == 2) then
                    self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - 2.0d0 * self%solution(i, self%N_y-1) * self%coeffY - &
                        (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end if
            end do
            !$OMP end do nowait

            ! Red Left/right columns
            !$OMP do 
            do j = 3, self%N_y-2, 2
                j_finest = j * self%y_indx_step - self%y_indx_step + 1

                ! Left column
                oldSol = self%solution(1,j)
                if (self%world%boundary_conditions(1, j_finest) == 2) then
                    self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j-1) + self%solution(1, j+1)) * self%coeffY - &
                        2.0d0 * self%solution(2, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                else if (self%world%boundary_conditions(1,j_finest) == 3) then
                    self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                        (self%solution(2, j) + self%solution(self%N_x-1, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                    self%solution(self%N_x, j) = self%solution(1,j)
                end if

                ! Right column
                oldSol = self%solution(self%N_x,j)
                if (self%world%boundary_conditions(self%world%N_x, j_finest) == 2) then
                    self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j-1) + self%solution(self%N_x, j+1)) * self%coeffY - &
                        2.0d0 * self%solution(self%N_x-1, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end if
            end do
            !$OMP end do


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
                    end do
                end do
            end do
            !$OMP end do nowait

            ! Now go through black boundaries
            ! black Upper/lower rows
            !$OMP do 
            do i = 2, self%N_x-1, 2
                i_finest = i * self%x_indx_step - self%x_indx_step + 1

                ! Lower row
                oldSol = self%solution(i,1)
                if (self%world%boundary_conditions(i_finest, 1) == 2) then
                    self%solution(i,1) = (self%sourceTerm(i,1) - 2.0d0 * self%solution(i, 2) * self%coeffY - &
                        (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                else if (self%world%boundary_conditions(i_finest, 1) == 3) then
                    self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, self%N_y-1)) * self%coeffY - &
                        (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                    self%solution(i, self%N_y) = self%solution(i,1)
                end if

                ! Upper row
                oldSol = self%solution(i,self%N_y)
                if (self%world%boundary_conditions(i_finest, self%world%N_y) == 2) then
                    self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - 2.0d0 * self%solution(i, self%N_y-1) * self%coeffY - &
                        (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end if
            end do
            !$OMP end do nowait

            ! black Left/right columns
            !$OMP do 
            do j = 2, self%N_y-1, 2
                j_finest = j * self%y_indx_step - self%y_indx_step + 1

                ! Left column
                oldSol = self%solution(1,j)
                if (self%world%boundary_conditions(1, j_finest) == 2) then
                    self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j-1) + self%solution(1, j+1)) * self%coeffY - &
                        2.0d0 * self%solution(2, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                else if (self%world%boundary_conditions(1,j_finest) == 3) then
                    self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                        (self%solution(2, j) + self%solution(self%N_x-1, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                    self%solution(self%N_x, j) = self%solution(1,j)
                end if

                ! Right column
                oldSol = self%solution(self%N_x,j)
                if (self%world%boundary_conditions(self%world%N_x, j_finest) == 2) then
                    self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j-1) + self%solution(self%N_x, j+1)) * self%coeffY - &
                        2.0d0 * self%solution(self%N_x-1, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end if
            end do
            !$OMP end do
            
            !$OMP end parallel
        end do

    
    end subroutine smoothIterations_RedBlackEven

    function smoothWithRes_RedBlackEven(self) result(Res)
        ! Solve GS down to some tolerance
        class(RedBlackSolverEven), intent(in out) :: self
        real(real64) :: oldSol, Res, omega_inv, res_first
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, i_finest, j_finest
        omega_inv = 1.0d0 - self%omega
        res_first = 0.0d0

        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            oldSol = self%solution(1,1)
            self%solution(1,1) = (self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * self%coeffY - &
                2.0d0 * self%solution(2, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1,1)
            Res = Res + (self%solution(1,1) - oldSol)**2
        end if
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            oldSol = self%solution(self%N_x,1)
            self%solution(self%N_x, 1) = (self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, 1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x, 1)
            Res = Res + (self%solution(self%N_x,1) - oldSol)**2
        end if
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            oldSol = self%solution(1, self%N_y)
            self%solution(1, self%N_y) = (self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(2, self%N_y-1) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(1, self%N_y)
            Res = Res + (self%solution(1,self%N_y) - oldSol)**2   
        end if
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            oldSol = self%solution(self%N_x, self%N_y)
            self%solution(self%N_x, self%N_y) = (self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, self%N_y) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * self%solution(self%N_x, self%N_y)
            Res = Res + (self%solution(self%N_x, self%N_y) - oldSol)**2
        end if
        res_first = Res
        Res = 0.0d0
        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol, i_finest, j_finest) reduction(+:Res)

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
                    Res = Res + (self%solution(i, j) - oldSol)**2
                end do
            end do
        end do
        !$OMP end do nowait

        ! Now go through red boundaries
        ! Red Upper/lower rows
        !$OMP do 
        do i = 3, self%N_x-2, 2
            i_finest = i * self%x_indx_step - self%x_indx_step + 1

            ! Lower row
            oldSol = self%solution(i,1)
            if (self%world%boundary_conditions(i_finest, 1) == 2) then
                self%solution(i,1) = (self%sourceTerm(i,1) - 2.0d0 * self%solution(i, 2) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            else if (self%world%boundary_conditions(i_finest, 1) == 3) then
                self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, self%N_y-1)) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                self%solution(i, self%N_y) = self%solution(i,1)
            end if
            Res = Res + (self%solution(i, 1) - oldSol)**2
            ! Upper row
            oldSol = self%solution(i,self%N_y)
            if (self%world%boundary_conditions(i_finest, self%world%N_y) == 2) then
                self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - 2.0d0 * self%solution(i, self%N_y-1) * self%coeffY - &
                    (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end if
            Res = Res + (self%solution(i, self%N_y) - oldSol)**2
            
        end do
        !$OMP end do nowait

        ! Red Left/right columns
        !$OMP do 
        do j = 3, self%N_y-2, 2
            j_finest = j * self%y_indx_step - self%y_indx_step + 1

            ! Left column
            oldSol = self%solution(1,j)
            if (self%world%boundary_conditions(1, j_finest) == 2) then
                self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j-1) + self%solution(1, j+1)) * self%coeffY - &
                    2.0d0 * self%solution(2, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            else if (self%world%boundary_conditions(1,j_finest) == 3) then
                self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                    (self%solution(2, j) + self%solution(self%N_x-1, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                self%solution(self%N_x, j) = self%solution(1,j)
            end if
            Res = Res + (self%solution(1,j) - oldSol)**2

            ! Right column
            oldSol = self%solution(self%N_x,j)
            if (self%world%boundary_conditions(self%world%N_x, j_finest) == 2) then
                self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j-1) + self%solution(self%N_x, j+1)) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end if
            Res = Res + (self%solution(self%N_x,j) - oldSol)**2
        end do
        !$OMP end do


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

        ! Now go through black boundaries
        ! black Upper/lower rows
        !$OMP do 
        do i = 2, self%N_x-1, 2
            i_finest = i * self%x_indx_step - self%x_indx_step + 1

            ! Lower row
            oldSol = self%solution(i,1)
            if (self%world%boundary_conditions(i_finest, 1) == 2) then
                self%solution(i,1) = (self%sourceTerm(i,1) - 2.0d0 * self%solution(i, 2) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            else if (self%world%boundary_conditions(i_finest, 1) == 3) then
                self%solution(i,1) = (self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, self%N_y-1)) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                self%solution(i, self%N_y) = self%solution(i,1)
            end if
            Res = Res + (self%solution(i,1) - oldSol)**2

            ! Upper row
            oldSol = self%solution(i,self%N_y)
            if (self%world%boundary_conditions(i_finest, self%world%N_y) == 2) then
                self%solution(i,self%N_y) = (self%sourceTerm(i,self%N_y) - 2.0d0 * self%solution(i, self%N_y-1) * self%coeffY - &
                    (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end if
            Res = Res + (self%solution(i,self%N_y) - oldSol)**2
        end do
        !$OMP end do nowait

        ! black Left/right columns
        !$OMP do 
        do j = 2, self%N_y-1, 2
            j_finest = j * self%y_indx_step - self%y_indx_step + 1

            ! Left column
            oldSol = self%solution(1,j)
            if (self%world%boundary_conditions(1, j_finest) == 2) then
                self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j-1) + self%solution(1, j+1)) * self%coeffY - &
                    2.0d0 * self%solution(2, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            else if (self%world%boundary_conditions(1,j_finest) == 3) then
                self%solution(1,j) = (self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                    (self%solution(2, j) + self%solution(self%N_x-1, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                self%solution(self%N_x, j) = self%solution(1,j)
            end if
            Res = Res + (self%solution(1,j) - oldSol)**2

            ! Right column
            oldSol = self%solution(self%N_x,j)
            if (self%world%boundary_conditions(self%world%N_x, j_finest) == 2) then
                self%solution(self%N_x,j) = (self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j-1) + self%solution(self%N_x, j+1)) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, j) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end if
            Res = Res + (self%solution(self%N_x,j) - oldSol)**2
        end do
        !$OMP end do
        !$OMP end parallel
        
        Res = SQRT((Res+res_first)/ (self%number_solve_nodes))

    end function smoothWithRes_RedBlackEven


end module mod_RedBlackSolverEven