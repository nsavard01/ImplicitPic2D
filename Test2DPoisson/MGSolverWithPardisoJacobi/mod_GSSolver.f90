module mod_GSSolver
    use iso_fortran_env, only: int32, int64, real64
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
    public :: GSSolver

    type :: GSSolver
        ! store grid quantities
        real(real64), allocatable :: sourceTerm(:), solution(:), residual(:)
        integer(int32), allocatable :: black_NESW_indx(:,:), red_NESW_indx(:,:), restrictionIndx(:,:)
        integer(int32) :: numberBlackNodes, numberRedNodes, matDimension, iterNumber, numberRestrictionNodes, numberProlongation2PointNodes, numberProlongation4PointNodes
        real(real64) :: coeffX, coeffY, coeffSelf, omega
    contains
        procedure, public, pass(self) :: constructPoissonEven
        procedure, public, pass(self) :: constructRestrictionIndex
        ! procedure, public, pass(self) :: constructProlongationIndex
        procedure, public, pass(self) :: solveGS
        procedure, public, pass(self) :: calcResidual
        procedure, public, pass(self) :: smoothIterations
        procedure, public, pass(self) :: smoothIterationsWithRes
    end type

    interface GSSolver
        module procedure :: GSSolver_constructor
    end interface GSSolver

contains

    type(GSSolver) function GSSolver_constructor(omega) result(self)
        ! Construct object, set initial variables
        real(real64), intent(in) :: omega
        self%omega = omega
        self%iterNumber = 0
        self%Residual = 0
    end function GSSolver_constructor

    subroutine constructPoissonEven(self, N_x, N_y, delX, delY, NESW_wallBoundaries, boundaryConditions)
        class(GSSolver), intent(in out) :: self
        integer(int32), intent(in) :: N_x, N_y, NESW_wallBoundaries(4), boundaryConditions(N_x*N_y)
        real(real64), intent(in) :: delX, delY
        integer :: matDimension, upperBound, lowerBound, rightBound, leftBound
        integer :: i, j, k, tempNumRed, tempNumBlack
        matDimension = N_x * N_y
        self%matDimension = matDimension

        upperBound = NESW_wallBoundaries(1)
        rightBound = NESW_wallBoundaries(2)
        lowerBound = NESW_wallBoundaries(3)
        leftBound = NESW_wallBoundaries(4)


        self%coeffY = 1.0d0 / delY**2
        self%coeffX = 1.0d0 / delX**2
        self%coeffSelf = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
        self%coeffSelf = 1.0d0 / self%coeffSelf

        allocate(self%sourceTerm(matDimension), self%solution(matDimension), self%residual(matDimension))
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

        tempNumRed = 0
        tempNumBlack = 0
        ! Count number red and black nodes
        ! Also fill in any dirichlet boundaries solution vector
        !$OMP parallel reduction(+:tempNumRed, tempNumBlack)
        !$OMP do
        do k = 1, matDimension, 2
            if (boundaryConditions(k) /= 1) then
                tempNumBlack = tempNumBlack+1
            end if
        end do
        !$OMP end do

        !$OMP do
        do k = 2, matDimension-1, 2
            if (boundaryConditions(k) /= 1) then
                tempNumRed = tempNumRed+1
            end if
        end do
        !$OMP end do
        !$OMP end parallel
        self%numberBlackNodes = tempNumBlack
        self%numberRedNodes = tempNumRed
        allocate(self%black_NESW_indx(5, self%numberBlackNodes), self%red_NESW_indx(5, self%numberRedNodes))


        i = 0
        ! Fill black index array
        do k = 1, matDimension, 2
            if (boundaryConditions(k) /= 1) then
                i = i + 1
                self%black_NESW_indx(1,i) = k
                self%black_NESW_indx(2,i) = k + N_x !North
                self%black_NESW_indx(3,i) = k+1 ! East
                self%black_NESW_indx(4,i) = k-N_x ! South
                self%black_NESW_indx(5,i) = k-1 ! West
                if (boundaryConditions(k) /= 0) then
                    if (k  <= N_x) then
                        ! lower boundary, change south
                        if (boundaryConditions(k) == 2) then
                            self%black_NESW_indx(4,i) = k+N_x
                        else
                            self%black_NESW_indx(4,i) = matDimension - 2*N_x + k
                        end if
                    else if(k >= matDimension - N_x + 1) then
                        ! upper boundary, change north
                        if (boundaryConditions(k) == 2) then
                            self%black_NESW_indx(2,i) = k-N_x
                        else
                            self%black_NESW_indx(2,i) = k-matDimension + 2 * N_x
                        end if
                    end if
                    if (MOD(k, N_x) == 1) then
                        ! left boundary, change West
                        if (boundaryConditions(k) == 2) then
                            self%black_NESW_indx(5,i) = k+1   
                        else
                            self%black_NESW_indx(5,i) = k+N_x-2
                        end if
                    else if (MOD(k, N_x) == 0) then
                        ! right boundary, change East
                        if (boundaryConditions(k) == 2) then
                            self%black_NESW_indx(3,i) = k-1
                        else
                            self%black_NESW_indx(3,i) = k-N_x+2
                        end if
                    end if
                end if
            end if
        end do

        i = 0
        ! Fill red index array
        do k = 2, matDimension-1, 2
            if (boundaryConditions(k) /= 1) then
                i = i + 1
                self%red_NESW_indx(1,i) = k
                self%red_NESW_indx(2,i) = k + N_x !North
                self%red_NESW_indx(3,i) = k+1 ! East
                self%red_NESW_indx(4,i) = k-N_x ! South
                self%red_NESW_indx(5,i) = k-1 ! West
                if (boundaryConditions(k) /= 0) then
                    if (k  <= N_x) then
                        ! lower boundary, change south
                        if (boundaryConditions(k) == 2) then
                            self%red_NESW_indx(4,i) = k+N_x
                        else
                            self%red_NESW_indx(4,i) = matDimension - 2*N_x + k
                        end if
                    else if(k >= matDimension - N_x + 1) then
                        ! upper boundary, change north
                        if (boundaryConditions(k) == 2) then
                            self%red_NESW_indx(2,i) = k-N_x
                        else
                            self%red_NESW_indx(2,i) = k-matDimension + 2 * N_x
                        end if
                    end if
                    if (MOD(k, N_x) == 1) then
                        ! left boundary, change West
                        if (boundaryConditions(k) == 2) then
                            self%red_NESW_indx(5,i) = k+1   
                        else
                            self%red_NESW_indx(5,i) = k+N_x-2
                        end if
                    else if (MOD(k, N_x) == 0) then
                        ! right boundary, change East
                        if (boundaryConditions(k) == 2) then
                            self%red_NESW_indx(3,i) = k-1
                        else
                            self%red_NESW_indx(3,i) = k-N_x+2
                        end if
                    end if
                end if
            end if
        end do
        
        
    end subroutine constructPoissonEven

    subroutine constructRestrictionIndex(self, numResNodes, N_x, N_y)
        class(GSSolver), intent(in out) :: self
        integer(int32), intent(in) :: numResNodes, N_x, N_y
        integer :: k, indx_fine, column_fine, row_fine, column_coarse, row_coarse, indx_coarse, i
        ! Construct Restriction index
        self%numberRestrictionNodes = numResNodes
        allocate(self%restrictionIndx(2, self%numberRestrictionNodes))
        i = 0
        do k = 1, self%numberBlackNodes ! each coarse node is on fine grid black node
            indx_fine = self%black_NESW_indx(1, k) ! get index on fine grid
            row_fine = (indx_fine-1)/N_x + 1 ! row number
            column_fine = indx_fine - (row_fine - 1)  * N_x ! column number
            if (MOD(column_fine,2) == 1 .and. MOD(row_fine,2) == 1) then ! Determine if at coarse node
                i = i + 1
                row_coarse = (row_fine + 1)/2
                column_coarse = (column_fine + 1)/2
                self%restrictionIndx(1,i) = k ! index in black_NESW_indx to point to
                self%restrictionIndx(2,i) = column_coarse + (row_coarse-1) * (N_x/2 + 1) ! index in coarse array to interpolate to
            end if
        end do
         
    end subroutine constructRestrictionIndex

    subroutine solveGS(self, tol)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: tol
        real(real64) :: Res
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k
        Res = 1.0
        self%iterNumber = 0
        do while (Res > tol)
            Res = 0.0d0
            !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx) reduction(+:Res)
            !$OMP do
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                N_indx = self%black_NESW_indx(2,k)
                E_indx = self%black_NESW_indx(3,k)
                S_indx = self%black_NESW_indx(4,k)
                W_indx = self%black_NESW_indx(5,k)
                self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
                Res = Res + (self%residual(O_indx) - self%solution(O_indx))**2
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                N_indx = self%red_NESW_indx(2,k)
                E_indx = self%red_NESW_indx(3,k)
                S_indx = self%red_NESW_indx(4,k)
                W_indx = self%red_NESW_indx(5,k)
                self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
                Res = Res + (self%residual(O_indx) - self%solution(O_indx))**2
            end do
            !$OMP end do
            !$OMP do
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                self%solution(O_indx) = self%residual(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                self%solution(O_indx) = self%residual(O_indx)
            end do
            !$OMP end do
            !$OMP end parallel
            Res = SQRT(Res/self%matDimension)
            self%iterNumber = self%iterNumber + 1
        end do



    end subroutine solveGS

    subroutine smoothIterations(self, iterNum, resetBool)
        ! Solve GS down with a certain amount of iterations
        class(GSSolver), intent(in out) :: self
        integer, intent(in) :: iterNum
        logical :: resetBool
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i
        ! No real difference putting openmp within the iteration loop
        ! If lower stage need to reset solution to 0
        if (resetBool) then
            !$OMP parallel
            !$OMP workshare
            self%solution = 0.0d0
            !$OMP end workshare
            !$OMP end parallel
        end if
        do i = 1, iterNum
            !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx)
            !$OMP do
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                N_indx = self%black_NESW_indx(2,k)
                E_indx = self%black_NESW_indx(3,k)
                S_indx = self%black_NESW_indx(4,k)
                W_indx = self%black_NESW_indx(5,k)
                self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                N_indx = self%red_NESW_indx(2,k)
                E_indx = self%red_NESW_indx(3,k)
                S_indx = self%red_NESW_indx(4,k)
                W_indx = self%red_NESW_indx(5,k)
                self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do
            !$OMP do
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                self%solution(O_indx) = self%residual(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                self%solution(O_indx) = self%residual(O_indx)
            end do
            !$OMP end do
            !$OMP end parallel
        end do



    end subroutine smoothIterations

    function smoothIterationsWithRes(self, iterNum) result(Res)
        ! Solve GS down to some tolerance, whilst calculating residual step in solution for last iteration
        ! Used for final stage to determine 
        class(GSSolver), intent(in out) :: self
        integer, intent(in) :: iterNum
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i
        real(real64) :: Res
        ! No real difference putting openmp within the iteration loop
        ! If lower stage need to reset solution to 0
        do i = 1, iterNum-1
            !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx)
            !$OMP do
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                N_indx = self%black_NESW_indx(2,k)
                E_indx = self%black_NESW_indx(3,k)
                S_indx = self%black_NESW_indx(4,k)
                W_indx = self%black_NESW_indx(5,k)
                self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                N_indx = self%red_NESW_indx(2,k)
                E_indx = self%red_NESW_indx(3,k)
                S_indx = self%red_NESW_indx(4,k)
                W_indx = self%red_NESW_indx(5,k)
                self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do
            !$OMP do
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                self%solution(O_indx) = self%residual(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                self%solution(O_indx) = self%residual(O_indx)
            end do
            !$OMP end do
            !$OMP end parallel
        end do
        Res = 0.0d0
        ! Add 
        !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx) reduction(+:Res)
        !$OMP do
        do k = 1, self%numberBlackNodes
            O_indx = self%black_NESW_indx(1,k)
            N_indx = self%black_NESW_indx(2,k)
            E_indx = self%black_NESW_indx(3,k)
            S_indx = self%black_NESW_indx(4,k)
            W_indx = self%black_NESW_indx(5,k)
            self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            Res = Res + (self%residual(O_indx) - self%solution(O_indx))**2
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedNodes
            O_indx = self%red_NESW_indx(1,k)
            N_indx = self%red_NESW_indx(2,k)
            E_indx = self%red_NESW_indx(3,k)
            S_indx = self%red_NESW_indx(4,k)
            W_indx = self%red_NESW_indx(5,k)
            self%residual(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            Res = Res + (self%residual(O_indx) - self%solution(O_indx))**2
        end do
        !$OMP end do
        !$OMP do
        do k = 1, self%numberBlackNodes
            O_indx = self%black_NESW_indx(1,k)
            self%solution(O_indx) = self%residual(O_indx)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedNodes
            O_indx = self%red_NESW_indx(1,k)
            self%solution(O_indx) = self%residual(O_indx)
        end do
        !$OMP end do
        !$OMP end parallel
        Res = SQRT(Res/(self%numberRedNodes + self%numberBlackNodes))


    end function smoothIterationsWithRes

    subroutine calcResidual(self)
        ! For multigrid, caclulate b - Ax for current solution of x, store in sourceTerm
        ! For Dirichlet boundaries, sourceTerm already 0, so can ignore
        class(GSSolver), intent(in out) :: self
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k
        
        !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx)
        !$OMP do
        do k = 1, self%numberBlackNodes
            O_indx = self%black_NESW_indx(1,k)
            N_indx = self%black_NESW_indx(2,k)
            E_indx = self%black_NESW_indx(3,k)
            S_indx = self%black_NESW_indx(4,k)
            W_indx = self%black_NESW_indx(5,k)
            self%residual(O_indx) = self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY &
                - (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX - self%solution(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedNodes
            O_indx = self%red_NESW_indx(1,k)
            N_indx = self%red_NESW_indx(2,k)
            E_indx = self%red_NESW_indx(3,k)
            S_indx = self%red_NESW_indx(4,k)
            W_indx = self%red_NESW_indx(5,k)
            self%residual(O_indx) = self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY &
                - (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX - self%solution(O_indx)/self%coeffSelf
        end do
        !$OMP end do
        !$OMP end parallel


    end subroutine calcResidual


end module mod_GSSolver