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
        integer(int32), allocatable :: black_InnerIndx(:), red_InnerIndx(:), black_NESW_BoundIndx(:,:), red_NESW_BoundIndx(:,:), restrictionInnerIndx(:,:), restrictionBoundIndx(:,:)
        integer(int32) :: numberBlackInnerNodes, numberRedInnerNodes, numberBlackBoundNodes, numberRedBoundNodes, matDimension, iterNumber, numberRestrictionInnerNodes, numberRestrictionBoundNodes, N_x
        real(real64) :: coeffX, coeffY, coeffSelf, omega
    contains
        procedure, public, pass(self) :: constructPoissonEven
        procedure, public, pass(self) :: constructRestrictionIndex
        procedure, public, pass(self) :: solveGS
        procedure, public, pass(self) :: calcResidual
        procedure, public, pass(self) :: smoothIterations
        procedure, public, pass(self) :: matMult
        procedure, public, pass(self) :: XAX_Mult
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
        integer :: i, j, k, tempNumInnerRed, tempNumInnerBlack, tempNumBoundRed, tempNumBoundBlack
        matDimension = N_x * N_y
        self%matDimension = matDimension
        self%N_x = N_x

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

        tempNumInnerRed = 0
        tempNumInnerBlack = 0
        tempNumBoundRed = 0
        tempNumBoundBlack = 0
        ! Count number red and black nodes
        ! Also fill in any dirichlet boundaries solution vector
        !$OMP parallel reduction(+:tempNumInnerRed, tempNumInnerBlack, tempNumBoundRed, tempNumBoundBlack)
        !$OMP do
        do k = 1, matDimension, 2
            if (boundaryConditions(k) == 0) then
                tempNumInnerBlack = tempNumInnerBlack+1
            else if (boundaryConditions(k) /= 1) then
                tempNumBoundBlack = tempNumBoundBlack + 1
            end if
        end do
        !$OMP end do

        !$OMP do
        do k = 2, matDimension-1, 2
            if (boundaryConditions(k) == 0) then
                tempNumInnerRed = tempNumInnerRed+1
            else if (boundaryConditions(k) /= 1) then
                tempNumBoundRed = tempNumBoundRed + 1
            end if
        end do
        !$OMP end do
        !$OMP end parallel
        self%numberBlackInnerNodes = tempNumInnerBlack
        self%numberRedInnerNodes = tempNumInnerRed
        self%numberBlackBoundNodes = tempNumBoundBlack
        self%numberRedBoundNodes = tempNumBoundRed

        allocate(self%black_InnerIndx(self%numberBlackInnerNodes), self%red_InnerIndx(self%numberRedInnerNodes), self%red_NESW_BoundIndx(9, self%numberRedBoundNodes), &
        self%black_NESW_BoundIndx(9, self%numberBlackBoundNodes))
       

        i = 0
        j = 0
        ! Fill black index array
        do k = 1, matDimension, 2
            if (boundaryConditions(k) == 0) then
                i = i + 1
                self%black_InnerIndx(i) = k ! node index
            else if (boundaryConditions(k) /= 1) then
                j = j + 1
                ! index 1 is Origin, index 2-9 is NW to W in clockwise direction
                self%black_NESW_BoundIndx(1,j) = k
                self%black_NESW_BoundIndx(2,j) = k + N_x -1 !NW
                self%black_NESW_BoundIndx(3,j) = k + N_x !N
                self%black_NESW_BoundIndx(4,j) = k + N_x +1 !NE
                self%black_NESW_BoundIndx(5,j) = k +1 !E
                self%black_NESW_BoundIndx(6,j) = k - N_x + 1 !SE
                self%black_NESW_BoundIndx(7,j) = k - N_x !S
                self%black_NESW_BoundIndx(8,j) = k - N_x -1 !SW
                self%black_NESW_BoundIndx(9,j) = k  -1 !W
                if (k  <= N_x) then
                    ! lower boundary, change south
                    if (boundaryConditions(k) == 2) then
                        self%black_NESW_BoundIndx(6:8,j) = self%black_NESW_BoundIndx(6:8,j) + 2*N_x
                    else
                        self%black_NESW_BoundIndx(6:8,j) = self%black_NESW_BoundIndx(6:8,j) + matDimension - N_x
                    end if
                else if(k >= matDimension - N_x + 1) then
                    ! upper boundary, change north
                    if (boundaryConditions(k) == 2) then
                        self%black_NESW_BoundIndx(2:4,j) = self%black_NESW_BoundIndx(2:4,j) - 2*N_x
                    else
                        self%black_NESW_BoundIndx(2:4,j) = self%black_NESW_BoundIndx(2:4,j) - matDimension + N_x
                    end if
                end if
                if (MOD(k, N_x) == 1) then
                    ! left boundary, change West
                    if (boundaryConditions(k) == 2) then
                        self%black_NESW_BoundIndx(8:9,j) = self%black_NESW_BoundIndx(8:9,j) + 2 
                        self%black_NESW_BoundIndx(2,j) = self%black_NESW_BoundIndx(2,j) + 2 
                    else
                        self%black_NESW_BoundIndx(8:9,j) = self%black_NESW_BoundIndx(8:9,j) + N_x - 1
                        self%black_NESW_BoundIndx(2,j) = self%black_NESW_BoundIndx(2,j) + N_x - 1
                    end if
                else if (MOD(k, N_x) == 0) then
                    ! right boundary, change East
                    if (boundaryConditions(k) == 2) then
                        self%black_NESW_BoundIndx(4:6,j) = self%black_NESW_BoundIndx(4:6,j) - 2
                    else
                        self%black_NESW_BoundIndx(4:6,j) = self%black_NESW_BoundIndx(4:6,j) - N_x + 1
                    end if
                end if
            end if
        end do

        i = 0
        j = 0
        ! Fill red index array
        do k = 2, matDimension-1, 2
            if (boundaryConditions(k) == 0) then
                i = i + 1
                self%red_InnerIndx(i) = k ! node index
            else if (boundaryConditions(k) /= 1) then
                j = j + 1
                ! index 1 is Origin, index 2-9 is NW to W in clockwise direction
                self%red_NESW_BoundIndx(1,j) = k
                self%red_NESW_BoundIndx(2,j) = k + N_x -1 !NW
                self%red_NESW_BoundIndx(3,j) = k + N_x !N
                self%red_NESW_BoundIndx(4,j) = k + N_x +1 !NE
                self%red_NESW_BoundIndx(5,j) = k +1 !E
                self%red_NESW_BoundIndx(6,j) = k - N_x + 1 !SE
                self%red_NESW_BoundIndx(7,j) = k - N_x !S
                self%red_NESW_BoundIndx(8,j) = k - N_x -1 !SW
                self%red_NESW_BoundIndx(9,j) = k  -1 !W
                if (k  <= N_x) then
                    ! lower boundary, change south
                    if (boundaryConditions(k) == 2) then
                        self%red_NESW_BoundIndx(6:8,j) = self%red_NESW_BoundIndx(6:8,j) + 2*N_x
                    else
                        self%red_NESW_BoundIndx(6:8,j) = self%red_NESW_BoundIndx(6:8,j) + matDimension - N_x
                    end if
                else if(k >= matDimension - N_x + 1) then
                    ! upper boundary, change north
                    if (boundaryConditions(k) == 2) then
                        self%red_NESW_BoundIndx(2:4,j) = self%red_NESW_BoundIndx(2:4,j) - 2*N_x
                    else
                        self%red_NESW_BoundIndx(2:4,j) = self%red_NESW_BoundIndx(2:4,j) - matDimension + N_x
                    end if
                end if
                if (MOD(k, N_x) == 1) then
                    ! left boundary, change West
                    if (boundaryConditions(k) == 2) then
                        self%red_NESW_BoundIndx(8:9,j) = self%red_NESW_BoundIndx(8:9,j) + 2 
                        self%red_NESW_BoundIndx(2,j) = self%red_NESW_BoundIndx(2,j) + 2 
                    else
                        self%red_NESW_BoundIndx(8:9,j) = self%red_NESW_BoundIndx(8:9,j) + N_x - 1
                        self%red_NESW_BoundIndx(2,j) = self%red_NESW_BoundIndx(2,j) + N_x - 1
                    end if
                else if (MOD(k, N_x) == 0) then
                    ! right boundary, change East
                    if (boundaryConditions(k) == 2) then
                        self%red_NESW_BoundIndx(4:6,j) = self%red_NESW_BoundIndx(4:6,j) - 2
                    else
                        self%red_NESW_BoundIndx(4:6,j) = self%red_NESW_BoundIndx(4:6,j) - N_x + 1
                    end if
                end if
            end if
        end do
        
        
    end subroutine constructPoissonEven

    subroutine constructRestrictionIndex(self, numInnerResNodes, numBoundResNodes, N_x)
        class(GSSolver), intent(in out) :: self
        integer(int32), intent(in) :: numInnerResNodes, numBoundResNodes, N_x
        integer :: k, indx_fine, column_fine, row_fine, column_coarse, row_coarse, i
        ! Construct Restriction index
        self%numberRestrictionInnerNodes = numInnerResNodes
        self%numberRestrictionBoundNodes = numBoundResNodes
        allocate(self%restrictionInnerIndx(2,self%numberRestrictionInnerNodes), self%restrictionBoundIndx(2,self%numberRestrictionBoundNodes))
        i = 0
        do k = 1, self%numberBlackInnerNodes ! each coarse node is on fine grid black node
            indx_fine = self%black_InnerIndx(k) ! get index on fine grid
            row_fine = (indx_fine-1)/N_x + 1 ! row number
            column_fine = indx_fine - (row_fine - 1)  * N_x ! column number
            if (MOD(column_fine,2) == 1 .and. MOD(row_fine,2) == 1) then ! Determine if at coarse node
                i = i + 1
                row_coarse = (row_fine + 1)/2
                column_coarse = (column_fine + 1)/2
                self%restrictionInnerIndx(1,i) = k ! index in black_NESW_InnerIndx to point to
                self%restrictionInnerIndx(2,i) = column_coarse + (row_coarse-1) * (N_x/2 + 1) ! index in coarse array to interpolate to
            end if
        end do

        i = 0
        do k = 1, self%numberBlackBoundNodes ! each coarse node is on fine grid black node
            indx_fine = self%black_NESW_BoundIndx(1,k) ! get index on fine grid
            row_fine = (indx_fine-1)/N_x + 1 ! row number
            column_fine = indx_fine - (row_fine - 1)  * N_x ! column number
            if (MOD(column_fine,2) == 1 .and. MOD(row_fine,2) == 1) then ! Determine if at coarse node
                i = i + 1
                row_coarse = (row_fine + 1)/2
                column_coarse = (column_fine + 1)/2
                self%restrictionBoundIndx(1,i) = k ! index in black_NESW_InnerIndx to point to
                self%restrictionBoundIndx(2,i) = column_coarse + (row_coarse-1) * (N_x/2 + 1) ! index in coarse array to interpolate to
            end if
        end do
         
    end subroutine constructRestrictionIndex

    subroutine matMult(self, x, y)
        ! Use gauss-seidel to calculate A*x, store solution into y
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: x(self%matDimension)
        real(real64), intent(in out) :: y(self%matDimension)
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k
        
        !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            O_indx = self%black_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            y(O_indx) = (x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            O_indx = self%black_NESW_BoundIndx(1,k)
            N_indx = self%black_NESW_BoundIndx(3,k)
            E_indx = self%black_NESW_BoundIndx(5,k)
            S_indx = self%black_NESW_BoundIndx(7,k)
            W_indx = self%black_NESW_BoundIndx(9,k)
            y(O_indx) = (x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            O_indx = self%red_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            y(O_indx) = (x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            O_indx = self%red_NESW_BoundIndx(1,k)
            N_indx = self%red_NESW_BoundIndx(3,k)
            E_indx = self%red_NESW_BoundIndx(5,k)
            S_indx = self%red_NESW_BoundIndx(7,k)
            W_indx = self%red_NESW_BoundIndx(9,k)
            y(O_indx) = (x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf
        end do
        !$OMP end do
        !$OMP end parallel
         
    end subroutine matMult

    function XAX_Mult(self, x) result(res)
        ! Use gauss-seidel to calculate x^T * A * x
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: x(self%matDimension)
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k
        real(real64) :: res
        res = 0.0d0
        !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx) reduction(+:res)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            O_indx = self%black_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            res = res + x(O_indx) * ((x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            O_indx = self%black_NESW_BoundIndx(1,k)
            N_indx = self%black_NESW_BoundIndx(3,k)
            E_indx = self%black_NESW_BoundIndx(5,k)
            S_indx = self%black_NESW_BoundIndx(7,k)
            W_indx = self%black_NESW_BoundIndx(9,k)
            res = res + x(O_indx) * ((x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            O_indx = self%red_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            res = res + x(O_indx) * ((x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            O_indx = self%red_NESW_BoundIndx(1,k)
            N_indx = self%red_NESW_BoundIndx(3,k)
            E_indx = self%red_NESW_BoundIndx(5,k)
            S_indx = self%red_NESW_BoundIndx(7,k)
            W_indx = self%red_NESW_BoundIndx(9,k)
            res = res + x(O_indx) * ((x(N_indx) + x(S_indx)) * self%coeffY &
               + (x(E_indx) + x(W_indx)) * self%coeffX + x(O_indx)/self%coeffSelf)
        end do
        !$OMP end do
        !$OMP end parallel
         
    end function XAX_Mult

    subroutine solveGS(self, tol)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: tol
        real(real64) :: oldSol, Res
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k
        Res = 1.0
        self%iterNumber = 0
        do while (Res > tol)
            Res = 0.0d0
            !$OMP parallel private(oldSol, O_indx, N_indx, E_indx, S_indx, W_indx) reduction(+:Res)
            !$OMP do
            do k = 1, self%numberBlackInnerNodes
                O_indx = self%black_InnerIndx(k)
                N_indx = O_indx + self%N_x
                E_indx = O_indx + 1
                S_indx = O_indx - self%N_x
                W_indx = O_indx - 1
                oldSol = self%solution(O_indx)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(O_indx) - oldSol)**2
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberBlackBoundNodes
                O_indx = self%black_NESW_BoundIndx(1,k)
                N_indx = self%black_NESW_BoundIndx(3,k)
                E_indx = self%black_NESW_BoundIndx(5,k)
                S_indx = self%black_NESW_BoundIndx(7,k)
                W_indx = self%black_NESW_BoundIndx(9,k)
                oldSol = self%solution(O_indx)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(O_indx) - oldSol)**2
            end do
            !$OMP end do
            ! Now red nodes
            !$OMP do
            do k = 1, self%numberRedInnerNodes
                O_indx = self%red_InnerIndx(k)
                N_indx = O_indx + self%N_x
                E_indx = O_indx + 1
                S_indx = O_indx - self%N_x
                W_indx = O_indx - 1
                oldSol = self%solution(O_indx)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(O_indx) - oldSol)**2
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedBoundNodes
                O_indx = self%red_NESW_BoundIndx(1,k)
                N_indx = self%red_NESW_BoundIndx(3,k)
                E_indx = self%red_NESW_BoundIndx(5,k)
                S_indx = self%red_NESW_BoundIndx(7,k)
                W_indx = self%red_NESW_BoundIndx(9,k)
                oldSol = self%solution(O_indx)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(O_indx) - oldSol)**2
            end do
            !$OMP end do
            !$OMP end parallel
            Res = SQRT(Res/(self%numberBlackInnerNodes + self%numberRedInnerNodes + self%numberBlackBoundNodes + self%numberRedBoundNodes))
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
            do k = 1, self%numberBlackInnerNodes
                O_indx = self%black_InnerIndx(k)
                N_indx = O_indx + self%N_x
                E_indx = O_indx + 1
                S_indx = O_indx - self%N_x
                W_indx = O_indx - 1
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberBlackBoundNodes
                O_indx = self%black_NESW_BoundIndx(1,k)
                N_indx = self%black_NESW_BoundIndx(3,k)
                E_indx = self%black_NESW_BoundIndx(5,k)
                S_indx = self%black_NESW_BoundIndx(7,k)
                W_indx = self%black_NESW_BoundIndx(9,k)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do
            ! Now red nodes
            !$OMP do
            do k = 1, self%numberRedInnerNodes
                O_indx = self%red_InnerIndx(k)
                N_indx = O_indx + self%N_x
                E_indx = O_indx + 1
                S_indx = O_indx - self%N_x
                W_indx = O_indx - 1
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedBoundNodes
                O_indx = self%red_NESW_BoundIndx(1,k)
                N_indx = self%red_NESW_BoundIndx(3,k)
                E_indx = self%red_NESW_BoundIndx(5,k)
                S_indx = self%red_NESW_BoundIndx(7,k)
                W_indx = self%red_NESW_BoundIndx(9,k)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
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
        real(real64) :: Res, oldSol
        ! No real difference putting openmp within the iteration loop
        ! If lower stage need to reset solution to 0
        do i = 1, iterNum-1
            !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx)
            !$OMP do
            do k = 1, self%numberBlackInnerNodes
                O_indx = self%black_InnerIndx(k)
                N_indx = O_indx + self%N_x
                E_indx = O_indx + 1
                S_indx = O_indx - self%N_x
                W_indx = O_indx - 1
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberBlackBoundNodes
                O_indx = self%black_NESW_BoundIndx(1,k)
                N_indx = self%black_NESW_BoundIndx(3,k)
                E_indx = self%black_NESW_BoundIndx(5,k)
                S_indx = self%black_NESW_BoundIndx(7,k)
                W_indx = self%black_NESW_BoundIndx(9,k)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do
            ! Now red nodes
            !$OMP do
            do k = 1, self%numberRedInnerNodes
                O_indx = self%red_InnerIndx(k)
                N_indx = O_indx + self%N_x
                E_indx = O_indx + 1
                S_indx = O_indx - self%N_x
                W_indx = O_indx - 1
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedBoundNodes
                O_indx = self%red_NESW_BoundIndx(1,k)
                N_indx = self%red_NESW_BoundIndx(3,k)
                E_indx = self%red_NESW_BoundIndx(5,k)
                S_indx = self%red_NESW_BoundIndx(7,k)
                W_indx = self%red_NESW_BoundIndx(9,k)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * self%solution(O_indx)
            end do
            !$OMP end do
            !$OMP end parallel
        end do
        Res = 0.0d0
        ! Add 
        !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx, oldSol) reduction(+:Res)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            O_indx = self%black_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            oldSol = self%solution(O_indx)
            self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(O_indx) - oldSol)**2
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            O_indx = self%black_NESW_BoundIndx(1,k)
            N_indx = self%black_NESW_BoundIndx(3,k)
            E_indx = self%black_NESW_BoundIndx(5,k)
            S_indx = self%black_NESW_BoundIndx(7,k)
            W_indx = self%black_NESW_BoundIndx(9,k)
            oldSol = self%solution(O_indx)
            self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(O_indx) - oldSol)**2
        end do
        !$OMP end do
        ! Now red nodes
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            O_indx = self%red_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            oldSol = self%solution(O_indx)
            self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(O_indx) - oldSol)**2
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            O_indx = self%red_NESW_BoundIndx(1,k)
            N_indx = self%red_NESW_BoundIndx(3,k)
            E_indx = self%red_NESW_BoundIndx(5,k)
            S_indx = self%red_NESW_BoundIndx(7,k)
            W_indx = self%red_NESW_BoundIndx(9,k)
            oldSol = self%solution(O_indx)
            self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(O_indx) - oldSol)**2
        end do
        !$OMP end do
        !$OMP end parallel
        Res = SQRT(Res/(self%numberBlackInnerNodes + self%numberRedInnerNodes + self%numberBlackBoundNodes + self%numberRedBoundNodes))


    end function smoothIterationsWithRes

    subroutine calcResidual(self)
        ! For multigrid, caclulate b - Ax for current solution of x, store in sourceTerm
        ! For Dirichlet boundaries, sourceTerm already 0, so can ignore
        class(GSSolver), intent(in out) :: self
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k
        
        !$OMP parallel private(O_indx, N_indx, E_indx, S_indx, W_indx)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            O_indx = self%black_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            self%residual(O_indx) = self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY &
                - (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX - self%solution(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            O_indx = self%black_NESW_BoundIndx(1,k)
            N_indx = self%black_NESW_BoundIndx(3,k)
            E_indx = self%black_NESW_BoundIndx(5,k)
            S_indx = self%black_NESW_BoundIndx(7,k)
            W_indx = self%black_NESW_BoundIndx(9,k)
            self%residual(O_indx) = self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY &
                - (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX - self%solution(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            O_indx = self%red_InnerIndx(k)
            N_indx = O_indx + self%N_x
            E_indx = O_indx + 1
            S_indx = O_indx - self%N_x
            W_indx = O_indx - 1
            self%residual(O_indx) = self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY &
                - (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX - self%solution(O_indx)/self%coeffSelf
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            O_indx = self%red_NESW_BoundIndx(1,k)
            N_indx = self%red_NESW_BoundIndx(3,k)
            E_indx = self%red_NESW_BoundIndx(5,k)
            S_indx = self%red_NESW_BoundIndx(7,k)
            W_indx = self%red_NESW_BoundIndx(9,k)
            self%residual(O_indx) = self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY &
                - (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX - self%solution(O_indx)/self%coeffSelf
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine calcResidual


end module mod_GSSolver