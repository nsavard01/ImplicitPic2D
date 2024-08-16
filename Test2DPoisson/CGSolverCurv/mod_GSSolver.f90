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
        real(real64), allocatable :: sourceTerm(:,:), solution(:,:), residual(:,:)
        integer(int32), allocatable :: black_InnerIndx(:,:), red_InnerIndx(:,:), black_NESW_BoundIndx(:,:), red_NESW_BoundIndx(:,:), restrictionInnerIndx(:), restrictionBoundIndx(:)
        ! For Bound index, have [i, j, i_East, i_West, j_North, j_South], basically direction of each NESW
        integer(int32) :: numberBlackInnerNodes, numberRedInnerNodes, numberBlackBoundNodes, numberRedBoundNodes, iterNumber, N_x, N_y, numberRestrictionInnerNodes, numberRestrictionBoundNodes
        real(real64) :: omega
        real(real64), allocatable :: matCoeffsInnerBlack(:,:), matCoeffsInnerRed(:,:), matCoeffsBoundBlack(:,:), matCoeffsBoundRed(:,:)! Coefficients in order of O, N, E, S, W
        ! Note that saving coefficients, while taking more memory, results in much faster (more than 50%) speedup. So assuming we are more limited by calculation time than memory, best to store
        ! everything beforehand
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
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

    subroutine constructPoissonOrthogonal(self, N_x, N_y, diffX, diffY, NESW_wallBoundaries, boundaryConditions)
        ! Construct orthogonal grid solver
        class(GSSolver), intent(in out) :: self
        integer(int32), intent(in) :: N_x, N_y, NESW_wallBoundaries(4), boundaryConditions(N_x, N_y)
        real(real64), intent(in) :: diffX(N_x-1), diffY(N_y-1)
        integer :: upperBound, lowerBound, rightBound, leftBound
        integer :: i, j, k, tempNumInnerRed, tempNumInnerBlack, tempNumBoundRed, tempNumBoundBlack, N_indx, S_indx, E_indx, W_indx
        real(real64) :: delX_E, delX_W, delY_N, delY_S, tempReal, C_N, C_S, C_E, C_W, C_O
        self%N_x = N_x
        self%N_y = N_y

        upperBound = NESW_wallBoundaries(1)
        rightBound = NESW_wallBoundaries(2)
        lowerBound = NESW_wallBoundaries(3)
        leftBound = NESW_wallBoundaries(4)


        allocate(self%sourceTerm(self%N_x, self%N_y), self%solution(self%N_x, self%N_y), self%residual(self%N_x, self%N_y))
        
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
        !$OMP do collapse(2)
        do j = 1, self%N_y
            do i = 1, self%N_x
                if (MOD(i,2) == MOD(j,2)) then
                    ! Black node
                    if (boundaryConditions(i,j) == 0) then
                        tempNumInnerBlack = tempNumInnerBlack+1
                    else if (boundaryConditions(i,j) /= 1) then
                        tempNumBoundBlack = tempNumBoundBlack + 1
                    end if
                else
                    !red node
                    if (boundaryConditions(i,j) == 0) then
                        tempNumInnerRed = tempNumInnerRed+1
                    else if (boundaryConditions(i,j) /= 1) then
                        tempNumBoundRed = tempNumBoundRed + 1
                    end if
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        self%numberBlackInnerNodes = tempNumInnerBlack
        self%numberRedInnerNodes = tempNumInnerRed
        self%numberBlackBoundNodes = tempNumBoundBlack
        self%numberRedBoundNodes = tempNumBoundRed
        


        allocate(self%black_InnerIndx(2, self%numberBlackInnerNodes), self%red_InnerIndx(2, self%numberRedInnerNodes), self%red_NESW_BoundIndx(6, self%numberRedBoundNodes), &
        self%black_NESW_BoundIndx(6, self%numberBlackBoundNodes), self%matCoeffsInnerBlack(5, self%numberBlackInnerNodes), self%matCoeffsInnerRed(5, self%numberRedInnerNodes), &
        self%matCoeffsBoundBlack(5, self%numberBlackBoundNodes), self%matCoeffsBoundRed(5, self%numberRedBoundNodes))
       

        tempNumInnerRed = 0
        tempNumInnerBlack = 0
        tempNumBoundRed = 0
        tempNumBoundBlack = 0
        ! Fill arrays 
        do j = 1, self%N_y
            do i = 1, self%N_x
                if (boundaryConditions(i,j) /= 1) then
                    if (boundaryConditions(i,j) == 0) then
                        delX_E = diffX(i)
                        delX_W = diffX(i-1)
                        delY_N = diffY(j)
                        delY_S = diffY(j-1)
                        ! Inner node, set coefficients
                        C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                        C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                        C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                        C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                        C_O = -1.0d0 / (C_N + C_E + C_W + C_S) ! O
                        if (MOD(i,2) == MOD(j,2)) then
                            ! Black node
                            tempNumInnerBlack = tempNumInnerBlack+1
                            self%black_InnerIndx(:, tempNumInnerBlack) = [i,j]
                            self%matCoeffsInnerBlack(:, tempNumInnerBlack) = [C_O, C_N, C_E, C_S, C_W]
                        else
                            !red node
                            tempNumInnerRed = tempNumInnerRed+1
                            self%red_InnerIndx(:, tempNumInnerRed) = [i,j]
                            self%matCoeffsInnerRed(:, tempNumInnerRed) = [C_O, C_N, C_E, C_S, C_W]
                        end if
                    else
                        ! Initial indices
                        N_indx = j+1
                        S_indx = j-1
                        E_indx = i+1
                        W_indx = i-1
                        if (boundaryConditions(i,j) == 2) then
                            !Neumann boundary
                            if (i == 1) then
                                ! left boundary, Switch west
                                W_indx = E_indx
                                delX_E = diffX(i)
                                delX_W = delX_E
                            else if (i == self%N_x) then
                                ! right boundary, Switch east
                                E_indx = W_indx
                                delX_W = diffX(i-1)
                                delX_E = delX_W
                            else
                                delX_E = diffX(i)
                                delX_W = diffX(i-1)
                            end if
                            if (j == 1) then
                                ! lower boundary, switch south
                                S_indx = N_indx
                                delY_N = diffY(j)
                                delY_S = delY_N
                            else if (j == self%N_y) then
                                ! upper boundary, switch north
                                N_indx = S_indx
                                delY_S = diffY(j-1)
                                delY_N = delY_S
                            else
                                delY_N = diffY(j)
                                delY_S = diffY(j-1)
                            end if
                        else
                            ! periodic boundary
                            if (i == 1) then
                                ! left boundary, Switch west
                                W_indx = self%N_x - 1
                            else if (i == self%N_x) then
                                ! right boundary, Switch east
                                E_indx = 2
                            end if
                            if (j == 1) then
                                ! lower boundary, switch south
                                S_indx = self%N_y - 1
                            else if (j == self%N_y) then
                                ! upper boundary, switch north
                                N_indx = 2
                            end if
                            delX_E = diffX(E_indx-1)
                            delX_W = diffX(W_indx)
                            delY_N = diffY(N_indx-1)
                            delY_S = diffY(S_indx)
                        end if
                        
                        ! Bound node, set coefficients
                        C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                        C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                        C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                        C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                        C_O = -1.0d0 / (C_N + C_E + C_W + C_S) ! O
                        if (MOD(i,2) == MOD(j,2)) then
                            ! Black node
                            tempNumBoundBlack = tempNumBoundBlack + 1
                            self%black_NESW_BoundIndx(:, tempNumBoundBlack) = [i,j,E_indx, W_indx, N_indx, S_indx]
                            self%matCoeffsBoundBlack(:, tempNumBoundBlack) = [C_O, C_N, C_E, C_S, C_W]
                        else
                            !red node
                            tempNumBoundRed = tempNumBoundRed + 1
                            self%red_NESW_BoundIndx(:, tempNumBoundRed) = [i,j,E_indx, W_indx, N_indx, S_indx]
                            self%matCoeffsBoundRed(:, tempNumBoundRed) = [C_O, C_N, C_E, C_S, C_W]
                        end if
                    end if
                end if
            end do
        end do
        
    end subroutine constructPoissonOrthogonal

    subroutine constructRestrictionIndex(self)
        class(GSSolver), intent(in out) :: self
        integer :: k, i_fine, j_fine, i_coarse, j_coarse, numInnerResNodes, numBoundResNodes
        ! Construct Restriction index
        ! coarse nodes overlap black nodes on odd rows (j index)
        numInnerResNodes = 0
        numBoundResNodes = 0
        !$OMP parallel private(j_fine) reduction(+:numInnerResNodes, numBoundResNodes)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            j_fine = self%black_InnerIndx(2, k)
            if (MOD(j_fine, 2) == 1) then
                numInnerResNodes = numInnerResNodes + 1
            end if
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            j_fine = self%black_NESW_BoundIndx(2, k)
            if (MOD(j_fine, 2) == 1) then
                numBoundResNodes = numBoundResNodes + 1
            end if
        end do
        !$OMP end do
        !$OMP end parallel
        self%numberRestrictionInnerNodes = numInnerResNodes 
        self%numberRestrictionBoundNodes = numBoundResNodes

        allocate(self%restrictionInnerIndx(self%numberRestrictionInnerNodes), self%restrictionBoundIndx(self%numberRestrictionBoundNodes))
        numInnerResNodes = 0
        do k = 1, self%numberBlackInnerNodes ! each coarse node is on fine grid black node
            j_fine = self%black_InnerIndx(2, k)
            if (MOD(j_fine, 2) == 1) then
                numInnerResNodes = numInnerResNodes + 1
                self%restrictionInnerIndx(numInnerResNodes) = k
            end if
        end do

        numBoundResNodes = 0
        do k = 1, self%numberBlackBoundNodes ! each coarse node is on fine grid black node
            j_fine = self%black_NESW_BoundIndx(2, k) ! get index on fine grid
            if (MOD(j_fine, 2) == 1) then
                numBoundResNodes = numBoundResNodes + 1
                self%restrictionBoundIndx(numBoundResNodes) = k
            end if
        end do
         
    end subroutine constructRestrictionIndex

    subroutine matMult(self, x, y)
        ! Use gauss-seidel to calculate A*x, store solution into y
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y)
        real(real64), intent(in out) :: y(self%N_x, self%N_y)
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j
        real(real64) :: C_N, C_E, C_O, C_W, C_S

        !$OMP parallel private(i, j, N_indx, E_indx, S_indx, W_indx, C_N, C_E, C_O, C_W, C_S)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            i = self%black_InnerIndx(1, k)
            j = self%black_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerBlack(1, k)
            C_N = self%matCoeffsInnerBlack(2, k)
            C_E = self%matCoeffsInnerBlack(3, k)
            C_S = self%matCoeffsInnerBlack(4, k)
            C_W = self%matCoeffsInnerBlack(5, k)
            y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i = self%black_NESW_BoundIndx(1, k)
            j = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundBlack(1, k)
            C_N = self%matCoeffsBoundBlack(2, k)
            C_E = self%matCoeffsBoundBlack(3, k)
            C_S = self%matCoeffsBoundBlack(4, k)
            C_W = self%matCoeffsBoundBlack(5, k)
            y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            i = self%red_InnerIndx(1, k)
            j = self%red_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerRed(1, k)
            C_N = self%matCoeffsInnerRed(2, k)
            C_E = self%matCoeffsInnerRed(3, k)
            C_S = self%matCoeffsInnerRed(4, k)
            C_W = self%matCoeffsInnerRed(5, k)
            y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            i = self%red_NESW_BoundIndx(1, k)
            j = self%red_NESW_BoundIndx(2, k)
            E_indx = self%red_NESW_BoundIndx(3, k)
            W_indx = self%red_NESW_BoundIndx(4, k)
            N_indx = self%red_NESW_BoundIndx(5, k)
            S_indx = self%red_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundRed(1, k)
            C_N = self%matCoeffsBoundRed(2, k)
            C_E = self%matCoeffsBoundRed(3, k)
            C_S = self%matCoeffsBoundRed(4, k)
            C_W = self%matCoeffsBoundRed(5, k)
            y(i, j) = x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O
        end do
        !$OMP end do
        !$OMP end parallel
         
    end subroutine matMult

    function XAX_Mult(self, x) result(res)
        ! Use gauss-seidel to calculate x^T * A * x
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: x(self%N_x, self%N_y)
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j
        real(real64) :: res, C_N, C_E, C_O, C_W, C_S
        res = 0.0d0
        !$OMP parallel private(i, j, N_indx, E_indx, S_indx, W_indx, C_N, C_E, C_O, C_W, C_S)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            i = self%black_InnerIndx(1, k)
            j = self%black_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerBlack(1, k)
            C_N = self%matCoeffsInnerBlack(2, k)
            C_E = self%matCoeffsInnerBlack(3, k)
            C_S = self%matCoeffsInnerBlack(4, k)
            C_W = self%matCoeffsInnerBlack(5, k)
            res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i = self%black_NESW_BoundIndx(1, k)
            j = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundBlack(1, k)
            C_N = self%matCoeffsBoundBlack(2, k)
            C_E = self%matCoeffsBoundBlack(3, k)
            C_S = self%matCoeffsBoundBlack(4, k)
            C_W = self%matCoeffsBoundBlack(5, k)
            res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            i = self%red_InnerIndx(1, k)
            j = self%red_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerRed(1, k)
            C_N = self%matCoeffsInnerRed(2, k)
            C_E = self%matCoeffsInnerRed(3, k)
            C_S = self%matCoeffsInnerRed(4, k)
            C_W = self%matCoeffsInnerRed(5, k)
            res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            i = self%red_NESW_BoundIndx(1, k)
            j = self%red_NESW_BoundIndx(2, k)
            E_indx = self%red_NESW_BoundIndx(3, k)
            W_indx = self%red_NESW_BoundIndx(4, k)
            N_indx = self%red_NESW_BoundIndx(5, k)
            S_indx = self%red_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundRed(1, k)
            C_N = self%matCoeffsBoundRed(2, k)
            C_E = self%matCoeffsBoundRed(3, k)
            C_S = self%matCoeffsBoundRed(4, k)
            C_W = self%matCoeffsBoundRed(5, k)
            res = res + x(i,j)*(x(i,N_indx)*C_N + x(i,S_indx) * C_S &
                + x(E_indx, j) *C_E + x(W_indx, j) * C_W + x(i,j)/C_O)
        end do
        !$OMP end do
        !$OMP end parallel
         
    end function XAX_Mult

    subroutine solveGS(self, tol)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: tol
        real(real64) :: oldSol, Res, C_N, C_E, C_O, C_W, C_S
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k
        Res = 1.0
        self%iterNumber = 0
        do while (Res > tol)
            Res = 0.0d0
            !$OMP parallel private(oldSol, i, j, N_indx, E_indx, S_indx, &
            !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W) reduction(+:Res)
            !$OMP do
            do k = 1, self%numberBlackInnerNodes
                i = self%black_InnerIndx(1, k)
                j = self%black_InnerIndx(2, k)
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffsInnerBlack(1, k)
                C_N = self%matCoeffsInnerBlack(2, k)
                C_E = self%matCoeffsInnerBlack(3, k)
                C_S = self%matCoeffsInnerBlack(4, k)
                C_W = self%matCoeffsInnerBlack(5, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberBlackBoundNodes
                i = self%black_NESW_BoundIndx(1, k)
                j = self%black_NESW_BoundIndx(2, k)
                E_indx = self%black_NESW_BoundIndx(3, k)
                W_indx = self%black_NESW_BoundIndx(4, k)
                N_indx = self%black_NESW_BoundIndx(5, k)
                S_indx = self%black_NESW_BoundIndx(6, k)
                C_O = self%matCoeffsBoundBlack(1, k)
                C_N = self%matCoeffsBoundBlack(2, k)
                C_E = self%matCoeffsBoundBlack(3, k)
                C_S = self%matCoeffsBoundBlack(4, k)
                C_W = self%matCoeffsBoundBlack(5, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
            !$OMP end do
            ! Now red nodes
            !$OMP do
            do k = 1, self%numberRedInnerNodes
                i = self%red_InnerIndx(1, k)
                j = self%red_InnerIndx(2, k)
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffsInnerRed(1, k)
                C_N = self%matCoeffsInnerRed(2, k)
                C_E = self%matCoeffsInnerRed(3, k)
                C_S = self%matCoeffsInnerRed(4, k)
                C_W = self%matCoeffsInnerRed(5, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedBoundNodes
                i = self%red_NESW_BoundIndx(1, k)
                j = self%red_NESW_BoundIndx(2, k)
                E_indx = self%red_NESW_BoundIndx(3, k)
                W_indx = self%red_NESW_BoundIndx(4, k)
                N_indx = self%red_NESW_BoundIndx(5, k)
                S_indx = self%red_NESW_BoundIndx(6, k)
                C_O = self%matCoeffsBoundRed(1, k)
                C_N = self%matCoeffsBoundRed(2, k)
                C_E = self%matCoeffsBoundRed(3, k)
                C_S = self%matCoeffsBoundRed(4, k)
                C_W = self%matCoeffsBoundRed(5, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
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
        integer :: N_indx, E_indx, S_indx, W_indx, k, i, j, l
        real(real64) :: C_N, C_E, C_O, C_W, C_S
        ! No real difference putting openmp within the iteration loop
        ! If lower stage need to reset solution to 0
        if (resetBool) then
            !$OMP parallel
            !$OMP workshare
            self%solution = 0.0d0
            !$OMP end workshare
            !$OMP end parallel
        end if
        do l = 1, iterNum
            !$OMP parallel private(i,j,N_indx, E_indx, S_indx, W_indx, C_N, C_E, C_O, C_W, C_S)
            !$OMP do
            do k = 1, self%numberBlackInnerNodes
                i = self%black_InnerIndx(1, k)
                j = self%black_InnerIndx(2, k)
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffsInnerBlack(1, k)
                C_N = self%matCoeffsInnerBlack(2, k)
                C_E = self%matCoeffsInnerBlack(3, k)
                C_S = self%matCoeffsInnerBlack(4, k)
                C_W = self%matCoeffsInnerBlack(5, k)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * self%solution(i,j)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberBlackBoundNodes
                i = self%black_NESW_BoundIndx(1, k)
                j = self%black_NESW_BoundIndx(2, k)
                E_indx = self%black_NESW_BoundIndx(3, k)
                W_indx = self%black_NESW_BoundIndx(4, k)
                N_indx = self%black_NESW_BoundIndx(5, k)
                S_indx = self%black_NESW_BoundIndx(6, k)
                C_O = self%matCoeffsBoundBlack(1, k)
                C_N = self%matCoeffsBoundBlack(2, k)
                C_E = self%matCoeffsBoundBlack(3, k)
                C_S = self%matCoeffsBoundBlack(4, k)
                C_W = self%matCoeffsBoundBlack(5, k)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * self%solution(i,j)
            end do
            !$OMP end do
            ! Now red nodes
            !$OMP do
            do k = 1, self%numberRedInnerNodes
                i = self%red_InnerIndx(1, k)
                j = self%red_InnerIndx(2, k)
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffsInnerRed(1, k)
                C_N = self%matCoeffsInnerRed(2, k)
                C_E = self%matCoeffsInnerRed(3, k)
                C_S = self%matCoeffsInnerRed(4, k)
                C_W = self%matCoeffsInnerRed(5, k)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * self%solution(i,j)
            end do
            !$OMP end do nowait
            !$OMP do
            do k = 1, self%numberRedBoundNodes
                i = self%red_NESW_BoundIndx(1, k)
                j = self%red_NESW_BoundIndx(2, k)
                E_indx = self%red_NESW_BoundIndx(3, k)
                W_indx = self%red_NESW_BoundIndx(4, k)
                N_indx = self%red_NESW_BoundIndx(5, k)
                S_indx = self%red_NESW_BoundIndx(6, k)
                C_O = self%matCoeffsBoundRed(1, k)
                C_N = self%matCoeffsBoundRed(2, k)
                C_E = self%matCoeffsBoundRed(3, k)
                C_S = self%matCoeffsBoundRed(4, k)
                C_W = self%matCoeffsBoundRed(5, k)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * self%solution(i,j)
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
        integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j
        real(real64) :: Res, oldSol, C_O, C_N, C_E, C_S, C_W
        ! No real difference putting openmp within the iteration loop
        ! If lower stage need to reset solution to 0
        ! Do a few iteratuibs
        call self%smoothIterations(iterNum-1, .false.)
        Res = 0.0d0
        !$OMP parallel private(oldSol, i, j, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W) reduction(+:Res)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            i = self%black_InnerIndx(1, k)
            j = self%black_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerBlack(1, k)
            C_N = self%matCoeffsInnerBlack(2, k)
            C_E = self%matCoeffsInnerBlack(3, k)
            C_S = self%matCoeffsInnerBlack(4, k)
            C_W = self%matCoeffsInnerBlack(5, k)
            oldSol = self%solution(i,j)
            self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(i,j) - oldSol)**2
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i = self%black_NESW_BoundIndx(1, k)
            j = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundBlack(1, k)
            C_N = self%matCoeffsBoundBlack(2, k)
            C_E = self%matCoeffsBoundBlack(3, k)
            C_S = self%matCoeffsBoundBlack(4, k)
            C_W = self%matCoeffsBoundBlack(5, k)
            oldSol = self%solution(i,j)
            self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(i,j) - oldSol)**2
        end do
        !$OMP end do
        ! Now red nodes
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            i = self%red_InnerIndx(1, k)
            j = self%red_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerRed(1, k)
            C_N = self%matCoeffsInnerRed(2, k)
            C_E = self%matCoeffsInnerRed(3, k)
            C_S = self%matCoeffsInnerRed(4, k)
            C_W = self%matCoeffsInnerRed(5, k)
            oldSol = self%solution(i,j)
            self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(i,j) - oldSol)**2
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            i = self%red_NESW_BoundIndx(1, k)
            j = self%red_NESW_BoundIndx(2, k)
            E_indx = self%red_NESW_BoundIndx(3, k)
            W_indx = self%red_NESW_BoundIndx(4, k)
            N_indx = self%red_NESW_BoundIndx(5, k)
            S_indx = self%red_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundRed(1, k)
            C_N = self%matCoeffsBoundRed(2, k)
            C_E = self%matCoeffsBoundRed(3, k)
            C_S = self%matCoeffsBoundRed(4, k)
            C_W = self%matCoeffsBoundRed(5, k)
            oldSol = self%solution(i,j)
            self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(i,j) - oldSol)**2
        end do
        !$OMP end do
        !$OMP end parallel
        Res = SQRT(Res/(self%numberBlackInnerNodes + self%numberRedInnerNodes + self%numberBlackBoundNodes + self%numberRedBoundNodes))


    end function smoothIterationsWithRes

    subroutine calcResidual(self)
        ! For multigrid, caclulate b - Ax for current solution of x, store in sourceTerm
        ! For Dirichlet boundaries, sourceTerm already 0, so can ignore
        class(GSSolver), intent(in out) :: self
        integer :: N_indx, E_indx, S_indx, W_indx, k, i, j
        real(real64) :: C_N, C_E, C_O, C_W, C_S
        
        !$OMP parallel private(i, j, N_indx, E_indx, S_indx, W_indx, C_N, C_E, C_O, C_W, C_S)
        !$OMP do
        do k = 1, self%numberBlackInnerNodes
            i = self%black_InnerIndx(1, k)
            j = self%black_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerBlack(1, k)
            C_N = self%matCoeffsInnerBlack(2, k)
            C_E = self%matCoeffsInnerBlack(3, k)
            C_S = self%matCoeffsInnerBlack(4, k)
            C_W = self%matCoeffsInnerBlack(5, k)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i = self%black_NESW_BoundIndx(1, k)
            j = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundBlack(1, k)
            C_N = self%matCoeffsBoundBlack(2, k)
            C_E = self%matCoeffsBoundBlack(3, k)
            C_S = self%matCoeffsBoundBlack(4, k)
            C_W = self%matCoeffsBoundBlack(5, k)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedInnerNodes
            i = self%red_InnerIndx(1, k)
            j = self%red_InnerIndx(2, k)
            N_indx = j + 1
            E_indx = i+1
            S_indx = j-1
            W_indx = i-1
            C_O = self%matCoeffsInnerRed(1, k)
            C_N = self%matCoeffsInnerRed(2, k)
            C_E = self%matCoeffsInnerRed(3, k)
            C_S = self%matCoeffsInnerRed(4, k)
            C_W = self%matCoeffsInnerRed(5, k)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            i = self%red_NESW_BoundIndx(1, k)
            j = self%red_NESW_BoundIndx(2, k)
            E_indx = self%red_NESW_BoundIndx(3, k)
            W_indx = self%red_NESW_BoundIndx(4, k)
            N_indx = self%red_NESW_BoundIndx(5, k)
            S_indx = self%red_NESW_BoundIndx(6, k)
            C_O = self%matCoeffsBoundRed(1, k)
            C_N = self%matCoeffsBoundRed(2, k)
            C_E = self%matCoeffsBoundRed(3, k)
            C_S = self%matCoeffsBoundRed(4, k)
            C_W = self%matCoeffsBoundRed(5, k)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine calcResidual


end module mod_GSSolver