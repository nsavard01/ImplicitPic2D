module mod_GS_Base
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
    public :: GS_Base

    type :: GS_Base
        ! store grid quantities
        real(real64), allocatable :: sourceTerm(:,:), solution(:,:), residual(:,:)
        integer(int32), allocatable :: black_NESW_BoundIndx(:,:), red_NESW_BoundIndx(:,:)
        ! For Bound index, have [i, j, i_East, i_West, j_North, j_South], basically direction of each NESW
        integer(int32) :: numberBlackBoundNodes, numberRedBoundNodes, iterNumber, N_x, N_y
        real(real64) :: omega
        real(real64), allocatable :: matCoeffs(:,:,:)! Coefficients in order of O, N, E, S, W
        ! Note that saving coefficients, while taking more memory, results in much faster (more than 50%) speedup. So assuming we are more limited by calculation time than memory, best to store
        ! everything beforehand
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: solveGS
        procedure, public, pass(self) :: smoothIterations
        procedure, public, pass(self) :: smoothWithRes
        procedure, public, pass(self) :: calcResidual
        procedure, public, pass(self) :: restriction
        procedure, public, pass(self) :: prolongation
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
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
        integer :: i, j, k, tempNumBoundRed, tempNumBoundBlack, N_indx, S_indx, E_indx, W_indx
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

        tempNumBoundRed = 0
        tempNumBoundBlack = 0
        ! Count number red and black nodes
        ! Also fill in any dirichlet boundaries solution vector
        !$OMP parallel reduction(+:tempNumBoundRed, tempNumBoundBlack)
        !$OMP do collapse(2)
        do j = 1, self%N_y
            do i = 1, self%N_x
                if (MOD(i,2) == MOD(j,2)) then
                    ! Black node
                    if (boundaryConditions(i,j) == 0) then
                        continue
                    else if (boundaryConditions(i,j) /= 1) then
                        tempNumBoundBlack = tempNumBoundBlack + 1
                    end if
                else
                    !red node
                    if (boundaryConditions(i,j) == 0) then
                        continue
                    else if (boundaryConditions(i,j) /= 1) then
                        tempNumBoundRed = tempNumBoundRed + 1
                    end if
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        self%numberBlackBoundNodes = tempNumBoundBlack
        self%numberRedBoundNodes = tempNumBoundRed
        


        allocate(self%red_NESW_BoundIndx(6, self%numberRedBoundNodes), &
        self%black_NESW_BoundIndx(6, self%numberBlackBoundNodes), self%matCoeffs(5, self%N_x, self%N_y))
       

        tempNumBoundRed = 0
        tempNumBoundBlack = 0
        ! Fill arrays 
        do j = 1, self%N_y
            do i = 1, self%N_x
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
                else
                    ! Initial indices
                    N_indx = j+1
                    S_indx = j-1
                    E_indx = i+1
                    W_indx = i-1
                    if (boundaryConditions(i,j) <= 2) then
                        !Neumann boundary or dirichlet, have reflected grid
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
                    if (boundaryConditions(i,j) /=1) then
                        if (MOD(i,2) == MOD(j,2)) then
                            ! Black node
                            tempNumBoundBlack = tempNumBoundBlack + 1
                            self%black_NESW_BoundIndx(:, tempNumBoundBlack) = [i,j,E_indx, W_indx, N_indx, S_indx]
                        else
                            !red node
                            tempNumBoundRed = tempNumBoundRed + 1
                            self%red_NESW_BoundIndx(:, tempNumBoundRed) = [i,j,E_indx, W_indx, N_indx, S_indx]
                        end if
                    end if
                end if
                self%matCoeffs(:, i, j) = [C_O, C_N, C_E, C_S, C_W]
            end do
        end do
        
    end subroutine constructPoissonOrthogonal


    ! subroutine matMult(self, x, y)
    !     ! Use gauss-seidel to calculate A*x, store solution into y
    !     class(GSSolver), intent(in out) :: self
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
    !     class(GSSolver), intent(in out) :: self
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

    subroutine solveGS(self, tol)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        real(real64), intent(in) :: tol
        real(real64) :: Res
        integer(int32) :: i
        Res = 1.0
        self%iterNumber = 0
        do while (Res > tol)
            ! single iterations slightly faster
            call self%smoothIterations(100)
            Res = self%smoothWithRes()
            self%iterNumber = self%iterNumber + 101
        end do



    end subroutine solveGS

    subroutine smoothIterations(self, iterNum)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        real(real64) :: oldSol, C_N, C_E, C_O, C_W, C_S
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, iter

        do iter = 1, iterNum
            !$OMP parallel private(oldSol, i, j, N_indx, E_indx, S_indx, &
            !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W)
            ! loop through inner black nodes
            !$OMP do collapse(2)
            do j = 3, self%N_y-2, 2
                do i = 3, self%N_x-2, 2
                    N_indx = j + 1
                    E_indx = i+1
                    S_indx = j-1
                    W_indx = i-1
                    C_O = self%matCoeffs(1, i, j)
                    C_N = self%matCoeffs(2, i, j)
                    C_E = self%matCoeffs(3, i, j)
                    C_S = self%matCoeffs(4, i, j)
                    C_W = self%matCoeffs(5, i, j)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                        self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                end do
            end do
            !$OMP end do nowait
            !$OMP do collapse(2)
            do j = 2, self%N_y-1, 2
                do i = 2, self%N_x-1, 2
                    N_indx = j + 1
                    E_indx = i+1
                    S_indx = j-1
                    W_indx = i-1
                    C_O = self%matCoeffs(1, i, j)
                    C_N = self%matCoeffs(2, i, j)
                    C_E = self%matCoeffs(3, i, j)
                    C_S = self%matCoeffs(4, i, j)
                    C_W = self%matCoeffs(5, i, j)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                        self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                end do
            end do
            !$OMP end do nowait
            ! Bound black nodes
            !$OMP do
            do k = 1, self%numberBlackBoundNodes
                i = self%black_NESW_BoundIndx(1, k)
                j = self%black_NESW_BoundIndx(2, k)
                E_indx = self%black_NESW_BoundIndx(3, k)
                W_indx = self%black_NESW_BoundIndx(4, k)
                N_indx = self%black_NESW_BoundIndx(5, k)
                S_indx = self%black_NESW_BoundIndx(6, k)
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            end do
            !$OMP end do
            ! loop through red nodes
            !$OMP do collapse(2)
            do j = 3, self%N_y-2, 2
                do i = 2, self%N_x-1, 2
                    N_indx = j + 1
                    E_indx = i+1
                    S_indx = j-1
                    W_indx = i-1
                    C_O = self%matCoeffs(1, i, j)
                    C_N = self%matCoeffs(2, i, j)
                    C_E = self%matCoeffs(3, i, j)
                    C_S = self%matCoeffs(4, i, j)
                    C_W = self%matCoeffs(5, i, j)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                        self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                end do
            end do
            !$OMP end do nowait
            !$OMP do collapse(2)
            do j = 2, self%N_y-1, 2
                do i = 3, self%N_x-2, 2
                    N_indx = j + 1
                    E_indx = i+1
                    S_indx = j-1
                    W_indx = i-1
                    C_O = self%matCoeffs(1, i, j)
                    C_N = self%matCoeffs(2, i, j)
                    C_E = self%matCoeffs(3, i, j)
                    C_S = self%matCoeffs(4, i, j)
                    C_W = self%matCoeffs(5, i, j)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                        self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                end do
            end do
            !$OMP end do nowait
            ! loop bound red nodes
            !$OMP do
            do k = 1, self%numberRedBoundNodes
                i = self%red_NESW_BoundIndx(1, k)
                j = self%red_NESW_BoundIndx(2, k)
                E_indx = self%red_NESW_BoundIndx(3, k)
                W_indx = self%red_NESW_BoundIndx(4, k)
                N_indx = self%red_NESW_BoundIndx(5, k)
                S_indx = self%red_NESW_BoundIndx(6, k)
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            end do
            !$OMP end do
            !$OMP end parallel
        end do
    end subroutine smoothIterations

    function smoothWithRes(self) result(Res)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        real(real64) :: oldSol, C_N, C_E, C_O, C_W, C_S, Res
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k
        Res = 0.0d0
        !$OMP parallel private(oldSol, i, j, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W) reduction(+:Res)
        ! loop through inner black nodes
        !$OMP do collapse(2)
        do j = 3, self%N_y-2, 2
            do i = 3, self%N_x-2, 2
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        do j = 2, self%N_y-1, 2
            do i = 2, self%N_x-1, 2
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
        end do
        !$OMP end do nowait
        ! Bound black nodes
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i = self%black_NESW_BoundIndx(1, k)
            j = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            C_O = self%matCoeffs(1, i, j)
            C_N = self%matCoeffs(2, i, j)
            C_E = self%matCoeffs(3, i, j)
            C_S = self%matCoeffs(4, i, j)
            C_W = self%matCoeffs(5, i, j)
            oldSol = self%solution(i,j)
            self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(i,j) - oldSol)**2
        end do
        !$OMP end do
        ! loop through red nodes
        !$OMP do collapse(2)
        do j = 3, self%N_y-2, 2
            do i = 2, self%N_x-1, 2
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        do j = 2, self%N_y-1, 2
            do i = 3, self%N_x-2, 2
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
        end do
        !$OMP end do nowait
        ! loop bound red nodes
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            i = self%red_NESW_BoundIndx(1, k)
            j = self%red_NESW_BoundIndx(2, k)
            E_indx = self%red_NESW_BoundIndx(3, k)
            W_indx = self%red_NESW_BoundIndx(4, k)
            N_indx = self%red_NESW_BoundIndx(5, k)
            S_indx = self%red_NESW_BoundIndx(6, k)
            C_O = self%matCoeffs(1, i, j)
            C_N = self%matCoeffs(2, i, j)
            C_E = self%matCoeffs(3, i, j)
            C_S = self%matCoeffs(4, i, j)
            C_W = self%matCoeffs(5, i, j)
            oldSol = self%solution(i,j)
            self%solution(i,j) = (self%sourceTerm(i,j) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S  - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * C_O * self%omega + (1.0d0 - self%omega) * oldSol
            Res = Res + (self%solution(i,j) - oldSol)**2
        end do
        !$OMP end do
        !$OMP end parallel
        Res = SQRT(Res/ (self%numberBlackBoundNodes + self%numberRedBoundNodes + (self%N_y-2) * (self%N_x-2)))

    end function smoothWithRes


    subroutine calcResidual(self)
        ! Solve GS down to some tolerance
        class(GSSolver), intent(in out) :: self
        real(real64) :: C_N, C_E, C_O, C_W, C_S
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, iter

      
        !$OMP parallel private(i, j, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W)
        ! loop through inner nodes
        !$OMP do collapse(2)
        do j = 2, self%N_y-1
            do i = 2, self%N_x-1
                N_indx = j + 1
                E_indx = i+1
                S_indx = j-1
                W_indx = i-1
                C_O = self%matCoeffs(1, i, j)
                C_N = self%matCoeffs(2, i, j)
                C_E = self%matCoeffs(3, i, j)
                C_S = self%matCoeffs(4, i, j)
                C_W = self%matCoeffs(5, i, j)
                self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                    - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
            end do
        end do
        !$OMP end do nowait
        ! Bound black nodes
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i = self%black_NESW_BoundIndx(1, k)
            j = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            C_O = self%matCoeffs(1, i, j)
            C_N = self%matCoeffs(2, i, j)
            C_E = self%matCoeffs(3, i, j)
            C_S = self%matCoeffs(4, i, j)
            C_W = self%matCoeffs(5, i, j)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                    - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
        !$OMP end do nowait
        ! ! bound red nodes
        !$OMP do
        do k = 1, self%numberRedBoundNodes
            i = self%red_NESW_BoundIndx(1, k)
            j = self%red_NESW_BoundIndx(2, k)
            E_indx = self%red_NESW_BoundIndx(3, k)
            W_indx = self%red_NESW_BoundIndx(4, k)
            N_indx = self%red_NESW_BoundIndx(5, k)
            S_indx = self%red_NESW_BoundIndx(6, k)
            C_O = self%matCoeffs(1, i, j)
            C_N = self%matCoeffs(2, i, j)
            C_E = self%matCoeffs(3, i, j)
            C_S = self%matCoeffs(4, i, j)
            C_W = self%matCoeffs(5, i, j)
            self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
        end do
        !$OMP end do
        !$OMP end parallel
    end subroutine calcResidual

    subroutine restriction(self, fineGrid, coarseGrid)
        ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
        class(GSSolver), intent(in) :: self
        real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
        real(real64), intent(out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
        real(real64) :: a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S
        integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k

        !$OMP parallel private(k, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, &
        !$OMP& a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S)
        !$OMP do collapse(2)
        ! loop over fine grid nodes which overlap with coarse grid nodes
        do j_fine = 3, self%N_y-2, 2
            do i_fine = 3, self%N_x-2, 2
                ! calculate fine indices around overlapping index
                i_coarse = (i_fine + 1)/2
                j_coarse = (j_fine + 1)/2
                N_indx = j_fine+1
                E_indx = i_fine+1
                S_indx = j_fine-1
                W_indx = i_fine-1

                ! Transpose of prolongation, seems to work best
                ! interpolation from north point
                C_N = self%matCoeffs(2, i_fine, N_indx)
                C_S = self%matCoeffs(4, i_fine, N_indx)
                a_N = C_S/(C_N + C_S)

                ! interpolation from south point
                C_N = self%matCoeffs(2, i_fine, S_indx)
                C_S = self%matCoeffs(4, i_fine, S_indx)
                a_S = C_N/(C_N + C_S)

                ! interpolation from east point
                C_E = self%matCoeffs(3, E_indx, j_fine)
                C_W = self%matCoeffs(5, E_indx, j_fine)
                a_E = C_W/(C_E + C_W)

                ! interpolation from west point
                C_E = self%matCoeffs(3, W_indx, j_fine)
                C_W = self%matCoeffs(5, W_indx, j_fine)
                a_W = C_E/(C_E + C_W)

                ! node based restriction
                ! C_N = self%matCoeffs(2, i_fine, j_fine)
                ! C_E = self%matCoeffs(3, i_fine, j_fine)
                ! C_S = self%matCoeffs(4, i_fine, j_fine)
                ! C_W = self%matCoeffs(5, i_fine, j_fine)
                ! a_N = C_N/(C_N + C_S)
                ! a_S = C_S/(C_N + C_S)
                ! a_E = C_E/(C_E + C_W)
                ! a_W = C_W/(C_E + C_W)


                ! Interpolate residual to coarse grid
                coarseGrid(i_coarse, j_coarse) = 0.25d0 * (fineGrid(i_fine, j_fine) + &
                fineGrid(E_indx, j_fine) * a_E + fineGrid(W_indx, j_fine) * a_W +  &
                fineGrid(i_fine, N_indx) * a_N + fineGrid(i_fine, S_indx) * a_S + &
                fineGrid(E_indx, N_indx) * a_E * a_N + fineGrid(E_indx, S_indx) * a_E * a_S + &
                fineGrid(W_indx, N_indx) * a_W * a_N + fineGrid(W_indx, S_indx) * a_W * a_S)

                ! Simple interpolation can be surprisingly good
                ! coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
                ! 0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
                ! fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
                ! 0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
                ! fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
            end do
        end do
        !$OMP end do nowait
        ! each black bound node is a coarse node
        !$OMP do
        do k = 1, self%numberBlackBoundNodes
            i_fine = self%black_NESW_BoundIndx(1, k)
            j_fine = self%black_NESW_BoundIndx(2, k)
            E_indx = self%black_NESW_BoundIndx(3, k)
            W_indx = self%black_NESW_BoundIndx(4, k)
            N_indx = self%black_NESW_BoundIndx(5, k)
            S_indx = self%black_NESW_BoundIndx(6, k)
            
            ! calculate fine indices around overlapping index
            i_coarse = (i_fine + 1)/2
            j_coarse = (j_fine + 1)/2
            ! From transpose of prolongation, works fairly well
            if (N_indx == S_indx) then
                ! neumann boundary
                if (j_fine == 1) then
                    ! bottom boundary
                    C_N = self%matCoeffs(2, i_fine, N_indx)
                    C_S = self%matCoeffs(4, i_fine, N_indx)
                    a_N = C_S/(C_N + C_S)
                    a_S = a_N
                else
                    ! top boundary
                    C_N = self%matCoeffs(2, i_fine, S_indx)
                    C_S = self%matCoeffs(4, i_fine, S_indx)
                    a_S = C_N/(C_N + C_S)
                    a_N = a_S
                end if
            else
                ! interpolation from north point
                C_N = self%matCoeffs(2, i_fine, N_indx)
                C_S = self%matCoeffs(4, i_fine, N_indx)
                a_N = C_S/(C_N + C_S)

                ! interpolation from south point
                C_N = self%matCoeffs(2, i_fine, S_indx)
                C_S = self%matCoeffs(4, i_fine, S_indx)
                a_S = C_N/(C_N + C_S)
            end if
            
            if (E_indx == W_indx) then
                ! neumann boundary
                if (i_fine == 1) then
                    ! left boundary
                    C_E = self%matCoeffs(3, E_indx, j_fine)
                    C_W = self%matCoeffs(5, E_indx, j_fine)
                    a_E = C_W/(C_E + C_W)
                    a_W = a_E
                else
                    C_E = self%matCoeffs(3, W_indx, j_fine)
                    C_W = self%matCoeffs(5, W_indx, j_fine)
                    a_W = C_E/(C_E + C_W)
                    a_E = a_W
                end if
            else
                ! interpolation from east point
                C_E = self%matCoeffs(3, E_indx, j_fine)
                C_W = self%matCoeffs(5, E_indx, j_fine)
                a_E = C_W/(C_E + C_W)

                ! interpolation from west point
                C_E = self%matCoeffs(3, W_indx, j_fine)
                C_W = self%matCoeffs(5, W_indx, j_fine)
                a_W = C_E/(C_E + C_W)
            end if

            ! node based restrition
            ! C_N = self%matCoeffs(2, i_fine, j_fine)
            ! C_E = self%matCoeffs(3, i_fine, j_fine)
            ! C_S = self%matCoeffs(4, i_fine, j_fine)
            ! C_W = self%matCoeffs(5, i_fine, j_fine)
            ! a_N = C_N/(C_N + C_S)
            ! a_S = C_S/(C_N + C_S)
            ! a_E = C_E/(C_E + C_W)
            ! a_W = C_W/(C_E + C_W)

            ! Interpolate residual to coarse grid
            coarseGrid(i_coarse, j_coarse) = 0.25d0 * (fineGrid(i_fine, j_fine) + &
            fineGrid(E_indx, j_fine) * a_E + fineGrid(W_indx, j_fine) * a_W +  &
            fineGrid(i_fine, N_indx) * a_N + fineGrid(i_fine, S_indx) * a_S + &
            fineGrid(E_indx, N_indx) * a_E * a_N + fineGrid(E_indx, S_indx) * a_E * a_S + &
            fineGrid(W_indx, N_indx) * a_W * a_N + fineGrid(W_indx, S_indx) * a_W * a_S)

            ! simple binomial filter, actually ok
            ! coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
            !     0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
            !     fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
            !     0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
            !     fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine restriction

    subroutine prolongation(self, fineGrid, coarseGrid)
        ! Prolongate operator from coarse to fine grid using gauss seidel data
        class(GSSolver), intent(in) :: self
        real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
        real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
        real(real64) :: w_1, w_2, w_3, w_4, w_tot
        integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse
        N_x_coarse = (self%N_x+1)/2
        N_y_coarse = (self%N_y+1)/2
        ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
        !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine, w_1, w_2, w_3, w_4, w_tot)
        ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
        !$OMP do collapse(2)
        ! set similar nodes
        do j_coarse = 1, N_y_coarse
            do i_coarse = 1, N_x_coarse
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2-1
                ! Add overlapping coarse nodes
                fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! go horizontal each row,
        do j_coarse = 1, N_y_coarse
            do i_coarse = 1, N_x_coarse-1
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2 - 1

                w_1 = self%matCoeffs(3, i_fine+1, j_fine) !E
                w_2 = self%matCoeffs(5, i_fine+1, j_fine) !W
                w_tot = w_1 + w_2
                w_1 = w_1 / w_tot
                w_2 = w_2 / w_tot
    
                ! Average fine nodes between coarse nodes horizontally
                fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
                ! simple interpolation
                ! fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! go vertical
        do j_coarse = 1, N_y_coarse-1
            do i_coarse = 1, N_x_coarse
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2 - 1
                w_1 = self%matCoeffs(2, i_fine, j_fine+1) !N
                w_2 = self%matCoeffs(4, i_fine, j_fine+1) !S
                w_tot = w_1 + w_2
                w_1 = w_1 / w_tot
                w_2 = w_2 / w_tot
                ! Average fine nodes between coarse nodes vertical direction
                fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
                ! Simple interpolation
                ! fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! add offset fine nodes which require 4 point interpolation
        ! first loop through SW corner coarse nodes
        do j_coarse = 1, N_y_coarse-1
            do i_coarse = 1, N_x_coarse-1
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2 - 1
                w_1 = self%matCoeffs(2, i_fine+1, j_fine+1) !N
                w_2 = self%matCoeffs(4, i_fine+1, j_fine+1) !S
                w_tot = w_1 + w_2
                w_1 = w_1/w_tot
                w_2 = w_2/w_tot

                w_3 = self%matCoeffs(3, i_fine+1, j_fine+1) !E
                w_4 = self%matCoeffs(5, i_fine+1, j_fine+1) !W
                w_tot = w_3 + w_4
                w_3 = w_3/w_tot
                w_4 = w_4/w_tot

                fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * w_2*w_4 + coarseGrid(i_coarse+1, j_coarse) * w_3 * w_2 &
                    + coarseGrid(i_coarse, j_coarse+1) * w_1 * w_4 + coarseGrid(i_coarse+1, j_coarse+1) * w_1 * w_3

                ! simple interpolation
                ! fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse) * 0.25d0 &
                !     + coarseGrid(i_coarse, j_coarse+1) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse+1) * 0.25d0
            end do
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine prolongation


end module mod_GSSolver