module mod_GSSolver
    use iso_fortran_env, only: int32, int64, real64
    use ieee_arithmetic
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
        real(real64), allocatable :: sourceTerm(:), solution(:)
        integer(int32), allocatable :: black_NESW_indx(:,:), red_NESW_indx(:,:)
        integer(int32) :: numberBlackNodes, numberRedNodes, matDimension, iterNumber
        real(real64) :: Residual, coeffX, coeffY, coeffSelf, omega
    contains
        procedure, public, pass(self) :: constructPoissonEven
        procedure, public, pass(self) :: solveGS
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

    subroutine constructPoissonEven(self, N_x, N_y, delX, delY, NESW_wallBoundaries, NESW_phiValues, boundaryConditions)
        class(GSSolver), intent(in out) :: self
        integer(int32), intent(in) :: N_x, N_y, NESW_wallBoundaries(4), boundaryConditions(N_x*N_y)
        real(real64), intent(in) :: NESW_phiValues(4), delX, delY
        integer :: matDimension, upperBound, lowerBound, rightBound, leftBound
        real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi
        integer :: i, j, k, tempNumRed, tempNumBlack
        matDimension = N_x * N_y
        self%matDimension = matDimension

        upperBound = NESW_wallBoundaries(1)
        rightBound = NESW_wallBoundaries(2)
        lowerBound = NESW_wallBoundaries(3)
        leftBound = NESW_wallBoundaries(4)

        upperPhi = NESW_phiValues(1)
        rightPhi = NESW_phiValues(2)
        lowerphi = NESW_phiValues(3)
        leftPhi = NESW_phiValues(4)

        self%coeffY = 1.0d0 / delY**2
        self%coeffX = 1.0d0 / delX**2
        self%coeffSelf = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
        self%coeffSelf = 1.0d0 / self%coeffSelf

        allocate(self%sourceTerm(matDimension), self%solution(matDimension))
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

        ! fill in Dirichlet solution
        self%solution = 0.0d0
        if (upperBound == 1)self%solution(matDimension-N_x+2:matDimension-1) = upperPhi
        if (rightBound == 1) self%solution(2*N_x:matDimension-N_x:N_x) = rightPhi
        if (lowerBound == 1) self%solution(2:N_x-1) = lowerPhi
        if (leftBound == 1) self%solution(N_x+1:matDimension-2*N_x + 1:N_x) = leftPhi
        if (boundaryConditions(1) == 1) then
            if (lowerBound == leftBound) then
                self%solution(1) = MIN(leftPhi, lowerPhi)
            else if (lowerBound == 1) then
                self%solution(1) = lowerPhi
            else 
                self%solution(1) = leftBound
            end if
        end if
        if (boundaryConditions(N_x) == 1) then
            if (lowerBound == rightBound) then
                self%solution(N_x) = MIN(rightPhi, lowerPhi)
            else if (lowerBound == 1) then
                self%solution(N_x) = lowerPhi
            else
                self%solution(N_x) = rightPhi
            end if
        end if
        if (boundaryConditions(matDimension-N_x+1) == 1) then
            if (upperBound == leftBound) then
                self%solution(matDimension-N_x+1) = MIN(leftPhi, upperPhi)
            else if (upperBound == 1) then
                self%solution(matDimension-N_x+1) = upperPhi
            else
                self%solution(matDimension-N_x+1) = leftPhi
            end if
        end if
        if (boundaryConditions(matDimension) == 1) then
            if (upperBound == rightBound) then
                self%solution(matDimension) = MIN(rightPhi, upperPhi)
            else if (upperBound == 1) then
                self%solution(matDimension) = upperPhi
            else
                self%solution(matDimension) = rightPhi
            end if
        end if

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
            do k = 1, self%numberBlackNodes
                O_indx = self%black_NESW_indx(1,k)
                N_indx = self%black_NESW_indx(2,k)
                E_indx = self%black_NESW_indx(3,k)
                S_indx = self%black_NESW_indx(4,k)
                W_indx = self%black_NESW_indx(5,k)
                oldSol = self%solution(O_indx)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(O_indx) - oldSol)**2
            end do
            !$OMP end do
            !$OMP do
            do k = 1, self%numberRedNodes
                O_indx = self%red_NESW_indx(1,k)
                N_indx = self%red_NESW_indx(2,k)
                E_indx = self%red_NESW_indx(3,k)
                S_indx = self%red_NESW_indx(4,k)
                W_indx = self%red_NESW_indx(5,k)
                oldSol = self%solution(O_indx)
                self%solution(O_indx) = (self%sourceTerm(O_indx) - (self%solution(N_indx) + self%solution(S_indx)) * self%coeffY - &
                    (self%solution(E_indx) + self%solution(W_indx)) * self%coeffX) * self%coeffSelf * self%omega + (1.0d0 - self%omega) * oldSol
                Res = Res + (self%solution(O_indx) - oldSol)**2
            end do
            !$OMP end do
            !$OMP end parallel
            Res = SQRT(Res/self%matDimension)
            self%iterNumber = self%iterNumber + 1
        end do
        print *, 'took', self%iterNumber, 'iterations'



    end subroutine solveGS


end module mod_GSSolver