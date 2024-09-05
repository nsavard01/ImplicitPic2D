module mod_ZebraSolverEven
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base_Even
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
    public :: ZebraSolverEven

    type, extends(GS_Base_Even) :: ZebraSolverEven
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal => constructPoissonOrthogonal_ZebraEven
        procedure, public, pass(self) :: solveGS => solveGS_ZebraEven
        procedure, public, pass(self) :: smoothIterations => smoothIterations_ZebraEven
        procedure, public, pass(self) :: smoothWithRes => smoothWithRes_ZebraEven
        procedure, public, pass(self) :: calcResidual => calcResidual_ZebraEven
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type
    interface ZebraSolverEven
        module procedure :: ZebraSolverEven_constructor
    end interface ZebraSolverEven

contains
    type(ZebraSolverEven) function ZebraSolverEven_constructor(omega, N_x, N_y) result(self)
        ! Construct object, set initial variables
        real(real64), intent(in) :: omega
        integer(int32), intent(in) :: N_x, N_y
        self%omega = omega
        self%iterNumber = 0
        self%Residual = 0
        self%N_x = N_x
        self%N_y = N_y
        self%numberRows = 0
        self%numberColumns = 0
        self%numberBoundNodes = 0
        self%startRow = 0
        self%endRow = 0
        self%startCol = 0
        self%endCol = 0
        self%iterNumber = 0
        allocate(self%sourceTerm(self%N_x, self%N_y), self%solution(self%N_x, self%N_y), self%residual(self%N_x, self%N_y))
    end function ZebraSolverEven_constructor

    subroutine constructPoissonOrthogonal_ZebraEven(self, diffX, diffY, NESW_wallBoundaries, boundaryConditions)
        ! Construct orthogonal grid solver
        class(ZebraSolverEven), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
        real(real64), intent(in) :: diffX(self%N_x-1), diffY(self%N_y-1)
        integer :: upperBound, lowerBound, rightBound, leftBound
        integer :: i, j, k, numberBound

        upperBound = NESW_wallBoundaries(1)
        rightBound = NESW_wallBoundaries(2)
        lowerBound = NESW_wallBoundaries(3)
        leftBound = NESW_wallBoundaries(4)
        
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

        numberBound = 0
        ! Count number red and black nodes
        ! Also fill in any dirichlet boundaries solution vector
        !$OMP parallel reduction(+:numberBound)
        !$OMP do collapse(2)
        do j = 1, self%N_y
            do i = 1, self%N_x
                if (boundaryConditions(i,j) == 0) then
                    continue
                else if (boundaryConditions(i,j) /= 1) then
                    numberBound = numberBound + 1
                end if
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        self%numberBoundNodes = numberBound
        self%startCol = 1
        self%startRow = 1
        self%endRow = self%N_y
        self%endCol = self%N_x
        if (lowerBound == 1) self%startRow = 2
        if (upperBound == 1) self%endRow = self%N_y-1
        self%numberRows = self%endRow - self%startRow + 1
        if (leftBound == 1) self%startCol = 2
        if (rightBound == 1) self%endCol = self%N_x-1
        self%numberColumns = self%endCol - self%startCol + 1  
        

        allocate(self%horzIndx(2,self%numberColumns), self%vertIndx(2,self%numberRows))

        self%coeffX = 1.0d0 / (diffX(1)**2)
        self%coeffY = 1.0d0 / (diffY(1)**2)
        self%centerCoeff = -1.0d0 / (2.0d0 * self%coeffX + 2.0d0 * self%coeffY)

        if (self%startCol == 1) then
            ! Start at neumann or periodic boundary
            if (leftBound == 2) then
                self%horzIndx(1,1) = 2 !E
                self%horzIndx(2,1) = 2 !W
            else
                self%horzIndx(1,1) = 2 !E
                self%horzIndx(2,1) = self%N_x-1 !W
            end if
        else
            ! start inner node
            self%horzIndx(1,1) = 3 !E
            self%horzIndx(2,1) = 1 !W
        end if
        ! inner node coeffs
        do k = 2, self%numberColumns-1
            i = self%startCol + k - 1
            self%horzIndx(1,k) = i+1 !E
            self%horzIndx(2,k) = i-1 !W
        end do

        if (self%endCol == self%N_x) then
            ! end at neumann or periodic boundary
            if (rightBound == 2) then
                self%horzIndx(1,self%numberColumns) = self%N_x-1 !E
                self%horzIndx(2,self%numberColumns) = self%N_x-1 !W
            else
                self%horzIndx(1,self%numberColumns) = 2 !E
                self%horzIndx(2,self%numberColumns) = self%N_x-1 !W
            end if
        else
            ! end inner node
            self%horzIndx(1,self%numberColumns) =  self%N_x!E
            self%horzIndx(2,self%numberColumns) = self%N_x-2 !W
        end if

        ! --------------- Initialize vertical components ------------------------ 
        if (self%startRow == 1) then
            ! Start at neumann or periodic boundary
            if (lowerBound == 2) then
                self%vertIndx(1,1) = 2 !N
                self%vertIndx(2,1) = 2 !S
            else
                self%vertIndx(1,1) = 2 !N
                self%vertIndx(2,1) = self%N_y-1 !S
            end if
        else
            ! start inner node
            self%vertIndx(1,1) = 3 !N
            self%vertIndx(2,1) = 1 !S
        end if
        ! inner node coeffs
        do k = 2, self%numberRows-1
            j = self%startRow + k - 1
            self%vertIndx(1,k) = j+1 !N
            self%vertIndx(2,k) = j-1 !S
        end do

        if (self%endRow == self%N_y) then
            ! end at neumann or periodic boundary
            if (upperBound == 2) then
                self%vertIndx(1,self%numberRows) = self%N_y-1 !N
                self%vertIndx(2,self%numberRows) = self%N_y-1 !S
            else
                self%vertIndx(1,self%numberRows) = 2 !N
                self%vertIndx(2,self%numberRows) = self%N_y-1 !S
            end if
        else
            self%vertIndx(1,self%numberRows) =  self%N_y !N
            self%vertIndx(2,self%numberRows) = self%N_y-2 !S
        end if

        ! Get starting/end points in GS smoother rows and cols for overlapping 
        if (self%startRow == 1) then
            self%startRowCoarse = 1
        else
            self%startRowCoarse = 2
        end if

        if (self%endRow == self%N_y) then
            self%endRowCoarse = self%numberRows
        else
            self%endRowCoarse = self%numberRows-1
        end if

        if (self%startCol == 1) then
            self%startColCoarse = 1
        else
            self%startColCoarse = 2
        end if

        if (self%endCol == self%N_x) then
            self%endColCoarse = self%numberColumns
        else
            self%endColCoarse = self%numberColumns-1
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

    subroutine solveGS_ZebraEven(self, tol)
        ! Solve GS down to some tolerance
        class(ZebraSolverEven), intent(in out) :: self
        real(real64), intent(in) :: tol
        real(real64) :: Res
        Res = 1.0
        self%iterNumber = 0
        do while (Res > tol)
            ! single iterations slightly faster
            call self%smoothIterations(100)
            Res = self%smoothWithRes()
            self%iterNumber = self%iterNumber + 101
        end do



    end subroutine solveGS_ZebraEven

    subroutine smoothIterations_ZebraEven(self, iterNum)
        ! Solve GS down to some tolerance
        class(ZebraSolverEven), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        real(real64) :: oldSol, omega_inv
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, iter
        omega_inv = 1.0d0 - self%omega
        ! p indexs in horizontal maps to i, k in vertical maps to j
        do iter = 1, iterNum 
            !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol)

            ! horizontal sweeps left to right
            !$OMP do
            ! Sweep odd numbers
            do k = 1, self%numberRows, 2
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                do p = 1, self%numberColumns
                    i = self%startCol + p - 1
                    E_indx = self%horzIndx(1, p)
                    W_indx = self%horzIndx(2, p)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do
            !$OMP do
            ! sweep even numbers
            do k = 2, self%numberRows, 2
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                do p = 1, self%numberColumns
                    i = self%startCol + p - 1
                    E_indx = self%horzIndx(1, p)
                    W_indx = self%horzIndx(2, p)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do

            ! vertical sweeps bottom to top
            !$OMP do
            ! Sweep odd numbers
            do p = 1, self%numberColumns, 2
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                do k = 1, self%numberRows
                    j = self%startRow + k - 1
                    N_indx = self%vertIndx(1, k)
                    S_indx = self%vertIndx(2, k)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do
            !$OMP do
            ! sweep even numbers
            do p = 2, self%numberColumns, 2
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                do k = 1, self%numberRows
                    j = self%startRow + k - 1
                    N_indx = self%vertIndx(1, k)
                    S_indx = self%vertIndx(2, k)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do

            ! horizontal sweeps right to left
            !$OMP do
            ! Sweep odd numbers
            do k = 1, self%numberRows, 2
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                do p = self%numberColumns, 1, -1
                    i = self%startCol + p - 1
                    E_indx = self%horzIndx(1, p)
                    W_indx = self%horzIndx(2, p)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do
            !$OMP do
            ! sweep even numbers
            do k = 2, self%numberRows, 2
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                do p = self%numberColumns, 1, -1
                    i = self%startCol + p - 1
                    E_indx = self%horzIndx(1, p)
                    W_indx = self%horzIndx(2, p)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do

            ! vertical sweeps top to bottom
            !$OMP do
            ! Sweep odd numbers
            do p = 1, self%numberColumns, 2
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                do k = self%numberRows, 1, -1
                    j = self%startRow + k - 1
                    N_indx = self%vertIndx(1, k)
                    S_indx = self%vertIndx(2, k)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do
            !$OMP do
            ! sweep even numbers
            do p = 2, self%numberColumns, 2
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                do k = self%numberRows, 1, -1
                    j = self%startRow + k - 1
                    N_indx = self%vertIndx(1, k)
                    S_indx = self%vertIndx(2, k)
                    oldSol = self%solution(i,j)
                    self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                end do
            end do
            !$OMP end do
            !$OMP end parallel
        end do

    
    end subroutine smoothIterations_ZebraEven

    function smoothWithRes_ZebraEven(self) result(Res)
        ! Solve GS down to some tolerance
        class(ZebraSolverEven), intent(in out) :: self
        real(real64) :: oldSol, Res, omega_inv
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p
        omega_inv = 1.0d0 - self%omega
        Res = 0.0d0
        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, oldSol) reduction(+:Res)

        ! horizontal sweeps left to right
        !$OMP do
        ! Sweep odd numbers
        do k = 1, self%numberRows, 2
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            do p = 1, self%numberColumns
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end do
        end do
        !$OMP end do
        !$OMP do
        ! sweep even numbers
        do k = 2, self%numberRows, 2
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            do p = 1, self%numberColumns
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end do
        end do
        !$OMP end do

        ! vertical sweeps bottom to top
        !$OMP do
        ! Sweep odd numbers
        do p = 1, self%numberColumns, 2
            i = self%startCol + p - 1
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)
            do k = 1, self%numberRows
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end do
        end do
        !$OMP end do
        !$OMP do
        ! sweep even numbers
        do p = 2, self%numberColumns, 2
            i = self%startCol + p - 1
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)
            do k = 1, self%numberRows
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end do
        end do
        !$OMP end do

        ! horizontal sweeps right to left
        !$OMP do
        ! Sweep odd numbers
        do k = 1, self%numberRows, 2
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            do p = self%numberColumns, 1, -1
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end do
        end do
        !$OMP end do
        !$OMP do
        ! sweep even numbers
        do k = 2, self%numberRows, 2
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            do p = self%numberColumns, 1, -1
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
            end do
        end do
        !$OMP end do

        ! vertical sweeps top to bottom
        !$OMP do
        ! Sweep odd numbers
        do p = 1, self%numberColumns, 2
            i = self%startCol + p - 1
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)
            do k = self%numberRows, 1, -1
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
        end do
        !$OMP end do
        !$OMP do
        ! sweep even numbers
        do p = 2, self%numberColumns, 2
            i = self%startCol + p - 1
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)
            do k = self%numberRows, 1, -1
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                oldSol = self%solution(i,j)
                self%solution(i,j) = (self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                    (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX) * self%centerCoeff * self%omega + omega_inv * oldSol
                Res = Res + (self%solution(i,j) - oldSol)**2
            end do
        end do
        !$OMP end do
        !$OMP end parallel
        
        Res = SQRT(Res/ (self%numberColumns * self%numberRows))

    end function smoothWithRes_ZebraEven


    subroutine calcResidual_ZebraEven(self)
        ! Solve GS down to some tolerance
        class(ZebraSolverEven), intent(in out) :: self
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p

      
        !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx)
        ! loop through inner nodes
        !$OMP do collapse(2)
        do k = 1, self%numberRows
            do p = 1, self%numberColumns
                i = self%startCol + p - 1
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1,k)
                E_indx = self%horzIndx(1,p)
                S_indx = self%vertIndx(2,k)
                W_indx = self%horzIndx(2,p)
                self%residual(i, j) = self%sourceTerm(i,j) - (self%solution(i,N_indx) + self%solution(i,S_indx)) * self%coeffY &
                    - (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX - self%solution(i,j)/self%centerCoeff
            end do
        end do
        !$OMP end do nowait
        !$OMP end parallel
    end subroutine calcResidual_ZebraEven


end module mod_ZebraSolverEven