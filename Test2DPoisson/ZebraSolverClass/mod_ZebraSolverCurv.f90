module mod_ZebraSolverCurv
    use iso_fortran_env, only: int32, int64, real64
    use mod_GS_Base_Curv
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
    public :: ZebraSolverCurv

    type, extends(GS_Base_Curv) :: ZebraSolverCurv
        ! store grid quantities
        integer(int32) :: leftBound, rightBound, lowerBound, upperBound
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal => constructPoissonOrthogonal_ZebraCurv
        procedure, public, pass(self) :: solveGS => solveGS_ZebraCurv
        procedure, public, pass(self) :: smoothIterations => smoothIterations_ZebraCurv
        procedure, public, pass(self) :: smoothWithRes => smoothWithRes_ZebraCurv
        procedure, public, pass(self) :: calcResidual => calcResidual_ZebraCurv
        ! procedure, public, pass(self) :: restriction
        ! procedure, public, pass(self) :: prolongation
        ! procedure, public, pass(self) :: matMult
        ! procedure, public, pass(self) :: XAX_Mult
    end type

    interface ZebraSolverCurv
        module procedure :: GSSolver_constructor
    end interface ZebraSolverCurv

contains

    type(ZebraSolverCurv) function GSSolver_constructor(omega, N_x, N_y) result(self)
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
    end function GSSolver_constructor

    subroutine constructPoissonOrthogonal_ZebraCurv(self, diffX, diffY, NESW_wallBoundaries, boundaryConditions)
        ! Construct orthogonal grid solver
        class(ZebraSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
        real(real64), intent(in) :: diffX(self%N_x-1), diffY(self%N_y-1)
        integer :: upperBound, lowerBound, rightBound, leftBound
        integer :: i, j, k, numberBound, N_indx, S_indx, E_indx, W_indx, i_2, j_2
        real(real64) :: delX_E, delX_W, delY_N, delY_S, tempReal, C_N, C_S, C_E, C_W, C_O

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
        self%upperBound = upperBound
        self%lowerBound = lowerBound
        self%rightBound = rightBound
        self%leftBound = leftBound
        if (self%numberColumns == self%N_x) then
            self%leftBound = 4
        end if
        if (self%numberRows == self%N_y) then
            self%lowerBound = 4
        end if

        allocate(self%horzIndx(2,self%numberColumns), self%vertIndx(2,self%numberRows))
        
        allocate(self%horzCoeffs(2,self%numberColumns), self%vertCoeffs(2,self%numberRows), &
        self%centerCoeffs(self%numberColumns,self%numberRows))
        !allocate(self%rowBoundCoeffs(2, 2, self%numberRows),self%colBoundCoeffs(2, 2, self%numberColumns))  
        

        ! --------------- Initialize horizontal components ------------------------ 
        if (self%startCol == 1) then
            ! Start at neumann or periodic boundary
            if (leftBound == 2) then
                delX_E = diffX(1)
                self%horzCoeffs(1, 1) = 1.0d0 / (delX_E**2) !E
                self%horzCoeffs(2, 1) = 1.0d0 / (delX_E**2) !W
                self%horzIndx(1,1) = 2 !E
                self%horzIndx(2,1) = 2 !W
            else
                delX_E = diffX(1)
                delX_W = diffX(self%N_x-1)
                self%horzCoeffs(1, 1) = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) !E
                self%horzCoeffs(2, 1) = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) !W
                self%horzIndx(1,1) = 2 !E
                self%horzIndx(2,1) = self%N_x-1 !W
            end if
        else
            ! start inner node
            delX_E = diffX(2)
            delX_W = diffX(1)
            self%horzCoeffs(1, 1) = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) !E
            self%horzCoeffs(2, 1) = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) !W
            self%horzIndx(1,1) = 3 !E
            self%horzIndx(2,1) = 1 !W
        end if
        ! inner node coeffs
        do k = 2, self%numberColumns-1
            i = self%startCol + k - 1
            delX_E = diffX(i)
            delX_W = diffX(i-1)
            self%horzCoeffs(1, k) = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) !E
            self%horzCoeffs(2, k) = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) !W 
            self%horzIndx(1,k) = i+1 !E
            self%horzIndx(2,k) = i-1 !W
        end do

        if (self%endCol == self%N_x) then
            ! end at neumann or periodic boundary
            if (rightBound == 2) then
                delX_W = diffX(self%N_x-1)
                self%horzCoeffs(1, self%numberColumns) = 1.0d0 / (delX_W**2) !E
                self%horzCoeffs(2, self%numberColumns) = 1.0d0 / (delX_W**2) !W
                self%horzIndx(1,self%numberColumns) = self%N_x-1 !E
                self%horzIndx(2,self%numberColumns) = self%N_x-1 !W
            else
                delX_E = diffX(1)
                delX_W = diffX(self%N_x-1)
                self%horzCoeffs(1, self%numberColumns) = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) !E
                self%horzCoeffs(2, self%numberColumns) = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) !W
                self%horzIndx(1,self%numberColumns) = 2 !E
                self%horzIndx(2,self%numberColumns) = self%N_x-1 !W
            end if
        else
            ! end inner node
            delX_E = diffX(self%N_x-1)
            delX_W = diffX(self%N_x-2)
            self%horzCoeffs(1, self%numberColumns) = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) !E
            self%horzCoeffs(2, self%numberColumns) = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) !W
            self%horzIndx(1,self%numberColumns) =  self%N_x!E
            self%horzIndx(2,self%numberColumns) = self%N_x-2 !W
        end if

        ! --------------- Initialize vertical components ------------------------ 
        if (self%startRow == 1) then
            ! Start at neumann or periodic boundary
            if (lowerBound == 2) then
                delY_N = diffY(1)
                self%vertCoeffs(1, 1) = 1.0d0 / (delY_N**2) !N
                self%vertCoeffs(2, 1) = 1.0d0 / (delY_N**2) !S
                self%vertIndx(1,1) = 2 !N
                self%vertIndx(2,1) = 2 !S
            else
                delY_N = diffY(1)
                delY_S = diffY(self%N_y-1)
                self%vertCoeffs(1, 1) = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) !N
                self%vertCoeffs(2, 1) = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) !S
                self%vertIndx(1,1) = 2 !N
                self%vertIndx(2,1) = self%N_y-1 !S
            end if
        else
            ! start inner node
            delY_N = diffY(2)
            delY_S = diffY(1)
            self%vertCoeffs(1, 1) = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) !N
            self%vertCoeffs(2, 1) = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) !S
            self%vertIndx(1,1) = 3 !N
            self%vertIndx(2,1) = 1 !S
        end if
        ! inner node coeffs
        do k = 2, self%numberRows-1
            j = self%startRow + k - 1
            delY_N = diffY(j)
            delY_S = diffY(j-1)
            self%vertCoeffs(1, k) = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) !N
            self%vertCoeffs(2, k) = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) !S
            self%vertIndx(1,k) = j+1 !N
            self%vertIndx(2,k) = j-1 !S
        end do

        if (self%endRow == self%N_y) then
            ! end at neumann or periodic boundary
            if (upperBound == 2) then
                delY_S = diffY(self%N_y-1)
                self%vertCoeffs(1, self%numberRows) = 1.0d0 / (delY_S**2) !N
                self%vertCoeffs(2, self%numberRows) = 1.0d0 / (delY_S**2) !S
                self%vertIndx(1,self%numberRows) = self%N_y-1 !N
                self%vertIndx(2,self%numberRows) = self%N_y-1 !S
            else
                delY_N = diffY(1)
                delY_S = diffY(self%N_y-1)
                self%vertCoeffs(1, self%numberRows) = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) !N
                self%vertCoeffs(2, self%numberRows) = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) !S
                self%vertIndx(1,self%numberRows) = 2 !N
                self%vertIndx(2,self%numberRows) = self%N_y-1 !S
            end if
        else
            ! end inner node
            delY_N = diffY(self%N_y-1)
            delY_S = diffY(self%N_y-2)
            self%vertCoeffs(1, self%numberRows) = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) !N
            self%vertCoeffs(2, self%numberRows) = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) !S
            self%vertIndx(1,self%numberRows) =  self%N_y !N
            self%vertIndx(2,self%numberRows) = self%N_y-2 !S
        end if

        
        ! save center coefficent in memory to avoid repeated calculations of it
        do j = 1, self%numberRows
            do i = 1, self%numberColumns
                ! convert to index solved nodes
                self%centerCoeffs(i, j) = -1.0d0 / (self%vertCoeffs(1,j) + self%vertCoeffs(2,j) &
                + self%horzCoeffs(1,i) + self%horzCoeffs(2,i))
            end do
        end do

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

       
        
    end subroutine constructPoissonOrthogonal_ZebraCurv


    ! subroutine matMult(self, x, y)
    !     ! Use gauss-seidel to calculate A*x, store solution into y
    !     class(ZebraSolver), intent(in out) :: self
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

    subroutine solveGS_ZebraCurv(self, tol)
        ! Solve GS down to some tolerance
        class(ZebraSolverCurv), intent(in out) :: self
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



    end subroutine solveGS_ZebraCurv

    subroutine smoothIterations_ZebraCurv(self, iterNum)
        ! Solve GS down to some tolerance
        class(ZebraSolverCurv), intent(in out) :: self
        integer(int32), intent(in) :: iterNum
        real(real64) :: omega_inv, C_N, C_E, C_S, C_W
        real(real64) :: cp_x(self%N_x-1), dp_x(self%N_x-1), cp_y(self%N_y-1), dp_y(self%N_y), sourceX(self%N_x), sourceY(self%N_y), m
        integer :: k,p, iter, E_indx, W_indx, N_indx, S_indx, i, j
        omega_inv = 1.0d0 - self%omega
        

        do iter = 1, iterNum
            !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, &
            !$OMP& C_N, C_S, C_E, C_W, cp_y, cp_x, dp_y, dp_x, sourceX, sourceY, m)

            ! horizontal sweeps left to right
            !$OMP do
            ! Sweep odd numbers forward
            do k = 1, self%numberRows, 2
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                C_N = self%vertCoeffs(1, k)
                C_S = self%vertCoeffs(2, k)
                ! initialize c-prime and d-prime
                select case (self%leftBound)
                case (1)
                    cp_x(1) = 0.0d0
                    sourceX(1) = self%solution(1,j)
                    dp_x(1) = sourceX(1)
                case (2)
                    C_E = 2.0d0 * self%horzCoeffs(1,1)
                    cp_x(1) = C_E * self%centerCoeffs(1,k)
                    sourceX(1) = self%sourceTerm(1, j) - C_N * self%solution(1, N_indx) - C_S * self%solution(1, S_indx)
                    dp_x(1) = (sourceX(1)) * self%centerCoeffs(1,k)
                case (3,4)
                    cp_x(1) = 0.0d0
                    C_E = self%horzCoeffs(1,1)
                    C_W = self%horzCoeffs(2,1)
                    E_indx = self%horzIndx(1,1)
                    W_indx = self%horzIndx(2,1)
                    sourceX(1) = (self%sourceTerm(1,j) - self%solution(1, N_indx) * C_N - self%solution(1, S_indx) * C_S - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * self%centerCoeffs(1,k)* self%omega + omega_inv * self%solution(1,j)
                    dp_x(1) = sourceX(1)
                end select
                ! solve for vectors c-prime and d-prime
                do i = 2,self%N_x-1
                    C_W = self%horzCoeffs(2, i+1-self%startCol)
                    C_E = self%horzCoeffs(1, i+1-self%startCol)
                    m = 1.0d0/self%centerCoeffs(i+1-self%startCol, k) - cp_x(i-1) * C_W
                    cp_x(i) = C_E / m
                    sourceX(i) = self%sourceTerm(i,j) - C_N * self%solution(i, N_indx) - C_S * self%solution(i, S_indx)
                    dp_x(i) = (sourceX(i) - dp_x(i-1) * C_W) / m
                end do
                select case(self%rightBound)
                case (1)
                    sourceX(self%N_x) = self%solution(self%N_x, j)
                case (2)
                    sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - C_N * self%solution(self%N_x, N_indx) - C_S * self%solution(self%N_x, S_indx)
                    C_W = 2.0d0 * self%horzCoeffs(2, self%numberColumns)
                    sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * C_W)/&
                        (1.0/self%centerCoeffs(self%numberColumns, k) - cp_x(self%N_x-1)  * C_W)
                case (3)
                    sourceX(self%N_x) = sourceX(1)
                end select
                do i = self%N_x-1, 1, -1
                    sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
                end do
                self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
            end do
            !$OMP end do
            ! Sweep horizontal even rows
            !$OMP do
            do k = 2, self%numberRows, 2
                j = self%startRow + k - 1
                N_indx = self%vertIndx(1, k)
                S_indx = self%vertIndx(2, k)
                C_N = self%vertCoeffs(1, k)
                C_S = self%vertCoeffs(2, k)
                select case (self%leftBound)
                case (1)
                    cp_x(1) = 0.0d0
                    sourceX(1) = self%solution(1,j)
                    dp_x(1) = sourceX(1)
                case (2)
                    C_E = 2.0d0 * self%horzCoeffs(1,1)
                    cp_x(1) = C_E * self%centerCoeffs(1,k)
                    sourceX(1) = self%sourceTerm(1, j) - C_N * self%solution(1, N_indx) - C_S * self%solution(1, S_indx)
                    dp_x(1) = (sourceX(1)) * self%centerCoeffs(1,k)
                case (3,4)
                    cp_x(1) = 0.0d0
                    C_E = self%horzCoeffs(1,1)
                    C_W = self%horzCoeffs(2,1)
                    E_indx = self%horzIndx(1,1)
                    W_indx = self%horzIndx(2,1)
                    sourceX(1) = (self%sourceTerm(1,j) - self%solution(1, N_indx) * C_N - self%solution(1, S_indx) * C_S - &
                    self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * self%centerCoeffs(1,k)* self%omega + omega_inv * self%solution(1,j)
                    dp_x(1) = sourceX(1)
                end select
                ! solve for vectors c-prime and d-prime
                do i = 2,self%N_x-1
                    C_W = self%horzCoeffs(2, i+1-self%startCol)
                    C_E = self%horzCoeffs(1, i+1-self%startCol)
                    m = 1.0d0/self%centerCoeffs(i+1-self%startCol, k) - cp_x(i-1) * C_W
                    cp_x(i) = C_E / m
                    sourceX(i) = self%sourceTerm(i,j) - C_N * self%solution(i, N_indx) - C_S * self%solution(i, S_indx)
                    dp_x(i) = (sourceX(i) - dp_x(i-1) * C_W) / m
                end do
                select case(self%rightBound)
                case (1)
                    sourceX(self%N_x) = self%solution(self%N_x, j)
                case (2)
                    sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - C_N * self%solution(self%N_x, N_indx) - C_S * self%solution(self%N_x, S_indx)
                    C_W = 2.0d0 * self%horzCoeffs(2, self%numberColumns)
                    sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * C_W)/&
                        (1.0/self%centerCoeffs(self%numberColumns, k) - cp_x(self%N_x-1)  * C_W)
                case (3)
                    sourceX(self%N_x) = sourceX(1)
                end select
                do i = self%N_x-1, 1, -1
                    sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
                end do
                self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
            end do
            !$OMP end do
            

            ! vertical sweeps bottom to top
            !$OMP do
            ! Sweep odd numbers upwards
            do p = 1, self%numberColumns, 2
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                C_E = self%horzCoeffs(1, p)
                C_W = self%horzCoeffs(2, p)
                select case (self%lowerBound)
                case (1)
                    cp_y(1) = 0.0d0
                    sourceY(1) = self%solution(i,1)
                    dp_y(1) = sourceY(1)
                case (2)
                    C_N = 2.0d0 * self%vertCoeffs(1,1)
                    cp_y(1) = C_N * self%centerCoeffs(p,1)
                    sourceY(1) = self%sourceTerm(i, 1) - C_E * self%solution(E_indx, 1) - C_W * self%solution(W_indx, 1)
                    dp_y(1) = (sourceY(1)) * self%centerCoeffs(p,1)
                case (3,4)
                    cp_y(1) = 0.0d0
                    C_N = self%vertCoeffs(1,1)
                    C_S = self%vertCoeffs(2,1)
                    N_indx = self%vertIndx(1,1)
                    S_indx = self%vertIndx(2,1)
                    sourceY(1) = (self%sourceTerm(i,1) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S - &
                    self%solution(E_indx, 1) * C_E - self%solution(W_indx, 1) * C_W) * self%centerCoeffs(p,1)* self%omega + omega_inv * self%solution(i,1)
                    dp_y(1) = sourceY(1)
                end select
                ! solve for vectors c-prime and d-prime
                do j = 2,self%N_y-1
                    C_S = self%vertCoeffs(2, j+1-self%startRow)
                    C_N = self%vertCoeffs(1, j+1-self%startRow)
                    m = 1.0d0/self%centerCoeffs(p, j+1-self%startRow) - cp_y(j-1) * C_S
                    cp_y(j) = C_N / m
                    sourceY(j) = self%sourceTerm(i, j) - C_E * self%solution(E_indx, j) - C_W * self%solution(W_indx, j)
                    dp_y(j) = (sourceY(j) - dp_y(j-1) * C_S) / m
                end do
                select case(self%upperBound)
                case (1)
                    sourceY(self%N_y) = self%solution(i, self%N_y)
                case (2)
                    sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - C_E * self%solution(E_indx, self%N_y) - C_W * self%solution(W_indx, self%N_y)
                    C_S = 2.0d0 * self%vertCoeffs(2, self%numberRows)
                    sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * C_S)/&
                        (1.0/self%centerCoeffs(p, self%numberRows) - cp_y(self%N_y-1)  * C_S)
                case (3)
                    sourceY(self%N_y) = sourceY(1)
                end select
                do j = self%N_y-1, 1, -1
                    sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
                end do
                self%solution(i,:) = sourceY *self%omega + omega_inv * self%solution(i,:)
            end do
            !$OMP end do
            !$OMP do
            ! sweep even numbered columns
            do p = 2, self%numberColumns, 2
                i = self%startCol + p - 1
                E_indx = self%horzIndx(1, p)
                W_indx = self%horzIndx(2, p)
                C_E = self%horzCoeffs(1, p)
                C_W = self%horzCoeffs(2, p)
                select case (self%lowerBound)
                case (1)
                    cp_y(1) = 0.0d0
                    sourceY(1) = self%solution(i,1)
                    dp_y(1) = sourceY(1)
                case (2)
                    C_N = 2.0d0 * self%vertCoeffs(1,1)
                    cp_y(1) = C_N * self%centerCoeffs(p,1)
                    sourceY(1) = self%sourceTerm(i, 1) - C_E * self%solution(E_indx, 1) - C_W * self%solution(W_indx, 1)
                    dp_y(1) = (sourceY(1)) * self%centerCoeffs(p,1)
                case (3,4)
                    cp_y(1) = 0.0d0
                    C_N = self%vertCoeffs(1,1)
                    C_S = self%vertCoeffs(2,1)
                    N_indx = self%vertIndx(1,1)
                    S_indx = self%vertIndx(2,1)
                    sourceY(1) = (self%sourceTerm(i,1) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S - &
                    self%solution(E_indx, 1) * C_E - self%solution(W_indx, 1) * C_W) * self%centerCoeffs(p,1)* self%omega + omega_inv * self%solution(i,1)
                    dp_y(1) = sourceY(1)
                end select
                ! solve for vectors c-prime and d-prime
                do j = 2,self%N_y-1
                    C_S = self%vertCoeffs(2, j+1-self%startRow)
                    C_N = self%vertCoeffs(1, j+1-self%startRow)
                    m = 1.0d0/self%centerCoeffs(p, j+1-self%startRow) - cp_y(j-1) * C_S
                    cp_y(j) = C_N / m
                    sourceY(j) = self%sourceTerm(i, j) - C_E * self%solution(E_indx, j) - C_W * self%solution(W_indx, j)
                    dp_y(j) = (sourceY(j) - dp_y(j-1) * C_S) / m
                end do
                select case(self%upperBound)
                case (1)
                    sourceY(self%N_y) = self%solution(i, self%N_y)
                case (2)
                    sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - C_E * self%solution(E_indx, self%N_y) - C_W * self%solution(W_indx, self%N_y)
                    C_S = 2.0d0 * self%vertCoeffs(2, self%numberRows)
                    sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * C_S)/&
                        (1.0/self%centerCoeffs(p, self%numberRows) - cp_y(self%N_y-1)  * C_S)
                case (3)
                    sourceY(self%N_y) = sourceY(1)
                end select
                do j = self%N_y-1, 1, -1
                    sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
                end do
                self%solution(i,:) = sourceY *self%omega + omega_inv * self%solution(i,:)
            end do
            !$OMP end do
            !$OMP end parallel
        end do
    end subroutine smoothIterations_ZebraCurv

    function smoothWithRes_ZebraCurv(self) result(Res)
        ! Solve GS down to some tolerance
        class(ZebraSolverCurv), intent(in out) :: self
        real(real64) :: omega_inv, C_N, C_E, C_S, C_W, Res
        real(real64) :: cp_x(self%N_x-1), dp_x(self%N_x-1), cp_y(self%N_y-1), dp_y(self%N_y), sourceX(self%N_x), sourceY(self%N_y), m
        integer :: k,p, iter, E_indx, W_indx, N_indx, S_indx, i, j
        omega_inv = 1.0d0 - self%omega
        Res = 0.0d0
       
        !$OMP parallel private(k, p, i, j, E_indx, N_indx, S_indx, W_indx, &
        !$OMP& C_N, C_S, C_E, C_W, cp_y, cp_x, dp_y, dp_x, sourceX, sourceY, m) reduction(+:Res)

        ! horizontal sweeps left to right
        !$OMP do
        ! Sweep odd numbers forward
        do k = 1, self%numberRows, 2
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            C_N = self%vertCoeffs(1, k)
            C_S = self%vertCoeffs(2, k)
            ! initialize c-prime and d-prime
            select case (self%leftBound)
            case (1)
                cp_x(1) = 0.0d0
                sourceX(1) = self%solution(1,j)
                dp_x(1) = sourceX(1)
            case (2)
                C_E = 2.0d0 * self%horzCoeffs(1,1)
                cp_x(1) = C_E * self%centerCoeffs(1,k)
                sourceX(1) = self%sourceTerm(1, j) - C_N * self%solution(1, N_indx) - C_S * self%solution(1, S_indx)
                dp_x(1) = (sourceX(1)) * self%centerCoeffs(1,k)
            case (3,4)
                cp_x(1) = 0.0d0
                C_E = self%horzCoeffs(1,1)
                C_W = self%horzCoeffs(2,1)
                E_indx = self%horzIndx(1,1)
                W_indx = self%horzIndx(2,1)
                sourceX(1) = (self%sourceTerm(1,j) - self%solution(1, N_indx) * C_N - self%solution(1, S_indx) * C_S - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * self%centerCoeffs(1,k)* self%omega + omega_inv * self%solution(1,j)
                dp_x(1) = sourceX(1)
            end select
            ! solve for vectors c-prime and d-prime
            do i = 2,self%N_x-1
                C_W = self%horzCoeffs(2, i+1-self%startCol)
                C_E = self%horzCoeffs(1, i+1-self%startCol)
                m = 1.0d0/self%centerCoeffs(i+1-self%startCol, k) - cp_x(i-1) * C_W
                cp_x(i) = C_E / m
                sourceX(i) = self%sourceTerm(i,j) - C_N * self%solution(i, N_indx) - C_S * self%solution(i, S_indx)
                dp_x(i) = (sourceX(i) - dp_x(i-1) * C_W) / m
            end do
            select case(self%rightBound)
            case (1)
                sourceX(self%N_x) = self%solution(self%N_x, j)
            case (2)
                sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - C_N * self%solution(self%N_x, N_indx) - C_S * self%solution(self%N_x, S_indx)
                C_W = 2.0d0 * self%horzCoeffs(2, self%numberColumns)
                sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * C_W)/&
                    (1.0/self%centerCoeffs(self%numberColumns, k) - cp_x(self%N_x-1)  * C_W)
            case (3)
                sourceX(self%N_x) = sourceX(1)
            end select
            do i = self%N_x-1, 1, -1
                sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
            end do
            self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
        end do
        !$OMP end do
        !$OMP do
        do k = 2, self%numberRows, 2
            j = self%startRow + k - 1
            N_indx = self%vertIndx(1, k)
            S_indx = self%vertIndx(2, k)
            C_N = self%vertCoeffs(1, k)
            C_S = self%vertCoeffs(2, k)
            select case (self%leftBound)
            case (1)
                cp_x(1) = 0.0d0
                sourceX(1) = self%solution(1,j)
                dp_x(1) = sourceX(1)
            case (2)
                C_E = 2.0d0 * self%horzCoeffs(1,1)
                cp_x(1) = C_E * self%centerCoeffs(1,k)
                sourceX(1) = self%sourceTerm(1, j) - C_N * self%solution(1, N_indx) - C_S * self%solution(1, S_indx)
                dp_x(1) = (sourceX(1)) * self%centerCoeffs(1,k)
            case (3,4)
                cp_x(1) = 0.0d0
                C_E = self%horzCoeffs(1,1)
                C_W = self%horzCoeffs(2,1)
                E_indx = self%horzIndx(1,1)
                W_indx = self%horzIndx(2,1)
                sourceX(1) = (self%sourceTerm(1,j) - self%solution(1, N_indx) * C_N - self%solution(1, S_indx) * C_S - &
                self%solution(E_indx, j) * C_E - self%solution(W_indx, j) * C_W) * self%centerCoeffs(1,k)* self%omega + omega_inv * self%solution(1,j)
                dp_x(1) = sourceX(1)
            end select
            ! solve for vectors c-prime and d-prime
            do i = 2,self%N_x-1
                C_W = self%horzCoeffs(2, i+1-self%startCol)
                C_E = self%horzCoeffs(1, i+1-self%startCol)
                m = 1.0d0/self%centerCoeffs(i+1-self%startCol, k) - cp_x(i-1) * C_W
                cp_x(i) = C_E / m
                sourceX(i) = self%sourceTerm(i,j) - C_N * self%solution(i, N_indx) - C_S * self%solution(i, S_indx)
                dp_x(i) = (sourceX(i) - dp_x(i-1) * C_W) / m
            end do
            select case(self%rightBound)
            case (1)
                sourceX(self%N_x) = self%solution(self%N_x, j)
            case (2)
                sourceX(self%N_x) = self%sourceTerm(self%N_x,j) - C_N * self%solution(self%N_x, N_indx) - C_S * self%solution(self%N_x, S_indx)
                C_W = 2.0d0 * self%horzCoeffs(2, self%numberColumns)
                sourceX(self%N_x) = (sourceX(self%N_x) - dp_x(self%N_x-1) * C_W)/&
                    (1.0/self%centerCoeffs(self%numberColumns, k) - cp_x(self%N_x-1)  * C_W)
            case (3)
                sourceX(self%N_x) = sourceX(1)
            end select
            do i = self%N_x-1, 1, -1
                sourceX(i) = dp_x(i)-cp_x(i)*sourceX(i+1)
            end do
            self%solution(:,j) = sourceX *self%omega + omega_inv * self%solution(:,j)
        end do
        !$OMP end do
        

        ! vertical sweeps bottom to top
        !$OMP do
        ! Sweep odd numbers upwards
        do p = 1, self%numberColumns, 2
            i = self%startCol + p - 1
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)
            C_E = self%horzCoeffs(1, p)
            C_W = self%horzCoeffs(2, p)
            select case (self%lowerBound)
            case (1)
                cp_y(1) = 0.0d0
                sourceY(1) = self%solution(i,1)
                dp_y(1) = sourceY(1)
            case (2)
                C_N = 2.0d0 * self%vertCoeffs(1,1)
                cp_y(1) = C_N * self%centerCoeffs(p,1)
                sourceY(1) = self%sourceTerm(i, 1) - C_E * self%solution(E_indx, 1) - C_W * self%solution(W_indx, 1)
                dp_y(1) = (sourceY(1)) * self%centerCoeffs(p,1)
            case (3,4)
                cp_y(1) = 0.0d0
                C_N = self%vertCoeffs(1,1)
                C_S = self%vertCoeffs(2,1)
                N_indx = self%vertIndx(1,1)
                S_indx = self%vertIndx(2,1)
                sourceY(1) = (self%sourceTerm(i,1) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S - &
                self%solution(E_indx, 1) * C_E - self%solution(W_indx, 1) * C_W) * self%centerCoeffs(p,1)* self%omega + omega_inv * self%solution(i,1)
                dp_y(1) = sourceY(1)
            end select
            ! solve for vectors c-prime and d-prime
            do j = 2,self%N_y-1
                C_S = self%vertCoeffs(2, j+1-self%startRow)
                C_N = self%vertCoeffs(1, j+1-self%startRow)
                m = 1.0d0/self%centerCoeffs(p, j+1-self%startRow) - cp_y(j-1) * C_S
                cp_y(j) = C_N / m
                sourceY(j) = self%sourceTerm(i, j) - C_E * self%solution(E_indx, j) - C_W * self%solution(W_indx, j)
                dp_y(j) = (sourceY(j) - dp_y(j-1) * C_S) / m
            end do
            select case(self%upperBound)
            case (1)
                sourceY(self%N_y) = self%solution(i, self%N_y)
            case (2)
                sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - C_E * self%solution(E_indx, self%N_y) - C_W * self%solution(W_indx, self%N_y)
                C_S = 2.0d0 * self%vertCoeffs(2, self%numberRows)
                sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * C_S)/&
                    (1.0/self%centerCoeffs(p, self%numberRows) - cp_y(self%N_y-1)  * C_S)
            case (3)
                sourceY(self%N_y) = sourceY(1)
            end select
            do j = self%N_y-1, 1, -1
                sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
            end do
            sourceY = sourceY *self%omega + omega_inv * self%solution(i,:)
            Res = Res + SUM((sourceY(self%startRow:self%endRow) - self%solution(i, self%startRow:self%endRow))**2)
            self%solution(i,:) = sourceY
        end do
        !$OMP end do
        !$OMP do
        ! sweep even numbers downwards
        do p = 2, self%numberColumns, 2
            i = self%startCol + p - 1
            E_indx = self%horzIndx(1, p)
            W_indx = self%horzIndx(2, p)
            C_E = self%horzCoeffs(1, p)
            C_W = self%horzCoeffs(2, p)
            select case (self%lowerBound)
            case (1)
                cp_y(1) = 0.0d0
                sourceY(1) = self%solution(i,1)
                dp_y(1) = sourceY(1)
            case (2)
                C_N = 2.0d0 * self%vertCoeffs(1,1)
                cp_y(1) = C_N * self%centerCoeffs(p,1)
                sourceY(1) = self%sourceTerm(i, 1) - C_E * self%solution(E_indx, 1) - C_W * self%solution(W_indx, 1)
                dp_y(1) = (sourceY(1)) * self%centerCoeffs(p,1)
            case (3,4)
                cp_y(1) = 0.0d0
                C_N = self%vertCoeffs(1,1)
                C_S = self%vertCoeffs(2,1)
                N_indx = self%vertIndx(1,1)
                S_indx = self%vertIndx(2,1)
                sourceY(1) = (self%sourceTerm(i,1) - self%solution(i, N_indx) * C_N - self%solution(i, S_indx) * C_S - &
                self%solution(E_indx, 1) * C_E - self%solution(W_indx, 1) * C_W) * self%centerCoeffs(p,1)* self%omega + omega_inv * self%solution(i,1)
                dp_y(1) = sourceY(1)
            end select
            ! solve for vectors c-prime and d-prime
            do j = 2,self%N_y-1
                C_S = self%vertCoeffs(2, j+1-self%startRow)
                C_N = self%vertCoeffs(1, j+1-self%startRow)
                m = 1.0d0/self%centerCoeffs(p, j+1-self%startRow) - cp_y(j-1) * C_S
                cp_y(j) = C_N / m
                sourceY(j) = self%sourceTerm(i, j) - C_E * self%solution(E_indx, j) - C_W * self%solution(W_indx, j)
                dp_y(j) = (sourceY(j) - dp_y(j-1) * C_S) / m
            end do
            select case(self%upperBound)
            case (1)
                sourceY(self%N_y) = self%solution(i, self%N_y)
            case (2)
                sourceY(self%N_y) = self%sourceTerm(i, self%N_y) - C_E * self%solution(E_indx, self%N_y) - C_W * self%solution(W_indx, self%N_y)
                C_S = 2.0d0 * self%vertCoeffs(2, self%numberRows)
                sourceY(self%N_y) = (sourceY(self%N_y) - dp_y(self%N_y-1) * C_S)/&
                    (1.0/self%centerCoeffs(p, self%numberRows) - cp_y(self%N_y-1)  * C_S)
            case (3)
                sourceY(self%N_y) = sourceY(1)
            end select
            do j = self%N_y-1, 1, -1
                sourceY(j) = dp_y(j)-cp_y(j)*sourceY(j+1)
            end do
            sourceY = sourceY*self%omega + omega_inv * self%solution(i,:)
            Res = Res + SUM((sourceY(self%startRow:self%endRow) - self%solution(i, self%startRow:self%endRow))**2)
            self%solution(i,:) = sourceY
        end do
        !$OMP end do
        !$OMP end parallel
        Res = SQRT(Res/ (self%numberColumns * self%numberRows))

    end function smoothWithRes_ZebraCurv


    subroutine calcResidual_ZebraCurv(self)
        ! Solve GS down to some tolerance
        class(ZebraSolverCurv), intent(in out) :: self
        real(real64) :: C_N, C_E, C_O, C_W, C_S
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, iter

        !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
        !$OMP&  W_indx, C_O, C_N, C_E, C_S, C_W)
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
                C_O = self%centerCoeffs(p,k)
                C_N = self%vertCoeffs(1,k)
                C_E = self%horzCoeffs(1,p)
                C_S = self%vertCoeffs(2,k)
                C_W = self%horzCoeffs(2,p)
                self%residual(i, j) = self%sourceTerm(i,j) - self%solution(i,N_indx)*C_N - self%solution(i,S_indx) * C_S &
                    - self%solution(E_indx, j) *C_E - self%solution(W_indx, j) * C_W - self%solution(i,j)/C_O
            end do
        end do
        !$OMP end do nowait
        !$OMP end parallel
    end subroutine calcResidual_ZebraCurv

    ! subroutine restriction(self, fineGrid, coarseGrid)
    !     ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
    !     class(ZebraSolver), intent(in) :: self
    !     real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
    !     real(real64), intent(out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    !     real(real64) :: a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S
    !     integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k
    !     if (self%evenGridBool) then
    !         !$OMP parallel private(k, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx)
    !         !$OMP do collapse(2)
    !         ! loop over fine grid nodes which overlap with coarse grid nodes
    !         do j_fine = 3, self%N_y-2, 2
    !             do i_fine = 3, self%N_x-2, 2
    !                 ! calculate fine indices around overlapping index
    !                 i_coarse = (i_fine + 1)/2
    !                 j_coarse = (j_fine + 1)/2
    !                 N_indx = j_fine+1
    !                 E_indx = i_fine+1
    !                 S_indx = j_fine-1
    !                 W_indx = i_fine-1

    !                 ! Simple interpolation can be surprisingly good
    !                 coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
    !                 0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
    !                 fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
    !                 0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
    !                 fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         ! each black bound node is a coarse node
    !         !$OMP do
    !         do k = 1, self%numberBlackBoundNodes
    !             i_fine = self%black_NESW_BoundIndx(1, k)
    !             j_fine = self%black_NESW_BoundIndx(2, k)
    !             E_indx = self%black_NESW_BoundIndx(3, k)
    !             W_indx = self%black_NESW_BoundIndx(4, k)
    !             N_indx = self%black_NESW_BoundIndx(5, k)
    !             S_indx = self%black_NESW_BoundIndx(6, k)
                
    !             ! calculate fine indices around overlapping index
    !             i_coarse = (i_fine + 1)/2
    !             j_coarse = (j_fine + 1)/2

    !             ! simple binomial filter
    !             coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
    !                 0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
    !                 fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
    !                 0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
    !                 fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
    !         end do
    !         !$OMP end do
    !         !$OMP end parallel

    !     else
    !         !$OMP parallel private(k, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, &
    !         !$OMP& a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S)
    !         !$OMP do collapse(2)
    !         ! loop over fine grid nodes which overlap with coarse grid nodes
    !         do j_fine = 3, self%N_y-2, 2
    !             do i_fine = 3, self%N_x-2, 2
    !                 ! calculate fine indices around overlapping index
    !                 i_coarse = (i_fine + 1)/2
    !                 j_coarse = (j_fine + 1)/2
    !                 N_indx = j_fine+1
    !                 E_indx = i_fine+1
    !                 S_indx = j_fine-1
    !                 W_indx = i_fine-1

    !                 ! Transpose of prolongation, seems to work best
    !                 ! interpolation from north point
    !                 C_N = self%matCoeffs(2, i_fine, N_indx)
    !                 C_S = self%matCoeffs(4, i_fine, N_indx)
    !                 a_N = C_S/(C_N + C_S)

    !                 ! interpolation from south point
    !                 C_N = self%matCoeffs(2, i_fine, S_indx)
    !                 C_S = self%matCoeffs(4, i_fine, S_indx)
    !                 a_S = C_N/(C_N + C_S)

    !                 ! interpolation from east point
    !                 C_E = self%matCoeffs(3, E_indx, j_fine)
    !                 C_W = self%matCoeffs(5, E_indx, j_fine)
    !                 a_E = C_W/(C_E + C_W)

    !                 ! interpolation from west point
    !                 C_E = self%matCoeffs(3, W_indx, j_fine)
    !                 C_W = self%matCoeffs(5, W_indx, j_fine)
    !                 a_W = C_E/(C_E + C_W)

    !                 ! node based restriction
    !                 ! C_N = self%matCoeffs(2, i_fine, j_fine)
    !                 ! C_E = self%matCoeffs(3, i_fine, j_fine)
    !                 ! C_S = self%matCoeffs(4, i_fine, j_fine)
    !                 ! C_W = self%matCoeffs(5, i_fine, j_fine)
    !                 ! a_N = C_N/(C_N + C_S)
    !                 ! a_S = C_S/(C_N + C_S)
    !                 ! a_E = C_E/(C_E + C_W)
    !                 ! a_W = C_W/(C_E + C_W)


    !                 ! Interpolate residual to coarse grid
    !                 coarseGrid(i_coarse, j_coarse) = 0.25d0 * (fineGrid(i_fine, j_fine) + &
    !                 fineGrid(E_indx, j_fine) * a_E + fineGrid(W_indx, j_fine) * a_W +  &
    !                 fineGrid(i_fine, N_indx) * a_N + fineGrid(i_fine, S_indx) * a_S + &
    !                 fineGrid(E_indx, N_indx) * a_E * a_N + fineGrid(E_indx, S_indx) * a_E * a_S + &
    !                 fineGrid(W_indx, N_indx) * a_W * a_N + fineGrid(W_indx, S_indx) * a_W * a_S)

    !                 ! Simple interpolation can be surprisingly good
    !                 ! coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
    !                 ! 0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
    !                 ! fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
    !                 ! 0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
    !                 ! fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         ! each black bound node is a coarse node
    !         !$OMP do
    !         do k = 1, self%numberBlackBoundNodes
    !             i_fine = self%black_NESW_BoundIndx(1, k)
    !             j_fine = self%black_NESW_BoundIndx(2, k)
    !             E_indx = self%black_NESW_BoundIndx(3, k)
    !             W_indx = self%black_NESW_BoundIndx(4, k)
    !             N_indx = self%black_NESW_BoundIndx(5, k)
    !             S_indx = self%black_NESW_BoundIndx(6, k)
                
    !             ! calculate fine indices around overlapping index
    !             i_coarse = (i_fine + 1)/2
    !             j_coarse = (j_fine + 1)/2
    !             ! From transpose of prolongation, works fairly well
    !             if (N_indx == S_indx) then
    !                 ! neumann boundary
    !                 if (j_fine == 1) then
    !                     ! bottom boundary
    !                     C_N = self%matCoeffs(2, i_fine, N_indx)
    !                     C_S = self%matCoeffs(4, i_fine, N_indx)
    !                     a_N = C_S/(C_N + C_S)
    !                     a_S = a_N
    !                 else
    !                     ! top boundary
    !                     C_N = self%matCoeffs(2, i_fine, S_indx)
    !                     C_S = self%matCoeffs(4, i_fine, S_indx)
    !                     a_S = C_N/(C_N + C_S)
    !                     a_N = a_S
    !                 end if
    !             else
    !                 ! interpolation from north point
    !                 C_N = self%matCoeffs(2, i_fine, N_indx)
    !                 C_S = self%matCoeffs(4, i_fine, N_indx)
    !                 a_N = C_S/(C_N + C_S)

    !                 ! interpolation from south point
    !                 C_N = self%matCoeffs(2, i_fine, S_indx)
    !                 C_S = self%matCoeffs(4, i_fine, S_indx)
    !                 a_S = C_N/(C_N + C_S)
    !             end if
                
    !             if (E_indx == W_indx) then
    !                 ! neumann boundary
    !                 if (i_fine == 1) then
    !                     ! left boundary
    !                     C_E = self%matCoeffs(3, E_indx, j_fine)
    !                     C_W = self%matCoeffs(5, E_indx, j_fine)
    !                     a_E = C_W/(C_E + C_W)
    !                     a_W = a_E
    !                 else
    !                     C_E = self%matCoeffs(3, W_indx, j_fine)
    !                     C_W = self%matCoeffs(5, W_indx, j_fine)
    !                     a_W = C_E/(C_E + C_W)
    !                     a_E = a_W
    !                 end if
    !             else
    !                 ! interpolation from east point
    !                 C_E = self%matCoeffs(3, E_indx, j_fine)
    !                 C_W = self%matCoeffs(5, E_indx, j_fine)
    !                 a_E = C_W/(C_E + C_W)

    !                 ! interpolation from west point
    !                 C_E = self%matCoeffs(3, W_indx, j_fine)
    !                 C_W = self%matCoeffs(5, W_indx, j_fine)
    !                 a_W = C_E/(C_E + C_W)
    !             end if

    !             ! node based restrition
    !             ! C_N = self%matCoeffs(2, i_fine, j_fine)
    !             ! C_E = self%matCoeffs(3, i_fine, j_fine)
    !             ! C_S = self%matCoeffs(4, i_fine, j_fine)
    !             ! C_W = self%matCoeffs(5, i_fine, j_fine)
    !             ! a_N = C_N/(C_N + C_S)
    !             ! a_S = C_S/(C_N + C_S)
    !             ! a_E = C_E/(C_E + C_W)
    !             ! a_W = C_W/(C_E + C_W)

    !             ! Interpolate residual to coarse grid
    !             coarseGrid(i_coarse, j_coarse) = 0.25d0 * (fineGrid(i_fine, j_fine) + &
    !             fineGrid(E_indx, j_fine) * a_E + fineGrid(W_indx, j_fine) * a_W +  &
    !             fineGrid(i_fine, N_indx) * a_N + fineGrid(i_fine, S_indx) * a_S + &
    !             fineGrid(E_indx, N_indx) * a_E * a_N + fineGrid(E_indx, S_indx) * a_E * a_S + &
    !             fineGrid(W_indx, N_indx) * a_W * a_N + fineGrid(W_indx, S_indx) * a_W * a_S)

    !             ! simple binomial filter, actually ok
    !             ! coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
    !             !     0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
    !             !     fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
    !             !     0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
    !             !     fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
    !         end do
    !         !$OMP end do
    !         !$OMP end parallel
    !     end if

    ! end subroutine restriction

    ! subroutine prolongation(self, fineGrid, coarseGrid)
    !     ! Prolongate operator from coarse to fine grid using gauss seidel data
    !     class(ZebraSolver), intent(in) :: self
    !     real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
    !     real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
    !     real(real64) :: w_1, w_2, w_3, w_4, w_tot
    !     integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse
    !     N_x_coarse = (self%N_x+1)/2
    !     N_y_coarse = (self%N_y+1)/2

    !     if (self%evenGridBool) then
    !         ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
    !         !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine)
    !         ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
    !         !$OMP do collapse(2)
    !         ! set similar nodes
    !         do j_coarse = 1, N_y_coarse
    !             do i_coarse = 1, N_x_coarse
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2-1
    !                 ! Add overlapping coarse nodes
    !                 fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         !$OMP do collapse(2)
    !         ! go horizontal each row,
    !         do j_coarse = 1, N_y_coarse
    !             do i_coarse = 1, N_x_coarse-1
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2 - 1
    !                 ! simple interpolation
    !                 fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         !$OMP do collapse(2)
    !         ! go vertical
    !         do j_coarse = 1, N_y_coarse-1
    !             do i_coarse = 1, N_x_coarse
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2 - 1
    !                 ! Simple interpolation
    !                 fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         !$OMP do collapse(2)
    !         ! add offset fine nodes which require 4 point interpolation
    !         ! first loop through SW corner coarse nodes
    !         do j_coarse = 1, N_y_coarse-1
    !             do i_coarse = 1, N_x_coarse-1
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2 - 1

    !                 ! simple interpolation
    !                 fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse) * 0.25d0 &
    !                     + coarseGrid(i_coarse, j_coarse+1) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse+1) * 0.25d0
    !             end do
    !         end do
    !         !$OMP end do
    !         !$OMP end parallel

    !     else
    !         ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
    !         !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine, w_1, w_2, w_3, w_4, w_tot)
    !         ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
    !         !$OMP do collapse(2)
    !         ! set similar nodes
    !         do j_coarse = 1, N_y_coarse
    !             do i_coarse = 1, N_x_coarse
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2-1
    !                 ! Add overlapping coarse nodes
    !                 fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         !$OMP do collapse(2)
    !         ! go horizontal each row,
    !         do j_coarse = 1, N_y_coarse
    !             do i_coarse = 1, N_x_coarse-1
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2 - 1

    !                 w_1 = self%matCoeffs(3, i_fine+1, j_fine) !E
    !                 w_2 = self%matCoeffs(5, i_fine+1, j_fine) !W
    !                 w_tot = w_1 + w_2
    !                 w_1 = w_1 / w_tot
    !                 w_2 = w_2 / w_tot
        
    !                 ! Average fine nodes between coarse nodes horizontally
    !                 fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
    !                 ! simple interpolation
    !                 ! fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         !$OMP do collapse(2)
    !         ! go vertical
    !         do j_coarse = 1, N_y_coarse-1
    !             do i_coarse = 1, N_x_coarse
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2 - 1
    !                 w_1 = self%matCoeffs(2, i_fine, j_fine+1) !N
    !                 w_2 = self%matCoeffs(4, i_fine, j_fine+1) !S
    !                 w_tot = w_1 + w_2
    !                 w_1 = w_1 / w_tot
    !                 w_2 = w_2 / w_tot
    !                 ! Average fine nodes between coarse nodes vertical direction
    !                 fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * w_1 + coarseGrid(i_coarse, j_coarse) * w_2
    !                 ! Simple interpolation
    !                 ! fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
    !             end do
    !         end do
    !         !$OMP end do nowait
    !         !$OMP do collapse(2)
    !         ! add offset fine nodes which require 4 point interpolation
    !         ! first loop through SW corner coarse nodes
    !         do j_coarse = 1, N_y_coarse-1
    !             do i_coarse = 1, N_x_coarse-1
    !                 i_fine = i_coarse*2 - 1
    !                 j_fine = j_coarse*2 - 1
    !                 w_1 = self%matCoeffs(2, i_fine+1, j_fine+1) !N
    !                 w_2 = self%matCoeffs(4, i_fine+1, j_fine+1) !S
    !                 w_tot = w_1 + w_2
    !                 w_1 = w_1/w_tot
    !                 w_2 = w_2/w_tot

    !                 w_3 = self%matCoeffs(3, i_fine+1, j_fine+1) !E
    !                 w_4 = self%matCoeffs(5, i_fine+1, j_fine+1) !W
    !                 w_tot = w_3 + w_4
    !                 w_3 = w_3/w_tot
    !                 w_4 = w_4/w_tot

    !                 fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * w_2*w_4 + coarseGrid(i_coarse+1, j_coarse) * w_3 * w_2 &
    !                     + coarseGrid(i_coarse, j_coarse+1) * w_1 * w_4 + coarseGrid(i_coarse+1, j_coarse+1) * w_1 * w_3

    !                 ! simple interpolation
    !                 ! fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse) * 0.25d0 &
    !                 !     + coarseGrid(i_coarse, j_coarse+1) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse+1) * 0.25d0
    !             end do
    !         end do
    !         !$OMP end do
    !         !$OMP end parallel
    !     end if

    ! end subroutine prolongation
    subroutine solve_tridiag(n, diagLower, diagUpper, diag, b, x)
        ! General tridiagonal solver
        integer(int32), intent(in) :: n
        real(real64), intent(in out) :: x(n)
        real(real64), intent(in) :: diagLower(n-1), diagUpper(n-1), diag(n), b(n)
        integer(int32) :: i !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, cp(n-1),dp(n)

    ! initialize c-prime and d-prime
        cp(1) = diagUpper(1)/diag(1)
        dp(1) = b(1)/diag(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n-1
            m = diag(i)-cp(i-1)*diagLower(i-1)
            cp(i) = diagUpper(i)/m
            dp(i) = (b(i)-dp(i-1)*diagLower(i-1))/m
        end do
        dp(n) = (b(n)-dp(n-1)*diagLower(n-1))/(diag(n)-cp(n-1)*diagLower(n-1))
        x(n) = dp(n)
        do i = n-1, 1, -1
            x(i) = dp(i)-cp(i)*x(i+1)
        end do
    end subroutine solve_tridiag

    subroutine solve_tridiag_replaceSource(n, diagLower, diagUpper, diag, b)
        ! General tridiagonal solver which replaces source term b with solution
        integer(int32), intent(in) :: n
        real(real64), intent(in out) :: b(n)
        real(real64), intent(in) :: diagLower(n-1), diagUpper(n-1), diag(n)
        integer(int32) :: i !n size dependent on how many points are inside (not boundary), so how many in rhs equation
        real(real64) :: m, cp(n-1),dp(n)

    ! initialize c-prime and d-prime
        cp(1) = diagUpper(1)/diag(1)
        dp(1) = b(1)/diag(1)
    ! solve for vectors c-prime and d-prime
        do i = 2,n-1
            m = diag(i)-cp(i-1)*diagLower(i-1)
            cp(i) = diagUpper(i)/m
            dp(i) = (b(i)-dp(i-1)*diagLower(i-1))/m
        end do
        dp(n) = (b(n)-dp(n-1)*diagLower(n-1))/(diag(n)-cp(n-1)*diagLower(n-1))
        b(n) = dp(n)
        do i = n-1, 1, -1
            b(i) = dp(i)-cp(i)*b(i+1)
        end do
    end subroutine solve_tridiag_replaceSource


end module mod_ZebraSolverCurv