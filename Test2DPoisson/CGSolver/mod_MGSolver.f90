module mod_MGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_GSSolver
    use mod_PardisoSolver
    use mod_CSRMAtrix
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: MGSolver

    type :: MGSolver
        type(GSSolver), allocatable :: GS_smoothers(:)
        type(pardisoSolver) :: directSolver
        integer :: numberStages, smoothNumber, maxIter, numberPreSmoothOper, numberPostSmoothOper, numIter ! number of stages is total, so includes finest and coarsest grid
        integer, allocatable :: N_x(:), N_y(:)
        real(real64) :: R2_current, stepResidual
    contains
        procedure, public, pass(self) :: makeSmootherStages
        procedure, public, pass(self) :: orthogonalGridRestriction
        procedure, public, pass(self) :: orthogonalGridProlongation
        procedure, public, pass(self) :: F_V_Cycle
        procedure, public, pass(self) :: F_Cycle
        procedure, public, pass(self) :: V_Cycle
        procedure, public, pass(self) :: MG_Cycle_Solve
    end type

    interface MGSolver
        module procedure :: MGSolver_constructor
    end interface MGSolver

contains

    type(MGSolver) function MGSolver_constructor(N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper) result(self)
        ! Construct object, set initial variables
        integer(int32), intent(in) :: N_x, N_y, numberStages, maxIter, numberPreSmoothOper, numberPostSmoothOper
        integer :: i
        self%numberStages = numberStages
        self%smoothNumber = numberStages-1 ! number smoothing GS stages (excluse coarsest grid)
        self%maxIter = maxIter
        self%numberPreSmoothOper = numberPreSmoothOper
        self%numberPostSmoothOper = numberPostSmoothOper
        self%numIter = 0
        allocate(self%GS_smoothers(self%smoothNumber), self%N_x(numberStages), self%N_y(numberStages))
        do i = 1, self%numberStages
            self%N_x(i) = (N_x + (2**(i-1) - 1))/(2**(i-1))
            self%N_y(i) = (N_y + (2**(i-1) - 1))/(2**(i-1))
            print *, 'Stage ', i, ':', self%N_x(i), self%N_y(i)
        end do
        self%directSolver = pardisoSolver(self%N_x(self%numberStages) * self%N_y(self%numberStages))
    end function MGSolver_constructor

    subroutine makeSmootherStages(self, delX, delY, NESW_wallBoundaries, boundaryConditions, omega)
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x(1)*self%N_y(1))
        real(real64), intent(in) :: delX, delY, omega
        real(real64) :: delY_temp, delX_temp, NESW_phiValTemp(4)
        integer :: stageIter, i, j, k_orig, k_temp, j_temp, i_temp
        integer, allocatable :: boundaryConditionsTemp(:)

        ! Construct first stage GS smoother (finest grid)
        self%GS_smoothers(1) = GSSolver(omega)
        call self%GS_smoothers(1)%constructPoissonEven(self%N_x(1), self%N_y(1), delX, delY, NESW_wallBoundaries, boundaryConditions)
        delY_temp = delY
        delX_temp = delX
    
        allocate(boundaryConditionsTemp(self%N_x(2) * self%N_y(2)))
        print *, 'Assembling new boundaries and generating new GS_smoothers'
        do stageIter = 2, self%smoothNumber
            delY_temp = delY_temp * 2.0d0
            delX_temp = delX_temp * 2.0d0
            do j = 1, self%N_y(1), 2**(stageIter-1)
                do i = 1, self%N_x(1), 2**(stageIter-1)
                    k_orig = (j-1) * self%N_x(1) + i
                    i_temp = (i-1)/(2**(stageIter-1)) + 1
                    j_temp = (j-1)/(2**(stageIter-1)) + 1
                    k_temp = (j_temp - 1) * self%N_x(stageIter) + i_temp
                    boundaryConditionsTemp(k_temp) = boundaryConditions(k_orig)
                end do
            end do
            self%GS_smoothers(stageIter) = GSSolver(omega)
            call self%GS_smoothers(stageIter)%constructPoissonEven(self%N_x(stageIter), self%N_y(stageIter), &
                delX_temp, delY_temp, NESW_wallBoundaries, boundaryConditionsTemp(1:(self%N_x(stageIter) * self%N_y(stageIter))))
        end do

        do stageIter = 1, self%smoothNumber-1
            ! get Restriction nodes for finer grids
            k_temp = self%GS_smoothers(stageIter+1)%numberBlackInnerNodes + self%GS_smoothers(stageIter+1)%numberRedInnerNodes ! number of total inner nodes on coarse grid
            j_temp = self%GS_smoothers(stageIter+1)%numberBlackBoundNodes + self%GS_smoothers(stageIter+1)%numberRedBoundNodes ! number of total non-dirichlet boundaries on coarse grid
            call self%GS_smoothers(stageIter)%constructRestrictionIndex(k_temp, j_temp, self%N_x(stageIter))
        end do
        
        ! Build direct pardiso solver
        delY_temp = delY_temp * 2.0d0
        delX_temp = delX_temp * 2.0d0
        NESW_phiValTemp = 0.0d0
        call buildEvenGridOrthogonalPoissonMatrix(self%N_x(self%numberStages), self%N_y(self%numberStages), delX_temp, delY_temp &
            , NESW_wallBoundaries, NESW_phiValTemp, self%directSolver%MatValues, self%directSolver%rowIndex, self%directSolver%columnIndex, self%directSolver%sourceTerm)
        call self%directSolver%initializePardiso(1, 11, 1, 0) ! use default pardiso initialization values
        ! Use directSolver to get final restriction data for final smoother
        k_temp = 0
        j_temp = 0
        do i = 1, self%directSolver%matDimension
            if (self%directSolver%rowIndex(i+1) - self%directSolver%rowIndex(i) /= 1) then
                if (MOD(i,self%N_x(self%numberStages)) <=1 .or. MOD((i-1)/self%N_x(self%numberStages) + 1,self%N_y(self%numberStages)) <=1 ) then
                    ! Add boundary nodes
                    j_temp = j_temp + 1
                else
                    ! Add inner nodes
                    k_temp = k_temp + 1
                end if
            end if
        end do
        call self%GS_smoothers(self%smoothNumber)%constructRestrictionIndex(k_temp, j_temp, self%N_x(self%smoothNumber))
    end subroutine makeSmootherStages

    subroutine orthogonalGridRestriction(self, stageInt, lowerSourceTerm, matDimension)
        ! Restriction of orthogonal grid to next coarse grid
        class(MGSolver), intent(in out) :: self
        real(real64), intent(in out) :: lowerSourceTerm(matDimension) ! source term of coarser stage
        integer(int32), intent(in) :: stageInt, matDimension !stage of finer grid to restrict
        integer(int32) :: nextStageInt, k, O_indx, N_indx, E_indx, S_indx, W_indx, O_indx_coarse, blackIdx, NE_indx, NW_indx, SE_indx, SW_indx
        real(real64) :: weightX, weightY, weightTotal, weightXY
        integer :: N_x
        N_x = self%N_x(stageInt)
        nextStageInt = stageInt+1 
        ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
        ! Weight interpolation does not seem to help, even with uneven delX and delY
        weightX = SQRT(self%GS_smoothers(stageInt)%coeffX) ! 1/delX
        weightY = SQRT(self%GS_smoothers(stageInt)%coeffY) ! 1/delY
        weightXY = 1.0d0/SQRT(1.0d0/self%GS_smoothers(stageInt)%coeffX + 1.0d0/self%GS_smoothers(stageInt)%coeffY) ! 1/sqrt(delX^2 + delY^2)
        weightTotal = 1.0d0 / (2.0d0 * weightX + 2.0d0 *weightY + 4.0d0 * weightXY)
        weightX = 0.75d0 * weightX * weightTotal ! Total weight for x dir fine nodes
        weightY = 0.75d0 * weightY * weightTotal
        weightXY = 0.75d0 * weightXY * weightTotal
        !$OMP parallel private(blackIdx, O_indx_coarse, O_indx, N_indx, E_indx, S_indx, W_indx, & 
        !$OMP NE_indx, NW_indx, SE_indx, SW_indx)
        !$OMP do
        do k = 1, self%GS_smoothers(stageInt)%numberRestrictionInnerNodes
            blackIdx = self%GS_smoothers(stageInt)%restrictionInnerIndx(1, k) ! index in black_Indx for overlapping fine grid node
            O_indx_coarse = self%GS_smoothers(stageInt)%restrictionInnerIndx(2, k) ! index in coarse grid
            O_indx = self%GS_smoothers(stageInt)%black_InnerIndx(blackIdx) ! overlapping index in fine grid
            ! calculate fine indices around overlapping index
            E_indx = O_indx + 1
            W_indx = O_indx - 1
            S_indx = O_indx - N_x
            N_indx = O_indx + N_x
            NW_indx = N_indx-1
            NE_indx = N_indx+1
            SW_indx = S_indx-1
            SE_indx = S_indx+1
            ! Interpolate residual to coarse grid
            lowerSourceTerm(O_indx_coarse) = 0.25d0 * self%GS_smoothers(stageInt)%residual(O_indx) + &
                0.125d0 * (self%GS_smoothers(stageInt)%residual(N_indx) + self%GS_smoothers(stageInt)%residual(S_indx) + &
                    self%GS_smoothers(stageInt)%residual(E_indx) + self%GS_smoothers(stageInt)%residual(W_indx)) + &
                    0.0625d0 * (self%GS_smoothers(stageInt)%residual(NE_indx) + self%GS_smoothers(stageInt)%residual(SE_indx) + &
                    self%GS_smoothers(stageInt)%residual(NW_indx) + self%GS_smoothers(stageInt)%residual(SW_indx))
            ! lowerSourceTerm(O_indx_coarse) = 0.25d0 * self%GS_smoothers(stageInt)%residual(O_indx) + &
            !     (self%GS_smoothers(stageInt)%residual(N_indx) + self%GS_smoothers(stageInt)%residual(S_indx)) * weightY + &
            !         (self%GS_smoothers(stageInt)%residual(E_indx) + self%GS_smoothers(stageInt)%residual(W_indx)) * weightX + &
            !         (self%GS_smoothers(stageInt)%residual(NE_indx) + self%GS_smoothers(stageInt)%residual(SE_indx) + &
            !         self%GS_smoothers(stageInt)%residual(NW_indx) + self%GS_smoothers(stageInt)%residual(SW_indx)) * weightXY
        end do
        !$OMP end do nowait
        !$OMP do
        do k = 1, self%GS_smoothers(stageInt)%numberRestrictionBoundNodes
            blackIdx = self%GS_smoothers(stageInt)%restrictionBoundIndx(1, k) ! index in black_Indx for overlapping fine grid node
            O_indx_coarse = self%GS_smoothers(stageInt)%restrictionBoundIndx(2, k) ! index in coarse grid
            O_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(1,blackIdx) ! overlapping index in fine grid
            ! calculate fine indices around overlapping index
            NW_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(2,blackIdx)
            N_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(3,blackIdx)
            NE_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(4,blackIdx)
            E_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(5,blackIdx)
            SE_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(6,blackIdx)
            S_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(7,blackIdx)
            SW_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(8,blackIdx)
            W_indx = self%GS_smoothers(stageInt)%black_NESW_BoundIndx(9,blackIdx)
            ! Interpolate residual to coarse grid
            lowerSourceTerm(O_indx_coarse) = 0.25d0 * self%GS_smoothers(stageInt)%residual(O_indx) + &
                0.125d0 * (self%GS_smoothers(stageInt)%residual(N_indx) + self%GS_smoothers(stageInt)%residual(S_indx) + &
                    self%GS_smoothers(stageInt)%residual(E_indx) + self%GS_smoothers(stageInt)%residual(W_indx)) + &
                    0.0625d0 * (self%GS_smoothers(stageInt)%residual(NE_indx) + self%GS_smoothers(stageInt)%residual(SE_indx) + &
                    self%GS_smoothers(stageInt)%residual(NW_indx) + self%GS_smoothers(stageInt)%residual(SW_indx))
            ! lowerSourceTerm(O_indx_coarse) = 0.25d0 * self%GS_smoothers(stageInt)%residual(O_indx) + &
            !     (self%GS_smoothers(stageInt)%residual(N_indx) + self%GS_smoothers(stageInt)%residual(S_indx)) * weightY + &
            !         (self%GS_smoothers(stageInt)%residual(E_indx) + self%GS_smoothers(stageInt)%residual(W_indx)) * weightX + &
            !         (self%GS_smoothers(stageInt)%residual(NE_indx) + self%GS_smoothers(stageInt)%residual(SE_indx) + &
            !         self%GS_smoothers(stageInt)%residual(NW_indx) + self%GS_smoothers(stageInt)%residual(SW_indx)) * weightXY
        end do
        !$OMP end do
        !$OMP end parallel
    end subroutine orthogonalGridRestriction

    subroutine orthogonalGridProlongation(self, stageInt, lowerSolution, matDimension)
        ! Restriction of orthogonal grid to next coarse grid
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: stageInt, matDimension !stage of finer grid to do correction to
        real(real64), intent(in out) :: lowerSolution(matDimension)
        real(real64) :: weightX, weightY, weightTotal, weightXY
        integer(int32) :: nextStageInt, O_indx,O_indx_coarse, i_coarse, j_coarse, N_x_fine, N_x_coarse
        weightX = SQRT(self%GS_smoothers(stageInt)%coeffX) ! 1/delX
        weightY = SQRT(self%GS_smoothers(stageInt)%coeffY) ! 1/delY
        weightXY = 1.0d0/SQRT(1.0d0/self%GS_smoothers(stageInt)%coeffX + 1.0d0/self%GS_smoothers(stageInt)%coeffY) ! 1/sqrt(delX^2 + delY^2)
        weightTotal = 1.0d0 / (2.0d0 * weightX + 2.0d0 *weightY + 4.0d0 * weightXY)
        weightX = 3.0d0 * weightX * weightTotal ! Total weight for x dir fine nodes, mult by 4 for normalization
        weightY = 3.0d0 * weightY * weightTotal
        weightXY = 3.0d0 * weightXY * weightTotal
        nextStageInt = stageInt+1 
        N_x_fine = self%N_x(stageInt)
        N_x_coarse = self%N_x(nextStageInt)
        ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
        !$OMP parallel private(i_coarse, j_coarse, O_indx_coarse, O_indx)
        ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
        !$OMP do collapse(2)
        ! set similar nodes
        do i_coarse = 1, N_x_coarse
            do j_coarse = 1, self%N_y(nextStageInt)
                O_indx_coarse = (j_coarse-1) * N_x_coarse + i_coarse
                O_indx = (2*j_coarse-2) * N_x_fine + 2*i_coarse-1
                ! Add overlapping coarse nodes
                self%GS_smoothers(stageInt)%solution(O_indx) = self%GS_smoothers(stageInt)%solution(O_indx) &
                    + lowerSolution(O_indx_coarse)
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! go horizontal each row,
        do i_coarse = 1, N_x_coarse-1
            do j_coarse = 1, self%N_y(nextStageInt)
                O_indx_coarse = (j_coarse-1) * N_x_coarse + i_coarse
                O_indx = (2*j_coarse-2) * N_x_fine + 2*i_coarse-1
                ! Average fine nodes between coarse nodes horizontally
                self%GS_smoothers(stageInt)%solution(O_indx+1) = self%GS_smoothers(stageInt)%solution(O_indx+1) + &
                    0.5d0 * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+1))
                ! self%GS_smoothers(stageInt)%solution(O_indx+1) = self%GS_smoothers(stageInt)%solution(O_indx+1) + &
                !     weightX * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+1))
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! go vertical
        do i_coarse = 1, N_x_coarse
            do j_coarse = 1, self%N_y(nextStageInt)-1
                O_indx_coarse = (j_coarse-1) * N_x_coarse + i_coarse
                O_indx = (2*j_coarse-2) * N_x_fine + 2*i_coarse-1
                ! Average fine nodes between coarse nodes vertical direction
                self%GS_smoothers(stageInt)%solution(O_indx+N_x_fine) = self%GS_smoothers(stageInt)%solution(O_indx+N_x_fine) + &
                    0.5d0 * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+N_x_coarse))
                ! self%GS_smoothers(stageInt)%solution(O_indx+N_x_fine) = self%GS_smoothers(stageInt)%solution(O_indx+N_x_fine) + &
                !     weightY * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+N_x_coarse))
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! add offset fine nodes which require 4 point interpolation
        ! first loop through SW corner coarse nodes
        do i_coarse = 1, N_x_coarse-1
            do j_coarse = 1, self%N_y(nextStageInt)-1
                O_indx_coarse = (j_coarse-1) * N_x_coarse + i_coarse
                O_indx = (2*j_coarse-1) * N_x_fine + 2*i_coarse ! index to uper right of coarse node
                self%GS_smoothers(stageInt)%solution(O_indx) =  self%GS_smoothers(stageInt)%solution(O_indx) + &
                    0.25d0 * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+1) &
                        + lowerSolution(O_indx_coarse+N_x_coarse) + lowerSolution(O_indx_coarse+N_x_coarse+1))
                ! self%GS_smoothers(stageInt)%solution(O_indx) =  self%GS_smoothers(stageInt)%solution(O_indx) + &
                !     weightXY * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+1) &
                !         + lowerSolution(O_indx_coarse+N_x_coarse) + lowerSolution(O_indx_coarse+N_x_coarse+1))
            end do
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine orthogonalGridProlongation


    subroutine MG_Cycle_Solve(self, stepTol, relTol, intSelect)
        ! V cycle for multigrid solver
        class(MGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        integer(int32), intent(in) :: intSelect
        real(real64) :: initRes
        integer(int32) :: i
        
        call self%GS_smoothers(1)%calcResidual()
        !$OMP parallel workshare
        initRes = SUM(self%GS_smoothers(1)%residual**2)
        !$OMP end parallel workshare
        self%stepResidual = self%GS_smoothers(1)%smoothIterationsWithRes(self%numberPreSmoothOper)
        
        call self%GS_smoothers(1)%calcResidual()
        !$OMP parallel workshare
        self%R2_current = SUM(self%GS_smoothers(1)%residual**2)
        !$OMP end parallel workshare
        
        i = 0
        if ((self%stepResidual > stepTol .and. self%R2_current/initRes > relTol)) then
            do i = 1, self%maxIter
                if (intSelect == 0) then
                    ! V cycle
                    call self%V_Cycle()
                else
                    ! F cycle
                    call self%F_Cycle()
                end if

                self%stepResidual = self%GS_smoothers(1)%smoothIterationsWithRes(self%numberPreSmoothOper)
                if (self%stepResidual < stepTol) then
                    print *, 'step residual lowered'
                    exit
                end if
                ! Final Residual calculation
                call self%GS_smoothers(1)%calcResidual()
                !$OMP parallel workshare
                self%R2_current = SUM(self%GS_smoothers(1)%residual**2)
                !$OMP end parallel workshare
                if (self%R2_current/initRes < relTol) then
                    print *, 'Starting residual lowered'
                    exit
                end if
            end do
        end if
        self%numIter = i-1
    end subroutine MG_Cycle_Solve

    subroutine V_Cycle(self)
        ! V cycle for multigrid solver
        class(MGSolver), intent(in out) :: self
        integer(int32) :: stageInt
        
        do stageInt = 1, self%smoothNumber-1
            ! Restrict to next smoother
            call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
            call self%GS_smoothers(stageInt+1)%smoothIterations(self%numberPreSmoothOper, .true.)
            call self%GS_smoothers(stageInt+1)%calcResidual() 
        end do
        ! Final Solver restriction and solution, then prolongation
        call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
        call self%directSolver%runPardiso()
        call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
        ! prolongation and smoothing to each additional grid
        do stageInt = self%smoothNumber, 2, -1
            ! Prolongation from each
            call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
            call self%orthogonalGridProlongation(stageInt-1, self%GS_smoothers(stageInt)%solution, self%GS_smoothers(stageInt)%matDimension)    
        end do

    end subroutine V_Cycle

    subroutine F_Cycle(self)
        ! first start with full multi-grid grid for good initial guess, then proceed with V-cycles
        class(MGSolver), intent(in out) :: self
        integer(int32) :: stageInt, j
        
        ! initialize residual to solve on fine grid (b - Ax) without smoothing
        ! then restrict all the way to coarsest grid
        ! For initial solution 0 this is same 
        ! do stageInt = 1, self%smoothNumber-1
        !     ! Reduce the sourceTerm to the lowest grid through series of restrictions, no smoothing
        !     call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
        !     !$OMP parallel
        !     !$OMP workshare
        !     self%GS_smoothers(stageInt+1)%residual = self%GS_smoothers(stageInt+1)%sourceTerm
        !     !$OMP end workshare
        !     !$OMP end parallel
        ! end do
        ! call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
        do stageInt = 1, self%smoothNumber-1
            ! Restrict to next smoother
            call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
            call self%GS_smoothers(stageInt+1)%smoothIterations(self%numberPreSmoothOper, .true.)
            call self%GS_smoothers(stageInt+1)%calcResidual() 
        end do
        ! Final Solver restriction and solution, then prolongation
        call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
        ! First solve on coarsest grid
        call self%directSolver%runPardiso()
        
        do j = self%smoothNumber, 2, -1
            ! j is each stage between coarsest and finest grid to to inverted V-cycle
            call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
            call self%GS_smoothers(self%smoothNumber)%smoothIterations(self%numberPostSmoothOper, .false.)
            do stageInt = self%smoothNumber-1, j, -1
                ! prolongate and smooth up to jth stage
                call self%orthogonalGridProlongation(stageInt, self%GS_smoothers(stageInt+1)%solution, self%GS_smoothers(stageInt+1)%matDimension)
                call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
            end do
            ! when prolongated (with) to j-th stage, calculate residual, restrict to next stage
            call self%GS_smoothers(j)%calcResidual()
            do stageInt = j, self%smoothNumber-1
                ! For each lower stage, restrict and recalculate residual
                call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
                call self%GS_smoothers(stageInt+1)%smoothIterations(self%numberPreSmoothOper, .true.)  
                call self%GS_smoothers(stageInt+1)%calcResidual()
            end do
            ! Restrict lowest grid
            call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
            ! Solve at lowest grid
            call self%directSolver%runPardiso()
        end do
        ! Now prolongate all the way to finest grid
        call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
        do stageInt = self%smoothNumber, 2, -1
            ! Prolongation from each
            call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
            call self%orthogonalGridProlongation(stageInt-1, self%GS_smoothers(stageInt)%solution, self%GS_smoothers(stageInt)%matDimension)
        end do

    end subroutine F_Cycle

    subroutine F_V_Cycle(self, stepTol, relTol)
        ! first start with full multi-grid grid for good initial guess, then proceed with V-cycles
        class(MGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        integer(int32) :: stageInt, j
        
        do stageInt = 1, self%smoothNumber-1
            ! Reduce the sourceTerm to the lowest grid through series of restrictions, no smoothing
            !$OMP parallel
            !$OMP workshare
            self%GS_smoothers(stageInt)%residual = self%GS_smoothers(stageInt)%sourceTerm
            !$OMP end workshare
            !$OMP end parallel
            call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
        end do
        !$OMP parallel
        !$OMP workshare
        self%GS_smoothers(self%smoothNumber)%residual = self%GS_smoothers(self%smoothNumber)%sourceTerm
        !$OMP end workshare
        !$OMP end parallel
        call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
        ! First solve on coarsest grid
        call self%directSolver%runPardiso()
        
        do j = self%smoothNumber, 2, -1
            ! j is each stage between coarsest and finest grid to to inverted V-cycle
            call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
            call self%GS_smoothers(self%smoothNumber)%smoothIterations(self%numberPostSmoothOper, .false.)
            do stageInt = self%smoothNumber-1, j, -1
                ! prolongate and smooth up to jth stage
                call self%orthogonalGridProlongation(stageInt, self%GS_smoothers(stageInt+1)%solution, self%GS_smoothers(stageInt+1)%matDimension)
                call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
            end do
            ! when prolongated (with) to j-th stage, calculate residual, restrict to next stage
            call self%GS_smoothers(j)%calcResidual()
            do stageInt = j, self%smoothNumber-1
                ! For each lower stage, restrict and recalculate residual
                call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
                call self%GS_smoothers(stageInt+1)%smoothIterations(self%numberPreSmoothOper, .true.)  
                call self%GS_smoothers(stageInt+1)%calcResidual()
            end do
            ! Restrict lowest grid
            call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
            ! Solve at lowest grid
            call self%directSolver%runPardiso()
        end do
        ! Now prolongate all the way to finest grid
        call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
        do stageInt = self%smoothNumber, 2, -1
            ! Prolongation from each
            call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
            call self%orthogonalGridProlongation(stageInt-1, self%GS_smoothers(stageInt)%solution, self%GS_smoothers(stageInt)%matDimension)
        end do
        
        ! No start V_Cycle
        call self%MG_Cycle_Solve(stepTol, relTol, 0)
        self%numIter = self%numIter + 1

    end subroutine F_V_Cycle




end module mod_MGSolver