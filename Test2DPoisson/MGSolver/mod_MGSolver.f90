module mod_MGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_GSSolver
    use mod_PardisoSolver
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: MGSolver

    type :: MGSolver
        type(GSSolver), allocatable :: GS_smoothers(:)
        integer :: numberStages, smoothNumber, maxIter, numberSmoothOper, numIter ! number of stages is total, so includes finest and coarsest grid
        integer, allocatable :: matDimensions(:), N_x(:), N_y(:)
        real(real64) :: residual_0, residualCurrent, stepResidual
    contains
        procedure, public, pass(self) :: makeSmootherStages
        procedure, public, pass(self) :: orthogonalGridRestriction
        procedure, public, pass(self) :: orthogonalGridProlongation
        procedure, public, pass(self) :: F_V_Cycle
        procedure, public, pass(self) :: V_Cycle
    end type

    interface MGSolver
        module procedure :: MGSolver_constructor
    end interface MGSolver

contains

    type(MGSolver) function MGSolver_constructor(N_x, N_y, numberStages, maxIter, numberSmoothOper) result(self)
        ! Construct object, set initial variables
        integer(int32), intent(in) :: N_x, N_y, numberStages, maxIter, numberSmoothOper
        integer :: i
        self%numberStages = numberStages
        self%smoothNumber = numberStages
        self%maxIter = maxIter
        self%numberSmoothOper = numberSmoothOper
        self%numIter = 0
        allocate(self%GS_smoothers(self%smoothNumber), self%N_x(numberStages), self%N_y(numberStages))
        do i = 1, self%numberStages
            self%N_x(i) = (N_x + (2**(i-1) - 1))/(2**(i-1))
            self%N_y(i) = (N_y + (2**(i-1) - 1))/(2**(i-1))
            print *, 'Stage ', i, ':', self%N_x(i), self%N_y(i)
        end do
    end function MGSolver_constructor

    subroutine makeSmootherStages(self, delX, delY, NESW_wallBoundaries, boundaryConditions, omega)
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x(1)*self%N_y(1))
        real(real64), intent(in) :: delX, delY, omega
        real(real64) :: delY_temp, delX_temp
        integer :: stageIter, i, j, k_orig, k_temp, j_temp, i_temp, matDimension
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
            k_temp = self%GS_smoothers(stageIter+1)%numberBlackNodes + self%GS_smoothers(stageIter+1)%numberRedNodes ! number of total non-dirichlet boundaries on coarse grid
            call self%GS_smoothers(stageIter)%constructRestrictionIndex(k_temp, self%N_x(stageIter), self%N_y(stageIter))
        end do


    end subroutine makeSmootherStages

    subroutine orthogonalGridRestriction(self, stageInt)
        ! Restriction of orthogonal grid to next coarse grid
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: stageInt !stage of finer grid to restrict
        integer(int32) :: nextStageInt, k, O_indx, N_indx, E_indx, S_indx, W_indx, O_indx_coarse, blackIdx
        nextStageInt = stageInt+1 
        ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes

        !$OMP parallel private(blackIdx, O_indx_coarse, O_indx, N_indx, E_indx, S_indx, W_indx)
        !$OMP do
        do k = 1, self%GS_smoothers(stageInt)%numberRestrictionNodes
            blackIdx = self%GS_smoothers(stageInt)%restrictionIndx(1, k) ! index in black_NESW for overlapping fine grid node
            O_indx_coarse = self%GS_smoothers(stageInt)%restrictionIndx(2, k) ! index in coarse grid
            O_indx = self%GS_smoothers(stageInt)%black_NESW_indx(1,blackIdx)
            N_indx = self%GS_smoothers(stageInt)%black_NESW_indx(2,blackIdx)
            E_indx = self%GS_smoothers(stageInt)%black_NESW_indx(3,blackIdx)
            S_indx = self%GS_smoothers(stageInt)%black_NESW_indx(4,blackIdx)
            W_indx = self%GS_smoothers(stageInt)%black_NESW_indx(5,blackIdx)
            self%GS_smoothers(stageInt+1)%sourceTerm(O_indx_coarse) = 0.5d0 * self%GS_smoothers(stageInt)%residual(O_indx) + &
                0.125d0 * (self%GS_smoothers(stageInt)%residual(N_indx) + self%GS_smoothers(stageInt)%residual(S_indx) + &
                self%GS_smoothers(stageInt)%residual(E_indx) + self%GS_smoothers(stageInt)%residual(W_indx))
            !self%GS_smoothers(stageInt+1)%sourceTerm(O_indx_coarse) = self%GS_smoothers(stageInt)%residual(O_indx)
        end do
        !$OMP end do
        !$OMP end parallel
    end subroutine orthogonalGridRestriction

    subroutine orthogonalGridProlongation(self, stageInt)
        ! Restriction of orthogonal grid to next coarse grid
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: stageInt !stage of finer grid to do correction to
        integer(int32) :: nextStageInt, O_indx,O_indx_coarse, i_coarse, j_coarse, N_x_fine, N_x_coarse
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
                    +self%GS_smoothers(nextStageInt)%solution(O_indx_coarse)
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
                    0.5d0 * (self%GS_smoothers(nextStageInt)%solution(O_indx_coarse) + self%GS_smoothers(nextStageInt)%solution(O_indx_coarse+1))
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
                    0.5d0 * (self%GS_smoothers(nextStageInt)%solution(O_indx_coarse) + self%GS_smoothers(nextStageInt)%solution(O_indx_coarse+N_x_coarse))
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
                    0.25d0 * (self%GS_smoothers(nextStageInt)%solution(O_indx_coarse) + self%GS_smoothers(nextStageInt)%solution(O_indx_coarse+1) &
                        + self%GS_smoothers(nextStageInt)%solution(O_indx_coarse+N_x_coarse) + self%GS_smoothers(nextStageInt)%solution(O_indx_coarse+N_x_coarse+1))
            end do
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine orthogonalGridProlongation


    subroutine V_Cycle(self, stepTol, relTol)
        ! V cycle for multigrid solver
        class(MGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        integer(int32) :: stageInt, i
        
        self%stepResidual = self%GS_smoothers(1)%smoothIterationsWithRes(self%numberSmoothOper)
        call self%GS_smoothers(1)%calcResidual()
        !$OMP parallel workshare
        self%residualCurrent = SUM(self%GS_smoothers(1)%residual**2)
        !$OMP end parallel workshare
        
        i = 0
        if ((self%stepResidual > stepTol .and. self%residualCurrent > relTol)) then
            do i = 1, self%maxIter
                call self%orthogonalGridRestriction(1)
                do stageInt = 2, self%numberStages-1
                    ! Series of smoothers and restriction for each stage
                    call self%GS_smoothers(stageInt)%smoothIterations(self%numberSmoothOper, .true.)  
                    call self%GS_smoothers(stageInt)%calcResidual()
                    call self%orthogonalGridRestriction(stageInt)
                end do
                ! Final Solver
                call self%GS_smoothers(self%numberStages)%solveGS(stepTol)
                do stageInt = self%numberStages-1, 2, -1
                    ! Prolongation from each
                    call self%orthogonalGridProlongation(stageInt)
                    call self%GS_smoothers(stageInt)%smoothIterations(self%numberSmoothOper, .false.)
                end do
                call self%orthogonalGridProlongation(1)
                self%stepResidual = self%GS_smoothers(1)%smoothIterationsWithRes(self%numberSmoothOper)
                ! Final Residual calculation
                call self%GS_smoothers(1)%calcResidual()
                !$OMP parallel workshare
                self%residualCurrent = SUM(self%GS_smoothers(1)%residual**2)
                !$OMP end parallel workshare
                if (self%stepResidual < stepTol .or. self%residualCurrent < relTol) exit
            end do
        end if
        self%numIter = i
    end subroutine V_Cycle

    subroutine F_V_Cycle(self, stepTol, relTol)
        ! first start with full multi-grid grid for good initial guess, then proceed with V-cycles
        class(MGSolver), intent(in out) :: self
        real(real64), intent(in) :: stepTol, relTol
        integer(int32) :: stageInt, j
        
        do stageInt = 1, self%numberStages-1
            ! Reduce the sourceTerm to the lowest grid through series of restrictions, no smoothing
            !$OMP parallel
            !$OMP workshare
            self%GS_smoothers(stageInt)%residual = self%GS_smoothers(stageInt)%sourceTerm
            !$OMP end workshare
            !$OMP end parallel
            call self%orthogonalGridRestriction(stageInt)
        end do

        ! First solve on coarsest grid
        call self%GS_smoothers(self%numberStages)%solveGS(stepTol)
        do j = self%numberStages-1, 2, -1
            ! j is each stage between coarsest and finest grid to to inverted V-cycle
            do stageInt = self%numberStages-1, j, -1
                ! prolongate and smooth up to jth stage
                call self%orthogonalGridProlongation(stageInt)
                call self%GS_smoothers(stageInt)%smoothIterations(self%numberSmoothOper, .false.)
            end do
            ! when prolongated (with) to j-th stage, calculate residual, restrict to next stage
            call self%GS_smoothers(j)%calcResidual()
            call self%orthogonalGridRestriction(j)
            do stageInt = j+1, self%numberStages-1
                ! For each lower stage, restrict and recalculate residual
                call self%GS_smoothers(stageInt)%smoothIterations(self%numberSmoothOper, .true.)  
                call self%GS_smoothers(stageInt)%calcResidual()
                call self%orthogonalGridRestriction(stageInt)
            end do
            ! Solve at lowest grid
            call self%GS_smoothers(self%numberStages)%solveGS(stepTol)
        end do
        ! Now prolongate all the way to finest grid
        do j = self%numberStages-1, 2, -1
            ! Prolongation from each
            call self%orthogonalGridProlongation(j)
            call self%GS_smoothers(j)%smoothIterations(2, .false.)
        end do
        call self%orthogonalGridProlongation(1)
        
        ! No start V_Cycle
        call self%V_Cycle(stepTol, relTol)
        self%numIter = self%numIter + 1

    end subroutine F_V_Cycle




end module mod_MGSolver