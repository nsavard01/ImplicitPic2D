module mod_MGSolver
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_MG_Stage
    use mod_PardisoSolver
    use mod_CSRMAtrix
    implicit none

    ! Stage for each multigrid, will define restriction and prolongation operations

    private
    public :: MGSolver

    type :: MGSolver
        type(MG_Stage), allocatable :: MG_smoothers(:)
        type(pardisoSolver) :: directSolver
        integer :: numberStages, smoothNumber, maxIter, numberPreSmoothOper, numberPostSmoothOper, numIter, N_x, N_y ! number of stages is total, so includes finest and coarsest grid
        real(real64) :: R2_current, stepResidual
    contains
        procedure, public, pass(self) :: makeSmootherStages
        ! procedure, public, pass(self) :: F_V_Cycle
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
        integer :: i, N_y_coarse, N_x_coarse
        self%numberStages = numberStages
        self%smoothNumber = numberStages-1 ! number smoothing GS stages (excluse coarsest grid)
        self%maxIter = maxIter
        self%numberPreSmoothOper = numberPreSmoothOper
        self%numberPostSmoothOper = numberPostSmoothOper
        self%numIter = 0
        self%N_x = N_x
        self%N_y = N_y
        allocate(self%MG_smoothers(self%smoothNumber))
        N_x_coarse = (self%N_x + (2**(self%numberStages-1) - 1))/(2**(self%numberStages-1))
        N_y_coarse = (self%N_y + (2**(self%numberStages-1) - 1))/(2**(self%numberStages-1))
        self%directSolver = pardisoSolver(N_x_coarse * N_y_coarse)
    end function MGSolver_constructor

    subroutine makeSmootherStages(self, diffX, diffY, NESW_wallBoundaries, boundaryConditions, omega, evenGridBool, redBlackBool)
        class(MGSolver), intent(in out) :: self
        integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x,self%N_y)
        real(real64), intent(in) :: diffX(self%N_x-1), diffY(self%N_y-1), omega
        logical, intent(in) :: evenGridBool, redBlackBool
        real(real64) :: diffX_temp(self%N_x-1), diffY_temp(self%N_y-1), NESW_phiValTemp(4), delX, delY
        integer :: stageIter, i, j, j_temp, i_temp, N_x_coarse, N_y_coarse
        integer, allocatable :: boundaryConditionsTemp(:, :)

        ! Construct first stage GS smoother (finest grid)
        self%MG_smoothers(1) = MG_Stage(omega, self%N_x, self%N_y, evenGridBool, redBlackBool)
        call self%MG_smoothers(1)%GS_Smoother%constructPoissonOrthogonal(diffX, diffY, NESW_wallBoundaries, boundaryConditions)
        N_x_coarse = self%N_x
        N_y_coarse = self%N_y
        diffX_temp = diffX
        diffY_temp = diffY
        print *, 'Assembling new boundaries and generating new GS_smoothers'
        do stageIter = 2, self%smoothNumber
            print *, 'Smoother number:', stageIter
            N_x_coarse = (N_x_coarse+1)/2
            N_y_coarse = (N_y_coarse+1)/2
            allocate(boundaryConditionsTemp(N_x_coarse, N_y_coarse))
            do i = 1, N_x_coarse-1
                diffX_temp(i) = diffX_temp(2*i-1) + diffX_temp(2*i)
            end do
            do i = 1, N_y_coarse-1
                diffY_temp(i) = diffY_temp(2*i-1) + diffY_temp(2*i)
            end do
            do j = 1, N_y_coarse
                do i = 1, N_x_coarse
                    i_temp = (2**(stageIter-1)) * i - (2**(stageIter-1)) +1
                    j_temp = (2**(stageIter-1)) * j - (2**(stageIter-1)) +1
                    boundaryConditionsTemp(i, j) = boundaryConditions(i_temp,j_temp)
                end do
            end do
            self%MG_smoothers(stageIter) = MG_Stage(omega, N_x_coarse, N_y_coarse, evenGridBool, redBlackBool)
            call self%MG_smoothers(stageIter)%GS_Smoother%constructPoissonOrthogonal(diffX_temp(1:N_x_coarse-1), diffY_temp(1:N_y_coarse-1), NESW_wallBoundaries, boundaryConditionsTemp)
            deallocate(boundaryConditionsTemp)
        end do
        
        ! Build direct pardiso solver
        NESW_phiValTemp = 0.0d0
        N_x_coarse = (N_x_coarse+1)/2
        N_y_coarse = (N_y_coarse+1)/2
        do i = 1, N_x_coarse-1
            diffX_temp(i) = diffX_temp(2*i-1) + diffX_temp(2*i)
        end do
        do i = 1, N_y_coarse-1
            diffY_temp(i) = diffY_temp(2*i-1) + diffY_temp(2*i)
        end do
        if (evenGridBool) then
            delX = diffX_temp(1)
            delY = diffY_temp(1)
            call buildEvenGridOrthogonalPoissonMatrix(N_x_coarse, N_y_coarse, delX, delY, NESW_wallBoundaries, &
            NESW_phiValTemp, self%directSolver%MatValues, self%directSolver%rowIndex, self%directSolver%columnIndex, self%directSolver%sourceTerm)
        else
            call buildCurvGridOrthogonalPoissonMatrix(N_x_coarse, N_y_coarse, diffX_temp(1:N_x_coarse-1), diffY_temp(1:N_y_coarse-1) &
                , NESW_wallBoundaries, NESW_phiValTemp, self%directSolver%MatValues, self%directSolver%rowIndex, self%directSolver%columnIndex, self%directSolver%sourceTerm)
        end if
        call self%directSolver%initializePardiso(1, 11, 1, 0) ! use default pardiso initialization values
        
    end subroutine makeSmootherStages



    subroutine MG_Cycle_Solve(self, stepTol, relTol, intSelect)
        ! V cycle for multigrid solver
        class(MGSolver), intent(in out) :: self
        integer, intent(in) :: intSelect
        real(real64), intent(in) :: stepTol, relTol
        real(real64) :: initRes
        integer(int32) :: i
        
        call self%MG_smoothers(1)%GS_Smoother%calcResidual()
        !$OMP parallel workshare
        initRes = SUM(self%MG_smoothers(1)%GS_Smoother%residual**2)
        !$OMP end parallel workshare
        call self%MG_smoothers(1)%GS_Smoother%smoothIterations(self%numberPreSmoothOper-1)
        self%stepResidual = self%MG_smoothers(1)%GS_Smoother%smoothWithRes()
        
        call self%MG_smoothers(1)%GS_Smoother%calcResidual()
        !$OMP parallel workshare
        self%R2_current = SUM(self%MG_smoothers(1)%GS_Smoother%residual**2)
        !$OMP end parallel workshare
        
        i = 0
        if ((self%stepResidual > stepTol .and. self%R2_current/initRes > relTol)) then
            do i = 1, self%maxIter
                if (intSelect == 0) then
                    call self%V_Cycle()
                else
                    call self%F_Cycle()
                end if
                call self%MG_smoothers(1)%GS_Smoother%smoothIterations(self%numberPreSmoothOper-1)
                self%stepResidual = self%MG_smoothers(1)%GS_Smoother%smoothWithRes()
                if (self%stepResidual < stepTol) then
                    print *, 'step residual lowered'
                    exit
                end if
                ! Final Residual calculation
                call self%MG_smoothers(1)%GS_Smoother%calcResidual()
                !$OMP parallel workshare
                self%R2_current = SUM(self%MG_smoothers(1)%GS_Smoother%residual**2)
                !$OMP end parallel workshare
                print *, 'MG iter:', i, 'res:', self%R2_current
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
            associate(smoother_fine => self%MG_smoothers(stageInt)%GS_smoother, smoother_coarse => self%MG_smoothers(stageInt+1)%GS_smoother)
            ! Restrict to next smoother
            call smoother_fine%restriction(smoother_fine%residual, smoother_coarse%sourceTerm)
            ! Reset solution to zero as first error guess
            !$OMP parallel workshare
            smoother_coarse%solution = 0.0d0
            !$OMP end parallel workshare
            call smoother_coarse%smoothIterations(self%numberPreSmoothOper)
            call smoother_coarse%calcResidual()
            end associate 
        end do
        ! Final Solver restriction and solution, then prolongation
        associate(smoother => self%MG_smoothers(self%smoothNumber)%GS_smoother)
        call smoother%restriction(smoother%residual, self%directSolver%sourceTerm)
        call self%directSolver%runPardiso()
        call smoother%prolongation(smoother%solution, self%directSolver%solution)
        end associate
        ! prolongation and smoothing to each additional grid
        do stageInt = self%smoothNumber, 2, -1
            ! Prolongation from each
            associate(smoother_coarse => self%MG_smoothers(stageInt)%GS_smoother, smoother_fine => self%MG_smoothers(stageInt-1)%GS_smoother)
            call smoother_coarse%smoothIterations(self%numberPostSmoothOper)
            call smoother_fine%prolongation(smoother_fine%solution, smoother_coarse%solution)
            end associate
        end do

    end subroutine V_Cycle

    subroutine F_Cycle(self)
        ! first start with full multi-grid grid for good initial guess
        class(MGSolver), intent(in out) :: self
        integer(int32) :: stageInt, j
        
        ! initialize residual to solve on fine grid (b - Ax)
        ! then restrict all the way to coarsest grid
        ! works better if restrict residual after smoothing instead of just rho
        do stageInt = 1, self%smoothNumber-1
            associate(smoother_fine => self%MG_smoothers(stageInt)%GS_smoother, smoother_coarse => self%MG_smoothers(stageInt+1)%GS_smoother)
            ! Restrict to next smoother
            call smoother_fine%restriction(smoother_fine%residual, smoother_coarse%sourceTerm)
            ! Reset solution to zero as first error guess
            !$OMP parallel workshare
            smoother_coarse%solution = 0.0d0
            !$OMP end parallel workshare
            call smoother_coarse%smoothIterations(self%numberPreSmoothOper)
            call smoother_coarse%calcResidual()
            end associate 
        end do
        ! Final Solver restriction and solution, then prolongation
        call self%MG_smoothers(self%smoothNumber)%GS_smoother%restriction(self%MG_smoothers(self%smoothNumber)%GS_smoother%residual, self%directSolver%sourceTerm)
        ! First solve on coarsest grid
        call self%directSolver%runPardiso()
        
        do j = self%smoothNumber, 2, -1
            ! j is each stage between coarsest and finest grid to to inverted V-cycle
            call self%MG_smoothers(self%smoothNumber)%GS_smoother%prolongation(self%MG_smoothers(self%smoothNumber)%GS_smoother%solution, self%directSolver%solution)
            call self%MG_smoothers(self%smoothNumber)%GS_smoother%smoothIterations(self%numberPostSmoothOper)
            do stageInt = self%smoothNumber-1, j, -1
                associate(smoother_fine => self%MG_smoothers(stageInt)%GS_smoother, smoother_coarse => self%MG_smoothers(stageInt+1)%GS_smoother)
                ! prolongate and smooth up to jth stage
                call smoother_fine%prolongation(smoother_fine%solution, smoother_coarse%solution)
                call smoother_fine%smoothIterations(self%numberPostSmoothOper)
                end associate
            end do
            ! when prolongated (with) to j-th stage, calculate residual, restrict to next stage
            call self%MG_smoothers(j)%GS_smoother%calcResidual()
            do stageInt = j, self%smoothNumber-1
                associate(smoother_fine => self%MG_smoothers(stageInt)%GS_smoother, smoother_coarse => self%MG_smoothers(stageInt+1)%GS_smoother)
                ! For each lower stage, restrict and recalculate residual
                call smoother_fine%restriction(smoother_fine%residual, smoother_coarse%sourceTerm)
                ! Reset solution to zero as first error guess
                !$OMP parallel workshare
                smoother_coarse%solution = 0.0d0
                !$OMP end parallel workshare
                call smoother_coarse%smoothIterations(self%numberPreSmoothOper)
                call smoother_coarse%calcResidual() 
                end associate
            end do
            ! Restrict lowest grid
            call self%MG_smoothers(self%smoothNumber)%GS_smoother%restriction(self%MG_smoothers(self%smoothNumber)%GS_smoother%residual, self%directSolver%sourceTerm)
            ! Solve at lowest grid
            call self%directSolver%runPardiso()
        end do
        ! Now prolongate all the way to finest grid
        call self%MG_smoothers(self%smoothNumber)%GS_smoother%prolongation(self%MG_smoothers(self%smoothNumber)%GS_smoother%solution, self%directSolver%solution)
        ! prolongation and smoothing to each additional grid
        do stageInt = self%smoothNumber, 2, -1
            ! Prolongation from each
            associate(smoother_coarse => self%MG_smoothers(stageInt)%GS_smoother, smoother_fine => self%MG_smoothers(stageInt-1)%GS_smoother)
            call smoother_coarse%smoothIterations(self%numberPostSmoothOper)
            call smoother_fine%prolongation(smoother_fine%solution, smoother_coarse%solution) 
            end associate  
        end do

    end subroutine F_Cycle

    ! subroutine F_V_Cycle(self, stepTol, relTol)
    !     ! first start with full multi-grid grid for good initial guess, then proceed with V-cycles
    !     class(MGSolver), intent(in out) :: self
    !     real(real64), intent(in) :: stepTol, relTol
    !     integer(int32) :: stageInt, j
        
    !     do stageInt = 1, self%smoothNumber-1
    !         ! Reduce the sourceTerm to the lowest grid through series of restrictions, no smoothing
    !         !$OMP parallel
    !         !$OMP workshare
    !         self%GS_smoothers(stageInt)%residual = self%GS_smoothers(stageInt)%sourceTerm
    !         !$OMP end workshare
    !         !$OMP end parallel
    !         call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
    !     end do
    !     !$OMP parallel
    !     !$OMP workshare
    !     self%GS_smoothers(self%smoothNumber)%residual = self%GS_smoothers(self%smoothNumber)%sourceTerm
    !     !$OMP end workshare
    !     !$OMP end parallel
    !     call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
    !     ! First solve on coarsest grid
    !     call self%directSolver%runPardiso()
        
    !     do j = self%smoothNumber, 2, -1
    !         ! j is each stage between coarsest and finest grid to to inverted V-cycle
    !         call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
    !         call self%GS_smoothers(self%smoothNumber)%smoothIterations(self%numberPostSmoothOper, .false.)
    !         do stageInt = self%smoothNumber-1, j, -1
    !             ! prolongate and smooth up to jth stage
    !             call self%orthogonalGridProlongation(stageInt, self%GS_smoothers(stageInt+1)%solution, self%GS_smoothers(stageInt+1)%matDimension)
    !             call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
    !         end do
    !         ! when prolongated (with) to j-th stage, calculate residual, restrict to next stage
    !         call self%GS_smoothers(j)%calcResidual()
    !         do stageInt = j, self%smoothNumber-1
    !             ! For each lower stage, restrict and recalculate residual
    !             call self%orthogonalGridRestriction(stageInt, self%GS_smoothers(stageInt+1)%sourceTerm, self%GS_smoothers(stageInt+1)%matDimension)
    !             call self%GS_smoothers(stageInt+1)%smoothIterations(self%numberPreSmoothOper, .true.)  
    !             call self%GS_smoothers(stageInt+1)%calcResidual()
    !         end do
    !         ! Restrict lowest grid
    !         call self%orthogonalGridRestriction(self%smoothNumber, self%directSolver%sourceTerm, self%directSolver%matDimension)
    !         ! Solve at lowest grid
    !         call self%directSolver%runPardiso()
    !     end do
    !     ! Now prolongate all the way to finest grid
    !     call self%orthogonalGridProlongation(self%smoothNumber, self%directSolver%solution, self%directSolver%matDimension)
    !     do stageInt = self%smoothNumber, 2, -1
    !         ! Prolongation from each
    !         call self%GS_smoothers(stageInt)%smoothIterations(self%numberPostSmoothOper, .false.)
    !         call self%orthogonalGridProlongation(stageInt-1, self%GS_smoothers(stageInt)%solution, self%GS_smoothers(stageInt)%matDimension)
    !     end do
        
    !     ! No start V_Cycle
    !     call self%MG_Cycle_Solve(stepTol, relTol, 0)
    !     self%numIter = self%numIter + 1

    ! end subroutine F_V_Cycle




end module mod_MGSolver