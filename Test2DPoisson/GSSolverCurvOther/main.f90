

program main
    use iso_fortran_env, only: int32, real64
    use mod_GSSolver
    use mod_pardisoSolver
    use mod_MGSolver
    use mod_CSRMAtrix
    use omp_lib
    implicit none

    real(real64), parameter :: e_const = 1.602176634d-19, eps_0 = 8.8541878188d-12, pi = 4.0d0*atan(1.0d0)
    integer(int32) :: N_x = 1001, N_y = 1001, numThreads = 7
    ! type(GSSolver) :: solver
    ! type(pardisoSolver) :: directSolver
    type(MGSolver) :: solver
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate, numberPreSmoothOper, numberPostSmoothOper, numberIter
    integer :: upperBound, lowerBound, rightBound, leftBound, stageInt
    integer, allocatable :: boundaryConditions(:, :)
    real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY
    real(real64) :: alpha, beta, R2_future, R2_init, resProduct_old, resProduct_new, solutionRes, relTol, stepTol
    real(real64), allocatable :: diffX(:), diffY(:)
    logical :: makeX

    
    
    numberStages = 5
    ! More skewed delX and delY, more smoothing operations needed
    numberPreSmoothOper = 10
    numberPostSmoothOper = 10
    numberIter = 200
    omega = 1.5d0
    relTol = 1.d-8
    stepTol = 1.d-6
    rho = e_const * 1d15
    NESW_wallBoundaries(1) = 1 ! North
    NESW_wallBoundaries(2) = 1 ! East
    NESW_wallBoundaries(3) = 2 ! South
    NESW_wallBoundaries(4) = 2 ! West

    NESW_phiValues(1) = 1000.0d0
    NESW_phiValues(2) = 0.0d0
    NESW_phiValues(3) = 0.0d0
    NESW_phiValues(4) = 0.0d0

    upperPhi = NESW_phiValues(1)
    rightPhi = NESW_phiValues(2)
    lowerphi = NESW_phiValues(3)
    leftPhi = NESW_phiValues(4)

    upperBound = NESW_wallBoundaries(1)
    rightBound = NESW_wallBoundaries(2)
    lowerBound = NESW_wallBoundaries(3)
    leftBound = NESW_wallBoundaries(4)
    
    
    call mkl_set_num_threads(numThreads)
    call omp_set_num_threads(numThreads)
    call checkNodeDivisionMG(N_x, N_y, numberStages)
    allocate(diffX(N_x-1), diffY(N_y-1))
    delX = 0.01d0 * Length/real(N_x-1)
    delY = 0.01d0 * Width/real(N_y-1)
    makeX = .false.
    call createCurvGrid(N_x, Length, diffX, delX)
    makeX = .true.
    call createCurvGrid(N_y, Width, diffY, delY)
    

    allocate(boundaryConditions(N_x, N_y))

    !$OMP parallel
        !$OMP workshare
    boundaryConditions = 0
    !$OMP end workshare
    !$OMP end parallel
    boundaryConditions(2:N_x-1, N_y) = NESW_wallBoundaries(1)
    boundaryConditions(N_x, 2:N_y-1) = NESW_wallBoundaries(2)
    boundaryConditions(2:N_x-1, 1) = NESW_wallBoundaries(3)
    boundaryConditions(1, 2:N_y-1) = NESW_wallBoundaries(4)


    ! Upper right corner
    boundaryConditions(N_x, N_y) = MIN(NESW_wallBoundaries(1), NESW_wallBoundaries(2))

    ! lower right corner
    boundaryConditions(N_x, 1) = MIN(NESW_wallBoundaries(3), NESW_wallBoundaries(2))

    ! lower left corner
    boundaryConditions(1, 1) = MIN(NESW_wallBoundaries(3), NESW_wallBoundaries(4))

    ! upper left corner
    boundaryConditions(1, N_y) = MIN(NESW_wallBoundaries(1), NESW_wallBoundaries(4))

    ! solver = GSSolver(omega)
    ! call solver%constructPoissonOrthogonal(N_x, N_y, diffX, diffY, NESW_wallBoundaries, boundaryConditions)

    solver = MGSolver(N_x, N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    call solver%makeSmootherStages(diffX, diffY, NESW_wallBoundaries, boundaryConditions, omega)
    
    ! Set phi values finer grid
    if (upperBound == 1) then
        solver%GS_smoothers(1)%solution(2:N_x-1, N_y) = upperPhi
    end if
    if (rightBound == 1) then
        solver%GS_smoothers(1)%solution(N_x, 2:N_y-1) = rightPhi
    end if
    if (lowerBound == 1) then 
        solver%GS_smoothers(1)%solution(2:N_x-1, 1) = lowerPhi
    end if
    if (leftBound == 1) then 
        solver%GS_smoothers(1)%solution(1, 2:N_y-1) = leftPhi
    end if
    if (boundaryConditions(1,1) == 1) then
        if (lowerBound == leftBound) then
            solver%GS_smoothers(1)%solution(1,1) = MIN(leftPhi, lowerPhi)
        else if (lowerBound == 1) then
            solver%GS_smoothers(1)%solution(1,1) = lowerPhi
        else 
            solver%GS_smoothers(1)%solution(1,1) = leftPhi
        end if
    end if
    if (boundaryConditions(N_x,1) == 1) then
        if (lowerBound == rightBound) then
            solver%GS_smoothers(1)%solution(N_x,1) = MIN(rightPhi, lowerPhi)
        else if (lowerBound == 1) then
            solver%GS_smoothers(1)%solution(N_x,1) = lowerPhi
        else
            solver%GS_smoothers(1)%solution(N_x,1) = rightPhi
        end if
    end if
    if (boundaryConditions(1, N_y) == 1) then
        if (upperBound == leftBound) then
            solver%GS_smoothers(1)%solution(1, N_y) = MIN(leftPhi, upperPhi)
        else if (upperBound == 1) then
            solver%GS_smoothers(1)%solution(1, N_y) = upperPhi
        else
            solver%GS_smoothers(1)%solution(1, N_y) = leftPhi
        end if
    end if
    if (boundaryConditions(N_x, N_y) == 1) then
        if (upperBound == rightBound) then
            solver%GS_smoothers(1)%solution(N_x, N_y) = MIN(rightPhi, upperPhi)
        else if (upperBound == 1) then
            solver%GS_smoothers(1)%solution(N_x, N_y) = upperPhi
        else
            solver%GS_smoothers(1)%solution(N_x, N_y) = rightPhi
        end if
    end if
    
    ! directSolver = pardisoSolver(N_x * N_y)
    ! call buildCurvGridOrthogonalPoissonMatrix(N_x, N_y, diffX, diffY &
    !         , NESW_wallBoundaries, NESW_phiValues, directSolver%MatValues, directSolver%rowIndex, directSolver%columnIndex, directSolver%sourceTerm)
    ! call directSolver%initializePardiso(1, 11, 1, 0) ! use default pardiso initialization values
    ! ! set source Term first stage
    !$OMP parallel
    !$OMP do collapse(2)
    do j = 1, N_y
        do i = 1, N_x
            if (boundaryConditions(i,j) /= 1) then
                solver%GS_smoothers(1)%sourceTerm(i,j) = -rho/eps_0
                ! k = (j-1) * N_x + i
                ! directSolver%sourceTerm(k) = -rho/eps_0
            end if
        end do
    end do
    !$OMP end do
    !$OMP end parallel


    ! !$OMP parallel workshare
    ! CG_Solver%solution = CG_Solver%GS_smoothers(1)%solution
    ! !$OMP end parallel workshare
    
    ! do stageInt = 1, solver%smoothNumber-1
    !     ! Restrict to next smoother
    !     call solver%GS_smoothers(stageInt)%restriction(solver%GS_smoothers(stageInt)%sourceTerm, solver%GS_smoothers(stageInt+1)%sourceTerm)
    ! end do
    ! call solver%GS_smoothers(solver%smoothNumber)%restriction(solver%GS_smoothers(solver%smoothNumber)%sourceTerm, solver%directSolver%sourceTerm)
    ! solver%directSolver%sourceTerm(solver%directSolver%matDimension - (solver%GS_smoothers(solver%smoothNumber)%N_x+1)/2 + 1:solver%directSolver%matDimension) = 1000.0d0
    

    call system_clock(count_rate = timingRate)
    call system_clock(startTime)
    call solver%MG_Cycle_Solve(stepTol, relTol)
    call system_clock(endTime)
    ! open(41,file='test.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    ! write(41) solver%directSolver%solution
    ! close(41)
    ! stop
    open(41,file='finalSol.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) solver%GS_smoothers(1)%solution
    close(41)

    
    print *, 'took', solver%numIter, 'MG cycles'
    print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'

  
    ! call restriction(solver, test)
    
    ! call prolongation(solver, test)

    

    
    
    ! call solver%restriction(directSolver%solution, test)
    

    !call solver%prolongation(solver%GS_smoothers(1)%solution, test)
    
    ! print *, 'Took', CG_Solver%numIter, 'CG iterations'
    ! print *, 'Took', i-1, 'main iterations'





contains
    subroutine checkNodeDivisionMG(N_x, N_y, numberStages)
        ! Check to make sure N_x and N_y is divisible by however many stages in multigrid we want
        integer, intent(in out) :: N_x, N_y, numberStages
        integer :: N_min, temp_min, i
        real(real64) :: temp_min_real
        N_min = 3
        ! what is minimum grid size in X
        temp_min_real = real(N_x + (2**(numberStages-1) - 1))/real(2**(numberStages-1))
        if (temp_min_real < N_min) then
            print *, 'Number of course nodes in X direction at coarsest grid is less than the minimum!'
            stop
        else if (ABS(temp_min_real - INT(temp_min_real)) > 1.d-8) then
            print *, 'N_x will not divide that many times, will find closest highest grid'
            temp_min = INT(temp_min_real) + 1
        else
            temp_min = (N_x + (2**(numberStages-1) - 1))/(2**(numberStages-1))
        end if
        if (MOD(temp_min, 2) /= 1) then
            temp_min = temp_min + 1 
        end if
        N_x = (2**(numberStages-1)) * temp_min - (2**(numberStages-1) - 1)
        print *, 'N_x stages are:'
        do i = 1, numberStages
            print *, i, ':', (N_x + (2**(i-1) - 1))/(2**(i-1))
        end do
        

        temp_min_real = real(N_y + (2**(numberStages-1) - 1))/real(2**(numberStages-1))
        if (temp_min < N_min) then
            print *, 'Number of course nodes in Y direction at coarsest grid is less than the minimum!'
            stop
        else if (ABS(temp_min_real - INT(temp_min_real)) > 1.d-8) then
            print *, 'N_y will not divide that many times, will find closest highest grid'
            temp_min = INT(temp_min_real) + 1
        else
            temp_min = (N_y + (2**(numberStages-1) - 1))/(2**(numberStages-1))
        end if
        if (MOD(temp_min, 2) /= 1) then
            temp_min = temp_min + 1 
        end if
        N_y = (2**(numberStages-1)) * temp_min - (2**(numberStages-1) - 1)
        print *, 'N_y stages are:'
        do i = 1, numberStages
            print *, i, ':', (N_y + (2**(i-1) - 1))/(2**(i-1))
        end do

        open(41,file='NumNodes.dat', form='UNFORMATTED', access = 'stream', status = 'new')
        write(41) N_x, N_y
        close(41)

    end subroutine checkNodeDivisionMG

    subroutine createCurvGrid(N_x, Length, diffX, delX)
        ! Check to make sure N_x and N_y is divisible by however many stages in multigrid we want
        integer, intent(in) :: N_x
        real(real64), intent(in) :: Length, delX
        real(real64), intent(in out) :: diffX(N_x-1)
        real(real64) :: grid(N_x)
        grid(1) = 0.0d0
        grid(N_x) = Length
        do i = 2,N_x-1
            grid(i) = Length * ((real(i)-1.0d0)/(real(N_x) - 1.0d0) - (1.0d0/(real(N_x) - 1.0d0) - delX/Length) &
            * SIN(2.0d0 * pi * (i-1) / real(N_x - 1)) / SIN(2.0d0 * pi / real(N_x - 1)) )
        end do
        do i = 1, N_x-1
            diffX(i) = grid(i+1) - grid(i)
        end do

        if (.not. makeX) then
            open(41,file='gridX.dat', form='UNFORMATTED', access = 'stream', status = 'new')
            write(41) grid
            close(41)
        else
            open(41,file='gridY.dat', form='UNFORMATTED', access = 'stream', status = 'new')
            write(41) grid
            close(41)
        end if

    end subroutine createCurvGrid    

end program main