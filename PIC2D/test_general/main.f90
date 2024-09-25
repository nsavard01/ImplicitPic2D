program main
    use iso_fortran_env, only: int32, real64
    use mod_CSRMAtrix
    use mod_pardisoSolver
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_GS_Base_Even
    use mod_RedBlackSolverEven
    use mod_ZebraSolverEven
    use omp_lib
    implicit none

    real(real64), parameter :: e_const = 1.602176634d-19, eps_0 = 8.8541878188d-12, pi = 4.0d0*atan(1.0d0)
    integer(int32) :: N_x = 401, N_y = 201, numThreads = 32
    type(pardisoSolver) :: directSolver
    class(domain_base), allocatable, target :: world
    class(GS_Base_Even), allocatable :: solver
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate, numberPreSmoothOper, numberPostSmoothOper, numberIter
    integer :: upperBound, lowerBound, rightBound, leftBound, stageInt, curv_grid_type_x, curv_grid_type_y, mat_dimension
    integer :: inner_box_last_y = 125, inner_box_first_y = 75, inner_box_first_x = 150, inner_box_last_x = 250
    real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi, innerPhi
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY
    real(real64) :: alpha, beta, R2_future, R2_init, resProduct_old, resProduct_new, solutionRes, relTol, stepTol
    logical :: evenGridBool, redBlackBool, PCG_bool, Krylov_bool, center_box_bool

    call execute_command_line("rm -r *.dat")
    call mkl_set_num_threads(numThreads)
    call omp_set_num_threads(numThreads)
    call omp_set_max_active_levels(numThreads)
    
    evenGridBool = .true.
    redBlackBool = .true.
    PCG_bool = .false.
    Krylov_bool = .false.
    center_box_bool = .true.
    curv_grid_type_x = 0
    curv_grid_type_y = 0
    
    numberStages = 1
    call checkNodeDivisionMG(N_x, N_y, numberStages)

    ! More skewed delX and delY, more smoothing operations needed
    numberPreSmoothOper = 4
    numberPostSmoothOper = 4
    numberIter = 50000
    omega = 1.5d0
    relTol = 1.d-12
    stepTol = 1.d-6
    rho = e_const * 1d15
    NESW_wallBoundaries(1) = 3 ! North
    NESW_wallBoundaries(2) = 1 ! East
    NESW_wallBoundaries(3) = 3 ! South
    NESW_wallBoundaries(4) = 1 ! West

    NESW_phiValues(1) = 0.0d0
    NESW_phiValues(2) = 0.0d0
    NESW_phiValues(3) = 0.0d0
    NESW_phiValues(4) = 0.0d0
    innerPhi = 1000.0d0
    
    
    upperPhi = NESW_phiValues(1)
    rightPhi = NESW_phiValues(2)
    lowerphi = NESW_phiValues(3)
    leftPhi = NESW_phiValues(4)

    upperBound = NESW_wallBoundaries(1)
    rightBound = NESW_wallBoundaries(2)
    lowerBound = NESW_wallBoundaries(3)
    leftBound = NESW_wallBoundaries(4)
    
    delX = 0.01d0 * Length/real(N_x-1)
    delY = 0.01d0 * Width/real(N_y-1)
    directSolver = pardisoSolver(N_x * N_y)
    if (evenGridBool) then
        world = domain_uniform(N_x, N_y, Length, Width)
    else
        world = domain_curv(N_x, N_y, Length, &
        Width, curv_grid_type_x, curv_grid_type_y, delX, delY)
    end if
    mat_dimension = world%N_x * world%N_y
    world%boundary_conditions(2:world%N_x-1, world%N_y) = upperBound
    world%boundary_conditions(world%N_x, 2:world%N_y-1) = rightBound
    world%boundary_conditions(2:world%N_x-1, 1) = lowerBound
    world%boundary_conditions(1, 2:world%N_y-1) = leftBound


    ! Upper right corner
    world%boundary_conditions(world%N_x, world%N_y) = MIN(upperBound, rightBound)

    ! lower right corner
    world%boundary_conditions(world%N_x, 1) = MIN(lowerBound, rightBound)

    ! lower left corner
    world%boundary_conditions(1, 1) = MIN(lowerBound, leftBound)

    ! upper left corner
    world%boundary_conditions(1, world%N_y) = MIN(upperBound, leftBound)

    if (center_box_bool) then
        !$OMP parallel
        !$OMP workshare
        world%boundary_conditions(inner_box_first_x:inner_box_last_x, inner_box_first_y:inner_box_last_y) = 1
        !$OMP end workshare
        !$OMP end parallel
    end if


    select type (world)
    type is (domain_uniform)
        call buildEvenGridOrthogonalPoissonMatrix(world%N_x, world%N_y, world%del_x, world%del_y, &
        world%boundary_conditions, directSolver%MatValues, directSolver%rowIndex, directSolver%columnIndex, directSolver%sourceTerm)
        if (redBlackBool) then
            solver = RedBlackSolverEven(omega, world, world%N_x, world%N_y)
        else
            solver = ZebraSolverEven(omega, world, world%N_x, world%N_y)
        end if
    end select
    call directSolver%initializePardiso(1, 11, 1, 0)

    open(41,file='gridX.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) world%grid_X
    close(41)

    open(41,file='gridY.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) world%grid_Y
    close(41)
    
    if (upperBound == 1 .and. upperPhi /= 0) directSolver%solution(mat_dimension - world%N_x + 1: mat_dimension) = upperPhi
    if (rightBound ==1 .and. rightPhi /=0) directSolver%solution(world%N_x:mat_dimension:world%N_x) = rightPhi
    if (lowerBound == 1 .and. lowerPhi /= 0) directSolver%solution(1:world%N_x) = lowerPhi
    if (leftBound == 1 .and. leftPhi /= 0) directSolver%solution(1:mat_dimension - world%N_x + 1:world%N_x) = leftPhi

    if (upperBound == 1 .and. upperPhi == 0) directSolver%solution(mat_dimension - world%N_x + 1: mat_dimension) = upperPhi
    if (rightBound ==1 .and. rightPhi ==0) directSolver%solution(world%N_x:mat_dimension:world%N_x) = rightPhi
    if (lowerBound == 1 .and. lowerPhi == 0) directSolver%solution(1:world%N_x) = lowerPhi
    if (leftBound == 1 .and. leftPhi == 0) directSolver%solution(1:mat_dimension - world%N_x + 1:world%N_x) = leftPhi

    if (upperBound == 1 .and. upperPhi /= 0) solver%solution(:, solver%N_y) = upperPhi
    if (rightBound ==1 .and. rightPhi /=0) solver%solution(solver%N_x, :) = rightPhi
    if (lowerBound == 1 .and. lowerPhi /= 0) solver%solution(:, 1) = lowerPhi
    if (leftBound == 1 .and. leftPhi /= 0) solver%solution(1,:) = leftPhi

    if (upperBound == 1 .and. upperPhi == 0) solver%solution(:, solver%N_y)  = upperPhi
    if (rightBound ==1 .and. rightPhi ==0) solver%solution(solver%N_x, :) = rightPhi
    if (lowerBound == 1 .and. lowerPhi == 0) solver%solution(:, 1) = lowerPhi
    if (leftBound == 1 .and. leftPhi == 0) solver%solution(1,:)  = leftPhi

    if (center_box_bool) then
        !$OMP parallel private(k)
        !$OMP do collapse(2)
        do j = inner_box_first_y, inner_box_last_y
            do i = inner_box_first_x, inner_box_last_x
                k = (j-1) * world%N_x + i 
                directSolver%solution(k) = innerPhi
                solver%solution(i,j) = innerPhi
            end do
        end do
        !$OMP end do
        !$OMP end parallel
    end if
    ! if (Krylov_bool) then
    !     if (PCG_bool) then
    !         solver = PreCondCGSolver(N_x, N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    !     else
    !         solver = BiCGSTAB_Solver(N_x, N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    !     end if
    ! else
    !     solver = MGSolver(N_x, N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    ! end if
    ! select type (world)
    ! type is (domain_uniform)
    !     call solver%makeSmootherStages_even(world%del_x, world%del_y, world%NESW_wall_boundaries, world%boundary_conditions, omega, redBlackBool)
    ! type is (domain_curv)
    !     call solver%makeSmootherStages_curv(world%del_x, world%del_y, world%NESW_wall_boundaries, world%boundary_conditions, omega, redBlackBool)
    ! end select
    ! associate( stageOne => solver%MG_smoothers(1)%GS_smoother)
    ! Set phi values finer grid

    ! if (upperBound == 1) then
    !     stageOne%solution(2:N_x-1, N_y) = upperPhi
    ! end if
    ! if (rightBound == 1) then
    !     stageOne%solution(N_x, 2:N_y-1) = rightPhi
    ! end if
    ! if (lowerBound == 1) then 
    !     stageOne%solution(2:N_x-1, 1) = lowerPhi
    ! end if
    ! if (leftBound == 1) then 
    !     stageOne%solution(1, 2:N_y-1) = leftPhi
    ! end if
    ! if (world%boundary_conditions(1,1) == 1) then
    !     if (lowerBound == leftBound) then
    !         stageOne%solution(1,1) = MIN(leftPhi, lowerPhi)
    !     else if (lowerBound == 1) then
    !         stageOne%solution(1,1) = lowerPhi
    !     else 
    !         stageOne%solution(1,1) = leftPhi
    !     end if
    ! end if
    ! if (world%boundary_conditions(N_x,1) == 1) then
    !     if (lowerBound == rightBound) then
    !         stageOne%solution(N_x,1) = MIN(rightPhi, lowerPhi)
    !     else if (lowerBound == 1) then
    !         stageOne%solution(N_x,1) = lowerPhi
    !     else
    !         stageOne%solution(N_x,1) = rightPhi
    !     end if
    ! end if
    ! if (world%boundary_conditions(1, N_y) == 1) then
    !     if (upperBound == leftBound) then
    !         stageOne%solution(1, N_y) = MIN(leftPhi, upperPhi)
    !     else if (upperBound == 1) then
    !         stageOne%solution(1, N_y) = upperPhi
    !     else
    !         stageOne%solution(1, N_y) = leftPhi
    !     end if
    ! end if
    ! if (world%boundary_conditions(N_x, N_y) == 1) then
    !     if (upperBound == rightBound) then
    !         stageOne%solution(N_x, N_y) = MIN(rightPhi, upperPhi)
    !     else if (upperBound == 1) then
    !         stageOne%solution(N_x, N_y) = upperPhi
    !     else
    !         stageOne%solution(N_x, N_y) = rightPhi
    !     end if
    ! end if

    !$OMP parallel
    !$OMP do collapse(2)
    do j = 1, N_y
        do i = 1, N_x
            k = (j-1) * N_x + i
            if (world%boundary_conditions(i,j) /= 1) then
                ! stageOne%sourceTerm(i,j) = -rho/eps_0
                solver%sourceTerm(i,j) = -rho/eps_0
                directSolver%sourceTerm(k) = -rho/eps_0
            else
                directSolver%sourceTerm(k) = directSolver%solution(k)
            end if
        end do
    end do
    !$OMP end do
    !$OMP end parallel

    
    
    call system_clock(count_rate = timingRate)
    call system_clock(startTime)
    call solver%solveGS(stepTol)
    call system_clock(endTime)

    ! print *, 'Took', solver%numIter, 'iterations'
    print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'
    print *, 'Took', solver%iterNumber, 'iterations'
  
    open(41,file='finalSol.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) solver%solution
    close(41)
   
    ! end associate
    ! ! !$OMP parallel workshare
    ! ! CG_Solver%solution = CG_Solver%GS_smoothers(1)%solution
    ! ! !$OMP end parallel workshare
    
    ! ! do stageInt = 1, solver%smoothNumber-1
    ! !     ! Restrict to next smoother
    ! !     call solver%GS_smoothers(stageInt)%restriction(solver%GS_smoothers(stageInt)%sourceTerm, solver%GS_smoothers(stageInt+1)%sourceTerm)
    ! ! end do
    ! ! call solver%GS_smoothers(solver%smoothNumber)%restriction(solver%GS_smoothers(solver%smoothNumber)%sourceTerm, solver%directSolver%sourceTerm)
    ! ! solver%directSolver%sourceTerm(solver%directSolver%matDimension - (solver%GS_smoothers(solver%smoothNumber)%N_x+1)/2 + 1:solver%directSolver%matDimension) = 1000.0d0
    

    
    ! ! open(41,file='test.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    ! ! write(41) solver%GS_smoothers(1)%directSolver%solution
    ! ! close(41)
    ! ! stop

    

    
    ! print *, 'took', solver%iterNumber
    ! print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'

  
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
        call execute_command_line("rm -r NumNodes.dat")
        open(41,file='NumNodes.dat', form='UNFORMATTED', access = 'stream', status = 'new')
        write(41) N_x, N_y
        close(41)

    end subroutine checkNodeDivisionMG

    ! subroutine createCurvGrid(N_x, Length, diffX, delX, evenGridBool)
    !     ! Check to make sure N_x and N_y is divisible by however many stages in multigrid we want
    !     integer, intent(in) :: N_x
    !     real(real64), intent(in) :: Length, delX
    !     real(real64), intent(in out) :: diffX(N_x-1)
    !     logical, intent(in) :: evenGridBool
    !     real(real64) :: grid(N_x), tempReal
    !     grid(1) = 0.0d0
    !     grid(N_x) = Length
    !     if (evenGridBool) then
    !         tempReal = Length/real(N_x-1)
    !         do i = 2,N_x-1
    !             grid(i) = grid(i-1) + tempReal
    !         end do
    !     else
    !         do i = 2,N_x-1
    !             grid(i) = Length * ((real(i)-1.0d0)/(real(N_x) - 1.0d0) - (1.0d0/(real(N_x) - 1.0d0) - delX/Length) &
    !             * SIN(2.0d0 * pi * (i-1) / real(N_x - 1)) / SIN(2.0d0 * pi / real(N_x - 1)) )
    !         end do
    !     end if
    !     do i = 1, N_x-1
    !         diffX(i) = grid(i+1) - grid(i)
    !     end do

    !     if (.not. makeX) then
    !         call execute_command_line("rm -r gridX.dat")
    !         open(41,file='gridX.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    !         write(41) grid
    !         close(41)
    !     else
    !         call execute_command_line("rm -r gridY.dat")
    !         open(41,file='gridY.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    !         write(41) grid
    !         close(41)
    !     end if

    ! end subroutine createCurvGrid    


end program main