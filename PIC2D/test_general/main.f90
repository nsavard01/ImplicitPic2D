program main
    use iso_fortran_env, only: int32, real64
    use constants
    use mod_CSRMAtrix
    use mod_pardisoSolver
    use mod_MGSolver
    use mod_PreCondCGSolver
    use mod_PreCondBiCGSTABSolver
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_particle
    use mod_rand_generator
    use omp_lib
    implicit none

<<<<<<< HEAD
    real(real64), parameter :: e_const = 1.602176634d-19, eps_0 = 8.8541878188d-12, pi = 4.0d0*atan(1.0d0)
    integer(int32) :: N_x = 1001, N_y = 1001, numThreads = 6
=======
    integer(int32) :: N_x = 1001, N_y = 2001, numThreads = 6
    type(Particle), allocatable :: particle_list(:)
>>>>>>> test1
    class(domain_base), allocatable, target :: world
    class(MGSolver), allocatable :: mg_solver
    type(rand_gen), allocatable :: random_gen(:)
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate, numberPreSmoothOper, numberPostSmoothOper, numberIter
    integer :: upperBound, lowerBound, rightBound, leftBound, stageInt, curv_grid_type_x, curv_grid_type_y, mat_dimension
    integer :: inner_box_first_y, inner_box_last_y, inner_box_first_x, inner_box_last_x, i_thread
    real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi, innerPhi
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY
    real(real64) :: relTol, stepTol, temp_real
    logical :: evenGridBool, redBlackBool, Krylov_bool, center_box_bool

    call execute_command_line("rm -r *.dat")
    call mkl_set_num_threads(numThreads)
    call omp_set_num_threads(numThreads)
    call omp_set_max_active_levels(numThreads)

    ! call change_global_thread(numThreads)
    
    
    evenGridBool = .false.
    redBlackBool = .false.
    Krylov_bool = .false.
    center_box_bool = .true.
    curv_grid_type_x = 0
    curv_grid_type_y = 0
    
    numberStages = 5
    call checkNodeDivisionMG(N_x, N_y, numberStages)
    ! call change_global_N(N_x, N_y)

    allocate(random_gen(numThreads))
    call random_seed()
    do i = 1, numThreads
        call random_number(temp_real)
        random_gen(i) = rand_gen(INT(temp_real * (huge(i)-1)) + 1)
    end do

    ! More skewed delX and delY, more smoothing operations needed
    numberPreSmoothOper = 8
    numberPostSmoothOper = 8
    numberIter = 50000
    omega = 1.5d0
    relTol = 1.d-8
    stepTol = 1.d-6
    rho = e_const * 1d15
    NESW_wallBoundaries(1) = 2 ! North
    NESW_wallBoundaries(2) = 2 ! East
    NESW_wallBoundaries(3) = 1 ! South
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
        inner_box_first_x = 3* N_x/8
        inner_box_last_x = 5 * N_x / 8
        inner_box_first_y = 3 * N_y / 8
        inner_box_last_y = 5 * N_y / 8
        !$OMP parallel
        !$OMP workshare
        world%boundary_conditions(inner_box_first_x:inner_box_last_x, inner_box_first_y:inner_box_last_y) = 1
        !$OMP end workshare
        !$OMP end parallel
    end if

    if (Krylov_bool) then
        if (evenGridBool) then
            mg_solver = PreCondCGSolver(world%N_x, world%N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
        else
            mg_solver = BiCGSTAB_Solver(world%N_x, world%N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
        end if
    else
        mg_solver = MGSolver(world%N_x, world%N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    end if
    call mg_solver%makeSmootherStages(world, omega, redBlackBool)

    associate(solver => mg_solver%MG_smoothers(1)%GS_smoother)

    open(41,file='gridX.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) world%grid_X
    close(41)

    open(41,file='gridY.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) world%grid_Y
    close(41)

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

    !$OMP parallel private(k, i, j)
    !$OMP do collapse(2)
    do j = 1, solver%N_y
        do i = 1, solver%N_x
            k = (j-1) * N_x + i
            if (world%boundary_conditions(i,j) /= 1) then
                ! stageOne%sourceTerm(i,j) = -rho/eps_0
                solver%sourceTerm(i,j) = -rho/epsilon_0
            end if
        end do
    end do
    !$OMP end do
    !$OMP end parallel


    call system_clock(count_rate = timingRate)
    call system_clock(startTime)
    call mg_solver%solve(stepTol, relTol, 1)
    call system_clock(endTime)
    

    print *, 'Took', mg_solver%numIter, 'iterations'
    print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'
    
    open(41,file='finalSol.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) solver%solution
    close(41)
    end associate
 





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

    ! subroutine readChargedParticleInputs(filename, irand, T_e, T_i, numThread, world, particleList)
    !     ! Read input file for particles
    !     type(Particle), allocatable, intent(out) :: particleList(:)
    !     type(Domain) :: world
    !     character(len=*), intent(in) :: filename
    !     integer(int32), intent(in) :: numThread
    !     integer(int32), intent(in out) :: irand(numThread)
    !     real(real64), intent(in) :: T_e, T_i
    !     integer(int32) :: j, io, numSpecies = 0, numParticles(100), particleIdxFactor(100), i, tempInt, distType(100), iThread
    !     character(len=15) :: name
    !     character(len=8) :: particleNames(100), char_i
    !     real(real64) :: mass(100), charge(100), Ti(100), tempReal, alpha(100), v_drift(100)
    !     logical :: boolVal

    !     print *, "Reading particle inputs:"
    !     if (.not. restartBool) then
    !         open(10,file='../InputData/'//filename, action = 'read')
    !     else
    !         open(10,file=restartDirectory//'/InputData/'//filename, action = 'read')
    !     end if
    !     do j=1, 10000
    !         read(10,*) name

    !         if( name(1:9).eq.'ELECTRONS') then
    !             read(10,*) name
    !             read(10,*) name
    !             read(10,'(A4)', ADVANCE = 'NO') name(1:2)
    !             numSpecies = numSpecies + 1
    !             read(10,*) numParticles(numSpecies), particleIdxFactor(numSpecies), alpha(numSpecies), distType(numSpecies), v_drift(numSpecies)
    !             Ti(numSpecies) = T_e
    !             mass(numSpecies) = m_e
    !             charge(numSpecies) = -1.0
    !             particleNames(numSpecies) = 'e'
    !             read(10,*) name
    !             read(10,*) name
    !         endif


    !         if(name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
    !             do while(name(1:4).ne.'----')
    !                 read(10,*) name
    !             end do
    !             read(10,'(A6)', ADVANCE = 'NO') name
    !             do while (name(1:4).ne.'----')
    !                 numSpecies = numSpecies + 1
    !                 read(10,*) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies), alpha(numSpecies), distType(numSpecies), v_drift(numSpecies)
    !                 Ti(numSpecies) = T_i
    !                 mass(numSpecies) = mass(numSpecies) * m_amu - charge(numSpecies) * m_e
    !                 particleNames(numSpecies) = trim(name)
    !                 read(10,'(A6)', ADVANCE = 'NO') name
    !             end do
    !         endif      

    !         if (name(1:7) == 'ENDFILE') then
    !             exit
    !         end if

    !     end do
    !     close(10)

    !     numberChargedParticles = numSpecies
    !     print *, 'Amount charged particles:', numberChargedParticles
    !     if (numberChargedParticles > 0) then
    !         ! Initialize and generate particles
    !         allocate(particleList(numberChargedParticles))
    !         do j=1, numberChargedParticles
    !             particleList(j) = Particle(mass(j), e * charge(j), 1.0d0, numParticles(j), numParticles(j) * particleIdxFactor(j), trim(particleNames(j)), numThread)
    !             call particleList(j) % initialize_n_ave(n_ave, world%L_domain)
    !             if (.not. restartBool) then
    !                 if (j==2 .and. numberChargedParticles == 2 .and. charge(2) == -charge(1) .and. numParticles(1) == numParticles(2)) then
    !                     ! If only ions and electrons (electrons come first) then set ion positions same as electrons for neutral start
    !                     print *, 'Neutral charge start!'
    !                     particleList(2)%phaseSpace(1, :, :) = particleList(1)%phaseSpace(1,:,:)
    !                 else
    !                     SELECT CASE(distType(j))
    !                     CASE(0)
    !                         call particleList(j)% initializeRandUniform(world, irand)
    !                     CASE(1)
    !                         call particleList(j)% initializeRandCosine(world, irand, alpha(j))
    !                     CASE(2)
    !                         call particleList(j)% initializeRandSine(world, irand, alpha(j))
    !                     CASE default
    !                         print *, 'Distribution type should be between 0 and 2!'
    !                         stop
    !                     END SELECT
    !                 end if
    !                 if (j == 1) then
    !                     call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, alpha(j), distType(j), v_drift(j))
    !                 else
    !                     call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, 1.236d0, distType(j), v_drift(j)) ! Used for IASW case, will likely change to general later
    !                 end if
    !             else
    !                 !$OMP parallel private(iThread, boolVal, char_i, io, i)
    !                 iThread = omp_get_thread_num() + 1
    !                 write(char_i, '(I3)'), iThread
    !                 INQUIRE(file=restartDirectory//'/PhaseSpace/phaseSpace_'//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", exist = boolVal)
    !                 if (.not. boolVal) then
    !                     print *, restartDirectory//'/PhaseSpace/phaseSpace_'//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", 'Does not exist'
    !                     stop
    !                 end if
    !                 open(iThread,file=restartDirectory//"/PhaseSpace/phaseSpace_"//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", form = 'UNFORMATTED', access = 'stream', status = 'old', IOSTAT=io)
    !                 i = 0
    !                 read(iThread,  IOSTAT = io) particleList(j)%phaseSpace(:, i+1, iThread)
    !                 do while (io == 0)
    !                     i = i + 1
    !                     read(iThread,  IOSTAT = io) particleList(j)%phaseSpace(:, i+1, iThread)
    !                 end do
    !                 particleList(j)%N_p(iThread) = i
    !                 close(iThread)
    !                 !$OMP end parallel
    !             end if
    !             print *, 'Initializing ', particleList(j) % name
    !             print *, 'Amount of macroparticles is:', SUM(particleList(j) % N_p)
    !             print *, "Particle mass is:", particleList(j)%mass
    !             print *, "Particle charge is:", particleList(j)%q
    !             print *, "Particle weight is:", particleList(j)%w_p
    !             print *, "Particle mean KE is:", particleList(j)%getKEAve()
    !             print *, 'Distribution type:', distType(j)
    !             print *, 'Drift velocity:', v_drift(j)
    !         end do
    !     end if

        
    !     print *, "---------------"
    !     print *, ""


    ! end subroutine readChargedParticleInputs



end program main