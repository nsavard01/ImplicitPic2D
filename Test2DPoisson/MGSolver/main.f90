

program main
    use iso_fortran_env, only: int32, real64
    use mod_MGSolver
    use mod_GSSolver
    use omp_lib
    implicit none

    real(real64), parameter :: e_const = 1.602176634d-19, eps_0 = 8.8541878188d-12
    integer(int32) :: N_x = 141, N_y = 141, numThreads = 6
    type(MGSolver) :: MG_Solver
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate
    integer, allocatable :: boundaryConditions(:)
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY

    numberStages = 2
    omega = 1.6d0
    rho = e_const * 1d15
    NESW_wallBoundaries(1) = 1 ! North
    NESW_wallBoundaries(2) = 3 ! East
    NESW_wallBoundaries(3) = 1 ! South
    NESW_wallBoundaries(4) = 3 ! West

    NESW_phiValues(1) = 0.0d0
    NESW_phiValues(2) = 0.0d0
    NESW_phiValues(3) = 1000.0d0
    NESW_phiValues(4) = 0.0d0
    
    
    call mkl_set_num_threads(numThreads)
    call omp_set_num_threads(numThreads)
    call checkNodeDivisionMG(N_x, N_y, numberStages)
    delX = Length/(N_x-1)
    delY = Width/(N_y-1)
    matDimension = N_x * N_y
    allocate(boundaryConditions(matDimension))

    !$OMP parallel
        !$OMP workshare
    boundaryConditions = 0
    !$OMP end workshare
    !$OMP end parallel
    boundaryConditions(matDimension-N_x+2:matDimension-1) = NESW_wallBoundaries(1)
    boundaryConditions(2*N_x:matDimension-N_x:N_x) = NESW_wallBoundaries(2)
    boundaryConditions(2:N_x-1) = NESW_wallBoundaries(3)
    boundaryConditions(N_x+1:matDimension-2*N_x + 1:N_x) = NESW_wallBoundaries(4)


    ! Upper right corner
    boundaryConditions(matDimension) = MIN(NESW_wallBoundaries(1), NESW_wallBoundaries(2))

    ! lower right corner
    boundaryConditions(N_x) = MIN(NESW_wallBoundaries(3), NESW_wallBoundaries(2))

    ! lower left corner
    boundaryConditions(1) = MIN(NESW_wallBoundaries(3), NESW_wallBoundaries(4))

    ! upper left corner
    boundaryConditions(matDimension-N_x+1) = MIN(NESW_wallBoundaries(1), NESW_wallBoundaries(4))


    MG_Solver = MGSolver(N_x, N_y, numberStages)
    call MG_Solver%makeSmootherStages(delX, delY, NESW_wallBoundaries, NESW_phiValues, boundaryConditions, omega)
    !$OMP parallel
    !$OMP do
    do i = 1, matDimension
        if (boundaryConditions(i) /= 1) then
            MG_Solver%GS_smoothers(1)%sourceTerm(i) = -rho/eps_0
        end if
    end do
    !$OMP end do
    !$OMP end parallel
    call system_clock(count_rate = timingRate)
    call system_clock(startTime)
    call MG_Solver%GS_smoothers(1)%smoothIterations(2)
    call system_clock(endTime)
    print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'
    call MG_Solver%GS_smoothers(1)%calcResidual()
    open(41,file='GSRes.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) MG_Solver%GS_smoothers(1)%sourceTerm
    close(41)
    call MG_Solver%orthogonalGridRestriction(1)
    open(41,file='GSRes_1.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) MG_Solver%GS_smoothers(2)%sourceTerm
    close(41)



contains
    subroutine checkNodeDivisionMG(N_x, N_y, numberStages)
        ! Check to make sure N_x and N_y is divisible by however many stages in multigrid we want
        integer, intent(in out) :: N_x, N_y, numberStages
        integer :: N_min, temp_min
        real(real64) :: temp_min_real
        N_min = 3
        ! what is minimum grid size in X
        temp_min_real = real(N_x + (2**(numberStages-1) - 1))/real(2**(numberStages-1))
        if (temp_min < N_min) then
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
        print *, 'N_x is:', N_x

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
        print *, 'N_y is:', N_y

    end subroutine checkNodeDivisionMG



    

end program main