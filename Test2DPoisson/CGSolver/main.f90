

program main
    use iso_fortran_env, only: int32, real64
    use mod_MGSolver
    use mod_GSSolver
    use omp_lib
    implicit none

    real(real64), parameter :: e_const = 1.602176634d-19, eps_0 = 8.8541878188d-12
    integer(int32) :: N_x = 1001, N_y = 10001, numThreads = 6
    type(MGSolver) :: MG_Solver
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate, numberPreSmoothOper, numberPostSmoothOper, numberIter
    integer :: upperBound, lowerBound, rightBound, leftBound
    integer, allocatable :: boundaryConditions(:)
    real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY
    real(real64) :: alpha, beta, R2_future, R2_init, resProduct_old, resProduct_new, solutionRes, relTol, stepTol
    real(real64), allocatable :: resFuture(:), D_vector(:), work(:), solutionCG(:), residualCG(:)

    numberStages = 8
    ! More skewed delX and delY, more smoothing operations needed
    numberPreSmoothOper = 10
    numberPostSmoothOper = 5
    numberIter = 200
    omega = 1.0d0
    relTol = 1.d-8
    stepTol = 1.d-6
    rho = e_const * 1d15
    NESW_wallBoundaries(1) = 1 ! North
    NESW_wallBoundaries(2) = 1 ! East
    NESW_wallBoundaries(3) = 2 ! South
    NESW_wallBoundaries(4) = 2 ! West

    NESW_phiValues(1) = 0.0d0
    NESW_phiValues(2) = 1000.0d0
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


    MG_Solver = MGSolver(N_x, N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    call MG_Solver%makeSmootherStages(delX, delY, NESW_wallBoundaries, boundaryConditions, omega)
    
    ! Set phi values finer grid
    if (upperBound == 1) then
        MG_Solver%GS_smoothers(1)%solution(matDimension-N_x+2:matDimension-1) = upperPhi
    end if
    if (rightBound == 1) then
        MG_Solver%GS_smoothers(1)%solution(2*N_x:matDimension-N_x:N_x) = rightPhi
    end if
    if (lowerBound == 1) then 
        MG_Solver%GS_smoothers(1)%solution(2:N_x-1) = lowerPhi
    end if
    if (leftBound == 1) then 
        MG_Solver%GS_smoothers(1)%solution(N_x+1:matDimension-2*N_x + 1:N_x) = leftPhi
    end if
    if (boundaryConditions(1) == 1) then
        if (lowerBound == leftBound) then
            MG_Solver%GS_smoothers(1)%solution(1) = MIN(leftPhi, lowerPhi)
        else if (lowerBound == 1) then
            MG_Solver%GS_smoothers(1)%solution(1) = lowerPhi
        else 
            MG_Solver%GS_smoothers(1)%solution(1) = leftPhi
        end if
    end if
    if (boundaryConditions(N_x) == 1) then
        if (lowerBound == rightBound) then
            MG_Solver%GS_smoothers(1)%solution(N_x) = MIN(rightPhi, lowerPhi)
        else if (lowerBound == 1) then
            MG_Solver%GS_smoothers(1)%solution(N_x) = lowerPhi
        else
            MG_Solver%GS_smoothers(1)%solution(N_x) = rightPhi
        end if
    end if
    if (boundaryConditions(matDimension-N_x+1) == 1) then
        if (upperBound == leftBound) then
            MG_Solver%GS_smoothers(1)%solution(matDimension-N_x+1) = MIN(leftPhi, upperPhi)
        else if (upperBound == 1) then
            MG_Solver%GS_smoothers(1)%solution(matDimension-N_x+1) = upperPhi
        else
            MG_Solver%GS_smoothers(1)%solution(matDimension-N_x+1) = leftPhi
        end if
    end if
    if (boundaryConditions(matDimension) == 1) then
        if (upperBound == rightBound) then
            MG_Solver%GS_smoothers(1)%solution(matDimension) = MIN(rightPhi, upperPhi)
        else if (upperBound == 1) then
            MG_Solver%GS_smoothers(1)%solution(matDimension) = upperPhi
        else
            MG_Solver%GS_smoothers(1)%solution(matDimension) = rightPhi
        end if
    end if

    ! set source Term first stage
    !$OMP parallel
    !$OMP do
    do i = 1, matDimension
        if (boundaryConditions(i) /= 1) then
            MG_Solver%GS_smoothers(1)%sourceTerm(i) = -rho/eps_0
        end if
    end do
    !$OMP end do
    !$OMP end parallel

    allocate(resFuture(matDimension), D_vector(matDimension), work(matDimension), solutionCG(matDimension), residualCG(matDimension))
    i = 0
    !$OMP parallel 
    !$OMP workshare
    solutionCG = MG_Solver%GS_smoothers(1)%solution
    !$OMP end workshare
    !$OMP workshare
    work = 0.0d0
    !$OMP end workshare
    !$OMP end parallel
    call system_clock(count_rate = timingRate)
    call system_clock(startTime)

    !call MG_Solver%MG_Cycle_Solve(stepTol, relTol, 1)
   
    ! ----------- MGCG ---------------------------------
    call MG_Solver%GS_smoothers(1)%calcResidual()
    !$OMP parallel workshare
    residualCG = MG_Solver%GS_smoothers(1)%residual
    !$OMP end parallel workshare
    
    !call MG_Solver%GS_smoothers(1)%smoothIterations(10, .false.)
    call MG_Solver%GS_smoothers(1)%smoothIterations(MG_Solver%numberPreSmoothOper, .false.) 
    call MG_Solver%GS_smoothers(1)%calcResidual()
    call MG_Solver%F_Cycle()
    solutionRes = MG_Solver%GS_smoothers(1)%smoothIterationsWithRes(MG_Solver%numberPostSmoothOper) 
    !$OMP parallel
    !$OMP workshare
    work = MG_Solver%GS_smoothers(1)%solution - solutionCG
    !$OMP end workshare
    !$OMP workshare
    D_vector = work
    !$OMP end workshare
    !$OMP end parallel
    


    !$OMP parallel
    !$OMP workshare
    R2_init = SUM(residualCG**2)
    !$OMP end workshare
    !$OMP workshare
    resProduct_old = SUM(residualCG * work)
    !$OMP end workshare
    !$OMP end parallel
    do i = 1, numberIter
        call MG_Solver%GS_smoothers(1)%matMult(D_vector, work)
        !$OMP parallel workshare
        alpha = SUM(D_vector * work)
        !$OMP end parallel workshare

        alpha = resProduct_old/alpha

        !$OMP parallel
        !$OMP workshare
        solutionCG = solutionCG + alpha * D_vector
        !$OMP end workshare
        !$OMP workshare
        MG_Solver%GS_smoothers(1)%solution = solutionCG
        !$OMP end workshare
        !$OMP end parallel


        call MG_Solver%GS_smoothers(1)%calcResidual
        !$OMP parallel workshare
        residualCG = MG_Solver%GS_smoothers(1)%residual
        !$OMP end parallel workshare

        !$OMP parallel workshare
        R2_future = SUM(residualCG**2)
        !$OMP end parallel workshare

        
        if (R2_future/R2_init < relTol) then
            print *, 'Exit relative error'
            exit
        end if

        !call MG_Solver%GS_smoothers(1)%smoothIterations(10, .false.)
        call MG_Solver%GS_smoothers(1)%smoothIterations(MG_Solver%numberPreSmoothOper, .false.) 
        call MG_Solver%GS_smoothers(1)%calcResidual()
        call MG_Solver%F_Cycle()
        solutionRes = MG_Solver%GS_smoothers(1)%smoothIterationsWithRes(MG_Solver%numberPostSmoothOper) 
        if (solutionRes < stepTol) then
            print *, 'Exit step error'
            exit
        end if
        !$OMP parallel workshare
        ! work contains preconditioned residual
        work = MG_Solver%GS_smoothers(1)%solution - solutionCG
        !$OMP end parallel workshare

        !$OMP parallel workshare
        resProduct_new = SUM(residualCG * work)
        !$OMP end parallel workshare
        beta = resProduct_new/resProduct_old

        !$OMP parallel workshare
        D_vector = work + beta * D_vector
        !$OMP end parallel workshare

        resProduct_old = resProduct_new


    end do


    !-------------------- mixed MG and CG ----------------------------
    ! call MG_Solver%GS_smoothers(1)%calcResidual()
    ! !$OMP parallel
    ! !$OMP workshare
    ! R2_init = SUM(MG_Solver%GS_smoothers(1)%residual**2)
    ! !$OMP end workshare
    ! !$OMP end parallel
    ! call MG_Solver%GS_smoothers(1)%smoothIterations(MG_Solver%numberPreSmoothOper, .false.) 
    ! call MG_Solver%GS_smoothers(1)%calcResidual()
    ! call MG_Solver%F_Cycle()
    ! solutionRes = MG_Solver%GS_smoothers(1)%smoothIterationsWithRes(MG_Solver%numberPostSmoothOper) 
    ! call MG_Solver%GS_smoothers(1)%calcResidual()
    ! !$OMP parallel
    ! !$OMP workshare
    ! D_vector = MG_Solver%GS_smoothers(1)%residual
    ! !$OMP end workshare
    ! !$OMP workshare
    ! resProduct_old = SUM(MG_Solver%GS_smoothers(1)%residual**2)
    ! !$OMP end workshare
    ! !$OMP end parallel
    
    ! do i = 1, numberIter
    !     call MG_Solver%GS_smoothers(1)%matMult(D_vector, work)
    !     !$OMP parallel workshare
    !     alpha = SUM(D_vector * work)
    !     !$OMP end parallel workshare

    !     alpha = resProduct_old/alpha

    !     !$OMP parallel
    !     !$OMP workshare
    !     MG_Solver%GS_smoothers(1)%solution = MG_Solver%GS_smoothers(1)%solution + alpha * D_vector
    !     !$OMP end workshare
    !     !$OMP end parallel
    !     call MG_Solver%GS_smoothers(1)%smoothIterations(MG_Solver%numberPreSmoothOper, .false.) 
    !     call MG_Solver%GS_smoothers(1)%calcResidual()
    !     call MG_Solver%F_Cycle()
    !     solutionRes = MG_Solver%GS_smoothers(1)%smoothIterationsWithRes(MG_Solver%numberPostSmoothOper) 
    !     if (solutionRes < stepTol) then
    !         print *, 'Exit step error'
    !         exit
    !     end if
    !     call MG_Solver%GS_smoothers(1)%calcResidual

    !     !$OMP parallel workshare
    !     R2_future = SUM(MG_Solver%GS_smoothers(1)%residual**2)
    !     !$OMP end parallel workshare

        
    !     if (R2_future/R2_init < relTol) then
    !         print *, 'Exit relative error'
    !         exit
    !     end if

    !     beta = R2_future/resProduct_old

    !     !$OMP parallel workshare
    !     D_vector = MG_Solver%GS_smoothers(1)%residual + beta * D_vector
    !     !$OMP end parallel workshare

    !     resProduct_old = R2_future


    ! end do



    call system_clock(endTime)
    print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'
    print *, 'Took', MG_Solver%numIter, 'MG iterations'
    print *, 'Took', i-1, 'CG iterations'
    print *, 'solutionRes:', solutionRes

    ! open(41,file='finalSol.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    ! write(41) MG_Solver%GS_smoothers(1)%solution
    ! close(41)

    ! open(41,file='finalRes.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    ! write(41) MG_Solver%GS_smoothers(1)%residual
    ! close(41)




contains
    subroutine checkNodeDivisionMG(N_x, N_y, numberStages)
        ! Check to make sure N_x and N_y is divisible by however many stages in multigrid we want
        integer, intent(in out) :: N_x, N_y, numberStages
        integer :: N_min, temp_min
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