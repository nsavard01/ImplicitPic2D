

program main
    use iso_fortran_env, only: int32, real64
    use mod_GSSolver
    use mod_pardisoSolver
    use mod_CSRMAtrix
    use omp_lib
    implicit none

    real(real64), parameter :: e_const = 1.602176634d-19, eps_0 = 8.8541878188d-12, pi = 4.0d0*atan(1.0d0)
    integer(int32) :: N_x = 101, N_y = 101, numThreads = 6
    type(GSSolver) :: solver
    type(pardisoSolver) :: directSolver
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate, numberPreSmoothOper, numberPostSmoothOper, numberIter
    integer :: upperBound, lowerBound, rightBound, leftBound
    integer, allocatable :: boundaryConditions(:, :)
    real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY
    real(real64) :: alpha, beta, R2_future, R2_init, resProduct_old, resProduct_new, solutionRes, relTol, stepTol
    real(real64), allocatable :: diffX(:), diffY(:), test(:,:)
    logical :: makeX

    
    
    numberStages = 2
    ! More skewed delX and delY, more smoothing operations needed
    numberPreSmoothOper = 10
    numberPostSmoothOper = 10
    numberIter = 200
    omega = 1.5d0
    relTol = 1.d-8
    stepTol = 1.d-6
    rho = e_const * 1d15
    NESW_wallBoundaries(1) = 1 ! North
    NESW_wallBoundaries(2) = 3 ! East
    NESW_wallBoundaries(3) = 2 ! South
    NESW_wallBoundaries(4) = 3 ! West

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
    allocate(diffX(N_x-1), diffY(N_y-1), test((N_x-1)/2 +1, (N_y-1)/2 + 1))
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

    solver = GSSolver(omega)
    call solver%constructPoissonOrthogonal(N_x, N_y, diffX, diffY, NESW_wallBoundaries, boundaryConditions)

    ! CG_Solver = PreCondCGSolver(N_x, N_y, numberStages, numberIter, numberPreSmoothOper, numberPostSmoothOper)
    ! call CG_Solver%makeSmootherStages(delX, delY, NESW_wallBoundaries, boundaryConditions, omega)
    
    ! Set phi values finer grid
    if (upperBound == 1) then
        solver%solution(2:N_x-1, N_y) = upperPhi
    end if
    if (rightBound == 1) then
        solver%solution(N_x, 2:N_y-1) = rightPhi
    end if
    if (lowerBound == 1) then 
        solver%solution(2:N_x-1, 1) = lowerPhi
    end if
    if (leftBound == 1) then 
        solver%solution(1, 2:N_y-1) = leftPhi
    end if
    if (boundaryConditions(1,1) == 1) then
        if (lowerBound == leftBound) then
            solver%solution(1,1) = MIN(leftPhi, lowerPhi)
        else if (lowerBound == 1) then
            solver%solution(1,1) = lowerPhi
        else 
            solver%solution(1,1) = leftPhi
        end if
    end if
    if (boundaryConditions(N_x,1) == 1) then
        if (lowerBound == rightBound) then
            solver%solution(N_x,1) = MIN(rightPhi, lowerPhi)
        else if (lowerBound == 1) then
            solver%solution(N_x,1) = lowerPhi
        else
            solver%solution(N_x,1) = rightPhi
        end if
    end if
    if (boundaryConditions(1, N_y) == 1) then
        if (upperBound == leftBound) then
            solver%solution(1, N_y) = MIN(leftPhi, upperPhi)
        else if (upperBound == 1) then
            solver%solution(1, N_y) = upperPhi
        else
            solver%solution(1, N_y) = leftPhi
        end if
    end if
    if (boundaryConditions(N_x, N_y) == 1) then
        if (upperBound == rightBound) then
            solver%solution(N_x, N_y) = MIN(rightPhi, upperPhi)
        else if (upperBound == 1) then
            solver%solution(N_x, N_y) = upperPhi
        else
            solver%solution(N_x, N_y) = rightPhi
        end if
    end if
    
    directSolver = pardisoSolver(N_x * N_y)
    call buildCurvGridOrthogonalPoissonMatrix(N_x, N_y, diffX, diffY &
            , NESW_wallBoundaries, NESW_phiValues, directSolver%MatValues, directSolver%rowIndex, directSolver%columnIndex, directSolver%sourceTerm)
    call directSolver%initializePardiso(1, 11, 1, 0) ! use default pardiso initialization values
    ! ! set source Term first stage
    !$OMP parallel
    !$OMP do collapse(2)
    do j = 1, N_y
        do i = 1, N_x
            if (boundaryConditions(i,j) /= 1) then
                solver%sourceTerm(i,j) = -rho/eps_0
                k = (j-1) * N_x + i
                directSolver%sourceTerm(k) = -rho/eps_0
            end if
        end do
    end do
    !$OMP end do
    !$OMP end parallel
    
    test = 0.0d0

    ! !$OMP parallel workshare
    ! CG_Solver%solution = CG_Solver%GS_smoothers(1)%solution
    ! !$OMP end parallel workshare
    

    call system_clock(count_rate = timingRate)
    call system_clock(startTime)
    call solver%solveGS(stepTol)
    call system_clock(endTime)
    

    
    print *, 'took', solver%iterNumber
    print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds'

  
    call restriction(solver, test)
    open(41,file='test.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) test
    close(41)
    call prolongation(solver, test)

    open(41,file='finalSol.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    write(41) solver%solution
    close(41)
    
    
    ! call solver%calcResidual()
    ! open(41,file='finalRes.dat', form='UNFORMATTED', access = 'stream', status = 'new')
    ! write(41) solver%solution
    ! close(41)
    ! print *, 'Took', CG_Solver%numIter, 'CG iterations'
    ! print *, 'Took', i-1, 'main iterations'





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

    subroutine restriction(solver, X)
        ! Check to make sure N_x and N_y is divisible by however many stages in multigrid we want
        ! weighted binomial filter type scheme
        type(GSSolver), intent(in) :: solver
        real(real64), intent(in out) :: X(:,:)
        real(real64) :: a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S
        integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, blackIdx
        !$OMP parallel private(blackIdx, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, &
        !$OMP& a_N, a_E, a_W, a_S, C_N, C_E, C_W, C_S)
        !$OMP do collapse(2)
        do j_fine = 3, solver%N_y-2, 2
            do i_fine = 3, solver%N_x-2, 2
                ! calculate fine indices around overlapping index
                i_coarse = (i_fine + 1)/2
                j_coarse = (j_fine + 1)/2
                N_indx = j_fine+1
                E_indx = i_fine+1
                S_indx = j_fine-1
                W_indx = i_fine-1

                ! ! interpolation from north point
                ! C_N = solver%matCoeffs(2, i_fine, N_indx)
                ! C_S = solver%matCoeffs(4, i_fine, N_indx)
                ! a_N = C_S/(C_N + C_S)

                ! ! interpolation from south point
                ! C_N = solver%matCoeffs(2, i_fine, S_indx)
                ! C_S = solver%matCoeffs(4, i_fine, S_indx)
                ! a_S = C_N/(C_N + C_S)

                ! ! interpolation from east point
                ! C_E = solver%matCoeffs(3, E_indx, j_fine)
                ! C_W = solver%matCoeffs(5, E_indx, j_fine)
                ! a_E = C_W/(C_E + C_W)

                ! ! interpolation from west point
                ! C_E = solver%matCoeffs(3, W_indx, j_fine)
                ! C_W = solver%matCoeffs(5, W_indx, j_fine)
                ! a_W = C_E/(C_E + C_W)

                C_N = solver%matCoeffs(2, i_fine, j_fine)
                C_E = solver%matCoeffs(3, i_fine, j_fine)
                C_S = solver%matCoeffs(4, i_fine, j_fine)
                C_W = solver%matCoeffs(5, i_fine, j_fine)
                a_N = C_N/(C_N + C_S)
                a_S = C_S/(C_N + C_S)
                a_E = C_E/(C_E + C_W)
                a_W = C_W/(C_E + C_W)


                ! Interpolate residual to coarse grid
                X(i_coarse, j_coarse) = 0.25d0 * (solver%solution(i_fine, j_fine) + &
                solver%solution(E_indx, j_fine) * a_E + solver%solution(W_indx, j_fine) * a_W +  &
                solver%solution(i_fine, N_indx) * a_N + solver%solution(i_fine, S_indx) * a_S + &
                solver%solution(E_indx, N_indx) * a_E * a_N + solver%solution(E_indx, S_indx) * a_E * a_S + &
                solver%solution(W_indx, N_indx) * a_W * a_N + solver%solution(W_indx, S_indx) * a_W * a_S)
            end do
        end do
        !$OMP end do nowait
        ! each black bound node is a coarse node
        !$OMP do
        do k = 1, solver%numberBlackBoundNodes
            i_fine = solver%black_NESW_BoundIndx(1, k)
            j_fine = solver%black_NESW_BoundIndx(2, k)
            E_indx = solver%black_NESW_BoundIndx(3, k)
            W_indx = solver%black_NESW_BoundIndx(4, k)
            N_indx = solver%black_NESW_BoundIndx(5, k)
            S_indx = solver%black_NESW_BoundIndx(6, k)
            
            ! calculate fine indices around overlapping index
            i_coarse = (i_fine + 1)/2
            j_coarse = (j_fine + 1)/2
            ! if (N_indx == S_indx) then
            !     ! neumann boundary
            !     if (j_fine == 1) then
            !         ! bottom boundary
            !         C_N = solver%matCoeffs(2, i_fine, N_indx)
            !         C_S = solver%matCoeffs(4, i_fine, N_indx)
            !         a_N = C_S/(C_N + C_S)
            !         a_S = a_N
            !     else
            !         ! top boundary
            !         C_N = solver%matCoeffs(2, i_fine, S_indx)
            !         C_S = solver%matCoeffs(4, i_fine, S_indx)
            !         a_S = C_N/(C_N + C_S)
            !         a_N = a_S
            !     end if
            ! else
            !     ! interpolation from north point
            !     C_N = solver%matCoeffs(2, i_fine, N_indx)
            !     C_S = solver%matCoeffs(4, i_fine, N_indx)
            !     a_N = C_S/(C_N + C_S)

            !     ! interpolation from south point
            !     C_N = solver%matCoeffs(2, i_fine, S_indx)
            !     C_S = solver%matCoeffs(4, i_fine, S_indx)
            !     a_S = C_N/(C_N + C_S)
            ! end if
            
            ! if (E_indx == W_indx) then
            !     ! neumann boundary
            !     if (i_fine == 1) then
            !         ! left boundary
            !         C_E = solver%matCoeffs(3, E_indx, j_fine)
            !         C_W = solver%matCoeffs(5, E_indx, j_fine)
            !         a_E = C_W/(C_E + C_W)
            !         a_W = a_E
            !     else
            !         C_E = solver%matCoeffs(3, W_indx, j_fine)
            !         C_W = solver%matCoeffs(5, W_indx, j_fine)
            !         a_W = C_E/(C_E + C_W)
            !         a_E = a_W
            !     end if
            ! else
            !     ! interpolation from east point
            !     C_E = solver%matCoeffs(3, E_indx, j_fine)
            !     C_W = solver%matCoeffs(5, E_indx, j_fine)
            !     a_E = C_W/(C_E + C_W)

            !     ! interpolation from west point
            !     C_E = solver%matCoeffs(3, W_indx, j_fine)
            !     C_W = solver%matCoeffs(5, W_indx, j_fine)
            !     a_W = C_E/(C_E + C_W)
            ! end if

            C_N = solver%matCoeffs(2, i_fine, j_fine)
            C_E = solver%matCoeffs(3, i_fine, j_fine)
            C_S = solver%matCoeffs(4, i_fine, j_fine)
            C_W = solver%matCoeffs(5, i_fine, j_fine)
            a_N = C_N/(C_N + C_S)
            a_S = C_S/(C_N + C_S)
            a_E = C_E/(C_E + C_W)
            a_W = C_W/(C_E + C_W)

            ! Interpolate residual to coarse grid
            X(i_coarse, j_coarse) = 0.25d0 * (solver%solution(i_fine, j_fine) + &
            solver%solution(E_indx, j_fine) * a_E + solver%solution(W_indx, j_fine) * a_W +  &
            solver%solution(i_fine, N_indx) * a_N + solver%solution(i_fine, S_indx) * a_S + &
            solver%solution(E_indx, N_indx) * a_E * a_N + solver%solution(E_indx, S_indx) * a_E * a_S + &
            solver%solution(W_indx, N_indx) * a_W * a_N + solver%solution(W_indx, S_indx) * a_W * a_S)
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine restriction

    subroutine prolongation(solver, X)
        ! Prolongate operator
        type(GSSolver), intent(in out) :: solver
        real(real64), intent(in) :: X(:,:)
        real(real64) :: w_1, w_2, w_3, w_4, w_tot
        integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse
        N_x_coarse = (solver%N_x+1)/2
        N_y_coarse = (solver%N_y+1)/2
        ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
        !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine, w_1, w_2, w_3, w_4, w_tot)
        ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
        !$OMP do collapse(2)
        ! set similar nodes
        do j_coarse = 1, N_y_coarse
            do i_coarse = 1, N_x_coarse
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2-1
                ! Add overlapping coarse nodes
                solver%solution(i_fine, j_fine) = X(i_coarse, j_coarse)
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! go horizontal each row,
        do j_coarse = 1, N_y_coarse
            do i_coarse = 1, N_x_coarse-1
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2 - 1

                w_1 = solver%matCoeffs(3, i_fine+1, j_fine) !E
                w_2 = solver%matCoeffs(5, i_fine+1, j_fine) !W
                w_tot = w_1 + w_2
                w_1 = w_1 / w_tot
                w_2 = w_2 / w_tot
                ! Average fine nodes between coarse nodes horizontally
                solver%solution(i_fine+1, j_fine) = X(i_coarse+1, j_coarse) * w_1 + X(i_coarse, j_coarse) * w_2
                ! self%GS_smoothers(stageInt)%solution(O_indx+1) = self%GS_smoothers(stageInt)%solution(O_indx+1) + &
                !     weightX * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+1))
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! go vertical
        do j_coarse = 1, N_y_coarse-1
            do i_coarse = 1, N_x_coarse
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2 - 1
                w_1 = solver%matCoeffs(2, i_fine, j_fine+1) !N
                w_2 = solver%matCoeffs(4, i_fine, j_fine+1) !S
                w_tot = w_1 + w_2
                w_1 = w_1 / w_tot
                w_2 = w_2 / w_tot
                ! Average fine nodes between coarse nodes vertical direction
                solver%solution(i_fine, j_fine+1) = X(i_coarse, j_coarse+1) * w_1 + X(i_coarse, j_coarse) * w_2
                ! self%GS_smoothers(stageInt)%solution(O_indx+N_x_fine) = self%GS_smoothers(stageInt)%solution(O_indx+N_x_fine) + &
                !     weightY * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+N_x_coarse))
            end do
        end do
        !$OMP end do nowait
        !$OMP do collapse(2)
        ! add offset fine nodes which require 4 point interpolation
        ! first loop through SW corner coarse nodes
        do j_coarse = 1, N_y_coarse-1
            do i_coarse = 1, N_x_coarse-1
                i_fine = i_coarse*2 - 1
                j_fine = j_coarse*2 - 1
                w_1 = solver%matCoeffs(2, i_fine+1, j_fine+1) !N
                w_2 = solver%matCoeffs(4, i_fine+1, j_fine+1) !S
                w_tot = w_1 + w_2
                w_1 = w_1/w_tot
                w_2 = w_2/w_tot

                w_3 = solver%matCoeffs(3, i_fine+1, j_fine+1) !E
                w_4 = solver%matCoeffs(5, i_fine+1, j_fine+1) !W
                w_tot = w_3 + w_4
                w_3 = w_3/w_tot
                w_4 = w_4/w_tot

                solver%solution(i_fine+1, j_fine+1) = X(i_coarse, j_coarse) * w_2*w_4 + X(i_coarse+1, j_coarse) * w_3 * w_2 &
                    + X(i_coarse, j_coarse+1) * w_1 * w_4 + X(i_coarse+1, j_coarse+1) * w_1 * w_3
                ! self%GS_smoothers(stageInt)%solution(O_indx) =  self%GS_smoothers(stageInt)%solution(O_indx) + &
                !     weightXY * (lowerSolution(O_indx_coarse) + lowerSolution(O_indx_coarse+1) &
                !         + lowerSolution(O_indx_coarse+N_x_coarse) + lowerSolution(O_indx_coarse+N_x_coarse+1))
            end do
        end do
        !$OMP end do
        !$OMP end parallel

    end subroutine prolongation



    

end program main