module mod_CSRMAtrix
    use iso_fortran_env, only: int32, real64
    use omp_lib
    implicit none

    ! Module to build CSR matrix, used for direct solver

contains


    subroutine buildEvenGridOrthogonalPoissonMatrix(N_x, N_y, delX, delY, boundaryConditions, MatValues, rowIndex, columnIndex, sourceTerm)
        ! Construct CSR matrix arrays used for Even grid orthogonal Gauss' law solver in phi
        ! For structured grid with fixed N_x, N_y
        integer, intent(in) :: N_x, N_y ! Node numbers and physical cell size in x and y
        integer, intent(in) :: boundaryConditions(N_x,N_y) ! NESW wall boundary value ! Phi values for each wall
        real(real64), intent(in) :: delX, delY
        real(real64), intent(in out) :: sourceTerm(N_x*N_y)
        integer(int32), intent(in out) :: rowIndex(N_x*N_y + 1) 
        real(real64), allocatable, intent(in out) :: MatValues(:) ! allocatable array for matrix values
        integer(int32), allocatable, intent(in out) :: columnIndex(:)
        integer(int32) :: matDimension, numberMatValues, i, j, k, matValIdx, upperBound, lowerBound, rightBound, leftBound, errorInt
        real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi
        
        matDimension = N_x * N_y

        

        numberMatValues = 0 ! inner nodes have five interdependencies
        ! For every boundary determine how many row values for matrix
        !$OMP parallel
        !$OMP workshare
        sourceTerm = 0.0d0
        !$OMP end workshare
        !$OMP end parallel

        !$OMP parallel private(j,i) reduction(+:numberMatValues)
        !$OMP do collapse(2)
        do j = 1, N_y
            do i = 1, N_x
                select case (boundaryConditions(i,j))
                case (0,3)
                    numberMatValues = numberMatValues + 5
                case (1)
                    numberMatValues = numberMatValues + 1
                case (2)
                    if ((i > 1 .and. i < N_x) .or. (j>1 .and. j < N_y))then
                        numberMatValues = numberMatValues + 4
                    else
                        numberMatValues = numberMatValues + 3
                    end if
                end select
            end do
        end do
        !$OMP end do
        !$OMP end parallel


        allocate(MatValues(numberMatValues), columnIndex(numberMatValues))
        errorInt = 0
        matValIdx = 0
        do j = 1, N_y
            do i = 1, N_x
                k = (j-1) * N_x + i
                rowIndex(k) = matValIdx+1
                if (boundaryConditions(i,j) == 0) then
                    ! Inner node
                    ! South
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k - N_x
                    MatValues(matValIdx) = 1.0d0 / delY**2
                    ! West
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k-1
                    MatValues(matValIdx) = 1.0d0 / delX**2
                    ! Origin
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k
                    MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                    ! East
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k+1
                    MatValues(matValIdx) = 1.0d0 / delX**2
                    ! North
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k + N_x
                    MatValues(matValIdx) = 1.0d0 / delY**2
                else
                    select case(boundaryConditions(i,j))
                    case(1)
                        columnIndex(matValIdx+1) = k
                        MatValues(matValIdx+1) = 1.0d0
                        matValIdx = matValIdx+1
                    case(2)
                        ! Neumann boundary conditions
                        if (k <= N_x) then
                            if (k /= 1 .and. k /= N_x) then  
                                ! lower boundary, have W, O, E, N
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = 2.0d0 / delY**2
                            else if (k == 1) then
                                ! lower left corner
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = 2.0d0 / delX**2
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = 2.0d0 / delY**2
                            else
                                ! lower right corner
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = 2.0d0 / delX**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = 2.0d0 / delY**2
                            end if
                        else if (k >= matDimension - N_x + 1) then
                            if (k/=matDimension-N_x + 1 .and. k/=matDimension) then
                                ! Upper Boundary
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = 2.0d0 / delY**2
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                            else if (k == matDimension-N_x+1) then
                                ! upper left boundary
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = 2.0d0 / delY**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = 2.0d0 / delX**2
                            else
                                ! upper right boundary
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = 2.0d0 / delY**2
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = 2.0d0 / delX**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                            end if
                        else if (MOD(k, N_x) == 1) then
                            ! left boundary
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                            ! East
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k+1
                            MatValues(matValIdx) = 2.0d0 / delX**2
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                        else 
                            ! Right boundary
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                            ! West
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k-1
                            MatValues(matValIdx) = 2.0d0 / delX**2
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                        end if
                    case(3)
                        ! Periodic boundary conditions
                        if (k <= N_x) then
                            if (k /= 1 .and. k /= N_x) then  
                                ! lower boundary, have W, O, E, N
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = 1.0d0 / delY**2
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = matDimension - 2*N_x + k
                                MatValues(matValIdx) = 1.0d0 / delY**2
                            else if (k == 1) then
                                ! lower left corner
                                errorInt = 1
                            else
                                ! lower right corner
                                errorInt = 1
                            end if
                        else if (k >= matDimension - N_x + 1) then
                            if (k/=matDimension-N_x + 1 .and. k/=matDimension) then
                                ! Upper Boundary
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-matDimension + 2 * N_x
                                MatValues(matValIdx) = 1.0d0 / delY**2
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = 1.0d0 / delY**2
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = 1.0d0 / delX**2
                            else if (k == matDimension-N_x+1) then
                                ! upper left boundary
                                errorInt = 1
                            else
                                ! upper right boundary
                                errorInt = 1
                            end if
                        else if (MOD(k, N_x) == 1) then
                            ! left boundary
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                            ! East
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k+1
                            MatValues(matValIdx) = 1.0d0 / delX**2
                            ! West
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k+N_x-2
                            MatValues(matValIdx) = 1.0d0 / delX**2
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                        else 
                            ! Right boundary
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                            ! East
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k-N_x+2
                            MatValues(matValIdx) = 1.0d0 / delX**2
                            ! West
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k-1
                            MatValues(matValIdx) = 1.0d0 / delX**2
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = -2.0d0 * (1.0d0 / delX**2 + 1.0d0/ delY**2)
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = 1.0d0 / delY**2
                        end if
                    end select    
                end if
            end do
        end do
        rowIndex(matDimension+1) = numberMatValues+1

    end subroutine buildEvenGridOrthogonalPoissonMatrix

    subroutine buildCurvGridOrthogonalPoissonMatrix(N_x, N_y, diffX, diffY, boundaryConditions, MatValues, rowIndex, columnIndex, sourceTerm)
        ! Construct CSR matrix arrays used for curvilinear grid orthogonal Gauss' law solver in phi
        ! For structured grid with fixed N_x, N_y
        integer, intent(in) :: N_x, N_y ! Node numbers and physical cell size in x and y
        integer, intent(in) :: boundaryConditions(N_x, N_y) ! boundary values node
        real(real64), intent(in) :: diffX(N_x-1), diffY(N_y-1) ! Phi values for each wall
        real(real64), intent(in out) :: sourceTerm(N_x*N_y)
        integer(int32), intent(in out) :: rowIndex(N_x*N_y + 1) 
        real(real64), allocatable, intent(in out) :: MatValues(:) ! allocatable array for matrix values
        integer(int32), allocatable, intent(in out) :: columnIndex(:)
        integer(int32) :: matDimension, numberMatValues, i, j, k, matValIdx, upperBound, lowerBound, rightBound, leftBound, errorInt
        real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi, delX_E, delX_W, delY_N, delY_S, C_N, C_S, C_E, C_W, C_O
        matDimension = N_x * N_y

        numberMatValues = 0 ! inner nodes have five interdependencies
        ! For every boundary determine how many row values for matrix
        !$OMP parallel
        !$OMP workshare
        sourceTerm = 0.0d0
        !$OMP end workshare
        !$OMP end parallel

        !$OMP parallel private(j,i) reduction(+:numberMatValues)
        !$OMP do collapse(2)
        do j = 1, N_y
            do i = 1, N_x
                select case (boundaryConditions(i,j))
                case (0,3)
                    numberMatValues = numberMatValues + 5
                case (1)
                    numberMatValues = numberMatValues + 1
                case (2)
                    if ((i > 1 .and. i < N_x) .or. (j>1 .and. j < N_y))then
                        numberMatValues = numberMatValues + 4
                    else
                        numberMatValues = numberMatValues + 3
                    end if
                end select
            end do
        end do
        !$OMP end do
        !$OMP end parallel


        allocate(MatValues(numberMatValues), columnIndex(numberMatValues))
        errorInt = 0
        matValIdx = 0
        do j = 1, N_y
            do i = 1, N_x
                k = (j-1) * N_x + i
                rowIndex(k) = matValIdx+1
                if (boundaryConditions(i,j) == 0) then
                    delX_E = diffX(i)
                    delX_W = diffX(i-1)
                    delY_N = diffY(j)
                    delY_S = diffY(j-1)
                    ! Inner node, set coefficients
                    C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                    C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                    C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                    C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                    C_O = -(C_N + C_E + C_W + C_S) ! O
                    ! Inner node
                    ! South
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k - N_x
                    MatValues(matValIdx) = C_S
                    ! West
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k-1
                    MatValues(matValIdx) = C_W
                    ! Origin
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k
                    MatValues(matValIdx) = C_O
                    ! East
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k+1
                    MatValues(matValIdx) = C_E
                    ! North
                    matValIdx = matValIdx + 1
                    columnIndex(matValIdx) = k + N_x
                    MatValues(matValIdx) = C_N
                else
                    select case(boundaryConditions(i,j))
                    case(1)
                        columnIndex(matValIdx+1) = k
                        MatValues(matValIdx+1) = 1.0d0
                        matValIdx = matValIdx+1
                    case(2)
                        ! Neumann boundary conditions
                        if (j == 1) then
                            delY_N = diffY(j)
                            C_N = 2.0d0 / (delY_N**2)
                            if (i /= 1 .and. i /= N_x) then 
                                delX_E = diffX(i)
                                delX_W = diffX(i-1)
                                C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                                C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                                C_O = -(C_N + C_E + C_W) ! O
                                ! lower boundary, have W, O, E, N
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = C_W
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = C_E
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = C_N
                            else if (i == 1) then
                                ! lower left corner
                                delX_E = diffX(i)
                                C_E = 2.0d0 / (delX_E**2) ! E
                                C_O = -(C_N + C_E) ! O
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = C_E
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = C_N
                            else
                                ! lower right corner
                                delX_W = diffX(i-1)
                                C_W = 2.0d0 / (delX_W**2) ! W
                                C_O = -(C_N + C_W) ! O
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = C_W
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = C_N
                            end if
                        else if (j == N_y) then
                            delY_S = diffY(j-1)
                            C_S = 2.0d0 / (delY_S**2)
                            if (i /= 1 .and. i /= N_x) then
                                ! Upper Boundary
                                delX_E = diffX(i)
                                delX_W = diffX(i-1)
                                C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                                C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                                C_O = -(C_S + C_E + C_W) ! O
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = C_S
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = C_W
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = C_E
                            else if (i == 1) then
                                ! upper left boundary
                                delX_E = diffX(i)
                                C_E = 2.0d0 / (delX_E**2) ! E
                                C_O = -(C_S + C_E) ! O
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = C_S
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = C_E
                            else
                                ! upper right boundary
                                delX_W = diffX(i-1)
                                C_W = 2.0d0 / (delX_W**2) ! W
                                C_O = -(C_S + C_W) ! O
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = C_S
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = C_W
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                            end if
                        else if (i == 1) then
                            ! left boundary
                            delX_E = diffX(i)
                            delY_N = diffY(j)
                            delY_S = diffY(j-1)
                            ! Inner node, set coefficients
                            C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                            C_E = 2.0d0 / (delX_E**2) ! E
                            C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                            C_O = -(C_N + C_E + C_S) ! O
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = C_S
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = C_O
                            ! East
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k+1
                            MatValues(matValIdx) = C_E
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = C_N
                        else 
                            ! Right boundary
                            delX_W = diffX(i-1)
                            delY_N = diffY(j)
                            delY_S = diffY(j-1)
                            ! Inner node, set coefficients
                            C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                            C_W = 2.0d0 / (delX_W**2)
                            C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                            C_O = -(C_N + C_W + C_S) ! O
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = C_S
                            ! West
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k-1
                            MatValues(matValIdx) = C_W
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = C_O
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = C_N
                        end if
                    case(3)
                        ! Periodic boundary conditions
                        if (j == 1) then
                            delY_N = diffY(j)
                            delY_S = diffY(N_y-1)
                            C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                            C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                            if (i /= 1 .and. i /= N_x) then  
                                ! lower boundary, have W, O, E, N
                                delX_E = diffX(i)
                                delX_W = diffX(i-1)
                                C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                                C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                                C_O = -(C_N + C_E + C_W + C_S) ! O
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = C_W
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = C_E
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k + N_x
                                MatValues(matValIdx) = C_N
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = matDimension - 2*N_x + k
                                MatValues(matValIdx) = C_S
                            else if (i == 1) then
                                ! lower left corner
                                errorInt = 1
                            else
                                ! lower right corner
                                errorInt = 1
                            end if
                        else if (j == N_y) then
                            delY_N = diffY(1)
                            delY_S = diffY(N_y-1)
                            C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                            C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                            if (i/=1 .and. i/=N_x) then
                                ! Upper Boundary
                                delX_E = diffX(i)
                                delX_W = diffX(i-1)
                                C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                                C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                                C_O = -(C_N + C_E + C_W + C_S) ! O
                                ! North
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-matDimension + 2 * N_x
                                MatValues(matValIdx) = C_N
                                ! South
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k - N_x
                                MatValues(matValIdx) = C_S
                                ! West
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k-1
                                MatValues(matValIdx) = C_W
                                ! Origin
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k
                                MatValues(matValIdx) = C_O
                                ! East
                                matValIdx = matValIdx + 1
                                columnIndex(matValIdx) = k+1
                                MatValues(matValIdx) = C_E
                            else if (i == 1) then
                                ! upper left boundary
                                errorInt = 1
                            else
                                ! upper right boundary
                                errorInt = 1
                            end if
                        else if (i == 1) then
                            ! left boundary
                            delY_N = diffY(j)
                            delY_S = diffY(j-1)
                            C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                            C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                            delX_E = diffX(1)
                            delX_W = diffX(N_x-1)
                            C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                            C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                            C_O = -(C_N + C_E + C_W + C_S) ! O
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = C_S
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = C_O
                            ! East
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k+1
                            MatValues(matValIdx) = C_E
                            ! West
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k+N_x-2
                            MatValues(matValIdx) = C_W
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = C_N
                        else 
                            ! Right boundary
                            delY_N = diffY(j)
                            delY_S = diffY(j-1)
                            C_N = 1.0d0 / (delY_N * 0.5d0 * (delY_N + delY_S)) ! N
                            C_S = 1.0d0 / (delY_S * 0.5d0 * (delY_N + delY_S)) ! S
                            delX_E = diffX(1)
                            delX_W = diffX(N_x-1)
                            C_E = 1.0d0 / (delX_E * 0.5d0 * (delX_E + delX_W)) ! E
                            C_W = 1.0d0 / (delX_W * 0.5d0 * (delX_E + delX_W)) ! W
                            C_O = -(C_N + C_E + C_W + C_S) ! O
                            ! South
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k - N_x
                            MatValues(matValIdx) = C_S
                            ! East
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k-N_x+2
                            MatValues(matValIdx) = C_E
                            ! West
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k-1
                            MatValues(matValIdx) = C_W
                            ! Origin
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k
                            MatValues(matValIdx) = C_O
                            ! North
                            matValIdx = matValIdx + 1
                            columnIndex(matValIdx) = k + N_x
                            MatValues(matValIdx) = C_N
                        end if
                    end select    
                end if
            end do
        end do
        rowIndex(matDimension+1) = numberMatValues+1

    end subroutine buildCurvGridOrthogonalPoissonMatrix

end module mod_CSRMAtrix