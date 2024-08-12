module mod_CSRMAtrix
    use iso_fortran_env, only: int32, real64
    use omp_lib
    implicit none

    ! Module to build CSR matrix, used for direct solver

contains


    subroutine buildEvenGridOrthogonalPoissonMatrix(N_x, N_y, delX, delY, NESW_wallBoundaries, NESW_phiValues, MatValues, rowIndex, columnIndex, sourceTerm)
        ! Construct CSR matrix arrays used for Even grid orthogonal Gauss' law solver in phi
        ! For structured grid with fixed N_x, N_y
        integer, intent(in) :: N_x, N_y ! Node numbers and physical cell size in x and y
        integer, intent(in) :: NESW_wallBoundaries(4) ! NESW wall boundary values
        real(real64), intent(in) :: NESW_phiValues(4) ! Phi values for each wall
        real(real64), intent(in out) :: sourceTerm(N_x*N_y)
        integer(int32), intent(in out) :: rowIndex(N_x*N_y + 1) 
        real(real64), allocatable, intent(in out) :: MatValues(:) ! allocatable array for matrix values
        integer(int32), allocatable, intent(in out) :: columnIndex(:)
        integer(int32) :: matDimension, numberMatValues, i, j, k, matValIdx, boundaryConditions(N_x*N_y), upperBound, lowerBound, rightBound, leftBound, errorInt
        real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi, delX, delY
        matDimension = N_x * N_y

        upperBound = NESW_wallBoundaries(1)
        rightBound = NESW_wallBoundaries(2)
        lowerBound = NESW_wallBoundaries(3)
        leftBound = NESW_wallBoundaries(4)

        upperPhi = NESW_phiValues(1)
        rightPhi = NESW_phiValues(2)
        lowerphi = NESW_phiValues(3)
        leftPhi = NESW_phiValues(4)

        numberMatValues = (N_x-2) * (N_y-2) * 5 ! inner nodes have five interdependencies
        ! For every boundary determine how many row values for matrix
        !$OMP parallel
        !$OMP workshare
        boundaryConditions = 0
        !$OMP end workshare
        !$OMP workshare
        sourceTerm = 0.0d0
        !$OMP end workshare
        !$OMP end parallel
        SELECT CASE (upperBound)
        CASE(1)
            numberMatValues = numberMatValues + (N_x-2)
            sourceTerm(matDimension-N_x+2:matDimension-1) = upperPhi
        CASE(2)
            numberMatValues = numberMatValues + (N_x-2) * 4
        CASE(3)
            numberMatValues = numberMatValues + (N_x-2) * 5
        END SELECT
        boundaryConditions(matDimension-N_x+2:matDimension-1) = upperBound

        SELECT CASE (rightBound)
        CASE(1)
            numberMatValues = numberMatValues + (N_y-2)
            sourceTerm(2*N_x:matDimension-N_x:N_x) = rightPhi
        CASE(2)
            numberMatValues = numberMatValues + (N_y-2) * 4
        CASE(3)
            numberMatValues = numberMatValues + (N_y-2) * 5
        END SELECT
        boundaryConditions(2*N_x:matDimension-N_x:N_x) = rightBound

        SELECT CASE (lowerBound)
        CASE(1)
            numberMatValues = numberMatValues + (N_x-2)
            sourceTerm(2:N_x-1) = lowerPhi
        case(2)
            numberMatValues = numberMatValues + (N_x-2) * 4
        case(3)
            numberMatValues = numberMatValues + (N_x-2) * 5
        end select
        boundaryConditions(2:N_x-1) = lowerBound

        select case (leftBound)
        case(1)
            numberMatValues = numberMatValues + (N_y-2)
            sourceTerm(N_x+1:matDimension-2*N_x + 1:N_x) = leftPhi
        case(2)
            numberMatValues = numberMatValues + (N_y-2) * 4
        case(3)
            numberMatValues = numberMatValues + (N_y-2) * 5
        end select
        boundaryConditions(N_x+1:matDimension-2*N_x + 1:N_x) = leftBound

        ! Upper right corner
        select case(MIN(upperBound, rightBound))
        case(1)
            numberMatValues = numberMatValues + 1
            if (upperBound == rightBound) then
                sourceTerm(matDimension) = MIN(upperPhi, rightPhi)
            else if (upperBound == 1) then
                sourceTerm(matDimension) = upperPhi
            else
                sourceTerm(matDimension) = rightPhi
            end if
        case(2)
            numberMatValues = numberMatValues + 3
        end select
        boundaryConditions(matDimension) = MIN(upperBound, rightBound)

        ! lower right corner
        select case(MIN(lowerBound, rightBound))
        case(1)
            numberMatValues = numberMatValues + 1
            if (lowerBound == rightBound) then
                sourceTerm(N_x) = MIN(lowerPhi, rightPhi)
            else if (lowerBound == 1) then
                sourceTerm(N_x) = lowerPhi
            else
                sourceTerm(N_x) = rightPhi
            end if
        case(2)
            numberMatValues = numberMatValues + 3
        end select
        boundaryConditions(N_x) = MIN(lowerBound, rightBound)

        ! lower left corner
        select case(MIN(lowerBound, leftBound))
        case(1)
            numberMatValues = numberMatValues + 1
            if (lowerBound == leftBound) then
                sourceTerm(1) = MIN(lowerPhi, leftPhi)
            else if (lowerBound == 1) then
                sourceTerm(1) = lowerPhi
            else
                sourceTerm(1) = leftPhi
            end if
        case(2)
            numberMatValues = numberMatValues + 3
        end select
        boundaryConditions(1) = MIN(lowerBound, leftBound)

        ! upper left corner
        select case(MIN(upperBound, leftBound))
        case(1)
            numberMatValues = numberMatValues + 1
            if (upperBound == leftBound) then
                sourceTerm(matDimension-N_x+1) = MIN(upperPhi, leftPhi)
            else if (upperBound == 1) then
                sourceTerm(matDimension-N_x+1) = upperPhi
            else
                sourceTerm(matDimension-N_x+1) = leftPhi
            end if
        case(2)
            numberMatValues = numberMatValues + 3
        end select
        boundaryConditions(matDimension-N_x+1) = MIN(upperBound, leftBound)

        allocate(MatValues(numberMatValues), columnIndex(numberMatValues))
        errorInt = 0
        matValIdx = 0
        do j = 1, N_y
            do i = 1, N_x
                k = (j-1) * N_x + i
                rowIndex(k) = matValIdx+1
                if (boundaryConditions(k) == 0) then
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
                    select case(boundaryConditions(k))
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

end module mod_CSRMAtrix