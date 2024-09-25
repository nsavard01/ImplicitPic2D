module mod_GS_Base_Even
    use iso_fortran_env, only: int32, int64, real64
    use omp_lib
    use mod_GS_Base
    use mod_domain_uniform
    implicit none


    type, extends(GS_Base) :: GS_Base_Even
        ! store grid quantities
        real(real64) :: coeffX, coeffY, centerCoeff
    contains
        procedure, public, pass(self) :: constructPoissonOrthogonal
        procedure, public, pass(self) :: initialize_GS_Even
        ! procedure, public, pass(self) :: restriction => restriction_even
        ! procedure, public, pass(self) :: prolongation => prolongation_even
        procedure, public, pass(self) :: calcResidual => calcResidual_even
        ! procedure, public, pass(self) :: AX_Mult => AX_Mult_even
        ! procedure, public, pass(self) :: YAX_Mult => YAX_Mult_even
    end type

contains

    subroutine initialize_GS_Even(self, world)
        ! Construct derived object for everything needed from GS_Base_Even class
        class(GS_Base_Even), intent(in out) :: self
        type(domain_uniform), intent(in), target :: world

        self%world => world
        self%coeffX = 1.0d0 / (world%del_x**2)
        self%coeffY = 1.0d0 / (world%del_y**2)
        self%centerCoeff = -1.0d0 / (2.0d0 * self%coeffX + 2.0d0 * self%coeffY)

    end subroutine initialize_GS_Even

    subroutine constructPoissonOrthogonal(self)
        ! Construct orthogonal grid solver
        class(GS_Base_Even), intent(in out) :: self
        ! integer(int32), intent(in) :: NESW_wallBoundaries(4), boundaryConditions(self%N_x, self%N_y)
        ! real(real64), intent(in) :: del_x, del_y
    end subroutine constructPoissonOrthogonal

!     subroutine AX_Mult_even(self, x, y)
!         ! Use gauss-seidel to calculate x^T * A * x
!         class(GS_Base_Even), intent(in out) :: self
!         real(real64), intent(in) :: x(self%N_x, self%N_y)
!         real(real64), intent(in out) :: y(self%N_x, self%N_y)
!         integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j
!         real(real64) :: inv_centerCoeff
!         inv_centerCoeff = 1.0d0/self%centerCoeff
!         !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
!         !$OMP&  W_indx)
!         ! loop through inner nodes
!         !$OMP do
!         do k = 1, self%number_inner_rows
!             do p = 1, self%number_row_sections(k)
!                 do i = self%start_inner_indx_x(p), self%end_inner_indx_x(p)
!                     j = self%start_row_indx + k - 1
!                     N_indx = j+1
!                     E_indx = i+1
!                     S_indx = j-1
!                     W_indx = i-1
!                     y(i,j) = self%coeffY * (x(i,N_indx) + x(i,S_indx)) + &
!                         self%coeffX * (x(E_indx, j) + x(W_indx, j)) + x(i,j)*inv_centerCoeff
!                 end do
!             end do
!         end do
!         !$OMP end do nowait
!         ! go through boundary conditions
!         !$OMP do
!         do k = 1, self%number_inner_rows
!             j = self%start_row_indx + k - 1
!             select case (self%start_boundary_x(1))
!             case(1)
!                 continue
!             case(2)
!                 if (j/=1 .and. j/= self%N_y) then
                    
!                 end if
!             case(3)

!             end select
!         end do
!         !$OMP end do nowait
!         !$OMP end parallel 
!     end subroutine AX_Mult_even

!     function YAX_Mult_even(self, x, y) result(res)
!         ! Use gauss-seidel to calculate x^T * A * x
!         class(GS_Base_Even), intent(in out) :: self
!         real(real64), intent(in) :: x(self%N_x, self%N_y), y(self%N_x, self%N_y)
!         integer :: O_indx, N_indx, E_indx, S_indx, W_indx, k, i, j, p
!         real(real64) :: res, inv_centerCoeff
!         inv_centerCoeff = 1.0d0/self%centerCoeff
!         res = 0.0d0
!         !$OMP parallel private(i, j, p, k, N_indx, E_indx, S_indx, &
!         !$OMP&  W_indx) reduction(+:res)
!         ! loop through inner nodes
!         !$OMP do collapse(2)
!         do k = 1, self%numberRows
!             do p = 1, self%numberColumns
!                 i = self%startCol + p - 1
!                 j = self%startRow + k - 1
!                 N_indx = self%vertIndx(1,k)
!                 E_indx = self%horzIndx(1,p)
!                 S_indx = self%vertIndx(2,k)
!                 W_indx = self%horzIndx(2,p)
!                 res = res + y(i,j) * (self%coeffY * (x(i,N_indx) + x(i,S_indx)) + &
!                     self%coeffX * (x(E_indx, j) + x(W_indx, j)) + x(i,j)*inv_centerCoeff)
!             end do
!         end do
!         !$OMP end do
!         !$OMP end parallel 
!     end function YAX_Mult_even

    subroutine calcResidual_even(self)
        ! Solve GS down to some tolerance
        class(GS_Base_Even), intent(in out) :: self
        integer :: N_indx, E_indx, S_indx, W_indx, i, j, k, p, i_finest, j_finest

        ! first do corners (due to inconvenience in vectorizing) which have to be neumann for changes
            ! these are all red points
        if (self%world%boundary_conditions(1,1) == 2) then
            ! lower left corner
            self%residual(1,1) = self%sourceTerm(1,1) - 2.0d0 * self%solution(1, 2) * self%coeffY - &
                2.0d0 * self%solution(2, 1) * self%coeffX  - self%solution(1,1)/self%centerCoeff
        end if
        if (self%world%boundary_conditions(self%world%N_x, 1) == 2) then
            ! lower right corner
            self%residual(self%N_x, 1) = self%sourceTerm(self%N_x, 1) - 2.0d0 * self%solution(self%N_x, 2) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, 1) * self%coeffX -  self%solution(self%N_x, 1)/self%centerCoeff
        end if
        if (self%world%boundary_conditions(1, self%world%N_y) == 2) then
            ! upper left corner
            self%residual(1, self%N_y) = self%sourceTerm(1, self%N_y) - 2.0d0 * self%solution(1, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(2, self%N_y-1) * self%coeffX - self%solution(1, self%N_y)/self%centerCoeff
        end if
        if (self%world%boundary_conditions(self%world%N_x, self%world%N_y) == 2) then
            ! upper right corner
            self%residual(self%N_x, self%N_y) = self%sourceTerm(self%N_x, self%N_y) - 2.0d0 * self%solution(self%N_x, self%N_y-1) * self%coeffY - &
                2.0d0 * self%solution(self%N_x-1, self%N_y) * self%coeffX - self%solution(self%N_x, self%N_y)/self%centerCoeff
        end if

        !$OMP parallel private(k, p, i, j, N_indx, W_indx, E_indx, S_indx, i_finest, j_finest)

        !$OMP do
        ! Sweep inner nodes
        do k = 1, self%number_inner_rows
            j = self%start_row_indx + k - 1
            N_indx = j + 1
            S_indx = j - 1
            do p = 1, self%number_row_sections(k)
                do i = self%start_inner_indx_x(p, k), self%end_inner_indx_x(p,k) 
                    E_indx = i+1
                    W_indx = i-1
                    self%residual(i,j) = self%sourceTerm(i,j) - (self%solution(i, N_indx) + self%solution(i, S_indx)) * self%coeffY - &
                        (self%solution(E_indx, j) + self%solution(W_indx, j)) * self%coeffX - self%solution(i,j) / self%centerCoeff
                end do
            end do
        end do
        !$OMP end do nowait

        ! Sweep upper/lower rows
        !$OMP do 
        do i = 2, self%N_x-1
            i_finest = i * self%x_indx_step - self%x_indx_step + 1

            ! Lower row
            if (self%world%boundary_conditions(i_finest, 1) == 2) then
                self%residual(i,1) = self%sourceTerm(i,1) - 2.0d0 * self%solution(i, 2) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX - self%solution(i,1)/self%centerCoeff
            else if (self%world%boundary_conditions(i_finest, 1) == 3) then
                self%residual(i,1) = self%sourceTerm(i,1) - (self%solution(i, 2) + self%solution(i, self%N_y-1)) * self%coeffY - &
                    (self%solution(i-1, 1) + self%solution(i+1, 1)) * self%coeffX -  self%solution(i,1)/self%centerCoeff
                self%solution(i, self%N_y) = self%solution(i,1)
            end if

            ! Upper row
            if (self%world%boundary_conditions(i_finest, self%world%N_y) == 2) then
                self%residual(i,self%N_y) = self%sourceTerm(i,self%N_y) - 2.0d0 * self%solution(i, self%N_y-1) * self%coeffY - &
                    (self%solution(i-1, self%N_y) + self%solution(i+1, self%N_y)) * self%coeffX -  self%solution(i, self%N_y)/self%centerCoeff
            end if
        end do
        !$OMP end do nowait

        ! Sweep left/right columns
        !$OMP do 
        do j = 2, self%N_y-1
            j_finest = j * self%y_indx_step - self%y_indx_step + 1

            ! Left column
            if (self%world%boundary_conditions(1, j_finest) == 2) then
                self%residual(1,j) = self%sourceTerm(1,j) - (self%solution(1, j-1) + self%solution(1, j+1)) * self%coeffY - &
                    2.0d0 * self%solution(2, j) * self%coeffX - self%solution(1,j)/self%centerCoeff
            else if (self%world%boundary_conditions(1,j_finest) == 3) then
                self%residual(1,j) = self%sourceTerm(1,j) - (self%solution(1, j+1) + self%solution(1, j-1)) * self%coeffY - &
                    (self%solution(2, j) + self%solution(self%N_x-1, j)) * self%coeffX - self%solution(1,j)/self%centerCoeff
                self%solution(self%N_x, j) = self%solution(1,j)
            end if

            ! Right column
            if (self%world%boundary_conditions(self%world%N_x, j_finest) == 2) then
                self%residual(self%N_x,j) = self%sourceTerm(self%N_x,j) - (self%solution(self%N_x, j-1) + self%solution(self%N_x, j+1)) * self%coeffY - &
                    2.0d0 * self%solution(self%N_x-1, j) * self%coeffX - self%solution(self%N_x,j)/self%centerCoeff
            end if
        end do
        !$OMP end do
        
        !$OMP end parallel
    end subroutine calcResidual_even

! subroutine restriction_even(self, fineGrid, coarseGrid)
!     ! Use gauss seidel information to interpolate array fineGrid to coarseGrid
!     class(GS_Base_Even), intent(in) :: self
!     real(real64), intent(in) :: fineGrid(self%N_x, self%N_y)
!     real(real64), intent(in out) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
!     integer :: i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx, k, p

!     !$OMP parallel private(k, p, i_fine, j_fine, i_coarse, j_coarse, N_indx, E_indx, S_indx, W_indx)
!     !$OMP do collapse(2)
!     ! loop over fine grid nodes which overlap with coarse grid nodes
!     do k = self%startRowCoarse, self%endRowCoarse, 2
!         do p = self%startColCoarse, self%endColCoarse, 2
!             ! calculate fine indices around overlapping index
!             j_fine = self%startRow + k - 1
!             j_coarse = (j_fine + 1)/2
!             N_indx = self%vertIndx(1, k)
!             S_indx = self%vertIndx(2, k)
!             i_fine = self%startCol + p - 1
!             i_coarse = (i_fine + 1)/2
!             E_indx = self%horzIndx(1, p)
!             W_indx = self%horzIndx(2, p)

!             ! Bilinear interpolation
!             coarseGrid(i_coarse, j_coarse) = 0.25d0 * fineGrid(i_fine, j_fine) + &
!             0.125d0 * (fineGrid(E_indx, j_fine) + fineGrid(W_indx, j_fine) +  &
!             fineGrid(i_fine, N_indx) + fineGrid(i_fine, S_indx)) + &
!             0.0625d0 * (fineGrid(E_indx, N_indx) + fineGrid(E_indx, S_indx) + &
!             fineGrid(W_indx, N_indx)+ fineGrid(W_indx, S_indx))
!         end do
!     end do
!     !$OMP end do
!     !$OMP end parallel
! end subroutine restriction_even

! subroutine prolongation_even(self, fineGrid, coarseGrid)
!     ! Prolongate operator from coarse to fine grid using gauss seidel data
!     class(GS_Base_Even), intent(in) :: self
!     real(real64), intent(in out) :: fineGrid(self%N_x, self%N_y)
!     real(real64), intent(in) :: coarseGrid((self%N_x+1)/2, (self%N_y+1)/2)
!     integer(int32) :: i_coarse, j_coarse, i_fine, j_fine, N_x_coarse, N_y_coarse
!     N_x_coarse = (self%N_x+1)/2
!     N_y_coarse = (self%N_y+1)/2
!     ! Each fine node which overlaps coarse node is black, take advantage of that by only going through black nodes
!     !$OMP parallel private(i_coarse, j_coarse, i_fine, j_fine)
!     ! do each do loop individually because can't think of fewer loops which do several calculations without doing same operation on each
!     !$OMP do collapse(2)
!     ! set similar nodes
!     do j_coarse = 1, N_y_coarse
!         do i_coarse = 1, N_x_coarse
!             i_fine = i_coarse*2 - 1
!             j_fine = j_coarse*2-1
!             ! Add overlapping coarse nodes
!             fineGrid(i_fine, j_fine) = fineGrid(i_fine, j_fine) + coarseGrid(i_coarse, j_coarse)
!         end do
!     end do
!     !$OMP end do nowait
!     !$OMP do collapse(2)
!     ! go horizontal each row,
!     do j_coarse = 1, N_y_coarse
!         do i_coarse = 1, N_x_coarse-1
!             i_fine = i_coarse*2 - 1
!             j_fine = j_coarse*2 - 1

!             ! simple interpolation
!             fineGrid(i_fine+1, j_fine) = fineGrid(i_fine+1, j_fine) + coarseGrid(i_coarse+1, j_coarse) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
!         end do
!     end do
!     !$OMP end do nowait
!     !$OMP do collapse(2)
!     ! go vertical
!     do j_coarse = 1, N_y_coarse-1
!         do i_coarse = 1, N_x_coarse
!             i_fine = i_coarse*2 - 1
!             j_fine = j_coarse*2 - 1
!             ! Simple interpolation
!             fineGrid(i_fine, j_fine+1) = fineGrid(i_fine, j_fine+1) + coarseGrid(i_coarse, j_coarse+1) * 0.5d0 + coarseGrid(i_coarse, j_coarse) * 0.5d0
!         end do
!     end do
!     !$OMP end do nowait
!     !$OMP do collapse(2)
!     ! add offset fine nodes which require 4 point interpolation
!     ! first loop through SW corner coarse nodes
!     do j_coarse = 1, N_y_coarse-1
!         do i_coarse = 1, N_x_coarse-1
!             i_fine = i_coarse*2 - 1
!             j_fine = j_coarse*2 - 1
!             ! simple interpolation
!             fineGrid(i_fine+1, j_fine+1) = fineGrid(i_fine+1, j_fine+1) + coarseGrid(i_coarse, j_coarse) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse) * 0.25d0 &
!                 + coarseGrid(i_coarse, j_coarse+1) * 0.25d0 + coarseGrid(i_coarse+1, j_coarse+1) * 0.25d0
!         end do
!     end do
!     !$OMP end do
!     !$OMP end parallel
! end subroutine prolongation_even


end module mod_GS_Base_Even