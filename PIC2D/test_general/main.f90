program main
    use iso_fortran_env, only: int32, real64, int64
    use constants
    use mod_CSRMAtrix
    use mod_pardisoSolver
    use mod_MGSolver
    use mod_PreCondCGSolver
    use mod_PreCondBiCGSTABSolver
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_GS_Base
    use mod_particle
    use mod_rand_generator
    use omp_lib
    implicit none

    integer(int32) :: N_x = 501, N_y = 501, numThreads = 6
    type(Particle), allocatable :: particle_list(:)
    class(domain_base), allocatable, target :: world
    class(MGSolver), allocatable :: mg_solver
    real(real64), allocatable :: E_Field(:, :, :)
    integer(int32) :: NESW_wallBoundaries(4), matDimension, i, j, k, numberStages, startTime, endTime, timingRate, numberPreSmoothOper, numberPostSmoothOper, numberIter
    integer :: upperBound, lowerBound, rightBound, leftBound, stageInt, curv_grid_type_x, curv_grid_type_y, mat_dimension
    integer :: inner_box_first_y, inner_box_last_y, inner_box_first_x, inner_box_last_x, i_thread
    real(real64) :: upperPhi, rightPhi, lowerPhi, leftPhi, innerPhi
    real(real64) :: NESW_phiValues(4), rho, omega
    real(real64) :: Length = 0.05, Width = 0.05, delX, delY
    real(real64) :: relTol, stepTol, temp_real, n_ave, del_t, T_e
    logical :: evenGridBool, redBlackBool, Krylov_bool, center_box_bool
    integer(int32) :: num_part_per_cell = 200
    integer(int64) :: num_part_total

    call execute_command_line("rm -r *.dat")
    call mkl_set_num_threads(numThreads)
    call omp_set_num_threads(numThreads)
    call omp_set_max_active_levels(numThreads)

    call change_global_thread(numThreads)
    
    
    evenGridBool = .true.
    redBlackBool = .true.
    Krylov_bool = .false.
    center_box_bool = .false.
    curv_grid_type_x = 0
    curv_grid_type_y = 0
    
    numberStages = 4
    call checkNodeDivisionMG(N_x, N_y, numberStages)
    ! call change_global_N(N_x, N_y)

    call initialize_rand_PCG(numThreads, .true.)
    

    ! More skewed delX and delY, more smoothing operations needed
    numberPreSmoothOper = 2
    numberPostSmoothOper = 2
    numberIter = 50000
    omega = 1.5d0
    relTol = 1.d-8
    stepTol = 1.d-6
    n_ave = 1.d12
    rho = e_charge * n_ave

    NESW_wallBoundaries(1) = 1 ! North
    NESW_wallBoundaries(2) = 2 ! East
    NESW_wallBoundaries(3) = 1 ! South
    NESW_wallBoundaries(4) = 2 ! West

    NESW_phiValues(1) = 0.0d0
    NESW_phiValues(2) = 0.0d0
    NESW_phiValues(3) = 0.0d0
    NESW_phiValues(4) = 0.0d0
    innerPhi = 0.0d0
    
    
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
        call world%form_boundary_conditions(upperBound, rightBound, lowerBound, leftBound, inner_box_first_x, inner_box_last_x, inner_box_first_y, inner_box_last_y, center_box_bool)
        call world%get_number_cells()
    else
        world = domain_curv(N_x, N_y, Length, &
        Width, curv_grid_type_x, curv_grid_type_y, delX, delY)
        call world%form_boundary_conditions(upperBound, rightBound, lowerBound, leftBound, inner_box_first_x, inner_box_last_x, inner_box_first_y, inner_box_last_y, center_box_bool)
        call world%get_number_cells()
        call world%generate_node_volume()
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


    num_part_total = num_part_per_cell * world%number_total_cells
    T_e = 2.0d0
    call change_global_numPart(1)
    allocate(particle_list(number_charged_particles))
    particle_list(1) = Particle(mass_electron, e_charge, 1.0d0, num_part_total, 2*num_part_total, 'e', world%N_x, world%N_y)
    call particle_list(1)%initialize_weight_from_n_ave(n_ave, world)
    call particle_list(1)%initialize_rand_uniform(world)
    call particle_list(1)%initialize_maxwellian_temperature(T_e)
    print *, 'Particle temp is:', particle_list(1)%getKEAve() * 2.0d0 / 3.0d0
    call system_clock(count_rate = timingRate)
    call system_clock(startTime)
    call particle_list(1)%interpolation_particle_to_nodes()
    call system_clock(endTime)
    print *, 'particle interpolation took', real(endTime - startTime)/real(timingRate), 'seconds'
    print *, ''



    
    select type (world)
    type is (domain_uniform)
        del_t = 0.5d0 * min(world%del_x, world%del_y)/sqrt(2.0d0 * T_e * e_charge/mass_electron)
    end select

    call system_clock(startTime)
    call get_poisson_source_term(solver, particle_list, 1, world)
    call system_clock(endTime)
    print *, 'particle collecting particle source term took', real(endTime - startTime)/real(timingRate), 'seconds'
    print *, ''
    

    ! ! !$OMP parallel private(k, i, j)
    ! ! !$OMP do collapse(2)
    ! ! do j = 1, solver%N_y
    ! !     do i = 1, solver%N_x
    ! !         k = (j-1) * N_x + i
    ! !         if (world%boundary_conditions(i,j) /= 1) then
    ! !             ! stageOne%sourceTerm(i,j) = -rho/eps_0
    ! !             solver%sourceTerm(i,j) = -rho/epsilon_0
    ! !         end if
    ! !     end do
    ! ! end do
    ! ! !$OMP end do
    ! ! !$OMP end parallel


    call system_clock(startTime)
    call mg_solver%solve(stepTol, relTol, 1)
    call system_clock(endTime)
    print *, 'Took', mg_solver%numIter, 'iterations for MG'
    print *, 'MG took', real(endTime - startTime)/real(timingRate), 'seconds'
    print *, ''

    ! get E-field
    allocate(E_Field(2,N_x, N_y))
    E_field = 0.0d0
    select type (world)
    type is (domain_uniform)
        call system_clock(startTime)
        call get_electric_field(E_Field, solver, world)
        call system_clock(endTime)
        print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds for making electric field'
        print *, ''

        open(41,file='E_x.dat', form='UNFORMATTED', access = 'stream', status = 'new')
        write(41) E_Field(1,:,:)
        close(41)

        open(41,file='E_y.dat', form='UNFORMATTED', access = 'stream', status = 'new')
        write(41) E_Field(2,:,:)
        close(41)

        ! print *, 'Doing particle push'
        ! print *, 'amount total particles', particle_list(1)%number_particles_thread
        ! call system_clock(startTime)
        ! call push_particles_uniform(particle_list, number_charged_particles, E_Field, world, del_t)
        ! ! call delete_particles_uniform(particle_list, number_charged_particles)
        ! call system_clock(endTime)
        ! print *, 'Took', real(endTime - startTime)/real(timingRate), 'seconds for particle push'
        ! print *, 'amount total particles', particle_list(1)%number_particles_thread
    end select

    

    

    
    
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

    subroutine get_electric_field(E_field, first_smoother, world)
        type(domain_uniform), intent(in) :: world
        real(real64), intent(in out) :: E_field(2,world%N_x, world%N_y)
        class(GS_Base), intent(in) :: first_smoother
        integer(int32) :: i, i_thread, j, p, k, part_num, start_indx, end_indx

        associate(phi => first_smoother%solution)
        !$OMP parallel private(k, p, i, j, part_num, i_thread, start_indx, end_indx)

        ! Inner Nodes
        !$OMP do
        do k = 1, first_smoother%number_inner_rows
            j = first_smoother%start_row_indx + k - 1
            end_indx = 1
            do p = 1, first_smoother%number_row_sections(k)
                start_indx = first_smoother%start_inner_indx_x(p, k)
                

                ! Left side boundary
                if (world%boundary_conditions(start_indx-1, j) == 1) then
                    ! Dirichlet to left
                    E_field(1, start_indx-1,j) = (phi(start_indx-1,j) - phi(start_indx, j))/world%del_x

                    ! Go backwards until latest end_indx and get E_y on dirichlet boundary
                    do i = start_indx - 1, end_indx, -1
                        if (world%boundary_conditions(i, j+1) /= 1) then
                            E_field(2, i, j) = (phi(i,j) - phi(i, j+1))/world%del_y
                        else
                            E_field(2, i, j) = (phi(i,j-1) - phi(i, j))/world%del_y
                        end if
                    end do
                else 
                    E_field(2, 1, j) = 0.5d0 * (phi(1,j-1) - phi(1, j+1))/world%del_y
                    ! Only do periodic in x since neumann stays at 0
                    if (world%boundary_conditions(start_indx-1, j) == 3) then
                        E_field(1, 1,j) = 0.5d0 * (phi(world%N_x-1,j) - phi(2, j))/world%del_x
                        E_field(1, world%N_x, j) = E_field(1, 1, j)
                    end if
                end if 

                end_indx = first_smoother%end_inner_indx_x(p,k)
                ! Inner Nodes
                do i = start_indx, end_indx 
                    E_field(1, i,j) = 0.5d0 * (phi(i-1,j) - phi(i+1, j))/world%del_x
                    E_field(2, i,j) = 0.5d0 * (phi(i,j-1) - phi(i, j+1))/world%del_y
                end do

                ! right-most side boundary
                if (world%boundary_conditions(end_indx+1, j) == 1) then
                    ! Dirichlet to right
                    E_field(1, end_indx+1,j) = (phi(end_indx,j) - phi(end_indx+1, j))/world%del_x
                    if (world%boundary_conditions(end_indx+1, j+1) /= 1) then
                        E_field(2,end_indx+1, j) = (phi(end_indx+1,j) - phi(end_indx+1, j+1))/world%del_y
                    else
                        E_field(2, end_indx+1, j) = (phi(end_indx+1,j-1) - phi(end_indx+1, j))/world%del_y
                    end if
                else
                    E_field(2, world%N_x, j) = 0.5d0 * (phi(world%N_x,j-1) - phi(world%N_x, j+1))/world%del_y
                end if 

                ! reset latest dirichlet end point for left E_y sweep on next loop
                end_indx = end_indx + 2
            end do

            ! if didn't get to end of row, need to still evaluate E_y along dirichlet edge to the right
            if (end_indx <= world%N_x) then
                do i = end_indx, world%N_x
                    if (world%boundary_conditions(i, j+1) /= 1) then
                        E_field(2,i, j) = (phi(i,j) - phi(i, j+1))/world%del_y
                    else
                        E_field(2, i, j) = (phi(i,j-1) - phi(i, j))/world%del_y
                    end if
                end do
            end if

            
        end do
        !$OMP end do nowait

        ! go through lowest and highest row without plasma nodes
        !$OMP sections
        !$OMP section
        ! Do lowest row
        if (first_smoother%start_row_indx > 2) then
            ! lowest row not along j = 1 boundary, so only dirichlet
            ! go through entire row and check if north node is non-dirichlet, in which case assign E_y
            j = first_smoother%start_row_indx - 1
            do i = 1, world%N_x
                if (world%boundary_conditions(i, j+1) /= 1) then
                    E_field(2,i, j) = (phi(i,j) - phi(i, j+1))/world%del_y
                end if
            end do
        else
            ! Go along bottom row
            end_indx = 1
            j = 1
            do p = 1, first_smoother%number_bottom_row_sections  
                start_indx = first_smoother%start_bottom_row_indx(p)
                
                ! Left side boundary
                if (world%boundary_conditions(start_indx-1, j) == 1) then
                    ! Dirichlet to left
                    E_field(1, start_indx-1,j) = (phi(start_indx-1,j) - phi(start_indx, j))/world%del_x
    
                    ! Go backwards until latest end_indx and get E_y on dirichlet boundary
                    do i = start_indx - 1, end_indx, -1
                        if (world%boundary_conditions(i, j+1) /= 1) then
                            E_field(2, i, j) = (phi(i,j) - phi(i, j+1))/world%del_y
                        end if
                    end do
                end if 
    
                end_indx = first_smoother%end_bottom_row_indx(p)
                ! Inner Nodes
                do i = start_indx, end_indx 
                    E_field(1, i,j) = 0.5d0 * (phi(i-1,j) - phi(i+1, j))/world%del_x
                    if (world%boundary_conditions(i, j) == 3) then
                        E_field(2, i,j) = 0.5d0 * (phi(i,world%N_y-1) - phi(i, 2))/world%del_y
                    end if
                end do
                ! right-most side boundary
                if (world%boundary_conditions(end_indx+1, j) == 1) then
                    ! Dirichlet to right
                    E_field(1, end_indx+1,j) = (phi(end_indx,j) - phi(end_indx+1, j))/world%del_x
                    if (world%boundary_conditions(end_indx+1, j+1) /= 1) then
                        E_field(2,end_indx+1, j) = (phi(end_indx+1,j) - phi(end_indx+1, j+1))/world%del_y
                    end if
                end if 

                end_indx = end_indx + 2
                
            end do

            ! if didn't get to end of row, need to still evaluate E_y along dirichlet edge to the right
            if (end_indx <= world%N_x) then
                do i = end_indx, world%N_x
                    if (world%boundary_conditions(i, j+1) /= 1) then
                        E_field(2,i, j) = (phi(i,j) - phi(i, j+1))/world%del_y
                    end if
                end do
            end if
            
        end if
        !$OMP section
        ! Do highest row
        if (first_smoother%end_row_indx < world%N_y-1) then
            ! highest row not along j = N_y-1 boundary, so only dirichlet
            ! go through entire row and check if south node is non-dirichlet, in which case assign E_y
            j = first_smoother%end_row_indx + 1
            do i = 1, world%N_x
                if (world%boundary_conditions(i, j-1) /= 1) then
                    E_field(2,i, j) = (phi(i,j-1) - phi(i, j))/world%del_y
                end if
            end do
        else
            ! Go along top row
            end_indx = 1
            j = world%N_y
            do p = 1, first_smoother%number_bottom_row_sections  
                start_indx = first_smoother%start_bottom_row_indx(p)
                
                ! Left side boundary
                if (world%boundary_conditions(start_indx-1, j) == 1) then
                    ! Dirichlet to left
                    E_field(1, start_indx-1,j) = (phi(start_indx-1,j) - phi(start_indx, j))/world%del_x
    
                    ! Go backwards until latest end_indx and get E_y on dirichlet boundary
                    do i = start_indx - 1, end_indx, -1
                        if (world%boundary_conditions(i, j-1) /= 1) then
                            E_field(2, i, j) = (phi(i,j-1) - phi(i, j))/world%del_y
                        end if
                    end do
                end if 
    
                end_indx = first_smoother%end_bottom_row_indx(p)
                ! Inner Nodes
                do i = start_indx, end_indx 
                    E_field(1, i,j) = 0.5d0 * (phi(i-1,j) - phi(i+1, j))/world%del_x
                    if (world%boundary_conditions(i, j) == 3) then
                        E_field(2, i,j) = 0.5d0 * (phi(i,world%N_y-1) - phi(i, 2))/world%del_y
                    end if
                end do
                ! right-most side boundary
                if (world%boundary_conditions(end_indx+1, j) == 1) then
                    ! Dirichlet to right
                    E_field(1, end_indx+1,j) = (phi(end_indx,j) - phi(end_indx+1, j))/world%del_x
                    if (world%boundary_conditions(end_indx+1, j-1) /= 1) then
                        E_field(2,end_indx+1, j) = (phi(end_indx+1,j-1) - phi(end_indx+1, j))/world%del_y
                    end if
                end if 

                end_indx = end_indx + 2
                
            end do

            ! if didn't get to end of row, need to still evaluate E_y along dirichlet edge to the right
            if (end_indx <= world%N_x) then
                do i = end_indx, world%N_x
                    if (world%boundary_conditions(i, j-1) /= 1) then
                        E_field(2,i, j) = (phi(i,j-1) - phi(i, j))/world%del_y
                    end if
                end do
            end if
            
        end if
        !$OMP end sections nowait
        !$OMP end parallel 
        end associate

    end subroutine get_electric_field

    subroutine get_poisson_source_term(first_smoother, particle_list, number_charged_particles, world)
        class(domain_base), intent(in) :: world
        class(GS_Base), intent(in out) :: first_smoother
        integer, intent(in) :: number_charged_particles
        type(Particle), intent(in) :: particle_list(number_charged_particles)
        integer(int32) :: i, i_thread, j, p, k, part_num
        real(real64) :: inv_epsilon_0
        inv_epsilon_0 = 1.0d0/epsilon_0

        ! can't do separation over threads since need to add to single array of (N_x,N_y), would cause race conditions
        ! Each thread goes over different rows in source_term, which then combines particles separated into different thread arrays

        select type (world)
        type is (domain_uniform)
            !$OMP parallel private(k, p, i, j, part_num, i_thread)
            !$OMP workshare
            first_smoother%sourceTerm = 0.0d0
            !$OMP end workshare

            ! Inner Nodes
            !$OMP do
            do k = 1, first_smoother%number_inner_rows
                j = first_smoother%start_row_indx + k - 1
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do p = 1, first_smoother%number_row_sections(k)
                            do i = first_smoother%start_inner_indx_x(p, k), first_smoother%end_inner_indx_x(p,k)  
                                first_smoother%sourceTerm(i,j) = first_smoother%sourceTerm(i,j) - particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(i,j,i_thread) * world%inv_node_volume*inv_epsilon_0
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            !first do corners
            !$OMP sections
            !$OMP section
            if (world%boundary_conditions(1,1) == 2) then
                ! lower left corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(1,1) = first_smoother%sourceTerm(1,1) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(1,1,i_thread) * world%inv_node_volume*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP section
            if (world%boundary_conditions(world%N_x, 1) == 2) then
                ! lower right corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(world%N_x,1) = first_smoother%sourceTerm(world%N_x,1) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(world%N_x,1,i_thread) * world%inv_node_volume*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP section
            if (world%boundary_conditions(1, world%N_y) == 2) then
                ! upper left corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(1,world%N_y) = first_smoother%sourceTerm(1,world%N_y) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(1,world%N_y,i_thread) * world%inv_node_volume*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP section
            if (world%boundary_conditions(world%N_x, world%N_y) == 2) then
                ! upper right corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(world%N_x,world%N_y) = first_smoother%sourceTerm(world%N_x,world%N_y) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(world%N_x,world%N_y,i_thread) * world%inv_node_volume*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP end sections nowait

            !lower boundary
            !$OMP do
            do p = 1, first_smoother%number_bottom_row_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do i = first_smoother%start_bottom_row_indx(p), first_smoother%end_bottom_row_indx(p)
                            first_smoother%sourceTerm(i,1) = first_smoother%sourceTerm(i,1) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(i,1,i_thread) * world%inv_node_volume*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            ! ! upper boundary
            !$OMP do
            do p = 1, first_smoother%number_top_row_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do i = first_smoother%start_top_row_indx(p), first_smoother%end_top_row_indx(p)
                            first_smoother%sourceTerm(i,world%N_y) = first_smoother%sourceTerm(i,world%N_y) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(i,world%N_y,i_thread) * world%inv_node_volume*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            ! left boundary
            !$OMP do
            do p = 1, first_smoother%number_left_column_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do j = first_smoother%start_left_column_indx(p), first_smoother%end_left_column_indx(p)
                            first_smoother%sourceTerm(1,j) = first_smoother%sourceTerm(1,j) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(1,j,i_thread) * world%inv_node_volume*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            ! right boundary
            !$OMP do
            do p = 1, first_smoother%number_right_column_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do j = first_smoother%start_right_column_indx(p), first_smoother%end_right_column_indx(p)
                            first_smoother%sourceTerm(world%N_x,j) = first_smoother%sourceTerm(world%N_x,j) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(world%N_x,j,i_thread) * world%inv_node_volume*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do
            !$OMP end parallel      
        type is (domain_curv)
            !$OMP parallel private(k, p, i, j, part_num, i_thread)
            !$OMP workshare
            first_smoother%sourceTerm = 0.0d0
            !$OMP end workshare

            ! Inner Nodes
            !$OMP do
            do k = 1, first_smoother%number_inner_rows
                j = first_smoother%start_row_indx + k - 1
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do p = 1, first_smoother%number_row_sections(k)
                            do i = first_smoother%start_inner_indx_x(p, k), first_smoother%end_inner_indx_x(p,k)  
                                first_smoother%sourceTerm(i,j) = first_smoother%sourceTerm(i,j) - particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(i,j,i_thread) * world%inv_node_volume(i,j)*inv_epsilon_0
                            end do
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            !first do corners
            !$OMP sections
            !$OMP section
            if (world%boundary_conditions(1,1) == 2) then
                ! lower left corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(1,1) = first_smoother%sourceTerm(1,1) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(1,1,i_thread) * world%inv_node_volume(1,1)*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP section
            if (world%boundary_conditions(world%N_x, 1) == 2) then
                ! lower right corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(world%N_x,1) = first_smoother%sourceTerm(world%N_x,1) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(world%N_x,1,i_thread) * world%inv_node_volume(world%N_x,1)*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP section
            if (world%boundary_conditions(1, world%N_y) == 2) then
                ! upper left corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(1,world%N_y) = first_smoother%sourceTerm(1,world%N_y) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(1,world%N_y,i_thread) * world%inv_node_volume(1,world%N_y)*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP section
            if (world%boundary_conditions(world%N_x, world%N_y) == 2) then
                ! upper right corner
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global   
                        first_smoother%sourceTerm(world%N_x,world%N_y) = first_smoother%sourceTerm(world%N_x,world%N_y) - 4.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(world%N_x,world%N_y,i_thread) * world%inv_node_volume(world%N_x,world%N_y)*inv_epsilon_0
                    end do
                end do
            end if
            !$OMP end sections nowait

            !lower boundary
            !$OMP do
            do p = 1, first_smoother%number_bottom_row_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do i = first_smoother%start_bottom_row_indx(p), first_smoother%end_bottom_row_indx(p)
                            first_smoother%sourceTerm(i,1) = first_smoother%sourceTerm(i,1) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(i,1,i_thread) * world%inv_node_volume(i,1)*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            ! ! upper boundary
            !$OMP do
            do p = 1, first_smoother%number_top_row_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do i = first_smoother%start_top_row_indx(p), first_smoother%end_top_row_indx(p)
                            first_smoother%sourceTerm(i,world%N_y) = first_smoother%sourceTerm(i,world%N_y) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(i,world%N_y,i_thread) * world%inv_node_volume(i,world%N_y)*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            ! left boundary
            !$OMP do
            do p = 1, first_smoother%number_left_column_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do j = first_smoother%start_left_column_indx(p), first_smoother%end_left_column_indx(p)
                            first_smoother%sourceTerm(1,j) = first_smoother%sourceTerm(1,j) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(1,j,i_thread) * world%inv_node_volume(1,j)*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do nowait

            ! right boundary
            !$OMP do
            do p = 1, first_smoother%number_right_column_sections
                do part_num = 1, number_charged_particles
                    do i_thread = 1, number_threads_global
                        do j = first_smoother%start_right_column_indx(p), first_smoother%end_right_column_indx(p)
                            first_smoother%sourceTerm(world%N_x,j) = first_smoother%sourceTerm(world%N_x,j) - 2.0d0 * particle_list(part_num)%q_times_weight * particle_list(part_num)%work_space(world%N_x,j,i_thread) * world%inv_node_volume(world%N_x,j)*inv_epsilon_0
                        end do
                    end do
                end do
            end do
            !$OMP end do
            !$OMP end parallel      
        end select

    end subroutine get_poisson_source_term





    ! subroutine push_particles_uniform(particle_list, number_charged_particles, E_Field, world, del_t)
    !     type(domain_uniform), intent(in) :: world
    !     real(real64), intent(in) :: E_field(2,world%N_x,world%N_y), del_t
    !     integer(int32), intent(in) :: number_charged_particles
    !     type(Particle), intent(in out) :: particle_list(number_charged_particles)
    !     real(real64) :: v_part(2), loc_i, loc_j, q_over_m, d_i, d_j, E_part(2),&
    !         E_SE(2), E_SW(2), E_NW(2), E_NE(2), inv_del_x, inv_del_y, &
    !         v_xi, v_eta, del_t_i, del_t_j, loc_i_new, loc_j_new
    !     integer(int32) :: i_thread, part_idx, corner_i, corner_j, wall_i, wall_j
    !     integer(int64) :: part_num, delete_idx
    !     logical :: delete_bool

    !     inv_del_x = 1.0d0 / world%del_x
    !     inv_del_y = 1.0d0 / world%del_y
    !     !$OMP parallel private(i_thread, q_over_m, d_i, d_j, E_part, part_idx, part_num, &
    !     !$OMP v_part, loc_i, loc_j, E_SE, E_SW, E_NW, E_NE, corner_i, corner_j, delete_idx, &
    !     !$OMP v_xi, v_eta, wall_i, wall_j, delete_bool)
    !     i_thread = omp_get_thread_num() + 1
    !     do part_idx = 1, number_charged_particles
    !         q_over_m = particle_list(part_idx)%q_over_m
    !         delete_idx = 0
    !         do part_num = 1, particle_list(part_idx)%number_particles_thread(i_thread)
    !             !get particle location and velocity in 2D
    !             loc_i = particle_list(part_idx)%logical_position(1, part_num, i_thread)
    !             loc_j = particle_list(part_idx)%logical_position(2, part_num, i_thread)
    !             v_part = particle_list(part_idx)%velocity(3:4, part_num, i_thread)
    !             corner_i = int(loc_i)
    !             corner_j = int(loc_j)
    !             d_i = loc_i - real(corner_i, kind = 8)
    !             d_j = loc_j - real(corner_j, kind = 8)

    !             ! first find E_x and E_y on the nodes
    !             E_SW = E_Field(:, corner_i, corner_j)
    !             E_SE = E_Field(:, corner_i+1, corner_j)
    !             E_NW = E_Field(:, corner_i, corner_j+1)
    !             E_NE = E_Field(:, corner_i+1, corner_j+1)

    !             ! interpolate to particle position
    !             E_part = E_SW * (1.0d0 - d_i) * (1.0d0 - d_j) + E_SE * (d_i) * (1.0d0-d_j) + &
    !                 E_NW * (1.0d0-d_i) * (d_j) + E_NE * (d_i) * (d_j)
                
    !             ! solve for new velocity and position
    !             v_part = v_part + q_over_m * E_part * del_t

    !             ! place on boundary if outside of boundary
    !             v_xi = v_part(1) * inv_del_x
    !             v_eta = v_part(2) * inv_del_y
    !             loc_i_new = loc_i + v_xi * del_t
    !             loc_j_new = loc_j + v_eta * del_t

    
    !             delete_bool = .false.

    !             ! Found in general that doing deletions after in new subroutine is slower, so do in place

    !             ! don't need to check boundaries if particle still in same cell
    !             if (corner_i /= int(loc_i_new) .or. corner_j /= int(loc_j_new)) then
    !                 if (loc_i_new > world%N_x) then
    !                     ! backtrack to wall position where it left
    !                     wall_i = world%N_x
    !                     del_t_i = (loc_i_new - real(wall_i, kind = 8))/v_xi
    
    !                     ! get boundary j indices for particle
    !                     wall_j = max(min(int(loc_j_new - v_eta * del_t_i), world%N_y-1), 1)
    !                     delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i, wall_j+1) == 1
    !                     if (.not. delete_bool) then
    !                         if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i, wall_j+1) == 2) then
    !                             ! Neumann boundary
    !                             loc_i_new = 2.0d0 * real(wall_i, kind = 8) - loc_i_new
    !                             v_xi = -v_xi
    !                             v_part(1) = -v_part(1)
    !                         else
    !                             loc_i_new = loc_i_new - real(world%N_x, kind = 8) + 1.0d0
    !                         end if
    !                     end if
    !                 else if (loc_i_new < 1) then
    !                     ! backtrack to wall position where it left
    !                     wall_i = 1
    !                     del_t_i = (loc_i_new - real(wall_i, kind = 8))/v_xi
    
    !                     ! get boundary j indices for particle
    !                     wall_j = max(min(int(loc_j_new - v_eta * del_t_i), world%N_y-1), 1)
    !                     delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i, wall_j+1) == 1
    !                     if (.not. delete_bool) then
    !                         if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i, wall_j+1) == 2) then
    !                             ! Neumann boundary
    !                             loc_i_new = 2.0d0 - loc_i_new
    !                             v_xi = -v_xi
    !                             v_part(1) = -v_part(1)
    !                         else
    !                             loc_i_new = real(world%N_x, kind = 8) - 1.0d0 + loc_i_new
    !                         end if
    !                     end if
    !                 end if
    
    
    !                 ! take care of any issue with particle outside eta boundary
    !                 ! it is possible for particle to be outside both x and y, so proceed if not deleted at left-right boundary
    !                 if (.not. delete_bool .and. loc_j_new > world%N_y) then
    !                     ! backtrack to wall position where it left
    !                     wall_j = world%N_y
    !                     del_t_j = (loc_j_new - real(wall_j, kind = 8))/v_eta
    
    !                     ! get boundary i indices for particle
    !                     wall_i = max(min(int(loc_i_new - v_xi * del_t_j), world%N_x-1), 1)
    !                     delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i+1, wall_j) == 1
    !                     if (.not. delete_bool) then
    !                         if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i+1, wall_j) == 2) then
    !                             ! Neumann boundary
    !                             loc_j_new = 2.0d0 * real(wall_j, kind = 8) - loc_j_new
    !                             v_eta = -v_eta
    !                             v_part(2) = -v_part(2)
    !                         else
    !                             loc_j_new = loc_j_new - real(world%N_y, kind = 8) + 1.0d0
    !                         end if
    !                     end if
    !                 else if (.not. delete_bool .and. loc_j_new < 1) then
    !                     ! backtrack to wall position where it left
    !                     wall_j = 1
    !                     del_t_j = (loc_j_new - real(wall_j, kind = 8))/v_eta
    
    !                     ! get boundary i indices for particle
    !                     wall_i = max(min(int(loc_i_new - v_xi * del_t_j), world%N_x-1), 1)
    !                     delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i+1, wall_j) == 1
    !                     if (.not. delete_bool) then
    !                         if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i+1, wall_j) == 2) then
    !                             ! Neumann boundary
    !                             loc_j_new = 2.0d0- loc_j_new
    !                             v_eta = -v_eta
    !                             v_part(2) = -v_part(2)
    !                         else
    !                             loc_j_new = real(world%N_y, kind = 8) - 1.0d0 + loc_j_new
    !                         end if
    !                     end if
    !                 end if

    !                 if (.not. delete_bool) then
    !                     ! particle pushed back into domain from boundary
    !                     wall_i = int(loc_i_new) + (1 - INT(SIGN(1.0d0, v_xi)))/2
    !                     wall_j = int(loc_j_new) + (1 - INT(SIGN(1.0d0, v_eta)))/2
                
    !                     ! We'll assume that particle doesn't go across many cells, so likelihood of passing dirichlet boundary and ending up in cell surrounded by plasma nodes is low
    !                     ! corner node (wall_i, wall_j) needs to be dirichlet for particle to have chance of passing dirichlet wall
    !                     if (world%boundary_conditions(wall_i, wall_j) == 1) then
    !                         ! Find time to each wall it could have passed through
    !                         del_t_i = (loc_i_new - wall_i)/v_xi
    !                         del_t_j = (loc_j_new - wall_j)/v_eta

    !                         if (del_t_i < del_t_j) then
    !                             ! hits along left-right wall, check other j-th index node
    !                             wall_j = int(loc_j_new)
    !                             delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i, wall_j+1) == 1
    !                         else
    !                             ! hits along up-down wall, check other i-th index node
    !                             wall_i = int(loc_i_new)
    !                             delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i+1, wall_j) == 1
    !                         end if
    !                     end if
    !                 end if
    !             end if


    !             if (.not. delete_bool) then
    !                 particle_list(part_idx)%logical_position(1, part_num - delete_idx, i_thread) = loc_i_new
    !                 particle_list(part_idx)%logical_position(2, part_num - delete_idx, i_thread) = loc_j_new
    !                 particle_list(part_idx)%velocity(3:4, part_num - delete_idx, i_thread) = v_part
    !                 particle_list(part_idx)%velocity(5, part_num - delete_idx, i_thread) = particle_list(part_idx)%velocity(5, part_num, i_thread)
    !             else
    !                 delete_idx = delete_idx + 1
    !             end if

                
                
    !         end do
    !         particle_list(part_idx)%number_particles_thread(i_thread) = particle_list(part_idx)%number_particles_thread(i_thread) - delete_idx
    !         particle_list(part_idx)%number_deletes_thread(i_thread) = delete_idx
    !     end do

    !     !$OMP end parallel

    ! end subroutine push_particles_uniform


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