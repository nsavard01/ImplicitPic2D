module mod_particle_per_cell

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_particle
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_rand_generator
    use omp_lib
    implicit none

    ! The following arrays will be used by all particles for threaded operations as temporaries, so only allocate a single time
    real(real64), allocatable, public, protected :: logical_position_overflow(:,:,:), velocity_overflow(:,:,:)
    integer(int64), allocatable, public, protected :: number_particles_overflow_thread(:)
    integer(int32), allocatable, public, protected :: number_particles_added_cell_thread(:,:,:)
    integer(int64), public, protected :: max_indx_overflow = 0, max_number_overflow = 0

    ! Particle contains particle properties and stored values in phase space
    type, extends(Particle) :: Particle_Per_Cell
        integer(int64), allocatable :: cell_starting_indx(:,:), cell_ending_indx(:,:)  ! starting/ending index (or particle number) for particle in cell per thread [i index, j index, thread #]
        ! We assume for the moment that particles will be equally spread through each thread, so starting/ending index is the same for all threads
        integer(int32), allocatable :: number_particles_cell_thread(:,:,:) ! total number of particles in a cell for a particular thread [i index, j index, thread #]

    contains
        procedure, public, pass(self) :: initialize_rand_uniform => initialize_rand_uniform_cell
        procedure, public, pass(self) :: interpolation_particle_to_nodes => interpolation_particle_to_nodes_cell
        procedure, public, pass(self) :: particle_mover_uniform => particle_mover_uniform_cell
        procedure, public, pass(self) :: initialize_maxwellian_temperature => initialize_maxwellian_temperature_cell
        procedure, public, pass(self) :: get_sum_totals => get_sum_totals_cell
    end type
    interface Particle_Per_Cell
        module procedure :: particle_per_cell_constructor
    end interface Particle_Per_Cell

contains

    type(Particle_Per_Cell) function particle_per_cell_constructor(mass, q, w_p, N_p, finalIdx, particleName, world) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, w_p
        integer(int64), intent(in) :: N_p, finalIdx
        class(domain_base), intent(in) :: world
        character(*), intent(in) :: particleName
        integer(int64) :: part_num_separation, k, leftOver_number
        integer(int32) :: i, j
        call self%initialize_base_variables(mass, q, w_p, N_p, finalIdx, particleName, world%N_x, world%N_y)
        max_indx_overflow = max(self%total_number_particles / number_threads_global, max_indx_overflow)
        if (.not. allocated(velocity_overflow)) allocate(velocity_overflow(3, max_indx_overflow, number_threads_global))
        if (.not. allocated(logical_position_overflow)) allocate(logical_position_overflow(2, max_indx_overflow, number_threads_global))
        if (.not. allocated(number_particles_overflow_thread)) allocate(number_particles_overflow_thread(number_threads_global))
        number_particles_overflow_thread = self%total_number_particles/number_threads_global

        if (.not. allocated(number_particles_added_cell_thread)) then
            allocate(number_particles_added_cell_thread(world%N_x-1, world%N_y-1, number_threads_global))
            number_particles_added_cell_thread = 0
        end if
        allocate(self%cell_starting_indx(world%N_x-1, world%N_y-1), self%cell_ending_indx(world%N_x-1, world%N_y-1), &
        self%number_particles_cell_thread(world%N_x-1, world%N_y-1, number_threads_global))
        part_num_separation = self%max_indx / world%number_total_cells
        leftOver_number = self%max_indx - part_num_separation * world%number_total_cells
        k = 1
        do j = 1, world%N_y-1
            do i = 1, world%N_x-1
                if (world%boundary_conditions(i, j) == 0 .or. world%boundary_conditions(i+1,j) == 0 &
                    .or. world%boundary_conditions(i,j+1) == 0 .or. world%boundary_conditions(i+1,j+1) == 0) then
                    if (leftOver_number > 0) then
                        ! while still left over numbers have total cell number + 1 what it would be otherwise
                        self%cell_starting_indx(i,j) = k
                        self%cell_ending_indx(i,j) = k + part_num_separation
                        k = k + part_num_separation + 1
                        leftOver_number = leftOver_number - 1
                    else
                        self%cell_starting_indx(i,j) = k
                        self%cell_ending_indx(i,j) = k + part_num_separation - 1
                        k = k + part_num_separation
                    end if
                end if
            end do
        end do
        self%number_particles_cell_thread = 0
        
    end function particle_per_cell_constructor


    subroutine initialize_rand_uniform_cell(self, world)
        ! distribute particles randomly over the domain
        class(Particle_Per_Cell), intent(in out) :: self
        class(domain_base), intent(in) :: world
        integer(int32) :: i_thread, int_xi, int_eta, last_xi, last_eta, number_cell
        integer(int64) :: part_num, number_particles, idx_num
        real(real64) :: particle_remainder, real_number_cell
        !$OMP parallel private(i_thread, part_num,int_xi, int_eta, idx_num, number_particles, last_xi, last_eta, number_cell, real_number_cell, particle_remainder)
        i_thread = omp_get_thread_num() + 1
        number_particles = number_particles_overflow_thread(i_thread)
        particle_remainder = 0.0
        idx_num = 0
        do int_eta = 1, world%N_y-1
            do int_xi = 1, world%N_x-1
                if (world%boundary_conditions(int_xi, int_eta) == 0 .or. world%boundary_conditions(int_xi+1,int_eta) == 0 &
                    .or. world%boundary_conditions(int_xi,int_eta+1) == 0 .or. world%boundary_conditions(int_xi+1,int_eta+1) == 0) then
                    select type (world)
                    type is (domain_uniform)
                        real_number_cell = real(number_particles, kind = 8) / real(world%number_total_cells, kind = 8)
                    type is (domain_curv)
                        real_number_cell = real(number_particles, kind = 8)  * world%del_x(int_xi) * world%del_y(int_eta)/ world%total_cell_area
                    end select
                    number_cell = int(real_number_cell + particle_remainder)
                    particle_remainder = particle_remainder + real_number_cell - real(number_cell, kind = 8)
                    last_xi = int_xi
                    last_eta = int_eta
                else
                    number_cell = 0
                end if
                self%number_particles_cell_thread(int_xi, int_eta, i_thread) = number_cell
                idx_num = self%cell_starting_indx(int_xi, int_eta)
                do part_num = 1, number_cell
                    self%logical_position(1, idx_num, i_thread) = pcg32_random_r(state_PCG) + real(int_xi, kind = 8)
                    self%logical_position(2, idx_num, i_thread) = pcg32_random_r(state_PCG) + real(int_eta, kind = 8)
                    idx_num = idx_num + 1
                end do
            end do
        end do

        !$OMP end parallel
    end subroutine initialize_rand_uniform_cell

    subroutine interpolation_particle_to_nodes_cell(self, i_thread, const)
        ! interpolate particles to work space array
        class(Particle_Per_Cell), intent(in) :: self
        integer(int32), intent(in) :: i_thread
        real(real64), intent(in) :: const
        integer(int32) :: i_cell, j_cell, N_x_cell, N_y_cell
        integer(int64) :: part_num, start_point, end_point
        real(real64) :: d_i, d_j, xi, eta, NW, NE, SW, SE, real_j_cell, real_i_cell
        N_x_cell = size(self%cell_starting_indx, DIM = 1)
        N_y_cell = size(self%cell_starting_indx, DIM = 2)
        do j_cell = 1, N_y_cell
            real_j_cell = real(j_cell, kind = 8)
            do i_cell = 1, N_x_cell
                real_i_cell = real(i_cell, kind = 8)
                start_point = self%cell_starting_indx(i_cell, j_cell)
                end_point = self%cell_starting_indx(i_cell, j_cell) + self%number_particles_cell_thread(i_cell, j_cell, i_thread) - 1
                SE = 0.0d0
                SW = 0.0d0
                NW = 0.0d0
                NE = 0.0d0
                do part_num = start_point, end_point
                    xi = self%logical_position(1,part_num,i_thread)
                    eta = self%logical_position(2,part_num,i_thread)
                    d_i = xi - real_i_cell
                    d_j = eta - real_j_cell
                    SW = SW + (1.0d0-d_i) * (1.0d0-d_j)
                    SE = SE + (d_i) * (1.0d0-d_j)
                    NW = NW + (1.0d0-d_i) * (d_j)
                    NE = NE + d_i * d_j
                end do
                particle_work_space(i_cell,j_cell, i_thread) = particle_work_space(i_cell,j_cell, i_thread) + SW * const
                particle_work_space(i_cell+1,j_cell, i_thread) = particle_work_space(i_cell+1,j_cell, i_thread) + SE * const
                particle_work_space(i_cell,j_cell+1, i_thread) = particle_work_space(i_cell,j_cell+1, i_thread) + NW * const
                particle_work_space(i_cell+1,j_cell+1, i_thread) = particle_work_space(i_cell+1,j_cell+1, i_thread) + NE * const
            end do
        end do
    end subroutine interpolation_particle_to_nodes_cell


    subroutine particle_mover_uniform_cell(self, i_thread, E_Field, world, del_t)
        ! place subroutine in particle, then have per thread subroutine rather than keeping track of private variables
        class(Particle_Per_Cell), intent(in out) :: self
        type(domain_uniform), intent(in) :: world
        real(real64), intent(in) :: E_Field(2,world%N_x,world%N_y), del_t
        integer(int32), intent(in) :: i_thread
        real(real64) :: v_part(2), loc_i, loc_j, d_i, d_j, E_part(2),&
        E_SE(2), E_SW(2), E_NW(2), E_NE(2), inv_del_x, inv_del_y, &
        loc_i_new, loc_j_new, i_cell_real, j_cell_real, v_z, del_t_i, del_t_j, v_xi, v_eta
        integer(int32) :: i_cell, j_cell, N_x_cell, N_y_cell, number_particles_cell, wall_i, wall_j, number_particles_outside_cell
        integer(int64) :: part_num, cell_end_indx, cell_start_indx, number_particles_overflow
        logical :: delete_bool
        
        number_particles_overflow = 0
        N_x_cell = world%N_x - 1
        N_y_cell = world%N_y - 1
        inv_del_x = 1.0d0 / world%del_x
        inv_del_y = 1.0d0 / world%del_y
        do j_cell = 1, N_y_cell
            j_cell_real = real(j_cell, kind = 8)
            do i_cell = 1, N_x_cell
                i_cell_real = real(i_cell, kind = 8)
                E_SW = E_Field(:, i_cell, j_cell)
                E_SE = E_Field(:, i_cell+1, j_cell)
                E_NW = E_Field(:, i_cell, j_cell+1)
                E_NE = E_Field(:, i_cell+1, j_cell+1)

                cell_start_indx = self%cell_starting_indx(i_cell, j_cell)
                number_particles_cell = 0
                do part_num = cell_start_indx, cell_start_indx + self%number_particles_cell_thread(i_cell, j_cell, i_thread) - 1
                    loc_i = self%logical_position(1, part_num, i_thread)
                    loc_j = self%logical_position(2, part_num, i_thread)
                    v_part = self%velocity(1:2, part_num, i_thread)
                    v_z = self%velocity(3, part_num, i_thread)
                    d_i = loc_i - i_cell_real
                    d_j = loc_j - j_cell_real
                    E_part = E_SW * (1.0d0 - d_i) * (1.0d0 - d_j) + E_SE * (d_i) * (1.0d0-d_j) + &
                    E_NW * (1.0d0-d_i) * (d_j) + E_NE * (d_i) * (d_j)
            
                    ! solve for new velocity and position
                    v_part = v_part + self%q_over_m * E_part * del_t
                    v_xi = v_part(1) * inv_del_x
                    v_eta = v_part(2) * inv_del_y
                    loc_i_new = loc_i + v_xi * del_t
                    loc_j_new = loc_j + v_eta * del_t

                    if (int(loc_i_new) == i_cell .and. int(loc_j_new) == j_cell) then
                        ! stays in cell, put at beginning of cell array, move in place
                        self%logical_position(1, cell_start_indx + number_particles_cell, i_thread) = loc_i_new
                        self%logical_position(2, cell_start_indx + number_particles_cell, i_thread) = loc_j_new
                        self%velocity(1:2, cell_start_indx + number_particles_cell, i_thread) = v_part
                        self%velocity(3, cell_start_indx + number_particles_cell, i_thread) = v_z
                        number_particles_cell = number_particles_cell + 1
                    else
                        delete_bool = .false.
                        do while ((loc_i_new > world%N_x .or. loc_i_new < 1) .and. .not. delete_bool) 
                            if (loc_i_new > world%N_x) then
                                ! backtrack to wall position where it left
                                wall_i = world%N_x
                                del_t_i = (loc_i_new - real(wall_i, kind = 8))/v_xi
            
                                ! get boundary j indices for particle
                                wall_j = max(min(int(loc_j_new - v_eta * del_t_i), world%N_y-1), 1)
                                delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i, wall_j+1) == 1
                                if (.not. delete_bool) then
                                    if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i, wall_j+1) == 2) then
                                        ! Neumann boundary
                                        loc_i_new = 2.0d0 * real(wall_i, kind = 8) - loc_i_new
                                        v_xi = -v_xi
                                        v_part(1) = -v_part(1)
                                    else
                                        loc_i_new = MODULO(loc_i_new, real(world%N_x-1, kind = 8))
                                    end if
                                end if
                            else
                                ! backtrack to wall position where it left
                                wall_i = 1
                                del_t_i = (loc_i_new - real(wall_i, kind = 8))/v_xi
            
                                ! get boundary j indices for particle
                                wall_j = max(min(int(loc_j_new - v_eta * del_t_i), world%N_y-1), 1)
                                delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i, wall_j+1) == 1
                                if (.not. delete_bool) then
                                    if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i, wall_j+1) == 2) then
                                        ! Neumann boundary
                                        loc_i_new = 2.0d0 - loc_i_new
                                        v_xi = -v_xi
                                        v_part(1) = -v_part(1)
                                    else
                                        loc_i_new = real(world%N_x, kind = 8) - MODULO(real(world%N_x, kind = 8) - loc_i_new, real(world%N_x-1, kind = 8))
                                    end if
                                end if
                            end if
                        end do

                        do while ((loc_j_new > world%N_y .or. loc_j_new < 1) .and. .not. delete_bool)
                            ! take care of any issue with particle outside eta boundary
                            ! it is possible for particle to be outside both x and y, so proceed if not deleted at left-right boundary
                            if (loc_j_new > world%N_y) then
                                ! backtrack to wall position where it left
                                wall_j = world%N_y
                                del_t_j = (loc_j_new - real(wall_j, kind = 8))/v_eta
            
                                ! get boundary i indices for particle
                                wall_i = max(min(int(loc_i_new - v_xi * del_t_j), world%N_x-1), 1)
                                delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i+1, wall_j) == 1
                                if (.not. delete_bool) then
                                    if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i+1, wall_j) == 2) then
                                        ! Neumann boundary
                                        loc_j_new = 2.0d0 * real(wall_j, kind = 8) - loc_j_new
                                        v_eta = -v_eta
                                        v_part(2) = -v_part(2)
                                    else
                                        loc_j_new = MODULO(loc_j_new, real(world%N_y-1, kind = 8))
                                    end if
                                end if
                            else
                                ! backtrack to wall position where it left
                                wall_j = 1
                                del_t_j = (loc_j_new - real(wall_j, kind = 8))/v_eta
            
                                ! get boundary i indices for particle
                                wall_i = max(min(int(loc_i_new - v_xi * del_t_j), world%N_x-1), 1)
                                delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i+1, wall_j) == 1
                                if (.not. delete_bool) then
                                    if (world%boundary_conditions(wall_i, wall_j) == 2 .or. world%boundary_conditions(wall_i+1, wall_j) == 2) then
                                        ! Neumann boundary
                                        loc_j_new = 2.0d0- loc_j_new
                                        v_eta = -v_eta
                                        v_part(2) = -v_part(2)
                                    else
                                        loc_j_new = real(world%N_y, kind = 8) - MODULO(real(world%N_y, kind = 8) - loc_j_new, real(world%N_y-1, kind = 8))
                                    end if
                                end if
                            end if

                        end do

                        if (.not. delete_bool) then

                            wall_i = int(loc_i_new) + (1 - INT(SIGN(1.0d0, v_xi)))/2
                            wall_j = int(loc_j_new) + (1 - INT(SIGN(1.0d0, v_eta)))/2
                    
                            ! We'll assume that particle doesn't go across many cells, so likelihood of passing dirichlet boundary and ending up in cell surrounded by plasma nodes is low
                            ! corner node (wall_i, wall_j) needs to be dirichlet for particle to have chance of passing dirichlet wall
                            if (world%boundary_conditions(wall_i, wall_j) == 1) then
                                ! Find time to each wall it could have passed through
                                del_t_i = (loc_i_new - wall_i)/v_xi
                                del_t_j = (loc_j_new - wall_j)/v_eta

                                if (del_t_i < del_t_j) then
                                    ! hits along left-right wall, check other j-th index node
                                    wall_j = int(loc_j_new)
                                    delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i, wall_j+1) == 1
                                else
                                    ! hits along up-down wall, check other i-th index node
                                    wall_i = int(loc_i_new)
                                    delete_bool = world%boundary_conditions(wall_i, wall_j) == 1 .and. world%boundary_conditions(wall_i+1, wall_j) == 1
                                end if
                            end if

                            if (.not. delete_bool) then
                                wall_i = int(loc_i_new)
                                wall_j = int(loc_j_new)
                                if ((wall_j < j_cell) .or. (wall_i < i_cell .and. wall_j == j_cell)) then
                                    ! if in previous cell can add directly to cell array and increase number of particle counter
                                    self%logical_position(1, self%cell_starting_indx(wall_i, wall_j) + self%number_particles_cell_thread(wall_i, wall_j, i_thread), i_thread) = loc_i_new
                                    self%logical_position(2, self%cell_starting_indx(wall_i, wall_j) + self%number_particles_cell_thread(wall_i, wall_j, i_thread), i_thread) = loc_j_new
                                    self%velocity(1:2, self%cell_starting_indx(wall_i, wall_j) + self%number_particles_cell_thread(wall_i, wall_j, i_thread), i_thread) = v_part
                                    self%velocity(3, self%cell_starting_indx(wall_i, wall_j) + self%number_particles_cell_thread(wall_i, wall_j, i_thread), i_thread) = v_z
                                    self%number_particles_cell_thread(wall_i, wall_j, i_thread) = self%number_particles_cell_thread(wall_i, wall_j, i_thread) + 1
                                else
                                    number_particles_outside_cell = self%number_particles_cell_thread(wall_i, wall_j, i_thread)
                                    cell_end_indx = self%cell_ending_indx(wall_i, wall_j)
                                    if (cell_end_indx - number_particles_added_cell_thread(wall_i, wall_j, i_thread) > &
                                        self%cell_starting_indx(wall_i, wall_j) + number_particles_outside_cell - 1) then
                                        ! Space within the cell array to place particles at end of the array
                                        self%logical_position(1, cell_end_indx - number_particles_added_cell_thread(wall_i, wall_j, i_thread), i_thread) = loc_i_new
                                        self%logical_position(2, cell_end_indx - number_particles_added_cell_thread(wall_i, wall_j, i_thread), i_thread) = loc_j_new
                                        self%velocity(1:2, cell_end_indx - number_particles_added_cell_thread(wall_i, wall_j, i_thread), i_thread) = v_part
                                        self%velocity(3, cell_end_indx - number_particles_added_cell_thread(wall_i, wall_j, i_thread), i_thread) = v_z
                                        number_particles_added_cell_thread(wall_i, wall_j, i_thread) = number_particles_added_cell_thread(wall_i, wall_j, i_thread) + 1
                                    else   
                                        ! no space left, place in overflow array
                                        number_particles_overflow = number_particles_overflow + 1
                                        logical_position_overflow(1, number_particles_overflow, i_thread) = loc_i_new
                                        logical_position_overflow(2, number_particles_overflow, i_thread) = loc_j_new
                                        velocity_overflow(1:2, number_particles_overflow, i_thread) = v_part
                                        velocity_overflow(3, number_particles_overflow, i_thread) = v_z
                                    end if
                                end if
                            end if

                        end if
                    end if

                end do  
                
                
                cell_end_indx = self%cell_ending_indx(i_cell, j_cell)
                ! go through any particles that have been added by earlier cells
                do part_num = cell_end_indx - number_particles_added_cell_thread(i_cell, j_cell, i_thread) + 1, cell_end_indx
                    ! place particle to end of current list
                    self%logical_position(:, cell_start_indx + number_particles_cell, i_thread) = self%logical_position(:, part_num, i_thread)
                    self%velocity(:, cell_start_indx + number_particles_cell, i_thread) = self%velocity(:, part_num, i_thread)
                    number_particles_cell = number_particles_cell + 1
                end do
                self%number_particles_cell_thread(i_cell, j_cell, i_thread) = number_particles_cell
                number_particles_added_cell_thread(i_cell,j_cell, i_thread) = 0
            end do
        end do
        
        
        ! ! put remaining particles in overflow array in proper cell
        do part_num = 1, number_particles_overflow
            loc_i = logical_position_overflow(1, part_num, i_thread)
            loc_j = logical_position_overflow(2, part_num, i_thread)
            wall_i = int(loc_i)
            wall_j = int(loc_j)
            cell_start_indx = self%cell_starting_indx(wall_i, wall_j)
            number_particles_cell = self%number_particles_cell_thread(wall_i, wall_j, i_thread)
            self%logical_position(1, cell_start_indx + number_particles_cell, i_thread) = loc_i
            self%logical_position(2, cell_start_indx + number_particles_cell, i_thread) = loc_j
            self%velocity(:, cell_start_indx + number_particles_cell, i_thread) = velocity_overflow(:, part_num, i_thread)
            self%number_particles_cell_thread(wall_i, wall_j, i_thread) = self%number_particles_cell_thread(wall_i, wall_j, i_thread) + 1
        end do
        
        number_particles_overflow_thread(i_thread) = number_particles_overflow
    end subroutine particle_mover_uniform_cell



    subroutine initialize_maxwellian_temperature_cell(self, T)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable
        class(Particle_Per_Cell), intent(in out) :: self
        real(real64), intent(in) :: T
        integer(int32) :: iThread, i_cell, j_cell, N_x_cell, N_y_cell
        integer(int64) :: number_particles_cell, part_num
        real(real64) :: U1, U2, U3, U4, v_therm
        v_therm = SQRT(T*e_charge/self%mass)
        N_x_cell = size(self%cell_starting_indx, DIM = 1)
        N_y_cell = size(self%cell_starting_indx, DIM = 2)
        !$OMP PARALLEL PRIVATE(iThread, part_num, number_particles_cell, U1, U2, U3, U4, j_cell, i_cell)
        iThread = omp_get_thread_num() + 1
        do j_cell = 1, N_y_cell
            do i_cell = 1, N_x_cell
                do part_num = self%cell_starting_indx(i_cell, j_cell), self%cell_starting_indx(i_cell, j_cell) + self%number_particles_cell_thread(i_cell, j_cell, iThread) - 1
                    U1 = pcg32_random_r(state_PCG)
                    U2 = pcg32_random_r(state_PCG)
                    U3 = pcg32_random_r(state_PCG)
                    U4 = pcg32_random_r(state_PCG)
                    self%velocity(1, part_num, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * COS(2.0d0 * pi_const * U2)
                    self%velocity(2, part_num, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * SIN(2.0d0 * pi_const * U2)
                    self%velocity(3, part_num, iThread) = v_therm * SQRT(-2.0d0 * LOG(U3)) * SIN(2.0d0 * pi_const * U4)
                end do
            end do
        end do
        !$OMP END PARALLEL
    end subroutine initialize_maxwellian_temperature_cell

    subroutine get_sum_totals_cell(self)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle_Per_Cell), intent(in out) :: self
        integer(int32) :: i_thread, i_cell, j_cell, N_x_cell, N_y_cell
        integer(int64) :: part_num
        real(real64) :: sum_v, sum_v_sqr, v(3)
        sum_v = 0.0d0
        sum_v_sqr = 0.0d0
        
        self%total_number_particles = SUM(self%number_particles_cell_thread)
        N_x_cell = size(self%cell_starting_indx, DIM = 1)
        N_y_cell = size(self%cell_starting_indx, DIM = 2)
        !$OMP parallel private(i_thread, j_cell, i_cell, part_num, v) reduction(+:sum_v, sum_v_sqr)
        i_thread = omp_get_thread_num() + 1
        do j_cell = 1, N_y_cell
            do i_cell = 1, N_x_cell
                do part_num = self%cell_starting_indx(i_cell, j_cell), self%cell_starting_indx(i_cell, j_cell) + self%number_particles_cell_thread(i_cell, j_cell, i_thread) - 1
                    v = self%velocity(:, part_num, i_thread)
                    sum_v = sum_v + v(1)
                    sum_v_sqr = sum_v_sqr + sum(v**2)
                end do
            end do
        end do
        !$OMP end parallel
        self%sum_v = sum_v
        self%sum_v_sqr = sum_v_sqr
    end subroutine get_sum_totals_cell


    ! subroutine resize_particle_arrays(particle_list)
    !     type(Particle), intent(in out) :: self
    !     integer(int32) :: i_thread
    !     integer(int64) :: max_indx_new
    !     real(real64), allocatable :: logical_position_copy(:,:,:), velocity_copy(:,:,:)
        
    !     max_indx_new = real(maxval(self%number_particles_thread))
    !     if (real(max_indx_new, kind = 8) < real(self%max_indx, kind = 8) * 0.8d0 .or. real(max_indx_new, kind = 8) > real(self%max_indx, kind = 8) * 0.95d0) then
    !         ! Try to keep max indx for particles within 20%
    !         self%max_indx = int(real(max_indx_new, kind = 8) / 0.875d0)
    !         ! copy arrays to allocated array
    !         allocate(logical_position_copy(2, self%max_indx, number_threads_global), velocity_copy(3, self%max_indx, number_threads_global))
    !         !$OMP parallel private(i_thread)
    !         i_thread = omp_get_thread_num() + 1
    !         logical_position_copy(:,1:self%number_particles_thread(i_thread), i_thread) = self%logical_position(:, 1:self%number_particles_thread(i_thread), i_thread)
    !         velocity_copy(:, 1:self%number_particles_thread(i_thread), i_thread) = self%velocity(:, 1:self%number_particles_thread(i_thread), i_thread)
    !         !$OMP end parallel

    !         deallocate(self%logical_position, self%velocity)
    !         call move_alloc(logical_position_copy, self%logical_position)
    !         call move_alloc(velocity_copy, self%velocity)
    !         ! allocate(self%logical_position(2, self%max_indx, number_threads_global), self%velocity(3, self%max_indx, number_threads_global))
    !         ! !$OMP parallel private(i_thread)
    !         ! i_thread = omp_get_thread_num() + 1
    !         ! self%logical_position(:,1:self%number_particles_thread(i_thread), i_thread) = logical_position_copy(:, 1:self%number_particles_thread(i_thread), i_thread)
    !         ! self%velocity(:, 1:self%number_particles_thread(i_thread), i_thread) = velocity_copy(:, 1:self%number_particles_thread(i_thread), i_thread)
    !         ! !$OMP end parallel
    !         ! deallocate(logical_position_copy, velocity_copy)
    !         if (allocated(logical_position_copy)) then
    !             print *, 'error with allocation'

    !         else if (allocated(velocity_copy)) then
    !             print *, 'error with allocation'
    !         end if

    !     end if

    ! end subroutine

    ! function getTotalMomentum(self) result(res)
    !     ! Get total momentum in domain
    !     class(Particle), intent(in) :: self
    !     real(real64) :: res(3), temp(3)
    !     integer(int32) :: iThread
    !     temp = 0
    !     !$OMP parallel private(iThread) reduction(+:temp)
    !     iThread = omp_get_thread_num() + 1
    !     temp = temp +  SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread), DIM = 2)
    !     !$OMP end parallel
    !     res = temp * self%w_p * self%mass
    ! end function getTotalMomentum

    ! function getTotalKE(self) result(res)
    !     ! calculate total KE in Joules/m^2
    !     class(Particle), intent(in) :: self
    !     real(real64) :: res
    !     integer(int32) :: iThread
    !     res = 0.0d0
    !     !$OMP parallel private(iThread) reduction(+:res)
    !     iThread = omp_get_thread_num() + 1
    !     res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread)**2) 
    !     !$OMP end parallel
    !     res = res * self % mass * 0.5d0 * self%w_p

    ! end function getTotalKE
    

    ! --------------------------- Writing Particle Data to File -----------------------------------

    ! subroutine writePhaseSpace(self, dirName)
    !     ! Writes particle phase space into binary file
    !     class(Particle), intent(in) :: self
    !     character(*), intent(in) :: dirName
    !     character(len=5) :: char_i
    !     integer(int32) :: iThread
    !     !write(char_i, '(I3)'), CurrentDiagStep
    !     !$OMP parallel private(iThread, char_i)
    !     iThread = omp_get_thread_num() + 1
    !     write(char_i, '(I3)'), iThread
    !     open(iThread,file=dirName//'/PhaseSpace/phaseSpace_'//self%name//"_thread"//trim(adjustl(char_i))//".dat", form='UNFORMATTED', access = 'STREAM', status = 'REPLACE')
    !     write(iThread) self%phaseSpace(:, 1:self%N_p(iThread), iThread)
    !     close(ithread)
    !     !$OMP end parallel
    ! end subroutine writePhaseSpace

    ! subroutine writeLocalTemperature(self, CurrentDiagStep, dirName, numCells)
    !     ! Write particle temperature averaged over local grid
    !     class(Particle), intent(in) :: self
    !     integer(int32), intent(in) :: CurrentDiagStep, numCells
    !     character(*), intent(in) :: dirName
    !     character(len=5) :: char_i
    !     integer(int32) :: j, index, counter(numCells, numThread), iThread
    !     real(real64) :: temp(numCells, numThread), EHist(numCells)
    !     temp = 0.0d0
    !     counter = 0
    !     !$OMP parallel private(iThread, j, index) 
    !     iThread = omp_get_thread_num() + 1
    !     do j = 1, self%N_p(iThread)
    !         index = INT(self%phaseSpace(1, j, iThread))
    !         if (index < NumberXNodes) then
    !             temp(index, iThread) = temp(index, iThread) + SUM(self%phaseSpace(2:4, j, iThread)**2)
    !             counter(index, iThread) = counter(index, iThread) + 1
    !         end if
    !     end do
    !     !$OMP end parallel
    !     do j = 1, numCells
    !         if (SUM(counter(j, :)) > 0) then
    !             EHist(j) = SUM(temp(j,:))*self%mass/SUM(counter(j, :))/3.0d0/e
    !         else
    !             EHist(j) = 0.0d0
    !         end if
    !     end do
    !     write(char_i, '(I3)'), CurrentDiagStep
    !     open(10,file=dirName//'/Temperature/Temp_'//self%name//"_"//trim(adjustl(char_i))//".dat", form='UNFORMATTED')
    !     write(10) EHist
    !     close(10)
        
    ! end subroutine writeLocalTemperature


    ! ------------------- read and initialize particles using world type ---------------------------------

    ! subroutine readChargedParticleInputs(filename, random_gen, T_e, T_i, n_ave, num_threads, world, particle_list)
    !     ! Read input file for particles
    !     type(Particle), allocatable, intent(out) :: particle_list(:)
    !     class(world), intent(in) :: world
    !     character(len=*), intent(in) :: filename
    !     integer(int32), intent(in) :: num_threads
    !     type(rand_gen), intent(in out) :: random_gen(num_threads)
    !     real(real64), intent(in) :: T_e, T_i, n_ave
    !     integer(int32) :: j, io, num_species = 0, num_particles(100), particle_idx_factor(100), i, tempInt, dist_type(100), i_thread, charge(100), i_thread
    !     character(len=15) :: name
    !     character(len=8) :: particleNames(100), char_i
    !     real(real64) :: mass(100), Ti(100), tempReal, alpha(100), v_drift(100)
    !     logical :: electrons_found

    !     print *, "Reading particle inputs:"
    !     open(10,file=filename, action = 'read')
    !     ! if (.not. restartBool) then
    !     !     open(10,file=filename, action = 'read')
    !     ! else
    !     !     open(10,file=restartDirectory//'/InputData/'//filename, action = 'read')
    !     ! end if
    !     ! first look for electrons
    !     electrons_found = .false.
    !     do j=1, 10000
    !         read(10,*) name

    !         if( name(1:9).eq.'ELECTRONS' .or. name(1:9).eq.'electrons' .or. name(1:9).eq.'electrons') then
    !             electrons_found = .true.
    !             read(10,*) name
    !             read(10,*) name
    !             read(10,'(A4)', ADVANCE = 'NO') name(1:2)
    !             numSpecies = numSpecies + 1
    !             read(10,*) numParticles(numSpecies), particleIdxFactor(numSpecies), alpha(numSpecies), distType(numSpecies), v_drift(numSpecies)
    !             Ti(numSpecies) = T_e
    !             mass(numSpecies) = mass_electron
    !             charge(numSpecies) = -1
    !             particleNames(numSpecies) = 'e'
    !             read(10,*) name
    !             read(10,*) name
    !             exit
    !         endif
   

    !         if (name(1:7) == 'ENDFILE') then
    !             exit
    !         end if

    !     end do

    !     if (.not. electrons_found) then
    !         print *, "WARNING: Electrons not found in particle list!"
    !     end if

    !     rewind(10)

    !     do j=1, 10000
    !         read(10,*) name


    !         if(name(1:4).eq.'IONS' .or. name(1:4).eq.'Ions' .or. name(1:4).eq.'ions' ) then
    !             do while(name(1:4).ne.'----')
    !                 read(10,*) name
    !             end do
    !             read(10,'(A6)', ADVANCE = 'NO') name
    !             do while (name(1:4).ne.'----')
    !                 numSpecies = numSpecies + 1
    !                 read(10,*) mass(numSpecies),charge(numSpecies), numParticles(numSpecies), particleIdxFactor(numSpecies), alpha(numSpecies), distType(numSpecies), v_drift(numSpecies)
    !                 Ti(numSpecies) = T_i
    !                 mass(numSpecies) = mass(numSpecies) * mass_amu - charge(numSpecies) * mass_electron
    !                 particleNames(numSpecies) = trim(name)
    !                 read(10,'(A6)', ADVANCE = 'NO') name
    !             end do
    !         endif      

    !         if (name(1:7) == 'ENDFILE') then
    !             exit
    !         end if

    !     end do
    !     close(10)

    !     number_charged_particles = numSpecies
    !     print *, 'Amount charged particles:', number_charged_particles
    !     if (number_charged_particles > 0) then
    !         ! Initialize and generate particles
    !         allocate(particleList(number_charged_particles))
    !         do j=1, number_charged_particles
    !             particle_list(j) = Particle(mass(j), e_charge * charge(j), 1.0d0, numParticles(j), numParticles(j) * particle_idx_factor(j), trim(particleNames(j)), num_threads)
    !             call particle_list(j) % initialize_weight_from_n_ave(n_ave, world)
    !             ! if (.not. restartBool) then
    !             call particleList(j)% initializeRandUniform(world, irand)
    !             ! if (j==2 .and. numberChargedParticles == 2 .and. charge(2) == -charge(1) .and. numParticles(1) == numParticles(2)) then
    !             !     ! If only ions and electrons (electrons come first) then set ion positions same as electrons for neutral start
    !             !     print *, 'Neutral charge start!'
    !             !     particleList(2)%phaseSpace(1, :, :) = particleList(1)%phaseSpace(1,:,:)
    !             ! else
    !             !     SELECT CASE(distType(j))
    !             !     CASE(0)
    !             !         call particleList(j)% initializeRandUniform(world, irand)
    !             !     CASE(1)
    !             !         call particleList(j)% initializeRandCosine(world, irand, alpha(j))
    !             !     CASE(2)
    !             !         call particleList(j)% initializeRandSine(world, irand, alpha(j))
    !             !     CASE default
    !             !         print *, 'Distribution type should be between 0 and 2!'
    !             !         stop
    !             !     END SELECT
    !             ! end if
    !             if (j == 1) then
    !                 call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, alpha(j), distType(j), v_drift(j))
    !             else
    !                 call particleList(j) % generate3DMaxwellian(Ti(j), world, irand, 1.236d0, distType(j), v_drift(j)) ! Used for IASW case, will likely change to general later
    !             end if
    !             ! else
    !             ! !$OMP parallel private(i_thread, boolVal, char_i, io, i)
    !             ! i_thread = omp_get_thread_num() + 1
    !             ! write(char_i, '(I3)'), iThread
    !             ! INQUIRE(file=restartDirectory//'/PhaseSpace/phaseSpace_'//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", exist = boolVal)
    !             ! if (.not. boolVal) then
    !             !     print *, restartDirectory//'/PhaseSpace/phaseSpace_'//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", 'Does not exist'
    !             !     stop
    !             ! end if
    !             ! open(iThread,file=restartDirectory//"/PhaseSpace/phaseSpace_"//particleList(j)%name//"_thread"//trim(adjustl(char_i))//".dat", form = 'UNFORMATTED', access = 'stream', status = 'old', IOSTAT=io)
    !             ! i = 0
    !             ! read(iThread,  IOSTAT = io) particleList(j)%phaseSpace(:, i+1, iThread)
    !             ! do while (io == 0)
    !             !     i = i + 1
    !             !     read(iThread,  IOSTAT = io) particleList(j)%phaseSpace(:, i+1, iThread)
    !             ! end do
    !             ! particleList(j)%N_p(iThread) = i
    !             ! close(iThread)
    !             !$OMP end parallel
    !             ! end if
    !             print *, 'Initializing ', particleList(j) % name
    !             print *, 'Amount of macroparticles is:', SUM(particleList(j) % N_p)
    !             print *, "Particle mass is:", particleList(j)%mass
    !             print *, "Particle charge is:", particleList(j)%q
    !             print *, "Particle weight is:", particleList(j)%w_p
    !             print *, "Particle mean KE is:", particleList(j)%getKEAve()
    !             print *, 'Distribution type:', distType(j)
    !             print *, 'Drift velocity:', v_drift(j)
    !         end do

    !     else
    !         print *, "ERROR: No particles input!"
    !         stop
    !     end if

        
    !     print *, "---------------"
    !     print *, ""


    ! end subroutine readChargedParticleInputs

end module mod_particle_per_cell