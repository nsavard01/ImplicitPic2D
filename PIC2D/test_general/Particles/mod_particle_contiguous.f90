module mod_particle_contiguous

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_particle
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_rand_generator
    use omp_lib
    implicit none

   

    ! Particle contains particle properties and stored values in phase space
    type, extends(Particle) :: Particle_Contiguous
        integer(int64), allocatable :: number_particles_thread(:)
        integer(int64), allocatable :: cell_indx_array(:,:,:)
        integer(int32), allocatable ::  cell_count(:,:,:)
        logical :: count_bool

    contains
        procedure, public, pass(self) :: initialize_rand_uniform => initialize_rand_uniform_contiguous
        procedure, public, pass(self) :: interpolation_particle_to_nodes => interpolation_particle_to_nodes_contiguous
        ! procedure, public, pass(self) :: interpolation_particle_to_nodes_and_sort
        ! procedure, public, pass(self) :: interpolation_particle_to_nodes_sorted
        procedure, public, pass(self) :: particle_sort
        procedure, public, pass(self) :: particle_mover_uniform => particle_mover_uniform_contiguous
        ! procedure, public, pass(self) :: initializeRandCosine
        ! procedure, public, pass(self) :: initializeRandSine
        procedure, public, pass(self) :: initialize_maxwellian_temperature => initialize_maxwellian_temperature_contiguous
        procedure, public, pass(self) :: get_sum_totals => get_sum_totals_contiguous
        ! procedure, public, pass(self) :: getTotalKE
        ! procedure, public, pass(self) :: getTotalMomentum
        ! procedure, public, pass(self) :: writePhaseSpace
        ! procedure, public, pass(self) :: writeLocalTemperature
    end type


    interface Particle_Contiguous
        module procedure :: particle_continguous_constructor
    end interface Particle_Contiguous

contains

    type(Particle_Contiguous) function particle_continguous_constructor(mass, q, w_p, N_p, finalIdx, particleName, world) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, w_p
        integer(int64), intent(in) :: N_p, finalIdx
        class(domain_base), intent(in) :: world
        character(*), intent(in) :: particleName
        call self%initialize_base_variables(mass, q, w_p, N_p, finalIdx, particleName, world%N_x, world%N_y)
        allocate(self%number_particles_thread(number_threads_global), self%cell_count(world%N_x-1, world%N_y-1, number_threads_global), &
            self%cell_indx_array(world%N_x-1, world%N_y-1, number_threads_global))
        self%count_bool = .false.
        self%number_particles_thread = self%total_number_particles/number_threads_global
       
    end function particle_continguous_constructor


    subroutine initialize_rand_uniform_contiguous(self, world)
        ! distribute particles randomly over the domain
        class(Particle_Contiguous), intent(in out) :: self
        class(domain_base), intent(in) :: world
        integer(int32) :: i_thread, int_xi, int_eta
        integer(int64) :: i
        real(real64) :: x_pos, y_pos, L_x, L_y, xi, eta
        L_x = world%end_X - world%start_X
        L_y = world%end_Y - world%start_Y
        !$OMP parallel private(i_thread, i, x_pos, y_pos, eta, xi, int_xi, int_eta)
        i_thread = omp_get_thread_num() + 1
        self%cell_count(:,:, i_thread) = 0
        do i = 1, self%number_particles_thread(i_thread)
            x_pos = pcg32_random_r(state_PCG) * L_x + world%start_X
            xi = world%get_xi_from_X(x_pos)
            int_xi = int(xi)

            y_pos = pcg32_random_r(state_PCG)  * L_y + world%start_Y
            eta = world%get_eta_from_Y(y_pos)
            int_eta = int(eta)
            do while (world%boundary_conditions(int_xi, int_eta) /= 0 .and. world%boundary_conditions(int_xi+1, int_eta) /= 0 &
                .and. world%boundary_conditions(int_xi, int_eta+1) /= 0 .and. world%boundary_conditions(int_xi+1, int_eta+1) /= 0)
                x_pos = pcg32_random_r(state_PCG)  * L_x + world%start_X
                xi = world%get_xi_from_X(x_pos)
                int_xi = int(xi)

                y_pos = pcg32_random_r(state_PCG) * L_y + world%start_Y
                eta = world%get_eta_from_Y(y_pos)
                int_eta = int(eta)
            end do
            
            self%logical_position(1,i, i_thread) = xi
            self%logical_position(2, i, i_thread) = eta
            self%cell_count(int_xi, int_eta, i_thread) = self%cell_count(int_xi, int_eta, i_thread) + 1
        end do
        !$OMP end parallel
    end subroutine initialize_rand_uniform_contiguous

    subroutine interpolation_particle_to_nodes_contiguous(self, i_thread, const)
        ! interpolate particles to work space array
        class(Particle_Contiguous), intent(in) :: self
        integer(int32), intent(in) :: i_thread
        real(real64), intent(in) :: const
        integer(int32) :: i_cell, j_cell
        integer(int64) :: part_num, number_particles
        real(real64) :: d_i, d_j, xi, eta

        number_particles = self%number_particles_thread(i_thread)
        do part_num = 1, number_particles
            xi = self%logical_position(1,part_num,i_thread)
            eta = self%logical_position(2,part_num,i_thread)
            i_cell = int(xi)
            j_cell = int(eta)
            d_i = xi - real(i_cell, kind = 8)
            d_j = eta - real(j_cell, kind = 8)

            particle_work_space(i_cell,j_cell, i_thread) = particle_work_space(i_cell,j_cell, i_thread) + (1.0d0-d_i) * (1.0d0-d_j) * const
            particle_work_space(i_cell+1,j_cell, i_thread) = particle_work_space(i_cell+1,j_cell, i_thread) + (d_i) * (1.0d0-d_j) * const
            particle_work_space(i_cell,j_cell+1, i_thread) = particle_work_space(i_cell,j_cell+1, i_thread) + (1.0d0-d_i) * (d_j) * const
            particle_work_space(i_cell+1,j_cell+1, i_thread) = particle_work_space(i_cell+1,j_cell+1, i_thread) + (d_i) * (d_j) * const
        end do

    end subroutine interpolation_particle_to_nodes_contiguous



    subroutine particle_sort(self, i_thread, N_x_cell, N_y_cell)
        ! sort particle by cell, with each cell going j = 1-> N_y-1, i = 1->N_x-1
        ! make sort in place so no need to 
        class(Particle_Contiguous), intent(in out) :: self
        integer(int32), intent(in) :: N_x_cell, N_y_cell
        integer(int64) :: part_num, cell_end_indx
        integer(int32) :: i_thread, eta, xi, i_cell, j_cell, cell_count
        real(real64) :: pos_curr(2), v_curr(3), pos_other(2), v_other(3)

        

        ! get final cell index of each bin
        self%cell_indx_array(1, 1, i_thread) = self%cell_count(1,1,i_thread)
        do j_cell = 1, N_y_cell-1
            do i_cell = 2, N_x_cell
                self%cell_indx_array(i_cell, j_cell, i_thread) = self%cell_count(i_cell, j_cell, i_thread) + self%cell_indx_array(i_cell-1, j_cell, i_thread)
            end do
            self%cell_indx_array(1, j_cell+1, i_thread) = self%cell_count(1, j_cell+1, i_thread) + self%cell_indx_array(N_x_cell, j_cell, i_thread)
        end do
        ! do last row
        j_cell = N_y_cell
        do i_cell = 2, N_x_cell
            self%cell_indx_array(i_cell, j_cell, i_thread) = self%cell_count(i_cell, j_cell, i_thread) + self%cell_indx_array(i_cell-1, j_cell, i_thread)
        end do
        
        ! order particles
        do j_cell = N_y_cell, 2, -1
            do i_cell = N_x_cell, 1, -1
                cell_end_indx = self%cell_indx_array(i_cell, j_cell, i_thread)
                cell_count = self%cell_count(i_cell, j_cell, i_thread)
                do part_num = cell_end_indx, cell_end_indx-cell_count+1, -1
                    pos_curr = self%logical_position(:,part_num, i_thread)
                    v_curr = self%velocity(:,part_num, i_thread)
                    xi = int(pos_curr(1))
                    eta = int(pos_curr(2))
                    do while (xi /= i_cell .or. eta /= j_cell)
                        pos_other = self%logical_position(:,self%cell_indx_array(xi, eta,i_thread), i_thread)
                        v_other = self%velocity(:,self%cell_indx_array(xi, eta,i_thread), i_thread)
                        ! put last index in current place and then reduce that section by 1
                        self%logical_position(:,self%cell_indx_array(xi,eta,i_thread), i_thread) = pos_curr
                        self%velocity(:,self%cell_indx_array(xi,eta,i_thread), i_thread) = v_curr
                        self%cell_indx_array(xi,eta,i_thread) = self%cell_indx_array(xi,eta,i_thread)-1
                        self%cell_count(xi,eta, i_thread) = self%cell_count(xi,eta,i_thread) - 1
                        pos_curr = pos_other
                        v_curr = v_other
                        xi = int(pos_curr(1))
                        eta = int(pos_curr(2))
                    end do
                    self%logical_position(:,part_num, i_thread) = pos_curr
                    self%velocity(:, part_num, i_thread) = v_curr
                end do
                self%cell_count(i_cell, j_cell, i_thread) = 0
                self%cell_indx_array(i_cell, j_cell, i_thread) = cell_end_indx-cell_count+1
            end do
        end do
        j_cell = 1
        do i_cell = N_x_cell, 2, -1
            cell_end_indx = self%cell_indx_array(i_cell, j_cell, i_thread)
            cell_count = self%cell_count(i_cell, j_cell, i_thread)
            do part_num = cell_end_indx, cell_end_indx-cell_count+1, -1
                pos_curr = self%logical_position(:,part_num, i_thread)
                v_curr = self%velocity(:,part_num, i_thread)
                xi = int(pos_curr(1))
                eta = int(pos_curr(2))
                do while (xi /= i_cell .or. eta /= j_cell)
                    pos_other = self%logical_position(:,self%cell_indx_array(xi, eta,i_thread), i_thread)
                    v_other = self%velocity(:,self%cell_indx_array(xi, eta,i_thread), i_thread)
                    ! put last index in current place and then reduce that section by 1
                    self%logical_position(:,self%cell_indx_array(xi,eta,i_thread), i_thread) = pos_curr
                    self%velocity(:,self%cell_indx_array(xi,eta,i_thread), i_thread) = v_curr
                    self%cell_indx_array(xi,eta,i_thread) = self%cell_indx_array(xi,eta,i_thread)-1
                    self%cell_count(xi,eta, i_thread) = self%cell_count(xi,eta,i_thread) - 1
                    pos_curr = pos_other
                    v_curr = v_other
                    xi = int(pos_curr(1))
                    eta = int(pos_curr(2))
                end do
                self%logical_position(:,part_num, i_thread) = pos_curr
                self%velocity(:, part_num, i_thread) = v_curr
            end do
            self%cell_count(i_cell, j_cell, i_thread) = 0
            self%cell_indx_array(i_cell, j_cell, i_thread) = cell_end_indx-cell_count+1
        end do
        self%cell_count(1,1,i_thread) = 0
        self%cell_indx_array(1,1,i_thread) = 1


    end subroutine particle_sort



    subroutine particle_mover_uniform_contiguous(self, i_thread, E_Field, world, del_t)
        class(Particle_Contiguous), intent(in out) :: self
        type(domain_uniform), intent(in) :: world
        real(real64), intent(in) :: E_Field(2,world%N_x,world%N_y), del_t
        integer(int32), intent(in) :: i_thread
        real(real64) :: v_part(2), loc_i, loc_j, d_i, d_j, E_part(2),&
            E_SE(2), E_SW(2), E_NW(2), E_NE(2), inv_del_x, inv_del_y, &
            v_xi, v_eta, del_t_i, del_t_j, loc_i_new, loc_j_new
        integer(int32) :: corner_i, corner_j, wall_i, wall_j
        integer(int64) :: part_num, delete_idx, number_particles
        logical :: delete_bool

        inv_del_x = 1.0d0/world%del_x
        inv_del_y = 1.0d0/world%del_y
        delete_idx = 0
        if (self%count_bool) self%cell_count(:,:,i_thread) = 0
        number_particles = self%number_particles_thread(i_thread)
        do part_num = 1, number_particles
            !get particle location and velocity in 2D
            loc_i = self%logical_position(1, part_num, i_thread)
            loc_j = self%logical_position(2, part_num, i_thread)
            v_part = self%velocity(1:2, part_num, i_thread)
            corner_i = int(loc_i)
            corner_j = int(loc_j)
            d_i = loc_i - real(corner_i, kind = 8)
            d_j = loc_j - real(corner_j, kind = 8)

            ! first find E_x and E_y on the nodes
            E_SW = E_Field(:, corner_i, corner_j)
            E_SE = E_Field(:, corner_i+1, corner_j)
            E_NW = E_Field(:, corner_i, corner_j+1)
            E_NE = E_Field(:, corner_i+1, corner_j+1)

            ! interpolate to particle position
            E_part = E_SW * (1.0d0 - d_i) * (1.0d0 - d_j) + E_SE * (d_i) * (1.0d0-d_j) + &
                E_NW * (1.0d0-d_i) * (d_j) + E_NE * (d_i) * (d_j)
            
            ! solve for new velocity and position
            v_part = v_part + self%q_over_m * E_part * del_t

            ! place on boundary if outside of boundary
            v_xi = v_part(1) * inv_del_x
            v_eta = v_part(2) * inv_del_y
            loc_i_new = loc_i + v_xi * del_t
            loc_j_new = loc_j + v_eta * del_t


            delete_bool = .false.

            ! Found in general that doing deletions after in new subroutine is slower, so do in place

            ! don't need to check boundaries if particle still in same cell
            if (corner_i /= int(loc_i_new) .or. corner_j /= int(loc_j_new)) then
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

                end if
            end if


            if (.not. delete_bool) then
                self%logical_position(1, part_num - delete_idx, i_thread) = loc_i_new
                self%logical_position(2, part_num - delete_idx, i_thread) = loc_j_new
                self%velocity(1:2, part_num - delete_idx, i_thread) = v_part
                self%velocity(3, part_num - delete_idx, i_thread) = self%velocity(3, part_num, i_thread)
                if (self%count_bool) self%cell_count(int(loc_i_new), int(loc_j_new), i_thread) = self%cell_count(int(loc_i_new), int(loc_j_new), i_thread) + 1
            else
                delete_idx = delete_idx + 1
            end if

            
            
        end do
        self%number_particles_thread(i_thread) = number_particles - delete_idx

        if (self%count_bool) call self%particle_sort(i_thread, world%N_x-1, world%N_y-1)

    end subroutine particle_mover_uniform_contiguous


    subroutine initialize_maxwellian_temperature_contiguous(self, T)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable
        class(Particle_Contiguous), intent(in out) :: self
        real(real64), intent(in) :: T
        integer(int32) :: iThread, i
        real(real64) :: U1, U2, U3, U4, v_therm
        v_therm = SQRT(T*e_charge/self%mass)
        !$OMP PARALLEL PRIVATE(iThread, i, U1, U2, U3, U4)
        iThread = omp_get_thread_num() + 1
        do i = 1, self%number_particles_thread(iThread)
            U1 = pcg32_random_r(state_PCG)
            U2 = pcg32_random_r(state_PCG)
            U3 = pcg32_random_r(state_PCG)
            U4 = pcg32_random_r(state_PCG)
            self%velocity(1, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * COS(2.0d0 * pi_const * U2)
            self%velocity(2, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * SIN(2.0d0 * pi_const * U2)
            self%velocity(3, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U3)) * SIN(2.0d0 * pi_const * U4)
            
        end do
        !$OMP END PARALLEL
    end subroutine initialize_maxwellian_temperature_contiguous

    subroutine get_sum_totals_contiguous(self)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle_Contiguous), intent(in out) :: self
        integer(int32) :: i_thread
        real(real64) :: sum_v, sum_v_sqr
        sum_v = 0.0d0
        sum_v_sqr = 0.0d0
        
        self%total_number_particles = SUM(self%number_particles_thread)
    
        !$OMP parallel private(i_thread) reduction(+:sum_v, sum_v_sqr)
        i_thread = omp_get_thread_num() + 1
        sum_v = sum_v + SUM(self%velocity(1, 1:self%number_particles_thread(i_thread), i_thread)) 
        sum_v_sqr = sum_v_sqr + SUM(self%velocity(:, 1:self%number_particles_thread(i_thread), i_thread)**2) 
        !$OMP end parallel
        self%sum_v = sum_v
        self%sum_v_sqr = sum_v_sqr
    end subroutine get_sum_totals_contiguous


    ! subroutine resize_particle_arrays(self)
    !     class(Particle), intent(in out) :: self
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

end module mod_particle_contiguous