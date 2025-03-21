module mod_particle_contiguous_decomp


    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_rand_generator
    use omp_lib
    use mod_particle_contiguous
    use mod_particle
    implicit none

    integer(int64), allocatable :: cell_indx_array_total(:,:)
    integer(int32), allocatable :: cell_count_total(:,:)
    real(real64), allocatable :: logical_position_total(:,:), velocity_total(:,:)

    ! Particle contains particle properties and stored values in phase space
    type, extends(Particle_Contiguous) :: Particle_Contiguous_Decomp

    contains
        procedure, public, pass(self) :: particle_sort => particle_sort_decomp
        procedure, public, pass(self) :: initialize_rand_uniform => initialize_rand_uniform_decomp
    end type


    interface Particle_Contiguous_Decomp
        module procedure :: particle_continguous_decomp_constructor
    end interface Particle_Contiguous_Decomp

contains

    type(Particle_Contiguous_Decomp) function particle_continguous_decomp_constructor(mass, q, w_p, N_p, finalIdx, particleName, world) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
        real(real64), intent(in) :: mass, q, w_p
        integer(int64), intent(in) :: N_p, finalIdx
        class(domain_base), intent(in) :: world
        character(*), intent(in) :: particleName
      
        call self%initialize_base_variables(mass, q, w_p, N_p, finalIdx, particleName, world%N_x, world%N_y)
        allocate(self%number_particles_thread(number_threads_global))
        if (.not. allocated(cell_count)) then
            allocate(cell_count(world%N_x-1, world%N_y-1, number_threads_global), &
            cell_indx_array(world%N_x-1, world%N_y-1, number_threads_global))
        end if
        if (.not. allocated(cell_indx_array_total)) then
            allocate(cell_indx_array_total(world%N_x-1, world%N_y-1), cell_count_total(world%N_x-1, world%N_y-1), &
                logical_position_total(2, self%total_number_particles), velocity_total(3, self%total_number_particles))
        end if
        self%count_bool = .false.
        self%number_particles_thread = self%total_number_particles/number_threads_global
       
    end function particle_continguous_decomp_constructor

    subroutine initialize_rand_uniform_decomp(self, world)
        ! distribute particles randomly over the domain
        class(Particle_Contiguous_Decomp), intent(in out) :: self
        class(domain_base), intent(in) :: world
        integer(int32) :: i_thread, int_xi, int_eta, last_xi, last_eta, number_cell
        integer(int64) :: part_num, number_particles, idx_num, start_point, end_point
        real(real64) :: particle_remainder, real_number_cell
        !$OMP parallel private(i_thread, part_num,int_xi, int_eta, idx_num, number_particles, last_xi, last_eta, number_cell, real_number_cell, particle_remainder, start_point, end_point)
        i_thread = omp_get_thread_num() + 1
        number_particles = self%number_particles_thread(i_thread)
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
                cell_count(int_xi, int_eta, i_thread) = number_cell
                do part_num = 1, number_cell
                    idx_num = idx_num + 1
                    self%logical_position(1, idx_num, i_thread) = pcg32_random_r(state_PCG) + real(int_xi, kind = 8)
                    self%logical_position(2, idx_num, i_thread) = pcg32_random_r(state_PCG) + real(int_eta, kind = 8)
                end do
            end do
        end do

        call self%particle_sort(i_thread, world)

        ! !$OMP barrier
        ! ! Get total across all threads
        ! start_point = world%omp_cell_thread_indices(1, i_thread)
        ! end_point = world%omp_cell_thread_indices(2, i_thread)
        ! cell_count_total(start_point:end_point,1) = 0
        ! do idx_num = 1, number_threads_global
        !     cell_count_total(start_point:end_point,1) = cell_count_total(start_point:end_point,1) + cell_count(start_point:end_point, 1, idx_num)
        ! end do

        ! !$OMP barrier
        ! ! Get end cell index for total array
        ! !$OMP single
        ! part_num = 0
        ! do int_eta = 1, world%N_y-1
        !     do int_xi = 1, world%N_x-1
        !         part_num = part_num + cell_count_total(int_xi, int_eta)
        !         cell_indx_array_total(int_xi, int_eta) = part_num
        !         ! wraps around so can use i_cell -1 even when i_cell = 0
        !     end do
        ! end do
        ! self%total_number_particles = part_num
        ! !$OMP end single

        ! !$OMP barrier
        ! do k_cell = start_point, end_point
        !     ! part_start is start of chunk total array
        !     part_start = cell_indx_array_total(k_cell,1) - cell_count_total(k_cell, 1) + 1
        !     ! part_end is last particle in k_cell chunk of total array
        !     do count = 1, number_threads_global
        !         thread_part_start = cell_indx_array(k_cell, 1, count) - cell_count(k_cell, 1, count) + 1
        !         thread_part_end = cell_indx_array(k_cell, 1,count)
        !         part_end = part_start + cell_count(k_cell,1, count) - 1
        !         logical_position_total(:, part_start:part_end) = self%logical_position(:, thread_part_start:thread_part_end, count)
        !         velocity_total(:, part_start:part_end) = self%velocity(:, thread_part_start:thread_part_end, count)
        !         part_start = part_end + 1
        !     end do
        ! end do
        ! !$OMP barrier
        ! ! Get indices to separate into different threads
    
        ! part_end = self%total_number_particles/number_threads_global - 1
        ! cell_end_indx = MOD(self%total_number_particles, number_threads_global)
        ! part_start = cell_end_indx + (cell_end_indx) * (part_end + 1) + 1
        ! if (i_thread <= cell_end_indx) then
        !     end_point = i_thread + (i_thread) * (part_end + 1)
        !     start_point = i_thread + (i_thread-1) * (part_end + 1)
        ! else
        !     start_point = part_start + (i_thread-cell_end_indx -1) * part_end + (i_thread-cell_end_indx -1)
        !     end_point = part_start + (i_thread-cell_end_indx) * part_end + (i_thread-cell_end_indx -1)
        ! end if
      
        
        ! self%number_particles_thread(i_thread) = end_point - start_point + 1
        ! self%logical_position(:, 1:self%number_particles_thread(i_thread), i_thread) = logical_position_total(:, start_point:end_point)
        ! self%velocity(:, 1:self%number_particles_thread(i_thread), i_thread) = velocity_total(:, start_point:end_point)
        !$OMP end parallel
        
    end subroutine initialize_rand_uniform_decomp


    subroutine particle_sort_decomp(self, i_thread, world)
        ! sort particle by cell, with each cell going j = 1-> N_y-1, i = 1->N_x-1
        ! make sort in place so no need to 
        class(Particle_Contiguous_Decomp), intent(in out) :: self
        class(domain_base), intent(in) :: world
        integer(int32), intent(in) :: i_thread
        integer(int64) :: part_num, cell_end_indx, start_point, end_point, total_cell_size, k_start, k_cell, part_start, part_end, thread_part_start, thread_part_end
        integer(int32) :: eta, xi, i_cell, j_cell, count, N_x_cell, N_y_cell
        real(real64) :: pos_curr(2), v_curr(3), pos_other(2), v_other(3)
        
        N_x_cell = world%N_x-1
        N_y_cell = world%N_y-1
        total_cell_size = N_x_cell * N_y_cell
        ! get final cell index of each bin

        part_num = 0
        do j_cell = 1, N_y_cell
            do i_cell = 1, N_x_cell
                part_num = part_num + cell_count(i_cell, j_cell,i_thread)
                cell_indx_array(i_cell, j_cell, i_thread) = part_num
                ! wraps around so can use i_cell -1 even when i_cell = 0
            end do
        end do

        !$OMP barrier
        ! collect cell counts into cell

        start_point = world%omp_cell_thread_indices(1, i_thread)
        end_point = world%omp_cell_thread_indices(2, i_thread)
        cell_count_total(start_point:end_point,1) = 0
        do count = 1, number_threads_global
            cell_count_total(start_point:end_point,1) = cell_count_total(start_point:end_point,1) + cell_count(start_point:end_point, 1, count)
        end do

        !$OMP barrier
        ! Get end cell index for total array
        !$OMP single
        part_num = 0
        do j_cell = 1, N_y_cell
            do i_cell = 1, N_x_cell
                part_num = part_num + cell_count_total(i_cell, j_cell)
                cell_indx_array_total(i_cell, j_cell) = part_num
                ! wraps around so can use i_cell -1 even when i_cell = 0
            end do
        end do
        self%total_number_particles = part_num
        !$OMP end single
    
        ! order particles in place for each cell
        do j_cell = N_y_cell, 1, -1
            do i_cell = N_x_cell, 1, -1
                cell_end_indx = cell_indx_array(i_cell, j_cell, i_thread)
                count = cell_count(i_cell, j_cell, i_thread)
                part_start = cell_end_indx
                part_end = cell_end_indx-count+1
                do part_num = part_start, part_end, -1
                    pos_curr = self%logical_position(:,part_num, i_thread)
                    v_curr = self%velocity(:,part_num, i_thread)
                    xi = int(pos_curr(1))
                    eta = int(pos_curr(2))
                    do while (xi /= i_cell .or. eta /= j_cell)
                        pos_other = self%logical_position(:,cell_indx_array(xi, eta,i_thread), i_thread)
                        v_other = self%velocity(:,cell_indx_array(xi, eta,i_thread), i_thread)
                        ! put last index in current place and then reduce that section by 1
                        self%logical_position(:,cell_indx_array(xi,eta,i_thread), i_thread) = pos_curr
                        self%velocity(:,cell_indx_array(xi,eta,i_thread), i_thread) = v_curr
                        cell_indx_array(xi,eta,i_thread) = cell_indx_array(xi,eta,i_thread)-1
                        cell_count(xi,eta, i_thread) = cell_count(xi,eta,i_thread) - 1
                        pos_curr = pos_other
                        v_curr = v_other
                        xi = int(pos_curr(1))
                        eta = int(pos_curr(2))
                    end do
                    self%logical_position(:,part_num, i_thread) = pos_curr
                    self%velocity(:, part_num, i_thread) = v_curr
                end do
                cell_count(i_cell, j_cell, i_thread) = 0
                cell_indx_array(i_cell, j_cell, i_thread) = cell_end_indx-count+1
            end do
        end do
        ! k_cell = 1
        ! do part_num = 1, self%number_particles_thread(i_thread)
        !     pos_curr = self%logical_position(:,part_num, i_thread)
        !     xi = int(pos_curr(1))
        !     eta = int(pos_curr(2))
        !     if (xi + (eta-1) * N_x_cell < k_cell) then
        !         print *, 'issue order'
        !         stop
        !     end if
        !     k_cell = xi + (eta-1) * N_x_cell
        ! end do

        ! Recalculate k_cell 
        do k_cell = 1, total_cell_size-1
            cell_count(k_cell, 1, i_thread) = cell_indx_array(k_cell+1, 1, i_thread) - cell_indx_array(k_cell, 1, i_thread)
            cell_indx_array(k_cell,1,i_thread) = cell_indx_array(k_cell+1, 1, i_thread)-1
        end do
        cell_count(total_cell_size, 1, i_thread) = self%number_particles_thread(i_thread) + 1 - cell_indx_array(total_cell_size, 1, i_thread)
        cell_indx_array(total_cell_size, 1, i_thread) = self%number_particles_thread(i_thread)

       
        !$OMP barrier
        do k_cell = start_point, end_point
            ! part_start is start of chunk total array
            part_start = cell_indx_array_total(k_cell,1) - cell_count_total(k_cell, 1) + 1
            ! part_end is last particle in k_cell chunk of total array
            do count = 1, number_threads_global
                thread_part_start = cell_indx_array(k_cell, 1, count) - cell_count(k_cell, 1, count) + 1
                thread_part_end = cell_indx_array(k_cell, 1,count)
                part_end = part_start + cell_count(k_cell,1, count) - 1
                logical_position_total(:, part_start:part_end) = self%logical_position(:, thread_part_start:thread_part_end, count)
                velocity_total(:, part_start:part_end) = self%velocity(:, thread_part_start:thread_part_end, count)
                part_start = part_end + 1
            end do
        end do
        !$OMP barrier
        ! Get indices to separate into different threads
    
        part_end = self%total_number_particles/number_threads_global - 1
        cell_end_indx = MOD(self%total_number_particles, number_threads_global)
        part_start = cell_end_indx + (cell_end_indx) * (part_end + 1) + 1
        if (i_thread <= cell_end_indx) then
            end_point = i_thread + (i_thread) * (part_end + 1)
            start_point = i_thread + (i_thread-1) * (part_end + 1)
        else
            start_point = part_start + (i_thread-cell_end_indx -1) * part_end + (i_thread-cell_end_indx -1)
            end_point = part_start + (i_thread-cell_end_indx) * part_end + (i_thread-cell_end_indx -1)
        end if
      
        
        self%number_particles_thread(i_thread) = end_point - start_point + 1
        self%logical_position(:, 1:self%number_particles_thread(i_thread), i_thread) = logical_position_total(:, start_point:end_point)
        self%velocity(:, 1:self%number_particles_thread(i_thread), i_thread) = velocity_total(:, start_point:end_point)

    
    end subroutine particle_sort_decomp




end module mod_particle_contiguous_decomp