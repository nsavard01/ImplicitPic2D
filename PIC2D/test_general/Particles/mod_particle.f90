module mod_particle

    use iso_fortran_env, only: int32, real64, output_unit
    use constants
    use mod_domain_base
    use mod_domain_uniform
    use mod_domain_curv
    use mod_rand_generator
    use omp_lib
    implicit none

    private
    public :: Particle, push_particles_uniform, interpolate_particle_charge_density
    real(real64), allocatable, public :: particle_work_space(:,:,:)

    ! Particle contains particle properties and stored values in phase space
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int64) :: max_indx, total_number_particles ! maximum particles per thread
        real(real64), allocatable :: densities(:,:)
        real(real64), allocatable :: logical_position(:,:, :), velocity(:,:,:) !particle phase space, represents [xi_x, xi_y, v_x, v_y, v_z]
        real(real64) :: mass, charge, weight, q_over_m, q_times_weight ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64) :: sum_v_sqr, sum_v

    contains
        procedure, public, pass(self) :: initialize_base_variables
        procedure, public, pass(self) :: initialize_weight_from_n_ave
        procedure, public, pass(self) :: initialize_rand_uniform
        procedure, public, pass(self) :: interpolation_particle_to_nodes
        ! procedure, public, pass(self) :: interpolation_particle_to_nodes_and_sort
        ! procedure, public, pass(self) :: interpolation_particle_to_nodes_sorted
        procedure, public, pass(self) :: resize_particle_arrays
        procedure, public, pass(self) :: particle_mover_uniform
        ! procedure, public, pass(self) :: initializeRandCosine
        ! procedure, public, pass(self) :: initializeRandSine
        procedure, public, pass(self) :: initialize_maxwellian_temperature
        procedure, public, pass(self) :: get_KE_Ave
        procedure, public, pass(self) :: get_sum_totals
        ! procedure, public, pass(self) :: getTotalKE
        ! procedure, public, pass(self) :: getTotalMomentum
        ! procedure, public, pass(self) :: writePhaseSpace
        ! procedure, public, pass(self) :: writeLocalTemperature
    end type Particle

contains

    subroutine initialize_base_variables(self, mass, q, w_p, N_p, finalIdx, particleName, N_x, N_y)
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: mass, q, w_p
        integer(int64), intent(in) :: N_p, finalIdx
        integer(int32), intent(in) :: N_x, N_y
        character(*), intent(in) :: particleName
        self % name = particleName
        self % mass = mass
        self % charge = q
        self % weight = w_p
        self%q_over_m = q/mass
        self%q_times_weight = q * w_p
        self%total_number_particles = (N_p/number_threads_global) * number_threads_global
        self %max_indx = finalIdx / number_threads_global
        allocate(self%logical_position(2,self%max_indx, number_threads_global), &
        self%velocity(3,self%max_indx, number_threads_global), self%densities(N_x, N_y))
        if (.not. allocated(particle_work_space)) then
            allocate(particle_work_space(N_x, N_y, number_threads_global))
            !$OMP parallel
            particle_work_space(:,:,omp_get_thread_num()+1) = 0.0d0
            !$OMP end parallel
        end if
    end subroutine 

    subroutine initialize_weight_from_n_ave(self, n_ave, world)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave
        class(domain_base), intent(in) :: world
        integer :: i, j
        real(real64) :: area

        select type (world)
        type is (domain_uniform)
            area = world%number_total_cells/world%inv_node_volume
        type is (domain_curv)
            area = 0
            !$OMP parallel private(j,i) reduction(+:area)
            !$OMP do
            do j = 1, world%N_y-1
                do i = 1, world%N_x-1
                    if (world%boundary_conditions(i,j) == 0 .or. world%boundary_conditions(i+1,j) == 0 &
                    .or. world%boundary_conditions(i,j+1) == 0 .or. world%boundary_conditions(i+1,j+1) == 0) then
                        area = area + world%del_x(i) * world%del_y(j)
                    end if
                end do
            end do
            !$OMP end do
            !$OMP end parallel
        end select 
        self%weight = n_ave * area / real(self%total_number_particles, kind = 8)
        self%q_times_weight = self%charge * self%weight
    end subroutine initialize_weight_from_n_ave

    subroutine initialize_rand_uniform(self, world)
        ! distribute particles randomly over the domain
        class(Particle), intent(in out) :: self
        class(domain_base), intent(in) :: world
        
    end subroutine initialize_rand_uniform

    subroutine interpolation_particle_to_nodes(self, i_thread, const)
        ! interpolate particles to work space array
        class(Particle), intent(in) :: self
        integer(int32), intent(in) :: i_thread
        real(real64), intent(in) :: const

    end subroutine interpolation_particle_to_nodes

    subroutine interpolate_particle_charge_density(particle_list)
        class(Particle), intent(in) :: particle_list(number_charged_particles)
        integer(int32) :: i_thread, part_idx

       !$OMP parallel private(i_thread, part_idx)
        i_thread = omp_get_thread_num() + 1
        particle_work_space(:,:, i_thread) = 0.0d0
        do part_idx = 1, number_charged_particles
            call particle_list(part_idx)%interpolation_particle_to_nodes(i_thread, particle_list(part_idx)%q_times_weight)
        end do
        !$OMP end parallel


    end subroutine interpolate_particle_charge_density


    subroutine particle_mover_uniform(self, i_thread, E_Field, world, del_t)
        class(Particle), intent(in out) :: self
        type(domain_uniform), intent(in) :: world
        integer(int32), intent(in) :: i_thread
        real(real64), intent(in) :: E_Field(2,world%N_x,world%N_y), del_t
    
    end subroutine particle_mover_uniform

    subroutine push_particles_uniform(particle_list, E_Field, world, del_t)
        class(Particle), intent(in out) :: particle_list(number_charged_particles)
        class(domain_uniform), intent(in) :: world
        real(real64), intent(in) :: E_Field(2,world%N_x,world%N_y), del_t
        integer(int32) :: i_thread, part_idx
        !$OMP parallel private(i_thread, part_idx)
        i_thread = omp_get_thread_num() + 1
        do part_idx = 1, number_charged_particles
            call particle_list(part_idx)%particle_mover_uniform(i_thread, E_field, world, del_t)
        end do
        !$OMP end parallel
    end subroutine

    subroutine get_sum_totals(self)
        class(Particle), intent(in out) :: self
    
    end subroutine get_sum_totals



    subroutine initialize_maxwellian_temperature(self, T)
        ! random velocity generator for the particle for temperature T (eV)
        ! Use box-muller method for random guassian variable
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: T
    end subroutine initialize_maxwellian_temperature

    function get_KE_Ave(self) result(res)
        ! calculate average kinetic energy (temperature) in eV
        class(Particle), intent(in) :: self
        real(real64) :: res
        res = 0.0d0
        res = self%sum_v_sqr * self % mass * 0.5d0 / e_charge / real(self%total_number_particles, kind = 8)
    end function get_KE_Ave

    subroutine resize_particle_arrays(self)
        class(Particle), intent(in out) :: self
    end subroutine

end module mod_particle