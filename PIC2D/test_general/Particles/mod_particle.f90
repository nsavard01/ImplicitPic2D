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
    public :: Particle
    integer(int32), protected, public :: number_charged_particles

    ! Particle contains particle properties and stored values in phase space
    type :: Particle
        character(:), allocatable :: name !name of the particle
        integer(int64), allocatable :: number_particles_thread(:)
        ! numToCollide saves number that can collided in nullCollision
        integer(int64) :: max_indx ! maximum particles per thread
        real(real64), allocatable :: densities(:,:)
        real(real64), allocatable :: phase_space(:,:, :) !particle phase space, represents [xi_x, xi_y, v_x, v_y, v_z]
        real(real64) :: mass, charge, weight, q_over_m, q_times_weight ! mass (kg), charge(C), and weight (N/m^2 in 1D) of particles. Assume constant weight for moment
        real(real64), allocatable :: work_space(:,:,:) ! Per thread densities and general workspace over domain

    contains
        procedure, public, pass(self) :: initialize_weight_from_n_ave
        procedure, public, pass(self) :: initialize_rand_uniform
        procedure, public, pass(self) :: interpolation_particle_to_nodes
        ! procedure, public, pass(self) :: initializeRandCosine
        ! procedure, public, pass(self) :: initializeRandSine
        ! procedure, public, pass(self) :: generate_3D_Maxwellian
        ! procedure, public, pass(self) :: getKEAve
        ! procedure, public, pass(self) :: getTotalKE
        ! procedure, public, pass(self) :: getTotalMomentum
        ! procedure, public, pass(self) :: writePhaseSpace
        ! procedure, public, pass(self) :: writeLocalTemperature
    end type Particle


    interface Particle
        module procedure :: particle_constructor
    end interface Particle

contains

    type(Particle) function particle_constructor(mass, q, w_p, N_p, finalIdx, particleName, N_x, N_y) result(self)
        ! Construct particle object, sizeIncrease is fraction larger stored array compared to initial amount of particles
        ! In future, use hash function for possible k = 1 .. Nx, m amount of boundaries, p = prime number  m < p < N_x. h(k) = (k%p)%m
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
        self %max_indx = finalIdx / number_threads_global
        allocate(self%number_particles_thread(number_threads_global), self%phase_space(5,self%max_indx, number_threads_global), &
        self%work_space(N_x, N_y, number_threads_global))
        !$OMP parallel
        self%work_space(:,:,omp_get_thread_num()+1) = 0.0d0
        !$OMP end parallel
        self%number_particles_thread = N_p/number_threads_global
        
    end function particle_constructor

    subroutine initialize_weight_from_n_ave(self, n_ave, world)
        ! initialize w_p based on initial average n (N/m^3) over domain
        class(Particle), intent(in out) :: self
        real(real64), intent(in) :: n_ave
        class(domain_base), intent(in) :: world
        integer :: i, j
        real(real64) :: area

        select type (world)
        type is (domain_uniform)
            area = world%number_total_cells/world%node_volume
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
        self%weight = n_ave * area / SUM(self % number_particles_thread)
        self%q_times_weight = self%charge * self%weight
    end subroutine initialize_weight_from_n_ave

    subroutine initialize_rand_uniform(self, world, random_gen)
        ! distribute particles randomly over the domain
        class(Particle), intent(in out) :: self
        class(domain_base), intent(in) :: world
        type(rand_gen), intent(in out) :: random_gen(number_threads_global)
        integer(int32) :: i_thread, int_xi, int_eta
        integer(int64) :: i
        real(real64) :: x_pos, y_pos, L_x, L_y, xi, eta
        L_x = world%end_X - world%start_X
        L_y = world%end_Y - world%start_Y
        !$OMP parallel private(i_thread, i, x_pos, y_pos, eta, xi, int_xi, int_eta)
        i_thread = omp_get_thread_num() + 1
        do i = 1, self%number_particles_thread(i_thread)
            x_pos = random_gen(i_thread)%get_rand_num() * L_x + world%start_X
            xi = world%get_xi_from_X(x_pos)
            int_xi = int(xi)

            y_pos = random_gen(i_thread)%get_rand_num() * L_y + world%start_Y
            eta = world%get_eta_from_Y(y_pos)
            int_eta = int(eta)
            do while (world%boundary_conditions(int_xi, int_eta) /= 0 .and. world%boundary_conditions(int_xi+1, int_eta) /= 0 &
                .and. world%boundary_conditions(int_xi, int_eta+1) /= 0 .and. world%boundary_conditions(int_xi+1, int_eta+1) /= 0)
                x_pos = random_gen(i_thread)%get_rand_num() * L_x + world%start_X
                xi = world%get_xi_from_X(x_pos)
                int_xi = int(xi)

                y_pos = random_gen(i_thread)%get_rand_num() * L_y + world%start_Y
                eta = world%get_eta_from_Y(y_pos)
                int_eta = int(eta)
            end do
            
            self%phase_space(1,i, i_thread) = xi
            self%phase_space(2, i, i_thread) = eta
        end do
        !$OMP end parallel
    end subroutine initialize_rand_uniform

    subroutine interpolation_particle_to_nodes(self)
        ! interpolate particles to work space array
        class(Particle), intent(in out) :: self
        integer(int32) :: i, j, i_thread, iter
        real(real64) :: d_i, d_j, xi, eta
        !$OMP parallel private(iter, i,j, d_i, d_j, xi, eta)
        i_thread = omp_get_thread_num() + 1
        do iter = 1, self%number_particles_thread(i_thread)
            xi = self%phase_space(1,iter,i_thread)
            eta = self%phase_space(2,iter,i_thread)
            i = int(xi)
            j = int(eta)
            d_i = xi - real(i)
            d_j = eta - real(j)

            self%work_space(i,j, i_thread) = self%work_space(i,j, i_thread) + (1.0d0-d_i) * (1.0d0-d_j)
            self%work_space(i+1,j, i_thread) = self%work_space(i+1,j, i_thread) + (d_i) * (1.0d0-d_j)
            self%work_space(i,j+1, i_thread) = self%work_space(i,j+1, i_thread) + (1.0d0-d_i) * (d_j)
            self%work_space(i+1,j+1, i_thread) = self%work_space(i+1,j+1, i_thread) + (d_i) * (d_j)
        end do
        !$OMP end parallel

    end subroutine interpolation_particle_to_nodes


    ! subroutine initializeRandCosine(self, world, irand, alpha)
    !     ! Use acceptance-rejection method to distribute particles
    !     ! f(x) = 1 + alpha * Cos(2 pi x/L)
    !     class(Particle), intent(in out) :: self
    !     type(Domain), intent(in) :: world
    !     real(real64), intent(in) :: alpha
    !     integer(int32), intent(in out) :: irand(numThread)
    !     integer(int32) :: iThread, i, randNum
    !     real(real64) :: Rand1, Rand2, x_pos, val
    !     !$OMP parallel private(iThread, i, randNum, Rand1, Rand2, x_pos, val)
    !     iThread = omp_get_thread_num() + 1
    !     randNum = irand(iThread)
    !     do i = 1, self%N_p(iThread)
    !         Rand1 = ran2(randNum)
    !         Rand2 = ran2(randNum)
    !         val = (1.0d0 + alpha * COS(2 * pi * Rand2))/(1.0d0 + alpha)
    !         do while (Rand1 > val)
    !             Rand1 = ran2(randNum)
    !             Rand2 = ran2(randNum)
    !             val = (1.0d0 + alpha * COS(2 * pi * Rand2))/(1.0d0 + alpha)
    !         end do
    !         x_pos = Rand2 * world%L_domain + world%startX
    !         self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
    !     end do
    !     irand(iThread) = randNum
    !     !$OMP end parallel
    ! end subroutine initializeRandCosine

    ! subroutine initializeRandSine(self, world, irand, alpha)
    !     ! Use acceptance-rejection method to distribute particles
    !     ! f(x) = 1 + alpha * Sin(2 pi x/L)
    !     class(Particle), intent(in out) :: self
    !     type(Domain), intent(in) :: world
    !     real(real64), intent(in) :: alpha
    !     integer(int32), intent(in out) :: irand(numThread)
    !     integer(int32) :: iThread, i, randNum
    !     real(real64) :: Rand1, Rand2, x_pos, val
    !     !$OMP parallel private(iThread, i, randNum, Rand1, Rand2, x_pos, val)
    !     iThread = omp_get_thread_num() + 1
    !     randNum = irand(iThread)
    !     do i = 1, self%N_p(iThread)
    !         Rand1 = ran2(randNum)
    !         Rand2 = ran2(randNum)
    !         val = (1.0d0 + alpha * SIN(2 * pi * Rand2))/(1.0d0 + alpha)
    !         do while (Rand1 > val)
    !             Rand1 = ran2(randNum)
    !             Rand2 = ran2(randNum)
    !             val = (1.0d0 + alpha * SIN(2 * pi * Rand2))/(1.0d0 + alpha)
    !         end do
    !         x_pos = Rand2 * world%L_domain + world%startX
    !         self%phaseSpace(1,i, iThread) = world%getLFromX(x_pos)
    !     end do
    !     irand(iThread) = randNum
    !     !$OMP end parallel
    ! end subroutine initializeRandSine

    ! ------------------------ Generating and initializing velocity with Maxwellian --------------------------


    ! subroutine generate_3D_Maxwellian(self, T, world, irand, alpha, distType, v_drift)
    !     ! random velocity generator for the particle for temperature T (eV)
    !     ! Use box-muller method for random guassian variable
    !     class(Particle), intent(in out) :: self
    !     class(Domain), intent(in) :: world
    !     real(real64), intent(in) :: T, alpha, v_drift
    !     integer(int32), intent(in) :: distType
    !     type(rand_gen), intent(in out) :: random_gen(:)
    !     integer(int32) :: iThread, i
    !     integer(int32) :: irand_thread
    !     real(real64) :: U1, U2, U3, U4, v_therm, v_extra, x_pos
    !     v_therm = SQRT(T*e/self%mass)
    !     !$OMP PARALLEL PRIVATE(iThread, i, U1, U2, U3, U4, irand_thread, v_extra, x_pos)
    !     iThread = omp_get_thread_num() + 1
    !     irand_thread = irand(iThread)
    !     do i = 1, self%N_p(iThread)
    !         U1 = ran2(irand_thread)
    !         U2 = ran2(irand_thread)
    !         U3 = ran2(irand_thread)
    !         U4 = ran2(irand_thread)
    !         SELECT CASE (distType)
    !         CASE(0)
    !             v_extra = v_drift
    !         CASE(1)
    !             x_pos = world%getXFromL(self%phaseSpace(1, i, iThread))/world%L_domain
    !             v_extra = v_drift * (1.0d0 + alpha * COS(2 * pi * x_pos))
    !         CASE(2)
    !             x_pos = world%getXFromL(self%phaseSpace(1, i, iThread))/world%L_domain
    !             v_extra = v_drift * (1.0d0 - alpha * SIN(2 * pi * x_pos)) 
    !         END SELECT
    !         self%phaseSpace(2, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * COS(2.0d0 * pi * U2) + v_extra
    !         self%phaseSpace(3, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U1)) * SIN(2.0d0 * pi * U2)
    !         self%phaseSpace(4, i, iThread) = v_therm * SQRT(-2.0d0 * LOG(U3)) * SIN(2.0d0 * pi * U4)
            
    !     end do
    !     irand(iThread) = irand_thread
    !     !$OMP END PARALLEL
    ! end subroutine generate_3D_Maxwellian

    ! function getKEAve(self) result(res)
    !     ! calculate average kinetic energy (temperature) in eV
    !     class(Particle), intent(in) :: self
    !     integer(int32) :: i_thread
    !     real(real64) :: res
    !     res = 0.0d0
    !     !$OMP parallel private(iThread) reduction(+:res)
    !     iThread = omp_get_thread_num() + 1
    !     res = res + SUM(self%phaseSpace(2:4, 1:self%N_p(iThread), iThread)**2) 
    !     !$OMP end parallel
    !     res = res * self % mass * 0.5d0 / e / SUM(self%N_p)
    ! end function getKEAve

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

end module mod_particle