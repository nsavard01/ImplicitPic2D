module mod_rand_generator

    use iso_fortran_env, only: int32, real64
    use iso_c_binding
    use omp_lib
    implicit none
    ! type which contains internal state random number generator, can easily be threaded
    ! Uses PCG

    public
    integer(c_int64_t), protected :: state_PCG
    ! thread private creates private state for each thread
    !$OMP threadprivate(state_PCG)

    interface 
        ! Interface for PCG generator run using C
        real(c_double) function pcg32_random_r(state) bind(c)
        use iso_c_binding
        integer(c_int64_t) :: state
        end function
    end interface

contains

    subroutine initialize_rand_PCG(num_thread, pre_determined_bool)
        integer(int32), intent(in) :: num_thread
        logical, intent(in) :: pre_determined_bool
        real(real64) :: rando(num_thread)
        integer(int32) :: i, i_thread

       
        print *, 'Initializing random number generator'
        if (.not. pre_determined_bool) call random_seed() ! Having seeding be randomly generated
        do i = 1, num_thread
            call random_number(rando(i))
        end do
        
        !$OMP parallel private(i_thread)
        i_thread = omp_get_thread_num() + 1
        state_PCG = INT((rando(i_thread)-0.5d0) * 2 * (huge(state_PCG-1)), kind = c_int64_t)
        !$OMP end parallel

    end subroutine initialize_rand_PCG

end module mod_rand_generator