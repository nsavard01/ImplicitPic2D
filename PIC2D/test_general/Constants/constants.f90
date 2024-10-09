module constants
    use iso_fortran_env, only: real64
    implicit none

    ! physical constants and global variables used among all modules
    real(real64), parameter :: speed_of_light = 299792458.0d0 ! m/s
    real(real64), parameter :: epsilon_0 = 8.8541878128d-12 !F/m
    real(real64), parameter :: mass_proton = 1.67262192369d-27 ! kg
    real(real64), parameter :: mass_amu = 1.66053906660d-27 ! kg
    real(real64), parameter :: mass_electron = 9.1093837015d-31 ! kg
    real(real64), parameter :: k_boltz = 1.380649d-23 ! m^2 kg s^-2 K^-1
    real(real64), parameter :: e_charge = 1.602176634d-19 ! C
    real(real64), parameter :: mu_0 = 1.25663706212d-6 ! m kg s^-2 A^-2
    real(real64), parameter :: pi_const = 4.0d0*atan(1.0d0) ! pi from atan

    ! global constants that will be read in once and used throughout the program
    integer, public, protected :: number_threads_global, N_x_global, N_y_global, number_charged_particles
contains
    ! subroutines to change global integers
    subroutine change_global_N(N_x, N_y)
        integer, intent(in) :: N_x, N_y
        N_x_global = N_x
        N_y_global = N_y

    end subroutine change_global_N

    subroutine change_global_thread(num_threads)
        integer, intent(in) :: num_threads
        number_threads_global = num_threads

    end subroutine change_global_thread

    subroutine change_global_numPart(num_part)
        integer, intent(in) :: num_part
        number_charged_particles = num_part

    end subroutine change_global_numPart

end module constants