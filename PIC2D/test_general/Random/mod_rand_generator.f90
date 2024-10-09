module mod_rand_generator

    use iso_fortran_env, only: int32, real64
    implicit none
    ! type which contains internal state random number generator, can easily be threaded
    ! based on ran from numberical recipes f90

    private
    public :: rand_gen
    integer(int32), parameter :: IA=16807,IM=2147483647,IQ=127773,IR=2836
    real(real64), parameter :: am = 1.0d0/IM

    type :: rand_gen
        integer :: states(2) ! contiguous is faster

    contains
        procedure, public, pass(self) :: get_rand_num
    end type

    interface rand_gen
        module procedure :: rand_gen_constructor
    end interface rand_gen

contains

    type(rand_gen) function rand_gen_constructor(int_random) result(self)
        ! input integer for initialization, presumably randomly generated
        integer(int32), intent(in) :: int_random
        self%states(2) = ior(ieor(888889999, abs(int_random)), 1)
        self%states(1) = ieor(777755555,abs(int_random))
    end function rand_gen_constructor

    function get_rand_num(self) result(res)
        class(rand_gen), intent(in out) :: self
        integer :: k
        real(real64) :: res
        self%states(1)=ieor(self%states(1),ishft(self%states(1),13))
        self%states(1)=ieor(self%states(1),ishft(self%states(1),-17))
        self%states(1)=ieor(self%states(1),ishft(self%states(1),5))
        k=self%states(2)/IQ 
        self%states(2)=IA*(self%states(2)-k*IQ)-IR*k
        if (self%states(2) < 0) self%states(2)=self%states(2)+IM
        res=am*ior(iand(IM,ieor(self%states(1),self%states(2))),1)
    end function get_rand_num

end module mod_rand_generator