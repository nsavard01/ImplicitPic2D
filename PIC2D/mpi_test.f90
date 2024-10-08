program mpi_test
    use iso_fortran_env, only: int32, real64
    use mpi
    implicit none

    integer :: ierr, mpi_rank, mpi_size
    call mpi_init(ierr)
    call mpi_comm_rank(mpi_comm_world, mpi_rank, ierr)
    call mpi_comm_size(mpi_comm_world, mpi_size, ierr)
    print *, 'Rank is:', mpi_rank, 'of total:', mpi_size
    call mpi_finalize(ierr)
end program mpi_test