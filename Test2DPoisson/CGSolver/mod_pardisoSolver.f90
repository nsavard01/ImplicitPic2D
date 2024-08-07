module mod_pardisoSolver
    use iso_fortran_env, only: int32, int64, real64
    implicit none

    ! Class for pardisoSolver, either as last portion of multigrid or stand alone solver

    private
    public :: pardisoSolver

    type :: pardisoSolver
        ! store grid quantities
        real(real64), allocatable :: sourceTerm(:), MatValues(:), solution(:)
        integer(int32), allocatable :: columnIndex(:), rowIndex(:)
        integer(int32) :: matDimension, maxfct, mnum, mtype, iparm(64), msglvl, error
        integer(int64) :: pt(64) ! pt needs to be int64, otherwise have issues with pardisoinit
    contains
        procedure, public, pass(self) :: initializePardiso
        procedure, public, pass(self) :: runPardiso
    end type

    interface pardisoSolver
        module procedure :: pardisoSolver_constructor
    end interface pardisoSolver

contains

    type(pardisoSolver) function pardisoSolver_constructor(matDimension) result(self)
        ! Construct object, set initial variables
        integer(int32), intent(in) :: matDimension
        self%matDimension = matDimension
        self%msglvl = 0
        self%error = 0
        self%mtype = 11
        self%mnum = 0
        self%maxfct = 1
        self%iparm = 0
        self%pt = 0
        allocate(self%sourceTerm(matDimension), self%rowIndex(matDimension+1), self%solution(matDimension))
    end function pardisoSolver_constructor

    subroutine initializePardiso(self, mnum, mtype, maxfct, msglvl)
        class(pardisoSolver), intent(in out) :: self
        integer(int32), intent(in) :: mnum, mtype, maxfct, msglvl
        self%mtype = mtype
        self%maxfct = maxfct
        self%msglvl = msglvl
        self%mnum = mnum
        call pardisoinit(self%pt,self%mtype,self%iparm)
        ! ignore perm, say nrhs = 1, setup phase 11
        call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, 11, self%matDimension, self%MatValues, &
            self%rowIndex, self%columnIndex, 1, 1, self%iparm, self%msglvl, self%sourceTerm, self%solution, self%error)
         
        ! phase 22 numerical factorization
        call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, 22, self%matDimension, self%MatValues, &
            self%rowIndex, self%columnIndex, 1, 1, self%iparm, self%msglvl, self%sourceTerm, self%solution, self%error)


    end subroutine initializePardiso

    subroutine runPardiso(self)
        class(pardisoSolver), intent(in out) :: self
        
        ! phase 33, run and produce output to the solution
        call pardiso(self%pt, self%maxfct, self%mnum, self%mtype, 33, self%matDimension, self%MatValues, &
            self%rowIndex, self%columnIndex, 1, 1, self%iparm, self%msglvl, self%sourceTerm, self%solution, self%error)


    end subroutine runPardiso

    


end module mod_pardisoSolver