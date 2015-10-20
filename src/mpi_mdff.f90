! mpi related module for mdff
MODULE mpimdff

#ifdef MPI
  INCLUDE "mpif.h"
#endif

  integer,           SAVE :: myrank           !< rank of the actual proc   
  integer,           SAVE :: numprocs         !< total number of procs    

  TYPE  :: decomposition 
    integer :: dim_data
    integer :: istart
    integer :: iend
    character(len=5) :: label
  END TYPE

CONTAINS

! *********************** SUBROUTINE do_split **********************************
!
!> \brief
!! this routine split the number of atoms into np procs
!! istart and iend are the atom index for myrank task
!! nex version mardi 15 oct
!
!> \param[in]  n     number of atoms
!> \param[in]  mrank current task index
!> \param[in]  np    number of mpi tasks
!> \param[out] dec decomposition type look at definition in module mpimdff 
!
!> \author
!! FMV
!
! ******************************************************************************
SUBROUTINE do_split ( n , mrank , np , dec , lab )

  USE io,                  ONLY :  ionode , stdout

  implicit none

  ! global
  integer :: n       !< number of atoms to be splitted
  integer :: np      !< number of mpi tasks
  integer :: mrank   !< local index proc
  character(len=5) :: lab
  TYPE (decomposition) :: dec

  ! local 
  integer :: imin,imax
  integer :: istartV(0:np-1),iendV(0:np-1)
  integer :: splitnumberV(0:np-1),isteps,x,y,me

#ifndef MPI
  ! quick return
  return
#endif        

  imin = 1
  imax = n

  isteps = (imax-imin)+1
  x = isteps/np
  y = MOD ( isteps , np )

  do me = 0,np-1
     if ((me.eq.0).or.(me.gt.y)) then
        splitnumberV(me) = x
     else if ((me.gt.0).or.(me.lt.(y+1))) then
        splitnumberV(me) = x+1
     endif
  enddo
  do me = 0,np-1
     if (me.eq.0) then
        istartV(0) = imin
        iendV(0)   = imin + ( x -1 )
     else if (me.gt.0) then
        istartV(me) = iendV(me-1) + 1
        iendV(me)   = istartV(me) + splitnumberV(me) - 1
     endif
  enddo

  dec%istart   = istartV(mrank)
  dec%iend     = iendV(mrank)
  dec%dim_data = ( iendV(mrank) - istartV(mrank) ) + 1
  dec%label    = lab

  if ( ionode ) then
    WRITE ( stdout ,'(a,a5,a)')      'paralelisation - ',dec%label,' decomposition'
    do me = 0,np-1
      WRITE ( stdout ,'(a,i4,a8,i8,a3,i8,a,i8)')      'rank =      ',me,dec%label,istartV(me),'to',iendV(me), ' load : ',( iendV(me) - istartV(me) ) + 1
    enddo     
  endif

  return

END SUBROUTINE do_split

! *********************** SUBROUTINE MPI_ALL_REDUCE_DOUBLE *********************
!
!
!
! ******************************************************************************
SUBROUTINE MPI_ALL_REDUCE_DOUBLE ( vec_result , ndim )

  USE constants, ONLY : dp
  implicit none

  ! global
  integer :: ndim
  real(kind=dp), dimension ( ndim ) :: vec_result
  ! local 
  integer :: ierr
  real(kind=dp), dimension ( : ) , allocatable :: vec_sum

  allocate ( vec_sum ( ndim ) )
  vec_sum=0.0_dp

#ifdef MPI
  CALL MPI_ALLREDUCE( vec_result , vec_sum , ndim , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , ierr )
#else
  vec_sum = vec_result
#endif
  vec_result = vec_sum

  deallocate ( vec_sum )

  return

END SUBROUTINE MPI_ALL_REDUCE_DOUBLE

! *********************** SUBROUTINE MPI_ALL_GATHER_DOUBLE *********************
!
!
!
! ******************************************************************************
SUBROUTINE MPI_ALL_GATHER_DOUBLE ( vec_result , ndim , dec )

  USE constants,        ONLY : dp
  implicit none

  ! global
  integer :: ndim
  real(kind=dp), dimension ( ndim ) :: vec_result
  TYPE(decomposition)               :: dec
  ! local 
  integer :: ierr

#ifdef MPI
  CALL MPI_ALLGATHER ( vec_result ( dec%istart : dec%iend ) , dec%dim_data , MPI_DOUBLE_PRECISION , vec_result , dec%dim_data , MPI_DOUBLE_PRECISION , MPI_COMM_WORLD,ierr)
#endif

  return

END SUBROUTINE MPI_ALL_GATHER_DOUBLE



! *********************** SUBROUTINE MPI_ALL_REDUCE_INTEGER ********************
!
!
!
! ******************************************************************************
SUBROUTINE MPI_ALL_REDUCE_INTEGER ( vec_result , ndim )

  implicit none

  ! global
  integer :: ndim
  integer , dimension ( ndim ) :: vec_result
  ! local 
  integer :: ierr
  integer , dimension ( : ) , allocatable :: vec_sum

  allocate ( vec_sum ( ndim ) )
  vec_sum=0

#ifdef MPI
  CALL MPI_ALLREDUCE( vec_result , vec_sum , ndim , MPI_INTEGER , MPI_SUM , MPI_COMM_WORLD , ierr )
#else
  vec_sum = vec_result
#endif
  vec_result = vec_sum

  deallocate ( vec_sum )

  return

END SUBROUTINE MPI_ALL_REDUCE_INTEGER

! *********************** SUBROUTINE MPI_ALL_REDUCE_DOUBLE_SCALAR **************
!
!
!
! ******************************************************************************
SUBROUTINE MPI_ALL_REDUCE_DOUBLE_SCALAR ( sresult )

  USE constants, ONLY : dp

  implicit none

  real(kind=dp), intent (inout) :: sresult
  ! local 
  integer          :: ierr
  real(kind=dp) :: ssum

  ssum=0.0_dp
#ifdef MPI
  CALL MPI_ALLREDUCE( sresult , ssum , 1 , MPI_DOUBLE_PRECISION , MPI_SUM , MPI_COMM_WORLD , ierr )
#else
  ssum = sresult
#endif
  sresult = ssum

  return

END SUBROUTINE


! *********************** SUBROUTINE MPI_ALL_REDUCE_DOUBLE_SCALAR **************
!
!
!
! ******************************************************************************
SUBROUTINE MPI_ALL_REDUCE_INTEGER_SCALAR ( sresult )

  USE constants, ONLY : dp

  implicit none

  integer, intent (inout) :: sresult
  ! local 
  integer          :: ierr
  integer          :: ssum

  ssum=0
#ifdef MPI
  CALL MPI_ALLREDUCE( sresult , ssum , 1 , MPI_INTEGER , MPI_SUM , MPI_COMM_WORLD , ierr )
#else
  ssum = sresult
#endif
  sresult = ssum

  return

END SUBROUTINE


#ifdef MPI
SUBROUTINE CLEAN_STOP

  implicit none

  integer :: rc , ierr

  if ( myrank .gt. 0 ) then
    CALL MPI_ABORT(MPI_COMM_WORLD, RC , ierr)
  else
    STOP
  endif

END SUBROUTINE CLEAN_STOP
#endif


END MODULE mpimdff
