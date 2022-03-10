! MDFF parallel Molecular Dynamics ... For Fun
! Copyright (C) 2011  F. Vasconcelos
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
! ===== fmV =====

! ======= Hardware =======
! ======= Hardware =======

! *********************** SUBROUTINE init_random_seed **************************
!> \brief
!> initialisation of seed for the fortran random number generator RANDOM_NUMBER
!> \note
!  modification for mpi execution to generate and propagates the same seed on
!  all processes.
! ******************************************************************************
SUBROUTINE init_random_seed(rank,nump)

  USE io,       ONLY :  stdout
  USE mpimdff

  implicit none
  integer :: rank,nump

  ! local
  integer :: i , n , clock , root, tag, proc
  INTEGER , dimension (:) , allocatable :: SEED
#ifdef MPI
  integer :: ierr
  integer status(MPI_STATUS_SIZE)
#else
  rank = 0
#endif
  if ( rank .eq. 0 ) then
    CALL RANDOM_SEED(SIZE = n)
    do proc = 1 , nump-1
      tag = 0
#ifdef MPI
      CALL MPI_SEND( n , 1, MPI_INTEGER, proc, tag, MPI_COMM_WORLD, ierr)
#endif
    enddo
  else
    root = 0
    tag = 0
#ifdef MPI
    call MPI_RECV ( n , 100, MPI_INTEGER, root, tag, MPI_COMM_WORLD, status, ierr)
#endif
  endif
  allocate(SEED(n))
  i=1
  CALL SYSTEM_CLOCK(COUNT = clock)

  seed = clock + 48 * (/ (i - 1, i = 1, n) /)

#ifdef MPI
  CALL MPI_ALL_REDUCE_INTEGER ( seed , n )
#endif
  CALL RANDOM_SEED(PUT = seed)
  !WRITE(stdout,'(3i12)') myrank,seed,n 
  !STOP
  deallocate(SEED)

  return

END SUBROUTINE


! *********************** SUBROUTINE knuth *************************************
!> \brief
!! Random variate from the standard normal distribution.
!! The distribution is gaussian with zero mean and unit variance         
!> \note
!! Knuth D, The Art of Computer Programming , (2nd edition Addison-Wesley), 1978                                     
!> \note
!! Adapted from F.24 (Allen-Tildesley)
!> \author
!! Allen-Tildesley
!> \param[in] mean mean value of the gaussian distribution
!> \param[in] sigma variance value of the gaussian distribution
!> \param[out] G random number from the gaussian distribution define by mean and sigma
! ******************************************************************************
SUBROUTINE knuth ( G , mean , sigma )

  USE constants, ONLY : dp
  implicit none

  ! global
  real(kind=dp), intent (out) :: G
  real(kind=dp), intent (in)  :: mean, sigma

  ! local
  real(kind=dp) :: A1, A3, A5, A7, A9
  PARAMETER ( A1 = 3.949846138_dp, A3 = 0.252408784_dp )
  PARAMETER ( A5 = 0.076542912_dp, A7 = 0.008355968_dp )
  PARAMETER ( A9 = 0.029899776_dp                   )
  real(kind=dp) :: x
  real(kind=dp) :: summ, r, r2
  integer :: i
  integer :: iseed

  summ = 0.0_dp
  do i = 1, 12
      CALL RANDOM_SEED(SIZE = iseed)
     CALL RANDOM_NUMBER(HARVEST = x)
     summ = summ + x
  enddo

  r  = ( summ - 6.0_dp ) / 4.0_dp
  r2 = r * r

  G = ((((A9*R2+A7)*R2+A5)*R2+A3)*R2+A1)*R

  G =  mean + G*sigma

  return

END SUBROUTINE knuth 

! *********************** SUBROUTINE boxmuller_polar ***************************
!> \brief
!! Box-Muller polar method
!> \note
!! http://en.wikipedia.org/wiki/Box-Muller_transform
!> \param[in] mean mean value of the gaussian distribution
!> \param[in] sigma variance value of the gaussian distribution
!> \param[out] G random number from the gaussian distribution define by mean and sigma
! ******************************************************************************

SUBROUTINE boxmuller_polar (G, mean, sigma)

  USE constants,        ONLY : dp 
  implicit none

  ! global
  real(kind=dp), intent (out) :: G
  real(kind=dp), intent (in)  :: mean, sigma

  ! local
  real(kind=dp) :: G1,U,V,S
  integer :: iseed
  
      CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = U)
      CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = V)
  U = (2.0_dp*U)-1.0_dp
  V = (2.0_dp*V)-1.0_dp
  S = U*U+V*V
  do while ( S .eq. 0 .or. S .ge. 1) 
      CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = U)
      CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = V)
    U = (2.0_dp*U)-1.0_dp
    V = (2.0_dp*V)-1.0_dp
    S = U * U + V * V 
  enddo
  G1 = -2.0_dp * LOG ( S )
  G = U * SQRT ( G1 / S ) 
  G = mean + G * sigma

  return

END SUBROUTINE boxmuller_polar

! *********************** SUBROUTINE boxmuller_basic ***************************
!> \brief
!! Box-Muller cartesian method
!> \note
!! http://en.wikipedia.org/wiki/Box-Muller_transform
!> \param[in] mean mean value of the gaussian distribution
!> \param[in] sigma variance value of the gaussian distribution
!> \param[out] G random number from the gaussian distribution define by mean and sigma
! ******************************************************************************
SUBROUTINE boxmuller_basic (G, mean, sigma)
  
  USE constants,        ONLY : dp , tpi
  implicit none

  ! global
  real(kind=dp), intent (out) :: G
  real(kind=dp), intent (in)  :: mean, sigma

  ! local
  real(kind=dp) :: C,U,V,R
  integer :: iseed

  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = U)
  CALL RANDOM_SEED(SIZE = iseed)
  CALL RANDOM_NUMBER(HARVEST = V)
  
  R = SQRT ( -2.0_dp * LOG ( U ) )
  C = COS ( tpi * V )
  G = R * C
  G = mean + G * sigma

  return

END SUBROUTINE boxmuller_basic

! *********************** SUBROUTINE gammadev **********************************
! adapted from resamplingkin.f90 giovannibussi
! https://sites.google.com/site/giovannibussi/Research/algorithms#TOC-Stochastic-velocity-rescaling
! gamma-distributed random number, implemented as described in numerical recipes
! ******************************************************************************
SUBROUTINE gammadev(G,n)

  USE constants,         ONLY : dp
  USE io,               ONLY :  stderr

  implicit none
  integer, intent(in) :: n
  integer j
  real(kind=dp) :: am,e,s,v1,v2,x,y,U,V,W,G
  integer :: iseed
  
  if(n.lt.1) then
   WRITE (stderr , '(a)') 'bad argument in gammadev'
   STOP
  endif
  if(n.lt.6)then
    x=1.
    do 11 j=1,n
      CALL RANDOM_SEED(SIZE = iseed)
      CALL RANDOM_NUMBER(HARVEST = U)
      x=x*U
11    continue
    x=-log(x)
  else
1   CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = U)
    CALL RANDOM_NUMBER(HARVEST = V)
    v1=2.*U-1.
    v2=2.*V-1.
    if(v1**2+v2**2.gt.1.)goto 1
    y=v2/v1
    am=n-1
    s=sqrt(2.*am+1.)
    x=s*y+am
    if(x.le.0.)goto 1
    e=(1.+y**2)*exp(am*log(x/am)-s*y)
    CALL RANDOM_NUMBER(HARVEST = W)
    if(W.gt.e)goto 1
  endif
  G=x
  return

END SUBROUTINE gammadev 


! ===== fmV =====
