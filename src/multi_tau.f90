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


! ======= Hardware =======
! ======= Hardware =======


! implementation of  the algorithm proposed by Ramirez et al. in 
! "Efficient on the fly calculation of time correlation functions in comuter
! simulations" J. CHem. Phys. 133 154103 (2010)
! 
! first we test it on mean square displacement to be compared to the Frenkel
! algorithm as implemented in msd.f90

! doesn't work ;)

MODULE multi_tau

  USE constants, ONLY : dp
  implicit none 

  integer ,PARAMETER :: S=12, p=8 , m=6

  real(kind=dp), dimension (:,:,:,:) , allocatable :: DD 
  real(kind=dp), dimension (:,:,:)   , allocatable :: AA 
  real(kind=dp), dimension (:,:,:)   , allocatable :: CC 
  integer*8    , dimension (:,:,:)   , allocatable :: MM
  integer*8    , dimension (:,:,:)   , allocatable :: NN
  integer :: lmax
  integer :: nprop

CONTAINS

! *********************** SUBROUTINE alloc *************************************
!
!
! ******************************************************************************
SUBROUTINE alloc

  USE config,   ONLY :  natm
  USE io,  ONLY :  stdout

  implicit none

  !tmp 
  integer*8 :: ii,lambda

  allocate( DD ( 3 , 0 : S , 0 : p - 1 , natm ) )
  allocate( CC ( 0 : S , 0 : p - 1 , natm ) )
  allocate( NN ( 0 : S , 0 : p - 1 , natm ) )
  allocate( AA ( 3 ,0 : S , natm ) )
  allocate( MM ( 3 , 0 : S , natm ) )


  DD = 0.0_dp
  CC = 0.0_dp
  AA = 0.0_dp

  NN = 0
  MM = 0

  lmax = 0

! tmp would be in check or print info
  do ii = 1 , S
    lambda = m ** (ii)
    WRITE ( stdout ,'(a,i16,a,i16)') 'average block every',lambda,' steps for level', ii 
  enddo

  return

END SUBROUTINE alloc

! *********************** SUBROUTINE dealloc ***********************************
!
!
! ******************************************************************************
SUBROUTINE dealloc

  implicit none

  deallocate( DD )
  deallocate( CC )
  deallocate( NN )
  deallocate( AA )
  deallocate( MM )

  return

END SUBROUTINE dealloc

! *********************** SUBROUTINE multi_tau_main ****************************
!
!
! ******************************************************************************
SUBROUTINE multi_tau_main ( wx , wy , wz , ncall)

  USE config,   ONLY :  natm
  USE md,       ONLY :  dt 

  implicit none

  ! global
  real(kind=dp) , dimension ( natm ) :: wx , wy , wz 
  integer                               :: ncall 
  ! local
  real(kind=dp) :: xtime, dtime
  integer :: j , lv , ia , lambda , l , ii  
  real(kind=dp), dimension ( : , : , : , :) , allocatable :: Dtmp

  allocate( Dtmp ( 3 , 0 : S , 0 : p - 1 , natm ) )

  dtime = nprop * dt
  ! ===================================================
  !  determine current maximum number of blocks: iblm
  ! ===================================================
  l = 1
  ii = ncall / m 
  do while  (ii.ne.0)
    l = l + 1
    ii = ii / m 
  ! ===========================================
  !  test maximu time not longer than tdifmax:
  ! ===========================================
    xtime = dtime * ( m ** ( l ) )
    if ( xtime .gt. 5.0_dp ) ii = 0 ! tmp
  enddo
!  print*,l,ncall
  if ( l.gt.lmax) lmax = l
  do lv = 0 , l
    lambda = m ** (lv)
    if ( mod ( ncall, lambda ) .eq. 0 ) then
      print*,'lv updated',lv
      Dtmp = DD
      do ia = 1 , natm
          do j = 1 , p - 1  
            DD ( 1 , lv , j , ia ) = Dtmp ( 1 , lv , j - 1 , ia ) 
            DD ( 2 , lv , j , ia ) = Dtmp ( 1 , lv , j - 1 , ia ) 
            DD ( 3 , lv , j , ia ) = Dtmp ( 1 , lv , j - 1 , ia ) 
          enddo
            DD ( 1 , lv , 0 , ia ) = wx ( ia )
            DD ( 2 , lv , 0 , ia ) = wy ( ia )
            DD ( 3 , lv , 0 , ia ) = wz ( ia )
!            print*,'DD',DD(1,lv,0,1)
         if ( lv .eq. 0 ) then
           do j = 0 , p - 1
             CC ( lv , j , ia ) = CC ( lv , j , ia ) + DD ( 1 , lv , 0 , ia ) * DD ( 1 , lv , j , ia ) + &
                                                       DD ( 2 , lv , 0 , ia ) * DD ( 2 , lv , j , ia ) + & 
                                                       DD ( 3 , lv , 0 , ia ) * DD ( 3 , lv , j , ia )
             NN ( lv , j , ia ) = NN ( lv , j , ia ) + 1
           enddo
         else
           do j = p / m , p - 1 
             CC ( lv , j , ia ) = CC ( lv , j , ia ) + DD ( 1 , lv , 0 , ia ) * DD ( 1 , lv , j , ia ) + &
                                                       DD ( 2 , lv , 0 , ia ) * DD ( 2 , lv , j , ia ) + &
                                                       DD ( 3 , lv , 0 , ia ) * DD ( 3 , lv , j , ia )
             NN ( lv , j , ia ) = NN ( lv , j , ia ) + 1
           enddo
         endif    
         AA ( 1 , lv , ia ) = AA ( 1 , lv , ia ) + wx ( ia )  
         AA ( 2 , lv , ia ) = AA ( 2 , lv , ia ) + wy ( ia )  
         AA ( 3 , lv , ia ) = AA ( 3 , lv , ia ) + wz ( ia )  
         MM ( 1 , lv , ia ) = MM ( 1 , lv , ia ) + 1   
         MM ( 2 , lv , ia ) = MM ( 2 , lv , ia ) + 1   
         MM ( 3 , lv , ia ) = MM ( 3 , lv , ia ) + 1 
         if ( MM ( 1 , lv , ia ) .eq. m ) then
           AA ( 1 , lv + 1 , ia ) = AA ( 1 , lv , ia ) / DBLE ( m )
           AA ( 1 , lv , ia ) = 0.0_dp   
           MM ( 1 , lv , ia ) = 0   
         endif  
         if ( MM ( 2 , lv , ia ) .eq. m ) then
           AA ( 2 , lv + 1 , ia ) = AA ( 2 , lv , ia ) / DBLE ( m )
           AA ( 2 , lv , ia ) = 0.0_dp   
           MM ( 2 , lv , ia ) = 0   
         endif  
         if ( MM ( 3 , lv , ia ) .eq. m ) then
           AA ( 3 , lv + 1 , ia ) = AA ( 3 , lv , ia ) / DBLE ( m )
           AA ( 3 , lv , ia ) = 0.0_dp   
           MM ( 3 , lv , ia ) = 0   
         endif  
       enddo
    endif
  enddo

  deallocate( Dtmp )

  return

END SUBROUTINE multi_tau_main


! *********************** SUBROUTINE multi_tau_write_output ********************
!
!
! ******************************************************************************

SUBROUTINE multi_tau_write_output

  USE io,  ONLY :  ionode , stdout 
  USE md,       ONLY :  dt 
  USE config,   ONLY :  natm

  implicit none

  ! local
  integer*8 :: lv , j ,ia
  real(kind=dp) :: dtime , tk ,fk

  dtime = nprop * dt
 
  print*,'lmax =',lmax
  print*,'natm =',natm

  do lv = 0 , lmax - 1
    if ( lv .eq. 0 ) then
      do j = 0 , p -1 
        tk = j * m ** ( lv ) * dtime
        do ia = 1 , natm 
          fk = fk + ( CC ( lv , j , ia ) / DBLE ( NN ( lv , j , ia ) ) )  
                      if ( NN ( lv , j , ia ) .eq. 0 ) WRITE ( stdout ,'(a,3i6)') 'N1 0', lv, j ,ia
                      if ( CC ( lv , j , ia ) .eq. 0 ) WRITE ( stdout ,'(a,3i6)') 'C1 0', lv, j ,ia
        enddo
        if ( ionode )  WRITE ( 1000 , * ) tk , fk  
      enddo
    else
      do j = p / m  , p -1 
        tk = j * m ** ( lv ) * dtime
        do ia = 1 , natm
          fk = fk + ( CC ( lv , j , ia ) / DBLE ( NN ( lv , j , ia ) )  )
                      if ( NN ( lv , j , ia ) .eq. 0 ) WRITE ( stdout ,'(a,3i6)') 'N1 p/m', lv, j ,ia
                      if ( CC ( lv , j , ia ) .eq. 0 ) WRITE ( stdout ,'(a,3i6)') 'C1 p/m', lv, j ,ia
        enddo
        if ( ionode )  WRITE ( 1000 , * ) tk , fk  
      enddo
  
    endif
  
  enddo

  return

END SUBROUTINE multi_tau_write_output

END MODULE multi_tau
