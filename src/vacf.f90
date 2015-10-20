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
#include "symbol.h"
! ======= Hardware =======

MODULE vacf 

! ==============================================================================================================
! note the commented statements are related to the convential scheme to compute mean square displacement 
! maybe add a tag in vacftag for that ... but it is really less efficient than the n-order scheme in msd.f90 
! removed 27/06/11
! ===============================================================================================================
! well I do not understand what the code is doing ;)
! tmax , t0max, tdifmax ????????????????? 
! ===============================================================================================================

  USE constants, ONLY : dp
  implicit none

  integer, PARAMETER  :: tmax=10000, t0max=10000 
  integer             :: it0  
  real(kind=dp)    :: tdifmax

  integer             :: tvacf , t0 
  real(kind=dp)    :: dtime 

  integer , dimension (:) , allocatable :: ttv0                           ! t0max
  integer , dimension (:) , allocatable :: nvacf                          ! tmax
  real(kind=dp) , dimension (:)   , allocatable :: vacff               ! tmax 
  real(kind=dp) , dimension (:,:) , allocatable :: vxt0 , vyt0 , vzt0  ! natm , t0max
!  integer :: nprop
!  logical :: lvacf

CONTAINS


! *********************** SUBROUTINE vacf_init *********************************
!
! init vacf calc
!
! ******************************************************************************
SUBROUTINE vacf_init

  USE control,          ONLY :  lvacf
  USE io,  ONLY :  stdin , stdout , ionode

  implicit none

  ! local
  integer            :: ioerr
  character(len=132) :: filename

  namelist /vacftag/   it0     ,  &
                       tdifmax 

  if ( .not. lvacf ) return

  CALL vacf_default_tag
  
  ! ============================
  ! reads vacftag tags
  ! ============================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , vacftag,iostat=ioerr)
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : vacftag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : vacftag wrong tag'
   STOP
  endif

  CLOSE ( stdin )

  CALL vacf_check_tag

  CALL vacf_print_info(stdout)

  return 
 
END SUBROUTINE vacf_init


! *********************** SUBROUTINE vacf_default_tag **************************
!
! set default values to vacf tag
!
! ******************************************************************************
SUBROUTINE vacf_default_tag

  implicit none

  tdifmax = 100.0_dp
  it0     = 1

  return 
 
END SUBROUTINE vacf_default_tag

! *********************** SUBROUTINE vacf_check_tag ****************************
!
! check vacf tag values
!
! ******************************************************************************
SUBROUTINE vacf_check_tag

  implicit none

  return 
 
END SUBROUTINE vacf_check_tag


! *********************** SUBROUTINE vacf_print_info ***************************
!
! print info to outputs (stdout)
!
! ******************************************************************************
SUBROUTINE vacf_print_info(kunit)

  USE io,  ONLY :  ionode 

  implicit none

  ! local
  integer :: kunit

  if ( ionode ) then
         blankline(kunit)
         WRITE ( kunit ,'(a)')           'velocity auto-correlation function   : '
         WRITE ( kunit ,'(a,f10.4)')     'tdifmax                              = ',tdifmax
         WRITE ( kunit ,'(a,i5)')        'it0                                  = ',it0
         WRITE ( kunit ,'(a)')           'output file                          : VACFFF'
  endif

  return 
 
END SUBROUTINE vacf_print_info

! *********************** SUBROUTINE vacf_alloc ********************************
! 
! allocate principal arrays for vacf
! 
! ******************************************************************************
SUBROUTINE vacf_alloc

  USE control,          ONLY :  lvacf
  USE md,               ONLY :  dt , npropr
  USE config,           ONLY :  natm

  implicit none

  ! local
  integer :: i

  if ( .not. lvacf ) return

  allocate (    vxt0 ( natm , t0max ) , vyt0 ( natm , t0max ) , vzt0 ( natm , t0max ) )
  allocate (    ttv0 ( t0max ) )                    
  allocate (   vacff ( tmax  ) )
  allocate (   nvacf ( tmax  ) )

  t0 = 0
  tvacf = 0
  dtime = npropr * dt 
  do i = 1, tmax
    nvacf(i) = 0
    vacff(i) = 0
  enddo 

  return 
 
END SUBROUTINE vacf_alloc

! *********************** SUBROUTINE vacf_dealloc ******************************
! 
! deallocate principal arrays for vacf
! 
! ******************************************************************************
SUBROUTINE vacf_dealloc

  USE control,          ONLY :  lvacf

  implicit none

  if ( .not. lvacf ) return 
 
  deallocate (   vacff ) 
  deallocate (   nvacf )
  deallocate (    ttv0 )                    
  deallocate (    vxt0 , vyt0 , vzt0 )

  return 
 
END SUBROUTINE vacf_dealloc

! *********************** SUBROUTINE vacf_main *********************************
!
! adapted from Frenkel and Smit :
!
! velocity autocorrelation function
!
! ******************************************************************************
SUBROUTINE vacf_main 

  USE config,           ONLY :  natm , rx , ry , rz , vx , vy , vz , rho 
  USE md,               ONLY :  dt
  USE io,               ONLY :  kunit_VACFFF
  USE constants,        ONLY :  pi
  USE time,             ONLY :  vacftimetot
  USE mpimdff

  implicit none 
  
  !INCLUDE 'mpif.h'

  ! local
  integer :: ia , dtt , t , ttel , ierr
  ! timeinfo
  real(kind=dp) :: ttt1 , ttt2 

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif

  tvacf = tvacf + 1
  ! ===============================================
  ! sample velocity auto correlation function and
  ! mean square displacement
  ! ===============================================
  if ( mod ( tvacf , it0 ) .eq. 0 ) then 
  ! ============
  !  new t=0
  ! ============
    t0 = t0 + 1
    ttel = mod ( t0 - 1 , t0max ) + 1
    ttv0 ( ttel ) = tvacf
    do ia = 1 , natm 
      vxt0 ( ia , ttel ) = vx ( ia )
      vyt0 ( ia , ttel ) = vy ( ia )
      vzt0 ( ia , ttel ) = vz ( ia )
    enddo 
  endif 
  do t = 1, MIN ( t0 , t0max )
    dtt = tvacf - ttv0( t ) + 1
    if ( dtt .lt. tmax .and. dtt * dtime .le. tdifmax) then 
      nvacf ( dtt ) = nvacf ( dtt ) + 1
      do ia = 1 , natm 
        vacff ( dtt ) = vacff ( dtt ) + &
        vx ( ia ) * vxt0 ( ia , t ) + vy ( ia ) * vyt0 ( ia , t ) + vz ( ia ) * vzt0 ( ia , t )
      enddo 
    endif 
  enddo 
 
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  vacftimetot = vacftimetot + ( ttt2 - ttt1 )
#endif

  return

END SUBROUTINE vacf_main

! *********************** SUBROUTINE vacf_write ********************************
!
! write vacf and fourier transform = > dos
!
! ******************************************************************************
SUBROUTINE vacf_write_output

  USE config,   ONLY :  natm 
  USE io,       ONLY :  ionode , stdout , kunit_VACFFF  
  USE time,     ONLY :  vacftimetot2
  USE mpimdff

  implicit none
  !INCLUDE 'mpif.h'

  ! local
  integer :: ihbmax , i , ierr
  real(kind=dp) :: vtime , dif
  real(kind=dp) :: thmax
  real(kind=dp) :: tauc , tau0 , errvacf 
  complex(kind=dp)   ,dimension (:), allocatable :: in , out 
  real(kind=dp) ,dimension (:), allocatable :: rout
  real(kind=dp) :: ttt1 , ttt2

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif
  
  allocate ( in ( tmax ) , out ( tmax ) , rout ( tmax ) ) 



  OPEN ( UNIT = kunit_VACFFF , FILE = 'VACFFF' ) 

  ! =======================================
  ! velocity auto-correlation function
  ! =======================================
  if ( tvacf .ne. 0 ) then 
    ! =======================================
    ! correlation time (for error estimate)
    ! =======================================
    tauc = 0
    do i = 1 , tmax
      if ( nvacf ( i ) * natm .ne. 0 ) then 
        tauc = tauc + ( vacff ( i ) / ( natm * nvacf ( i ) ) ) ** 2 * dtime
      endif
    enddo 
    ! =======================================
    ! normalisation
    ! =======================================
    if ( nvacf ( 1 ) * natm .ne. 0 ) tauc = tauc / ( vacff ( 1 ) / ( natm * nvacf ( 1 ) ) ** 2 )
    dif = 0
    ! =======================================
    ! total averaging time:
    ! =======================================
    tau0 = dtime*it0*t0
    errvacf = SQRT ( 2.0_dp * tauc * ( vacff ( 1 ) / ( natm * nvacf ( 1 ) ) ** 2 / tau0 ) )
    thmax = 0
    ihbmax = 0

    ! =======================
    !  normalisation
    ! =======================
    do i = 1 , tmax
     if ( nvacf ( i ) .ne. 0 ) then
       vacff ( i ) = vacff ( i ) / ( natm * nvacf ( i ) )
     endif
    enddo
    ! ========
    !   FFT
    ! ========
    in  = vacff
    CALL fft_1D_complex ( in , out , tmax )
    rout = DBLE ( out ) 

    do i = 1 , tmax
      vtime = dtime * ( i - 1 )
      if ( nvacf ( i ) .ne. 0 ) then 
        dif = dif + vacff ( i ) * dtime
        WRITE ( kunit_VACFFF , '(4e16.8)' ) vtime, vacff ( i ) , 1.0_dp / vtime , rout ( i ) 
        if ( vtime .gt. thmax ) then 
          ihbmax = nvacf ( i )
          thmax = vtime
        endif 
      endif 
    enddo 

    if ( ionode ) then
      WRITE ( stdout , 99002) tauc
      WRITE ( stdout , 99004 ) tvacf , t0 , dtime , dtime * it0 , dif/3.0_dp
      WRITE ( stdout , '(a)' ) 'Diffusion calculated with conventional scheme '
      WRITE ( stdout , 99001 ) 2 * dtime , nvacf(3) , thmax , ihbmax 
    endif
  endif 

  CLOSE (kunit_VACFFF)

  deallocate ( in , out , rout )

#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  vacftimetot2 = vacftimetot2 + ( ttt2 - ttt1 )
#endif

  RETURN

99001 FORMAT ('Number of samples for tmin = ', f8.3, ' is : ', i10, /, &
              'Number of samples for tmax = ', f8.3, ' is : ', i10)
99002 FORMAT ('Decorrelation time ', e12.4)
99004 FORMAT ('Velocity auto correlation function and mean square dis.:',&
             /, '   Number of samples       : ', i8, /, &
                '   Number of t=0           : ', i8, /, &
                '   Timestep between samples: ', f8.3, /, &
                '   Timestep between t=0    : ', f8.5, /, &
                '   Diffusion coef.         : ', f8.5)
 
 
END SUBROUTINE vacf_write_output


END MODULE vacf 
! ===== fmV =====
