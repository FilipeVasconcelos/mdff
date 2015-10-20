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
!#define debug
!#define debug_otf_frenkel
!#define debug_otf_multiwin
! ======= Hardware =======

! *********************** MODULE msd *******************************************
!
!> \brief
!! module related to mean square diplacement calculation
!
!> \note
!! see Dubbeldam et al. Molecular Simulation  v35 n12 (2009) p1084
!
! ******************************************************************************
MODULE msd

  USE constants,  ONLY : dp
  USE mpimdff 
  implicit none

  character(len=60), SAVE                         :: msdalgo           !< algorithm for sampling msd
  character(len=60), SAVE                         :: msdalgo_allowed(3)
  data msdalgo_allowed / 'otf_frenkel', 'otf_multwin' , 'pp_win' /

  integer                                         :: nblcel   ! max number of block elements
  integer                                         :: nblocks  ! max number of block

  TYPE :: block_correlation
    real(kind=dp) , dimension (:,:,:) , allocatable :: x, y , z  ! vector of data 
    real(kind=dp) , dimension (:,:,:)   , allocatable :: aver        ! averaging
    integer       , dimension (:,:,:)   , allocatable :: cnt         ! count
    integer       , dimension (:)     , allocatable :: ibl         ! block length for each block
    integer                                         :: n           ! number of blocks:
  END TYPE

  TYPE ( block_correlation ) :: msd_data 

!  real(kind=dp) , dimension (:,:,:) , allocatable ::  msd_data%x , msd_data%y , msd_data%z
!  real(kind=dp) , dimension (:,:) , allocatable   ::  msd_data%aver   ! msd average
!  integer ,dimension (:,:) , allocatable          :: msd_data%cnt  ! msd count
!  integer ,dimension (:)   , allocatable          :: ibl    ! block length

  real(kind=dp) :: tcormax  !< correlation time maximum
  real(kind=dp) :: dtime    !< correlation time minimum ( dt * npropr )

CONTAINS

! *********************** SUBROUTINE msd_init **********************************
!
!> \brief
!! mean square diplacement initialisation
!
! ******************************************************************************
SUBROUTINE msd_init

  USE control,          ONLY :  lmsd
  USE io,               ONLY :  stdin, stdout, ionode

  implicit none

  ! local
  integer            :: ioerr
  character(len=132) :: filename


  namelist /msdtag/  nblcel  , &
                     nblocks , &
                     msdalgo , &
                     tcormax 

  if ( .not. lmsd ) return

  CALL msd_default_tag
  
  ! ==================
  !  read msdtag tags
  ! ==================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , msdtag, iostat=ioerr)
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : msdtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : msdtag wrong tag'
   STOP
  endif

  CLOSE ( stdin )

  CALL msd_check_tag

  CALL msd_print_info(stdout)

  return 
 
END SUBROUTINE msd_init


! *********************** SUBROUTINE msd_default_tag ***************************
!
!> \brief
!! set default values to msd tags
!
! ******************************************************************************
SUBROUTINE msd_default_tag

  implicit none

  msdalgo = 'otf_frenkel'
  tcormax = 100.0_dp
  nblcel  = 10
  nblocks   = 20

  return 
 
END SUBROUTINE msd_default_tag



! *********************** SUBROUTINE msd_check_tag *****************************
!
!> \brief
!! check msd tag values
!
! ******************************************************************************
SUBROUTINE msd_check_tag

  USE constants,        ONLY :  time_unit 
  USE io,               ONLY :  ionode, stdout

  implicit none

  ! local
  logical :: allowed
  integer :: i


  ! ======
  ! msdalgo
  ! ======
  do i = 1 , size( msdalgo_allowed )
   if ( trim(msdalgo) .eq. msdalgo_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
      if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR msdtag: msdalgo should be ', msdalgo_allowed
      STOP
  endif
  allowed = .false.


  tcormax = tcormax * time_unit

  return 
 
END SUBROUTINE msd_check_tag


! *********************** SUBROUTINE msd_print_info ****************************
!
!> \brief
!! print msd information to standard output
!
! ******************************************************************************
SUBROUTINE msd_print_info(kunit)

  USE io,  ONLY :  ionode 

  implicit none
  
  !local
  integer :: kunit

  if ( ionode  ) then
    blankline(kunit)
    WRITE ( kunit ,'(a)')           'mean square displacement:'
    if ( msdalgo .eq. 'otf_frenkel' ) then
      WRITE ( kunit ,'(a)')         'using on-the-fly order-n frenkel-smit algorithm'
    endif
    if ( msdalgo .eq. 'otf_multwin' ) then
      WRITE ( kunit ,'(a)')         'using on-the-fly multi-window algorithm '
      WRITE ( kunit ,'(a)')         'Dubbeldam et al. Molecular Simulation  v35 n12 (2009) p1084'
    endif
    WRITE ( kunit ,'(a,f10.4)')     'maximum time correlation     (tcormax)  = ',tcormax
    WRITE ( kunit ,'(a,i10)')       'max number of blocks         (nblocks)  = ',nblocks
    WRITE ( kunit ,'(a,i10)')       'max number of block elements (nblcel) = ',nblcel
    WRITE ( kunit ,'(a)')           'output file                             : MSDFF'
  endif

  return 
 
END SUBROUTINE msd_print_info


! *********************** SUBROUTINE msd_alloc *********************************
!
!> \brief
!! allocate and initialize variables
!
! ******************************************************************************
SUBROUTINE msd_alloc

  USE control,          ONLY :  lmsd
  USE md,               ONLY :  dt , npropr
  USE config,           ONLY :  natm, ntype

  implicit none

  if ( .not. lmsd ) return

  allocate ( msd_data%x ( nblocks , nblcel , natm ) , msd_data%y ( nblocks , nblcel , natm )  , msd_data%z ( nblocks , nblcel , natm )  ) 
  allocate ( msd_data%aver ( nblocks , nblcel , 0:ntype ) )
  allocate ( msd_data%cnt ( nblocks , nblcel , 0:ntype ) , msd_data%ibl ( nblocks ) )

  dtime = npropr * dt
  
  msd_data % ibl  = 1 
  msd_data % cnt  = 0
  msd_data % aver = 0.0_dp
  msd_data % x    = 0.0_dp
  msd_data % y    = 0.0_dp
  msd_data % z    = 0.0_dp

  return 
 
END SUBROUTINE msd_alloc


! *********************** SUBROUTINE msd_dealloc *******************************
!
!> \brief
!! deallocate variables
!
! ******************************************************************************
SUBROUTINE msd_dealloc

  USE control,          ONLY :  lmsd

  implicit none

  if ( .not. lmsd ) return

    deallocate ( msd_data % x , msd_data % y  , msd_data % z  )
    deallocate ( msd_data % aver )
    deallocate ( msd_data % cnt , msd_data % ibl )

  return 
 
END SUBROUTINE msd_dealloc

SUBROUTINE msd_sample ( nmsd )

  implicit none

  ! global
  integer :: nmsd

  if ( msdalgo .eq. 'otf_frenkel' ) CALL msd_onthefly_frenkel     ( nmsd ) 
  if ( msdalgo .eq. 'otf_multwin' ) CALL msd_onthefly_multiwindow ( nmsd ) 

  return

END SUBROUTINE msd_sample


! *********************** SUBROUTINE msd_onthefly_frenkel *******************
!
!> \brief
!! Determine the mean square displacement using Algorithm 9
!
!> \note
!! order-n algorithm based on the velocities
!! Adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE msd_onthefly_frenkel ( nmsd )

  USE constants,        ONLY :  time_unit 
  USE config,           ONLY :  natm , vx , vy , vz , itype
  USE md,               ONLY :  npropr
  USE time,             ONLY :  msdtimetot
  USE io,               ONLY :  stdout , ionode , ioprintnode     

  implicit none

  ! global
  integer :: nmsd

  ! local
  integer :: ia , iblock, ib, ie,  inp, ii, current_block_length, ierr , iblm , it
  real(kind=dp) :: delx, dely, delz, xtime
#ifdef debug 
  real(kind=dp) :: r2asum
#endif
  ! timeinfo
  real(kind=dp) :: ttt1 , ttt2 

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif
 
  ! ===================================================
  !  determine current maximum number of blocks: iblm
  ! ===================================================
  iblm =  1
  ii   = nmsd / nblcel
  do while  (ii.ne.0)
    iblm = iblm+ 1
    ii   = ii / nblcel
  ! ===========================================
  !  test maximu time not longer than tcormax:
  ! ===========================================
    xtime = dtime * ( nblcel ** ( iblm ) )
    if ( xtime .gt. tcormax ) ii = 0
  enddo
  msd_data % n = iblm

#ifdef debug_otf_frenkel
  if ( ioprintnode ) then
    WRITE ( stdout ,'(a)')                  ' ------------------------------------------------------------------------------'
    WRITE ( stdout ,'(a)')                  '   msd on-the-fly : order -n frenkel and smit ' 
    WRITE ( stdout ,'(a)')                  ' ------------------------------------------------------------------------------'
    WRITE ( stdout ,'(a,i12)')              ' number of elements/block    : ', nblcel 
    WRITE ( stdout ,'(a,i12)')              ' property calc. period       : ', npropr
    WRITE ( stdout ,'(a,f12.4)')            ' dt min                      : ', dtime/time_unit
    WRITE ( stdout ,'(a,i12,a6,i12,a)')     ' number of blocks      (max) : ', msd_data % n,           '(',nblocks,')'
    WRITE ( stdout ,'(a,f12.3,a6,f12.3,a)') ' current time interval (max) : ', xtime/time_unit,'(',tcormax/time_unit,')'
  endif
#endif

  ! =====================================
  !  limit the maximum number of blocks ! well needed only if memory is an issue
  ! =====================================
  msd_data%n=MIN( msd_data%n , nblocks )
  blocks : do ib = 1, msd_data % n 
    iblock = nblcel ** ( ib - 1 )
    ! ==============================
    !  test for blocking operation
    ! ==============================
    if ( MOD ( nmsd , iblock ) .ne. 0 ) cycle

    msd_data%ibl ( ib ) = msd_data%ibl ( ib ) + 1
    ! ==============================
    !  limit to length n (=nblcel)
    ! ==============================
    current_block_length = MIN( msd_data%ibl ( ib ) , nblcel )
#ifdef debug
     r2asum = 0
#endif

#ifdef debug_otf_frenkel
    ! print info
    if ( ioprintnode ) then
      WRITE ( stdout ,'(a)')                         ' ------------------------------------------------------------------------------'
      WRITE ( stdout ,'(a,i5,a,i5)')                 ' block : ' , ib , ' /', msd_data%n 
      WRITE ( stdout ,'(a,i5,a,i5,a,f12.5)')         ' block : ' , ib , ' /', msd_data%n, ' time : ' , dtime / time_unit * ( nblcel ** ( ib - 1 ) )
    endif
#endif

    ions : do ia = 1 , natm 

      it = itype ( ia ) ! type speciation 

      if (ib.eq.1) then
        ! ================================
        !  zero block: ordinary velocity
        ! ================================
        delx = vx ( ia )
        dely = vy ( ia )
        delz = vz ( ia )
      ELSE
        ! ===============================================
        !  (ib)th block: coarsed velocity previous block
        ! ===============================================
        delx = msd_data%x ( ib-1 , 1 , ia )
        dely = msd_data%y ( ib-1 , 1 , ia )
        delz = msd_data%z ( ib-1 , 1 , ia )
      endif
      do ie = 1, current_block_length
        inp = ie
        if ( msd_data%ibl ( ib ) .gt. nblcel ) inp = ie + 1
        if ( ie .lt. current_block_length ) then
          msd_data%x ( ib , ie , ia ) = msd_data%x ( ib , inp , ia ) + delx
          msd_data%y ( ib , ie , ia ) = msd_data%y ( ib , inp , ia ) + dely
          msd_data%z ( ib , ie , ia ) = msd_data%z ( ib , inp , ia ) + delz
        else
          msd_data%x ( ib , ie , ia ) = delx
          msd_data%y ( ib , ie , ia ) = dely
          msd_data%z ( ib , ie , ia ) = delz
        endif
      enddo
      do ie = 1, current_block_length
        msd_data%cnt  ( ib , ie , 0  ) = msd_data%cnt  ( ib , ie , 0  ) + 1
        msd_data%cnt  ( ib , ie , it ) = msd_data%cnt  ( ib , ie , it ) + 1

        msd_data%aver ( ib , ie , 0 ) = msd_data%aver ( ib , ie , 0 ) &
                                               + ( msd_data%x ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp &
                                               + ( msd_data%y ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp &
                                               + ( msd_data%z ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp

        msd_data%aver ( ib , ie , it ) = msd_data%aver ( ib , ie , it ) &
                                               + ( msd_data%x ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp &
                                               + ( msd_data%y ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp &
                                               + ( msd_data%z ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp

#ifdef debug
        !  print*,'debug'
        if ( ie .eq. 1 ) r2asum = r2asum       + ( msd_data%x ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp &
                                               + ( msd_data%y ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp &
                                               + ( msd_data%z ( ib , current_block_length - ie + 1 , ia ) * dtime ) ** 2.0_dp
#endif
      enddo

    enddo ions

#ifdef debug
      r2asum = r2asum / natm
      ! ============================================================
      !  print mean square displacement to file for t=1,10,100,etc
      ! ============================================================
      WRITE ( 79 + ib , *) 'debug: ',dtime * ( nblcel ** ( ib - 1 ) ) / time_unit , r2asum
#endif
  enddo blocks
 
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  msdtimetot = msdtimetot + ( ttt2 - ttt1 )
#endif
 
  return

END SUBROUTINE msd_onthefly_frenkel 


! *********************** SUBROUTINE msd_postprocess_window **************************
!
!> \brief
!! write results to file MSDFF
!
! ******************************************************************************
!SUBROUTINE msd_postprocess_window ( nmsd )
!
!  implicit none
!
!  ! global
!  integer :: nmsd
!
!  return
!
!END SUBROUTINE msd_postprocess_window

! *********************** SUBROUTINE msd_onthefly_multiwindow **********************
!
!> \brief
!! write results to file MSDFF
!! use positions maybe is an issue with periodic boundary conditions
!
! ******************************************************************************
SUBROUTINE msd_onthefly_multiwindow ( nmsd ) 

  USE constants,        ONLY : time_unit
  USE config,           ONLY : natm , rx, ry ,rz , itype
  USE md,               ONLY :  npropr
  USE io,               ONLY : stdout, ioprintnode

  implicit none

  ! global
  integer :: nmsd

  ! local
  integer :: ia , ie , ib , ii , it , iblock, current_block_length , iblm
  integer :: index_origin , index_t , key

  real(kind=dp ) :: x0 , y0 , z0 ! origin 

  ! ===================================================
  !  determine current maximum number of blocks: iblm
  ! ===================================================
  iblm= 1
  ii = nmsd / nblcel
  do while  (ii.ne.0)
    iblm= iblm+1
    ii = ii / nblcel
  enddo
  msd_data%n=iblm

#ifdef debug_otf_multiwin
  if ( ioprintnode ) then
    WRITE ( stdout ,'(a)')                   ' ------------------------------------------------------------------------------'
    WRITE ( stdout ,'(a)')                   '                         msd on-the-fly : multi window              ' 
    WRITE ( stdout ,'(a)')                   ' ------------------------------------------------------------------------------'
    WRITE ( stdout ,'(a,i12)')               ' property calc. period             : ', npropr
    WRITE ( stdout ,'(a,f12.4,a6,f12.4,a)')  ' Dt correlation min          (max) : ', dtime/time_unit,'(',tcormax,')'
    WRITE ( stdout ,'(a,i12)')               ' maximum number of blocks          : ', nblocks
    WRITE ( stdout ,'(a,i12)')               ' maximum number of elements/blocks : ', nblcel
    do ib = 1 , msd_data%n
      WRITE ( stdout ,'(a,i5,a,i5)')         ' block : ',ib, '    filled elem. : ', msd_data%ibl(ib)
    enddo 
  endif
#endif

  ! ============================================================
  ! loop over all the blocks to test wich blocks need sampling
  ! ============================================================
  blocks : do ib = 1 , msd_data%n 

    key = 0 ! just for printing not fondamental 

    iblock = nblcel ** ( ib - 1 )
    ! ==============================
    !  test for blocking operation
    ! ==============================
    if ( MOD ( nmsd , iblock ) .ne. 0 ) cycle

    ! ==============================
    !  compute the current length of the block 
    !  current_block_length 
    !  limit to length n (=nblcel)
    ! ==============================
    current_block_length = MIN( msd_data%ibl ( ib ) , nblcel )

#ifdef debug_otf_multiwin
    ! print info
    if ( ioprintnode ) then
      WRITE ( stdout ,'(a)')                         ' ------------------------------------------------------------------------------'
      WRITE ( stdout ,'(a,i5,a,i5,a,f12.5)')         ' block : ' , ib , ' /', msd_data%n, ' time : ' , dtime / time_unit * ( nblcel ** ( ib - 1 ) )
    endif
#endif

    ! =================
    ! loop over ions 
    ! =================
    ions : do ia = 1 , natm 
      
      it = itype ( ia ) 

      index_origin = nblcel-current_block_length+1
#ifdef debug_otf_multiwin
      if ( ioprintnode .and. ia .eq. 1 ) write( stdout , '(a,i5)') ' index_origin     = ',index_origin
#endif
      ! shift to the left, and set last index to the correlation value
      do ie= index_origin+1 , nblcel
#ifdef debug_otf_multiwin
          if ( ioprintnode .and. ia .eq. 1 ) write( stdout , '(a,i5,a,i5,a,i5)') ' shift on ', current_block_length , ' elements : ',ie-1,' <= ',ie
#endif
          msd_data%x ( ib , ie-1 , ia ) = msd_data%x ( ib , ie , ia )
          msd_data%y ( ib , ie-1 , ia ) = msd_data%y ( ib , ie , ia )
          msd_data%z ( ib , ie-1 , ia ) = msd_data%z ( ib , ie , ia ) 
      enddo
      msd_data%x ( ib , nblcel , ia ) = rx ( ia )
      msd_data%y ( ib , nblcel , ia ) = ry ( ia )
      msd_data%z ( ib , nblcel , ia ) = rz ( ia ) 
#ifdef debug_otf_multiwin
      if ( ioprintnode .and. ia .eq. 1 ) write( stdout , '(a,i5)') ' store new positions to = ',nblcel 
      if ( ioprintnode .and. ia .eq. 1 ) write( stdout , '(a)')    ''
#endif
      ! get the origin, take into account that blocks can be partially filled
      x0 = msd_data%x ( ib , index_origin , ia )  
      y0 = msd_data%y ( ib , index_origin , ia )  
      z0 = msd_data%z ( ib , index_origin , ia )  
      ! sample msd using proper reference position
      block_element : do ie = 1 , current_block_length
        index_t = index_origin + ie - 1 

        msd_data%cnt  ( ib , ie , 0 )  = msd_data%cnt  ( ib , ie , 0 ) + 1
        msd_data%aver ( ib , ie , 0 )  = msd_data%aver ( ib , ie , 0 ) &
                                               + ( msd_data%x ( ib , index_t , ia ) - x0 ) ** 2.0_dp &
                                               + ( msd_data%y ( ib , index_t , ia ) - y0 ) ** 2.0_dp &
                                               + ( msd_data%z ( ib , index_t , ia ) - z0 ) ** 2.0_dp

        msd_data%cnt  ( ib , ie , it ) = msd_data%cnt  ( ib , ie , it ) + 1
        msd_data%aver ( ib , ie , it ) = msd_data%aver ( ib , ie , it ) &
                                               + ( msd_data%x ( ib , index_t , ia ) - x0 ) ** 2.0_dp &
                                               + ( msd_data%y ( ib , index_t , ia ) - y0 ) ** 2.0_dp &
                                               + ( msd_data%z ( ib , index_t , ia ) - z0 ) ** 2.0_dp

#ifdef debug_otf_multiwin
        if ( ioprintnode ) then
          if ( ia .eq. 1 .and. ie .eq. 1 ) write ( stdout , '(a)' ) '     ia      ib     ie_max   ie     r(0)    r(t)           aver        aver_tot'

          write(stdout,'(6i8,2e16.8)') ia,ib,current_block_length,ie,index_origin,index_t , &
                                                                        ( msd_data%x ( ib , index_t , ia ) - x0 ) ** 2.0_dp & 
                                                                      + ( msd_data%y ( ib , index_t , ia ) - y0 ) ** 2.0_dp &
                                                                      + ( msd_data%z ( ib , index_t , ia ) - z0 ) ** 2.0_dp , msd_data%aver ( ib , ie , 0 )
 
        endif
#endif
      enddo block_element

    enddo ions

#ifdef debug_otf_multiwin
      if ( ioprintnode ) then
        write(stdout,'(a)') '      ib      ie    ie_max    aver(ib,ie)    cnt(ib,ie)'
        do ie = 1, current_block_length
          write(stdout,'(3i8,f16.6,i14)') ib,ie,current_block_length,msd_data%aver ( ib , ie , 0 ) ,  msd_data%cnt ( ib , ie , 0 )
        enddo
      endif
#endif

    if ( msd_data%ibl(ib) .lt. nblcel ) msd_data%ibl(ib) = msd_data%ibl(ib)+1

  enddo blocks

  return

END SUBROUTINE msd_onthefly_multiwindow



! *********************** SUBROUTINE msd_write_output **************************
!
!> \brief
!! write results to file MSDFF
!
! ******************************************************************************
SUBROUTINE msd_write_output 

  USE config,           ONLY :  ntype, natm
  USE constants,        ONLY :  time_unit
  USE io,               ONLY :  ionode , kunit_MSDFF , stdout
  USE md,               ONLY :  npropr

  implicit none

  ! local
  integer :: j, ib, ihbmax , it , ie
  real(kind=dp) :: thmax

  if ( ionode ) then
    OPEN( kunit_MSDFF , file = 'MSDFF' ) 
    WRITE ( kunit_MSDFF, 99003 ) nblcel , nblocks , tcormax
    thmax = 0
    ihbmax = 0
    do ib = 1 , MIN( nblocks , msd_data % n  )
      if ( msdalgo .eq. 'otf_frenkel' ) then
        do j = 2, MIN( msd_data%ibl(ib), nblcel )
          if ( msd_data%cnt ( ib , j , 0 ) .lt. nblcel*natm ) cycle 
          WRITE ( kunit_MSDFF , 99002 )  &
          REAL ( j , kind = dp ) * dtime / time_unit * ( nblcel ** ( ib - 1 ) ), &
          ( msd_data%aver ( ib , j , it ) / REAL ( msd_data%cnt ( ib , j , it ) , kind = dp ), it = 0,ntype ) , msd_data%cnt ( ib , j , 0 )
        enddo
      endif
    
      if ( msdalgo .eq. 'otf_multwin' ) then
        do j = 2 , MIN( msd_data%ibl(ib), nblcel ) 
          if ( msd_data%cnt ( ib , j , 0 ) .lt. nblcel*natm ) cycle 
          WRITE ( kunit_MSDFF , 99002 )  &
          REAL ( j-1 , kind = dp ) * dtime / time_unit * ( nblcel ** ( ib - 1 ) ), &
          ( msd_data%aver ( ib , j , it ) / REAL ( msd_data%cnt ( ib , j , it ) , kind = dp ), it = 0,ntype ) , msd_data%cnt ( ib , j , 0 )
        enddo    
      endif
    enddo

#ifdef debug
     do ib=1, msd_data % n 
       WRITE(100000,'(a,i16)') 'ib : ',ib
       WRITE(100000,'(<nblcel>e16.8)') ( msd_data%aver ( ib , ie , 0 ) , ie = 1 , size ( msd_data%aver ( ib , : , 0 ) ) )
       WRITE(100000,'(<nblcel>i16)')   ( msd_data%cnt  ( ib , ie , 0 ) , ie = 1 , size ( msd_data%aver ( ib , : , 0 ) ) )
     enddo
#endif

    CLOSE ( kunit_MSDFF ) 
  endif

  return 
#ifdef GFORTRAN
99002 FORMAT ( 2x,e16.8,8(2x,e16.8),2x,i10 )
#else
99002 FORMAT ( 2x,e16.8,<ntype+1>(2x,e16.8),2x,i10 )
#endif
99003 FORMAT ('# nblcel =', i6, ' nblocks = ', i6 , ' tcormax = ', f10.4 )

END SUBROUTINE msd_write_output

END MODULE msd
! ===== fmV =====
