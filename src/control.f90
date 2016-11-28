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

! *********************** MODULE control ***************************************
!> \brief
!! Module related to main parameters of the code
! ******************************************************************************
MODULE control 

  USE constants , ONLY : dp 
  USE mpimdff
      
  implicit none

  integer                 :: itraj_start      !< write trajectory from step itraj_start
  integer                 :: itraj_period     !< write trajectory each itraj_period steps 
  integer                 :: itraj_format     !< format of TRAJFF file ( = 0 BINARY, = 1 FORMATED)
  integer                 :: iscff_format     !< format of ISCFF  file      0 = BINARY
  integer                 :: iefall_format    !< format of EFALL  file  0 = BINARY
  integer                 :: iefgall_format   !< format of EFGALL file  0 = BINARY
  integer                 :: idipall_format   !< format of EFGALL file  0 = BINARY
  character(len=8),  SAVE :: DATE             !< execution DATE
  character(len=10), SAVE :: HOUR             !< execution HOUR

  logical,           SAVE :: ltraj            !< save trajectory                                    
  logical,           SAVE :: lstatic          !< no MD                                                
  logical,           SAVE :: lvnlist          !< verlet list if .true.                            
  logical,           SAVE :: lwrite_restart   !< control RESTART file 
  logical,           SAVE :: full_restart     !< full restart from RESTART file
  logical,           SAVE :: lreduced         !< print reduced units ( see reduced_units subroutine in constants.f90 )
  logical,           SAVE :: lreducedN        !< print reduced thermo quantites by the number of atoms (natm)
  logical,           SAVE :: lnmlj            !< n-m lennard-jones potential
  logical,           SAVE :: lbmhft           !< born huggins Mayer potential
  logical,           SAVE :: lbmhftd          !< born huggins Mayer potential + damping
  logical,           SAVE :: lharm            !< harmonic oscilaltor ( to test integration )
  logical,           SAVE :: lmorse           !< morse potential 
  logical,           SAVE :: lcoulomb         !< coulombic potential
  logical,           SAVE :: lsurf            !< add surface contribution in electrostatic quantities  
  logical,           SAVE :: ltest            !< testing flag
  logical,           SAVE :: lcsvr            !< Stochastic velocity rescaling
  logical,           SAVE :: lmsd             !< mean square displacement switch
  logical,           SAVE :: lvacf            !< velocity auto-correlation function switch

  real(kind=dp),     SAVE :: cutlongrange     !< longrange cutoff
  real(kind=dp),     SAVE :: cutshortrange    !< shortrange cutoff
  real(kind=dp),     SAVE :: skindiff         !< verlet-list cutoff ( cutoff + skindiff )
  TYPE(decomposition), SAVE :: kpt_dec


  ! =====================================================
  !   type of calculation : md, opt, vib, efg ...              
  ! =====================================================
  character(len=60), SAVE :: calc             !< type of calculation : md, opt, vib, efg ...  
  character(len=60), SAVE :: calc_allowed(14)    
  data calc_allowed / 'md'       , 'opt'     , 'vib'        , 'vib+fvib'       , 'vib+gmod' , &
                      'vib+band' , 'vib+dos' , 'efg'        , 'efg+acf'        , 'gr'       , &
                      'vois1'    , 'rmc'     , 'dist'       , 'stochio'  /

  ! =====================================================
  ! algorithm for long-range calculation
  ! =====================================================
  character(len=60), SAVE :: longrange        !< algorithm for long-range interaction  
  character(len=60), SAVE :: longrange_allowed(2) 
  data longrange_allowed / 'ewald' , 'direct' /

  ! =====================================================
  !   algorithm for gaussian distribution 
  !   (no fondamental difference just for fun)  
  ! =====================================================
  character(len=60), SAVE :: dgauss           !< algorithm for gaussian distribution
  character(len=60), SAVE :: dgauss_allowed(3) 
  data dgauss_allowed / 'boxmuller_basic', 'boxmuller_polar' , 'knuth' /

  ! =====================================================
  !  format of TRAJFF allowed  
  ! =====================================================
  character(len=3), SAVE :: trajff_data
  character(len=3), SAVE :: restart_data
  character(len=3), SAVE :: iscff_data
  character(len=3), SAVE :: data_allowed(4)
  data data_allowed / 'rvf' , 'rnn' , 'rnf' , 'rvn' /

  ! ==============================================================
  !  non-bonded potential if one of lnmlj, lbmhftd, lbmhftd,lmorse is true 
  ! ==============================================================
  logical, SAVE           :: non_bonded 
   
CONTAINS

! *********************** SUBROUTINE control_init ******************************
!
!> \brief
!! Initialization of control parameters.
!! Set default values, read  and check consistenstency of control parameters
!
! ******************************************************************************
SUBROUTINE control_init ( MDFF )

  USE io,  ONLY :  ionode , stdin , stdout

  implicit none

  ! local
  integer            :: ioerr
  character(len=80)  :: MDFF
  character(len=132) :: filename

  namelist /controltag/  lnmlj          , & ! #1
                         lcoulomb       , &
                         lmorse         , &
                         lbmhft         , &
                         lbmhftd        , &
                         lsurf          , &
                         lcsvr          , &
                         lharm          , &
                         ltraj          , &
                         lvnlist        , &
                         lstatic        , &
                         lreduced       , & 
                         lreducedN      , & 
                         ltest          , &
                         lmsd           , &
                         lvacf          , &
                         lwrite_restart , &
                         full_restart   , &
                         cutlongrange   , &
                         cutshortrange  , &
                         calc           , &
                         dgauss         , &
                         longrange      , &
                         itraj_start    , & 
                         itraj_period   , & 
                         itraj_format   , & 
                         trajff_data    , & 
                         iscff_format   , & 
                         iscff_data     , & 
                         iefall_format  , & 
                         iefgall_format , & 
                         idipall_format , & 
                         restart_data   , & 
                         skindiff     
               
  ! ======================
  !  set default values
  ! ======================
  CALL control_default_tag
  ! =======================
  !  read control namelist
  ! =======================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , controltag, iostat=ioerr)
  if ( ioerr .lt. 0 )  then
    io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : controltag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : controltag wrong tag',ioerr
    STOP    
  endif
  CLOSE ( stdin )
  ! ======================
  !      check tags
  ! ======================
  CALL control_check_tag
  ! ======================
  !  print info to output
  ! ======================
  if ( .not. full_restart ) CALL control_print_info( stdout , MDFF )

  return

END SUBROUTINE control_init


! *********************** SUBROUTINE control_default_tag ***********************
!
!> \brief
!! set default values to control tag
!
! ******************************************************************************
SUBROUTINE control_default_tag

  implicit none

  ! ================
  !  default values
  ! ================
  lnmlj         = .false.
  lbmhft        = .false.
  lbmhftd       = .false.
  lmorse        = .false.
  lcoulomb      = .false.
  lsurf         = .false.
  lcsvr         = .false.
  lharm         = .false.
  ltraj         = .false.
  lvnlist       = .true.
  lstatic       = .false.
  lreduced      = .false.
  lreducedN     = .false.
  ltest         = .false.
  lmsd          = .false.
  lvacf         = .false.
  lwrite_restart= .false.
  full_restart  = .false.
  calc          = 'md'
  dgauss        = 'boxmuller_basic'
  longrange     = 'ewald'
  skindiff      = 0.15_dp
  cutshortrange = 0.0_dp
  cutlongrange  = 0.0_dp
  itraj_start   = 1          
  itraj_period  = 10000
  itraj_format  = 1
  iscff_format  = 1
  iefall_format = 1
  iefgall_format= 1
  idipall_format= 1
  trajff_data   = 'rnn'
  iscff_data    = 'rnn'
  restart_data  = 'rnn'

  return 
 
END SUBROUTINE control_default_tag

! *********************** SUBROUTINE control_check_tag *************************
!
!> \brief
!! check control tag values
!
! ******************************************************************************
SUBROUTINE control_check_tag

  USE io,               ONLY :  stdout , ionode
  USE constants,        ONLY : reduced_units 

  implicit none

  ! local
  logical :: restart_file_exists
  logical :: allowed
  integer :: i


  ! ======
  !  calc
  ! ======
  do i = 1 , size( calc_allowed ) 
   if ( trim(calc) .eq. calc_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
      if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR controltag: calc should be ', calc_allowed
      STOP 
  endif
  allowed = .false.
  ! ===========
  !  longrange
  ! ===========
  do i = 1 , size( longrange_allowed )
   if ( trim(longrange) .eq. longrange_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: longrange should be ', longrange_allowed
  endif
  allowed = .false.
  ! =========
  !  dgauss
  ! =========
  do i = 1 , size( dgauss_allowed )
   if ( trim(dgauss) .eq. dgauss_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: dgauss should be ', dgauss_allowed
  endif
  ! =========
  ! trajff_data  
  ! =========
  allowed = .false.
  do i = 1 , size( data_allowed )
   if ( trim(trajff_data) .eq. data_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: trajf_data should be ', data_allowed
  endif
  allowed = .false.
  do i = 1 , size( data_allowed )
   if ( trim(iscff_data) .eq. data_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: iscff_data should be ', data_allowed
  endif
  allowed = .false.
  do i = 1 , size( data_allowed )
   if ( trim(restart_data) .eq. data_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR controltag: restart_data should be ', data_allowed
  endif

  if ( lnmlj .or. lmorse .or. lbmhftd .or. lbmhft ) then
    non_bonded = .true.
  endif

  if ( non_bonded .and. cutshortrange .eq. 0.0_dp ) then
    io_node WRITE ( stdout , '(a)' ) 'controltag: cutshortrange is null', cutshortrange
    STOP
  endif

  if ( lcoulomb .and. cutlongrange .eq. 0.0_dp ) then
    io_node WRITE ( stdout , '(a)' ) 'controltag: cutlongrange is null', cutlongrange
    STOP  
  endif
  
  if ( lreduced ) then
    CALL reduced_units 
  endif
  
  if ( calc .ne. 'md' ) return 

  if ( .not. lnmlj .and. .not. lcoulomb .and. .not. lmorse .and. .not. lharm .and. .not. lbmhftd .and. .not. lbmhft ) then
   if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR controltag: nmlj, harm , morse, bmhftd or coulomb or all of them . Anyway make a choice !! '
   STOP
  endif

  ! full restart
  !to be removed nov 2016
 ! INQUIRE(FILE="RESTART", EXIST=restart_file_exists)
 ! if ( lwrite_restart .and. restart_file_exists ) then
 !   full_restart = .true.
 ! endif


  return

END SUBROUTINE control_check_tag

! *********************** SUBROUTINE control_print_info ************************
!
!> \brief
!! print information to standard output about general parameters
!
!! \note
!! should have more information about other control parameters
!
! ******************************************************************************
SUBROUTINE control_print_info( kunit , MDFF )

  USE io,  ONLY :  ionode 
  USE mpimdff,  ONLY :  numprocs

  implicit none
 
  !local 
  integer :: kunit
  integer :: status
#ifdef GFORTRAN
  character(len=80) :: host 
#endif
  character(len=80) :: MDFF
  character(len=30) :: user_name
#ifndef HOST
#define HOST "unknown"
#endif
 
  ! ===============
  !  date time info      
  ! ===============
  CALL DATE_AND_TIME( DATE, HOUR)
  ! ===============
  !  user name 
  ! ===============
  CALL GETENV ( 'USER', user_name )

#ifdef GFORTRAN
  ! ===============
  !  hostname info      
  ! ===============
  status = hostnm(host)
#else   
  status = 0
#endif 
  ! =================
  !  standard output
  ! =================

  if ( .not. full_restart ) then
    if ( ionode ) call dumb_guy(kunit)
  endif
  if ( ionode ) then
  if ( .not. full_restart ) then
     WRITE ( kunit ,'(a)')       "          ____    ____  ______   ________  ________  "
     WRITE ( kunit ,'(a)')       "         |_   \  /   _||_   _ `.|_   __  ||_   __  | "
     WRITE ( kunit ,'(a)')       "           |   \/   |    | | `. \ | |_ \_|  | |_ \_| "
     WRITE ( kunit ,'(a)')       "           | |\  /| |    | |  | | |  _|     |  _|    "
     WRITE ( kunit ,'(a)')       "          _| |_\/_| |_  _| |_.' /_| |_     _| |_     "
     WRITE ( kunit ,'(a)')       "         |_____||_____||______.'|_____|   |_____|    "
  endif
     blankline(kunit)
     separator(kunit)
     blankline(kunit)
     WRITE ( kunit ,'(a)')       'MOLECULAR DYNAMICS ...for fun                 '
     WRITE ( kunit ,'(a)')       MDFF
     WRITE ( kunit ,'(a)')       'filipe.manuel.vasconcelos@gmail.com  '
     WRITE ( kunit ,'(a,i4,a)')  'Running on  : ',numprocs,' nodes                  '
     WRITE ( kunit ,'(a,a)')     'by user     : ',user_name
     if ( status == 0 ) then
#ifdef GFORTRAN
     WRITE ( kunit ,'(a,a)')     'host        : ',trim(host)
#else
     WRITE ( kunit ,'(a,a)')     'host        : ',HOST
#endif
     endif
     WRITE ( kunit ,'(a,a4,a1,a2,a1,a2,a4,a2,a1,a2,a1,a2)') &
                                'date        : ',DATE(1:4),'/',DATE(5:6),'/',DATE(7:8),'   ',HOUR(1:2),&
                                ':',HOUR(3:4),':',HOUR(5:6)
     blankline(kunit)
     WRITE ( kunit ,'(a)'  )     'CONTROL MODULE ... WELCOME'
     blankline(kunit)
     WRITE ( kunit ,'(a,a)')     'calc           =  ', calc 
     WRITE ( kunit ,'(a,l2)')    'lnmlj          = ', lnmlj 
     WRITE ( kunit ,'(a,l2)')    'lbmhft         = ', lbmhft
     WRITE ( kunit ,'(a,l2)')    'lbmhftd        = ', lbmhftd
     WRITE ( kunit ,'(a,l2)')    'lmorse         = ', lmorse 
     WRITE ( kunit ,'(a,l2)')    'lcoulomb       = ', lcoulomb
     WRITE ( kunit ,'(a,l2)')    'lsurf          = ', lsurf
     WRITE ( kunit ,'(a,l2)')    'lvnlist        = ', lvnlist
     WRITE ( kunit ,'(a,l2)')    'lstatic        = ', lstatic
     WRITE ( kunit ,'(a,l2)')    'lreduced       = ', lreduced 
     WRITE ( kunit ,'(a,l2)')    'lreducedN      = ', lreducedN 
     WRITE ( kunit ,'(a,l2)')    'lwrite_restart = ', lwrite_restart 
     WRITE ( kunit ,'(a,l2)')    'full_restart   = ', full_restart 
     if ( full_restart ) then
       WRITE ( kunit ,'(a)')     'restarting from RESTART file (full restart) '
     else 
       WRITE ( kunit ,'(a)')     'restart from positions and velocities only ( partial restart )' 
     endif
     
  endif

  return 
 
END SUBROUTINE control_print_info




SUBROUTINE warning_print_info

! __          __     _____  _   _ _____ _   _  _____ 
! \ \        / /\   |  __ \| \ | |_   _| \ | |/ ____|
!  \ \  /\  / /  \  | |__) |  \| | | | |  \| | |  __ 
!   \ \/  \/ / /\ \ |  _  /| . ` | | | | . ` | | |_ |
!    \  /\  / ____ \| | \ \| |\  |_| |_| |\  | |__| |
!     \/  \/_/    \_\_|  \_\_| \_|_____|_| \_|\_____|


END SUBROUTINE warning_print_info 



END MODULE control 
! ===== fmV =====
