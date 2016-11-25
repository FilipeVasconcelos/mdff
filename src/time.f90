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

! *********************** MODULE time ******************************************
!> @brief
!> this MODULE implements a small timer utility to time
!> (see START_TIMING)
!> @author FMV
! ******************************************************************************
MODULE time

  USE constants, ONLY : dp

  implicit none
  
  real(kind=dp) :: timetot         !< total wall time        
  real(kind=dp) :: forcetimetot    !< engforce wall time       
  real(kind=dp) :: fcoultimetot1   !< engforce_coul direct space wall time     
  real(kind=dp) :: fcoultimetot2   !< engforce_coul fourier space wall time     
  real(kind=dp) :: fcoultimetot2_2 !< engforce_coul fourier space wall time     
  real(kind=dp) :: fcoultimetot3   !< engforce_coul wall time (comm only )
  real(kind=dp) :: efgtimetot1    !< EFG direct space wall time          
  real(kind=dp) :: efgtimetot2    !< EFG fourier space wall time        
  real(kind=dp) :: efgtimetot3    !< communication only wall time         
  real(kind=dp) :: rhoktimetot    !< facteur de structure wall time
  real(kind=dp) :: mdsteptimetot  !<                            
  real(kind=dp) :: vacftimetot    !<                            
  real(kind=dp) :: vacftimetot2   !<                            
  real(kind=dp) :: vnlisttimetot  !<                            
  real(kind=dp) :: msdtimetot     !<                            
  real(kind=dp) :: grtimetot      !<                            
  real(kind=dp) :: grtimetot_comm !<                            
  real(kind=dp) :: rmcgrtimetot_comm
  real(kind=dp) :: opttimetot     !<                          
  real(kind=dp) :: vibtimetot     !<                           
  real(kind=dp) :: hessiantimetot !<                           
  real(kind=dp) :: diaghessiantimetot !<                           
  real(kind=dp) :: bandtimetot    !<                           
  real(kind=dp) :: doskpttimetot  !<                           
  real(kind=dp) :: fvibtimetot    !<                           
  real(kind=dp) :: chisqvartimetot
  real(kind=dp) :: chisqvartimetotu
  real(kind=dp) :: chisqvartime2
  real(kind=dp) :: chisqvartime2u
  real(kind=dp) :: chisqvartime3
  real(kind=dp) :: chisqvartime3u
  real(kind=dp) :: chisqvartime4
  real(kind=dp) :: chisqvartime4u
  real(kind=dp) :: chisqvartime5
  real(kind=dp) :: kardirtottime 
  real(kind=dp) :: dirkartottime 
  real(kind=dp) :: time_moment_from_pola
  real(kind=dp) :: integratimetot
  real(kind=dp) :: timetmp        !< use in opt               

CONTAINS


!*********************** SUBROUTINE time_init *********************************
!> @brief
!> set timers to zero
!> @author FMV
!******************************************************************************
SUBROUTINE time_init

  implicit none

  timetot       = 0.0_dp
  efgtimetot1   = 0.0_dp
  efgtimetot2   = 0.0_dp
  efgtimetot3   = 0.0_dp
  rhoktimetot   = 0.0_dp 
  mdsteptimetot = 0.0_dp
  vnlisttimetot = 0.0_dp
  msdtimetot    = 0.0_dp
  grtimetot     = 0.0_dp 
  grtimetot_comm= 0.0_dp 
  forcetimetot  = 0.0_dp
  fcoultimetot1 = 0.0_dp
  fcoultimetot2 = 0.0_dp
  fcoultimetot2_2 = 0.0_dp
  fcoultimetot3 = 0.0_dp
  vibtimetot    = 0.0_dp
  hessiantimetot= 0.0_dp
  diaghessiantimetot = 0.0_dp
  bandtimetot   = 0.0_dp
  doskpttimetot = 0.0_dp
  opttimetot    = 0.0_dp
  vacftimetot   = 0.0_dp                             
  vacftimetot2  = 0.0_dp                             
  chisqvartimetot  = 0.0_dp
  chisqvartimetotu  = 0.0_dp
  chisqvartime2 = 0.0_dp
  chisqvartime2u = 0.0_dp
  chisqvartime3 = 0.0_dp
  chisqvartime3u = 0.0_dp
  chisqvartime4 = 0.0_dp
  chisqvartime4u = 0.0_dp
  chisqvartime5 = 0.0_dp
  kardirtottime = 0.0_dp
  dirkartottime = 0.0_dp
  time_moment_from_pola = 0.0_dp
  integratimetot = 0.0_dp

  return 
 
END SUBROUTINE time_init

! *********************** SUBROUTINE print_time_info ***************************
!> @author FMV
!> @brief
!> print final time information
!> @param[in] kunit unit file 
! ******************************************************************************
SUBROUTINE print_time_info ( kunit )  

  USE io,       ONLY :  ionode 
  USE control,  ONLY :  longrange , calc
  USE dumb

  implicit none

  !global
  integer, intent (in) :: kunit

  ! local
  real(kind=dp) :: dstime , estime
  real(kind=dp) :: fcoultime_es , fcoultime_ds

  dstime = efgtimetot1 + efgtimetot3               ! direct
  estime = rhoktimetot + efgtimetot1 + efgtimetot2 + efgtimetot3 ! ewald
  fcoultime_es = fcoultimetot1 + fcoultimetot2 + fcoultimetot3 + rhoktimetot
  fcoultime_ds = fcoultimetot1 + fcoultimetot3 

  !timing information at the end
  if ( ionode ) then
                                        separator(kunit)    
                                        blankline(kunit)
                                        WRITE ( kunit , '(a)' )         'TIME MODULE ... WELCOME'
                                        blankline(kunit)
                                        lseparator_noionode(kunit) 
                                        WRITE ( kunit ,110)             'TOTAL'  , timetot
    if ( mdsteptimetot .ne. 0.0_dp )    WRITE ( kunit ,110)             'MD'     , mdsteptimetot
    if ( opttimetot    .ne. 0.0_dp )    WRITE ( kunit ,110)             'OPT'    , opttimetot
    if ( vibtimetot    .ne. 0.0_dp )    WRITE ( kunit ,110)             'VIB'    , vibtimetot
    if ( fvibtimetot   .ne. 0.0_dp )    WRITE ( kunit ,110)             'FVIB'   , fvibtimetot
    if ( dstime .ne.0.0_dp .and. &
         efgtimetot2   .eq. 0.0_dp )    WRITE ( kunit ,110)             'EFG_DS' , dstime
    if ( efgtimetot2   .ne. 0.0_dp )    WRITE ( kunit ,110)             'EFG_ES' , estime
    if ( fcoultime_ds  .ne. 0.0_dp .and. &
         fcoultimetot2 .eq. 0.0_dp )    WRITE ( kunit ,110)             'FIELD_DS', fcoultime_ds
    if ( fcoultimetot2 .ne. 0.0_dp )    WRITE ( kunit ,110)             'FIELD_ES', fcoultime_es
    if ( kardirtottime .ne. 0.0_dp )    WRITE ( kunit ,110)             'kardir'  , kardirtottime
    if ( dirkartottime .ne. 0.0_dp )    WRITE ( kunit ,110)             'dirkar'  , dirkartottime
    if ( grtimetot     .ne. 0.0_dp )    WRITE ( kunit ,110)             'GR'                   , grtimetot 
    if ( grtimetot_comm.ne. 0.0_dp )    WRITE ( kunit ,110)             'GR(comm only)'        , grtimetot_comm 
    if ( chisqvartimetot.ne. 0.0_dp )    WRITE ( kunit ,110)            'eval_INVERT'        , chisqvartimetot 
    if ( chisqvartime2 .ne. 0.0_dp )    WRITE ( kunit ,110)             '       - dab '        , chisqvartime2
    if ( chisqvartime3 .ne. 0.0_dp )    WRITE ( kunit ,110)             '       - merge/shift' , chisqvartime3
    if ( chisqvartime4 .ne. 0.0_dp )    WRITE ( kunit ,110)             '       - average'     , chisqvartime4
    if ( chisqvartime5 .ne. 0.0_dp )    WRITE ( kunit ,110)             '       - chisq_var'   , chisqvartime5
    if ( chisqvartimetotu.ne. 0.0_dp )    WRITE ( kunit ,110)            'eval_INVERT_update'        , chisqvartimetotu 
    if ( chisqvartime2u .ne. 0.0_dp )    WRITE ( kunit ,110)             '      - dab_update '        , chisqvartime2u
    if ( chisqvartime3u .ne. 0.0_dp )    WRITE ( kunit ,110)             '       - merge/shift/average' , chisqvartime3u
    if ( chisqvartime4u .ne. 0.0_dp )    WRITE ( kunit ,110)             '       -chisq_var_update '     , chisqvartime4u
    if ( msdtimetot    .ne. 0.0_dp )    WRITE ( kunit ,110)             'MSD'    , msdtimetot
    if ( vacftimetot   .ne. 0.0_dp )    WRITE ( kunit ,110)             'VACF'   , vacftimetot
    if ( vacftimetot2  .ne. 0.0_dp )    WRITE ( kunit ,110)             'VACF2'  , vacftimetot2
                                        blankline(kunit)
                                        lseparator_noionode(kunit) 
                                        WRITE ( kunit ,'(a)')           'main subroutines:'
                                        lseparator_noionode(kunit) 
    if ( calc .eq. 'md' ) then
                                        WRITE ( kunit ,'(a)')           'MD:'
                                        lseparator_noionode(kunit) 
      if ( forcetimetot   .ne. 0.0_dp )  WRITE ( kunit ,110)             'engforce_nmlj      ' , forcetimetot        
      if ( vnlisttimetot  .ne. 0.0_dp )  WRITE ( kunit ,110)             'vnlistcheck        ' , vnlisttimetot
      if ( integratimetot .ne. 0.0_dp )  WRITE ( kunit ,110)             'integrator         ' , integratimetot
                                        lseparator_noionode(kunit) 
    endif
    
    if ( longrange .eq. 'direct' )   then
      if ( dstime        .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           'EFG:'
                                        lseparator_noionode(kunit) 
      if ( efgtimetot1   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_DS(efg  only)       ', efgtimetot1 
      if ( efgtimetot3   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_DS(comm only)       ', efgtimetot3
      if ( fcoultimetot1 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_DS            ', fcoultimetot1
      if ( fcoultimetot3 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_DS(comm only) ', fcoultimetot3
    endif
    if ( longrange .eq. 'ewald')   then
      if ( calc .eq. 'efg' ) then
        if ( estime        .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           'EFG:'
        if ( estime        .ne. 0.0_dp )  lseparator_noionode(kunit) 
        if ( rhoktimetot   .ne. 0.0_dp )  WRITE ( kunit ,110)             'struct fact             ', rhoktimetot
        if ( efgtimetot1   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_ES(real  part)      ', efgtimetot1
        if ( efgtimetot2   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_ES(recip part)      ', efgtimetot2 
        if ( efgtimetot3   .ne. 0.0_dp )  WRITE ( kunit ,110)             'efg_ES(comm  only)      ', efgtimetot3
      endif
      if ( fcoultime_es  .ne. 0.0_dp )  lseparator_noionode(kunit) 
      if ( fcoultime_es  .ne. 0.0_dp )  WRITE ( kunit ,'(a)')           'Coulomb:'
      if ( fcoultime_es  .ne. 0.0_dp )  lseparator_noionode(kunit) 
      if ( fcoultimetot1 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_ES(real  part)', fcoultimetot1
      if ( fcoultimetot2 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_ES(recip part)', fcoultimetot2 
      if ( fcoultimetot2_2 .ne. 0.0_dp )WRITE ( kunit ,110)             'multipole_ES(recip part rhonk)', fcoultimetot2_2 
      if ( fcoultimetot3 .ne. 0.0_dp )  WRITE ( kunit ,110)             'multipole_ES(comm only) ', fcoultimetot3
      if ( rhoktimetot .ne. 0.0_dp )    WRITE ( kunit ,110)             'charge_density_k        ',rhoktimetot
      if ( time_moment_from_pola .ne. 0.0_dp )  WRITE ( kunit ,110)     'moment_from_pola        ', time_moment_from_pola
    endif
    if ( vibtimetot.ne.0.0_dp )         lseparator_noionode(kunit) 
    if ( vibtimetot.ne.0.0_dp )         WRITE ( kunit ,'(a)')           'VIB:'
    if ( vibtimetot.ne.0.0_dp )         lseparator_noionode(kunit) 
    if ( hessiantimetot.ne.0.0_dp )     WRITE ( kunit ,110)             'hessian                ', hessiantimetot
    if ( diaghessiantimetot.ne.0.0_dp ) WRITE ( kunit ,110)             'full diag : hessian    ', diaghessiantimetot
    if ( bandtimetot.ne.0.0_dp )        WRITE ( kunit ,110)             'band                   ', bandtimetot
    if ( doskpttimetot.ne.0.0_dp )      WRITE ( kunit ,110)             'doskpttimetot          ', doskpttimetot

  !CALL print_dumb
  endif
  separator(kunit)

110   FORMAT(2X,A30,' :  cpu time',F9.2)


END SUBROUTINE print_time_info

END MODULE time
! ===== fmV =====
