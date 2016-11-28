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
!#define debug_read_pos
! ======= Hardware =======

! *********************** SUBROUTINE read_pos **********************************

!> \brief
!  this subroutine read configuration (pos,vel,force) from POSFF file
!
!> \todo
!! find a way to have weither only pos or pos/vel or pos/vel/for in POSFF and
!! other configurations file
!
! ******************************************************************************
SUBROUTINE read_pos 
 
  USE constants,                ONLY :  dp 
  USE control,                  ONLY :  calc , posff_data
  USE config,                   ONLY :  rx , ry , rz , vx , vy , vz , fx , fy , fz , atype , atypei , itype , &
                                        natmi , natm , dipia , qia , ipolar , rho , system , ntype , config_alloc , &
                                        simu_cell , config_print_info , coord_format_allowed , write_CONTFF
  USE field,                    ONLY :  qch , dip , quad, ldip_polar , field_init
  USE io,                       ONLY :  ionode , stdout , kunit_POSFF
  USE cell,                     ONLY :  lattice, periodicbc , dirkar, kardir

  implicit none

  ! local
  integer           :: it , ia , i
  logical           :: allowed
  character(len=60) :: cpos

  separator(stdout) 
  IF ( ionode ) then
    blankline(stdout)    
    WRITE ( stdout ,'(a)')   'reading configuration'
    WRITE ( stdout ,'(a)')   'config from file POSFF'
  endif

  OPEN ( kunit_POSFF , file = 'POSFF' ) 
  READ ( kunit_POSFF ,* ) natm 
  READ ( kunit_POSFF ,* ) system
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
  READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
  READ ( kunit_POSFF ,* ) ntype
  READ ( kunit_POSFF ,* ) ( atypei ( it ) , it = 1 , ntype )
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
  READ( kunit_POSFF ,*) ( natmi ( it ) , it = 1 , ntype )
  READ( kunit_POSFF ,*) cpos
  ! ======
  !  cpos
  ! ======
  do i = 1 , size( coord_format_allowed )
   if ( trim(cpos) .eq. coord_format_allowed(i))  allowed = .true.
  enddo
  if ( .not. allowed ) then
    if ( ionode )  WRITE ( stdout , '(a)' ) 'ERROR in POSFF at line 9 should be ', coord_format_allowed
    STOP
  endif
  if ( cpos .eq. 'Direct' .or. cpos .eq. 'D' ) then
    io_node WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in POSFF'
  else if ( cpos .eq. 'Cartesian' .or. cpos .eq. 'C' ) then
    io_node WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in cartesian coordinates in POSFF'
  endif

 
  CALL lattice ( simu_cell )
  rho = DBLE ( natm ) / simu_cell%omega

  CALL config_alloc 

  ! ===============================
  ! config info
  ! ===============================
  CALL config_print_info(stdout)

  ! =========================================      
  ! read positions, velocities, forces from disk 
  ! =========================================      
  if ( posff_data == "rnn" ) then
    READ  ( kunit_POSFF , * ) ( atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , ia = 1 , natm )
    
  else if ( posff_data == "rvn" ) then
    READ  ( kunit_POSFF , * ) ( atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , &
                                               vx ( ia ) , vy ( ia ) , vz ( ia ) , ia = 1 , natm ) 
  else if ( posff_data == "rvf" ) then  
    READ  ( kunit_POSFF , * ) ( atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , &
                                               vx ( ia ) , vy ( ia ) , vz ( ia ) , &
                                               fx ( ia ) , fy ( ia ) , fz ( ia ) , ia = 1 , natm )
  endif
  !27-11-16 a-t-on besoin de ca ?
  !if ( .not. lwrite_restart ) then  
  !  if ( ionode .and. (( posff_data == "rvn" ) .or. ( posff_data == "rvf" )) ) WRITE ( stdout ,'(A,20A3)' ) &
  !  'WARNING in non restart mode velocities and forces are not considered even present in the input file POSFF'
  !  vx=0.0_dp
  !  vy=0.0_dp
  !  vz=0.0_dp
  !  fx=0.0_dp
  !  fy=0.0_dp
  !  fz=0.0_dp
  !endif
  ! ===============================
  ! init force_field 
  ! ===============================
  CALL field_init

  if ( cpos .eq. 'Direct' .or. cpos .eq. 'D' ) then
    io_node WRITE ( stdout ,'(A,20A3)' ) 'input positions in Direct coordinates'
    ! ======================================
    !         direct to cartesian
    ! ======================================
    CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
  else if ( cpos .eq. 'Cartesian' .or. cpos .eq. 'C' ) then
    io_node WRITE ( stdout ,'(A,20A3)' ) 'input positions in Cartesian coordinates'
  endif 

  CLOSE ( kunit_POSFF )

  CALL kardir     ( natm , rx , ry , rz , simu_cell%B )
  CALL periodicbc ( natm , rx , ry , rz  )
  CALL dirkar     ( natm , rx , ry , rz , simu_cell%A )

  CALL typeinfo_init
#ifdef debug_read_pos
  CALL distance_tab 
#endif

  return

END SUBROUTINE read_pos

! *********************** SUBROUTINE typeinfo_init *****************************
!
!> \brief
!! this subroutine initialize main atomic quanitites defined from type
!! definitions  
!
!> \note
!! see config module for definitions of the main quantities
!
! ******************************************************************************
SUBROUTINE typeinfo_init

  USE constants, ONLY :  dp
  USE config ,  ONLY :  atype , atypei , itype , natmi , natm , ntype , massia, dipia , quadia, qia , &
                        quadia_nuc , ipolar , poldipia , polquadia, invpoldipia 
  USE field ,   ONLY :  mass, qch , quad_nuc , dip , quad , ldip_polar , poldip , polquad
  USE io,       ONLY :  stdout
  
  implicit none

  ! local  
  integer :: ccs , cc
  integer :: ia , it, i , j 
  integer, parameter :: LWORK=1000
  real(kind=dp) :: WORK ( LWORK )
  integer :: ipiv ( 3 )
  integer :: ierr


  !print*,'in typeinfo_init'
  ! ==========================
  !  set some type parameters
  ! ==========================      
  natmi ( 0 ) = 0
  cc = 0
  do it = 1 , ntype
    ccs = cc
    cc = cc + natmi ( it )
    do ia = ccs + 1 , cc
      atype ( ia )                = atypei ( it )
      itype ( ia )                = it
      qia   ( ia )                = qch ( it )
      massia( ia )                = mass( it )
      quadia_nuc ( ia )           = quad_nuc ( it )   
      dipia ( : , ia )            = dip ( it , : )
      quadia ( : , : , ia )       = quad (it , : , : )
      ipolar ( ia )               = ldip_polar ( it )
      poldipia ( : , :  , ia )    = poldip ( it , : , : )
      polquadia( : , : , : , ia ) = polquad ( it , : , : , : )
      if ( .not. ldip_polar ( it ) ) cycle
      invpoldipia ( : , : , ia )  = poldipia ( : , : , ia )
      CALL DGETRF( 3, 3, invpoldipia(:,:,ia), 3, ipiv, ierr )  
      if ( ierr.lt.0 ) then
        WRITE( stdout , '(a,i6)' ) 'ERROR call to DGETRF failed in induced_moment',ierr
        STOP
      endif
      CALL DGETRI( 3 , invpoldipia(:,:,ia) , 3 ,  ipiv , WORK, LWORK, ierr )
      if ( ierr.lt.0 ) then
        WRITE( stdout , '(a,i6)' ) 'ERROR call to DGETRI failed in induced_moment',ierr
        STOP
      endif
    enddo
  enddo
  natmi  ( 0 ) = natm
  atypei ( 0 ) = 'ALL'

#ifdef debug
do ia=1,natm
 WRITE ( stdout ,'(a,a2,a,f12.4)')'polarizability on ion ', atype(ia),' : '
 WRITE ( stdout ,'(3f12.4)')      ( poldipia ( 1 , j , ia ) , j = 1 , 3 )
 WRITE ( stdout ,'(3f12.4)')      ( poldipia ( 2 , j , ia ) , j = 1 , 3 )
 WRITE ( stdout ,'(3f12.4)')      ( poldipia ( 3 , j , ia ) , j = 1 , 3 )

 WRITE ( stdout ,'(a,a2,a,f12.4)')'inverse polarizability on ion ', atype(ia),' : '
 WRITE ( stdout ,'(3f12.4)')      ( invpoldipia ( 1 , j , ia ) , j = 1 , 3 )
 WRITE ( stdout ,'(3f12.4)')      ( invpoldipia ( 2 , j , ia ) , j = 1 , 3 )
 WRITE ( stdout ,'(3f12.4)')      ( invpoldipia ( 3 , j , ia ) , j = 1 , 3 )
enddo
#endif

  return

END SUBROUTINE typeinfo_init



