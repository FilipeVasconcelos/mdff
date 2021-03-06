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
! along with this program; if not, WRITE to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
! ===== fmV =====

! ======= Hardware =======
#include "symbol.h"
#define debug
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE field 

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE kspace,                           ONLY :  kmesh 
  USE rspace,                           ONLY :  rmesh
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  USE tensors_rk,                       ONLY :  interaction_dd
  USE mpimdff

  implicit none

  integer :: cccc=0
  logical, SAVE     :: lwrite_dip        !< write dipoles 
  logical, SAVE     :: lwrite_quad       !< write quadrupoles to QUADFF
  logical, SAVE     :: lwrite_efg        !< write electric field gradient to EFGALL
  logical, SAVE     :: lwrite_ef         !< write electric field s to EFALL


  ! ============================================================  
  !                    force field type info 
  ! ============================================================  

  ! type dependent properties
  real(kind=dp)    :: mass     ( ntypemax )            !< masses ( not yet tested everywhere )



CONTAINS

! *********************** SUBROUTINE field_default_tag *************************
!> \brief
!! set default values to field tags
! ******************************************************************************
SUBROUTINE field_default_tag

  USE coulomb,                  ONLY :  qch, dip, quad, quad_nuc 

  implicit none

  ! =================
  !  default values
  ! =================

  ! field
  mass          = 1.0_dp
  qch           = 0.0_dp  ! charge
  dip           = 0.0_dp  ! dipolar moment
  quad          = 0.0_dp  ! quadrupolar moment
  quad_nuc      = 0.0_dp  ! quadrupolar moment
  lwrite_dip    = .false.            
  lwrite_quad   = .false.            
  lwrite_ef     = .false.            
  lwrite_efg    = .false.            

  return

END SUBROUTINE field_default_tag


! *********************** SUBROUTINE field_check_tag ***************************
!> \brief
!! check field tag values
! ******************************************************************************
SUBROUTINE field_check_tag

  USE control,                  ONLY :  lbmhftd , lbmhft , lcoulomb
  USE config,                   ONLY :  ntype , natm , dipia
  USE tt_damp,                  ONLY :  maximum_of_TT_expansion , get_TT_damp

  implicit none

  ! local
  !integer :: i, it, it2
  !logical :: ldamp , ldip, lqua, lqch


  return 

END SUBROUTINE field_check_tag

! *********************** SUBROUTINE field_init ********************************
!> \brief
!! force field initialisation
! ******************************************************************************
SUBROUTINE field_init

  USE control,                  ONLY :  calc , lnmlj , lcoulomb , lmorse , longrange, lnon_bonded
  USE non_bonded,               ONLY :  non_bonded_init
  USE coulomb,                  ONLY :  qch, dip, quad, quad_nuc ,coulomb_init, coulomb_print_info 
  USE pim,                      ONLY :  pim_init,pim_print_info

  implicit none

  ! local
  character(len=132)    :: filename
  integer               :: ioerr

  namelist /fieldtag/    mass          , &
                         qch           , &
                         quad_nuc      , &
                         dip           , &
                         quad          
  
                 

  ! ================================
  ! defaults values for field tags 
  ! ================================
  CALL field_default_tag

  ! ================================
  ! read field tags
  ! ================================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
    READ ( stdin , fieldtag, iostat=ioerr)
    if ( ioerr .lt. 0 )  then
      io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : fieldtag section is absent'
      STOP
    elseif ( ioerr .gt. 0 )  then
      io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : fieldtag wrong tag',ioerr
      STOP
    endif
  CLOSE ( stdin )

  ! ================================
  ! check field tags values
  ! ================================
  CALL field_check_tag

  ! ===============================================
  !  this routines generates the ewald parameters
  ! ===============================================
!  if ( longrange .eq. 'ewald' .and. lcoulomb ) CALL ewald_param
 
  ! ================================
  ! initialize constant parameters
  ! ================================
  if ( lnon_bonded )    then
    CALL non_bonded_init
  endif

  if ( lcoulomb ) then
    CALL coulomb_init
    CALL pim_init
    ! ================================
    !  print coulomb info
    ! ================================
    CALL coulomb_print_info(stdout)
    CALL pim_print_info(stdout)
  endif

  ! ================================
  !  pint field info
  ! ================================
  CALL field_print_info(stdout,quiet=.false.)

  return

END SUBROUTINE field_init

! *********************** SUBROUTINE field_print_info **************************
!> \brief
!! print force field information to standard output
! ******************************************************************************
SUBROUTINE field_print_info ( kunit , quiet )

  USE config,           ONLY :  ntype, atype, atypei, natmi , rho , simu_cell, rhoN!natm , ntype , atype , atypei , natmi , simu_cell , rhoN, rho , massia 
!  USE control,          ONLY :  calc , cutshortrange , lnmlj , lmorse , lbmhft , lbmhftd , lcoulomb , longrange , lreducedN , cutlongrange
  USE constants,        ONLY :  rho_unit
  USE coulomb,          ONLY :  qch, dip, quad, quad_nuc 

  implicit none
  ! global
  logical , optional :: quiet
  integer            :: kunit

  !local 
  integer            :: i, it
  logical            :: lquiet
  real(kind=dp)      :: total_mass, qtot, qtot2,  mu_sum(3) , theta_sum(3,3)


  if ( ( present ( quiet ) .and. quiet ) .and. .not. lquiet ) then
    lquiet = .true.
  else if ( ( present ( quiet ) .and. quiet ) .and. lquiet ) then
    return
  endif

  total_mass = 0.0_dp
  do it = 1 , ntype
    total_mass = total_mass + mass(it) * natmi(it)
  enddo

  ! ==============
  !  mass density
  ! ==============
  rho = total_mass / simu_cell%omega

  qtot   = 0.0_dp
  qtot2  = 0.0_dp
  mu_sum = 0.0_dp
  theta_sum = 0.0_dp
  do it = 1 , ntype
      qtot  = qtot  +   qch(it) * natmi ( it )
      qtot2 = qtot2 + ( qch(it) * natmi ( it ) ) * ( qch(it) * natmi ( it ) )
      mu_sum = mu_sum + dip(it,:)
      theta_sum = theta_sum + quad(it,:,:)
  enddo

  if ( ionode ) then
    separator(kunit)
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'FIELD MODULE ... WELCOME'
    blankline(kunit)
    lseparator(kunit)
    WRITE ( kunit ,'(a,f12.4,a)')       'total mass            = ', total_mass ,' a.m '
    WRITE ( kunit ,'(a,2f12.4,a)')      'density               = ', rhoN , rho * rho_unit ,' g/cm^3 '
    print*,rho_unit,rho
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'point charges: '
    lseparator(kunit)
    do it = 1 , ntype
      WRITE ( kunit ,'(a,a,a,e10.3)')   'q_',atypei(it),'                   = ',qch(it)
    enddo
    WRITE ( kunit ,'(a,e10.3)')         'total charge            = ',  qtot
    WRITE ( kunit ,'(a,e10.3)')         'second moment of charge = ',  qtot2
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'quadrupolar nuclear moment: '
    lseparator(kunit)
    do it = 1 , ntype
      WRITE ( kunit ,'(a,a,a,e10.3,a)') 'Q_',atypei(it),'                 = ',quad_nuc(it),' mb'
    enddo
    blankline(kunit)
    lseparator(kunit)
    WRITE ( kunit ,'(a)')               'static dipoles: '
    lseparator(kunit)
    do it = 1 , ntype
      WRITE ( kunit ,'(a,a,a,3e12.3)')  'mu_',atypei(it),'      = ',dip(it,1),dip(it,2),dip(it,3)
    enddo
    WRITE ( kunit ,'(a,3e12.3)')        'sum         = ',mu_sum(1),mu_sum(2),mu_sum(3)
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'static quadrupoles: '
    lseparator(kunit)
    do it = 1 , ntype
      WRITE ( kunit ,'(a,a,a)')         'theta_',atypei(it),'      = '
      do i = 1 , 3
        WRITE ( kunit ,'(3e16.8)')      quad(it,i,1) , quad(it,i,2) , quad(it,i,3)
      enddo
      WRITE ( kunit ,'(a,e16.8)')       'iso = ',(quad(it,1,1) + quad(it,2,2) + quad(it,3,3))/3.0_dp
    enddo
      WRITE ( kunit ,'(a)')             'sum   = '
      do i = 1 , 3
        WRITE ( kunit ,'(3e16.8)')      theta_sum(i,1) , theta_sum(i,2) , theta_sum(i,3)
      enddo
      WRITE ( kunit ,'(a,e16.8)')       'iso = ',(theta_sum(1,1) + theta_sum(2,2) + theta_sum(3,3))/3.0_dp
    blankline(kunit)

  endif

  return

END SUBROUTINE field_print_info 


! *********************** SUBROUTINE write_DIPFF ******************************
!
!>\brief
! write dipoles at ions to DIPFF file
!
! ******************************************************************************
SUBROUTINE write_DIPFF

  USE io,                       ONLY :  kunit_DIPFF
  USE coulomb,                  ONLY :  mu_t
  USE cell,                     ONLY :  kardir , periodicbc , dirkar
  USE control,                  ONLY :  lstatic
  USE config,                   ONLY :  system , natm , ntype , atype , simu_cell, atypei, natmi

  implicit none

  ! local
  integer :: ia , it

  if ( ionode ) then

  write(stdout,'(a)') 'writing DIPFF'
  if ( lstatic ) OPEN ( kunit_DIPFF ,file = 'DIPFF',STATUS = 'UNKNOWN')

    WRITE ( kunit_DIPFF , * )  natm
    WRITE ( kunit_DIPFF , * )  system
    WRITE ( kunit_DIPFF , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
    WRITE ( kunit_DIPFF , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
    WRITE ( kunit_DIPFF , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
    WRITE ( kunit_DIPFF , * )  ntype
    WRITE ( kunit_DIPFF , * )  ( atypei ( it ) , it = 1 , ntype )
    WRITE ( kunit_DIPFF , * )  ( natmi  ( it ) , it = 1 , ntype )
    WRITE ( kunit_DIPFF ,'(a)') &
              '      ia type                   mux                  muy                 muz'
    do ia= 1 , natm
      WRITE ( kunit_DIPFF , '(i8,2x,a3,6e24.16)' ) ia , atype( ia ) , mu_t ( 1 , ia ) , mu_t ( 2 , ia ) , mu_t ( 3 , ia )
    enddo
  endif

  if ( lstatic ) CLOSE( kunit_DIPFF )

  return

END SUBROUTINE write_DIPFF

! *********************** SUBROUTINE write_QUADFF ******************************
!
!>\brief
! write dipoles at ions to QUADFF file
!
! ******************************************************************************
SUBROUTINE write_QUADFF 

  USE io,                       ONLY :  kunit_QUADFF
  USE coulomb,                  ONLY :  theta_t
  USE cell,                     ONLY :  kardir , periodicbc , dirkar
  USE control,                  ONLY :  lstatic
  USE config,                   ONLY :  system , natm , ntype , atype , simu_cell, atypei, natmi

  implicit none

  ! local
  integer :: ia , it

  if ( ionode ) then

  write(stdout,'(a)') 'writing QUADFF'
  if ( lstatic ) OPEN ( kunit_QUADFF ,file = 'QUADFF',STATUS = 'UNKNOWN')

    WRITE ( kunit_QUADFF , * )  natm
    WRITE ( kunit_QUADFF , * )  system
    WRITE ( kunit_QUADFF , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
    WRITE ( kunit_QUADFF , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
    WRITE ( kunit_QUADFF , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
    WRITE ( kunit_QUADFF , * )  ntype
    WRITE ( kunit_QUADFF , * )  ( atypei ( it ) , it = 1 , ntype )
    WRITE ( kunit_QUADFF , * )  ( natmi  ( it ) , it = 1 , ntype )
    WRITE ( kunit_QUADFF ,'(a)') &
      '      ia type                   thetaxx                     thetayy                     thetazz                     thetaxy                  thetaxz                     thetayz'
    do ia= 1 , natm
          WRITE ( kunit_QUADFF ,'(i8,2x,a3,6e24.16)') ia , atype ( ia ) , theta_t ( 1 , 1, ia) , theta_t ( 2 , 2, ia) , &
                                                                          theta_t ( 3 , 3, ia) , theta_t ( 1 , 2, ia) , &
                                                                          theta_t ( 1 , 3, ia) , theta_t ( 2 , 3, ia)
    enddo
  endif

  if ( lstatic ) CLOSE( kunit_QUADFF )

  return

END SUBROUTINE write_QUADFF

! *********************** SUBROUTINE write_EFGALL ******************************
!
!>\brief
! write Electric Field Gradient at ions to EFGALL file
!
! ******************************************************************************
SUBROUTINE write_EFGALL

  USE io,                       ONLY :  kunit_EFGALL
  USE coulomb,                  ONLY :  efg_t
  USE pim,                      ONLY :  lwfc
  USE cell,                     ONLY :  kardir , periodicbc , dirkar
  USE control,                  ONLY :  lstatic, iefgall_format
  USE config,                   ONLY :  system , natm , ntype , atype , itype , simu_cell, atypei, natmi

  implicit none

  ! local
  integer :: ia , it

  if ( ionode ) then

    if ( lstatic ) OPEN ( kunit_EFGALL ,file = 'EFGALL',STATUS = 'UNKNOWN')
    if ( iefgall_format .ne. 0 ) then
      WRITE ( kunit_EFGALL , * )  natm
      WRITE ( kunit_EFGALL , * )  system
      WRITE ( kunit_EFGALL , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      WRITE ( kunit_EFGALL , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      WRITE ( kunit_EFGALL , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      WRITE ( kunit_EFGALL , * )  ntype
      WRITE ( kunit_EFGALL , * )  ( atypei ( it ) , it = 1 , ntype )
      WRITE ( kunit_EFGALL , * )  ( natmi  ( it ) , it = 1 , ntype )
      WRITE ( kunit_EFGALL ,'(a)') &
      '      ia type                   vxx                     vyy                     vzz                     vxy                     vxz                     vyz'
      do ia = 1 , natm
        it = itype ( ia )
        if ( lwfc( it ) .ge. 0 ) then
          WRITE ( kunit_EFGALL ,'(i8,2x,a3,6e24.16)') ia , atype ( ia ) , efg_t ( 1 , 1 , ia ) , efg_t ( 2 , 2 , ia ) , &
                                                                          efg_t ( 3 , 3 , ia ) , efg_t ( 1 , 2 , ia ) , &
                                                                          efg_t ( 1 , 3 , ia ) , efg_t ( 2 , 3 , ia )
        endif
      enddo
    endif

    if ( iefgall_format .eq. 0 ) then
      WRITE ( kunit_EFGALL )  natm
      WRITE ( kunit_EFGALL )  system
      WRITE ( kunit_EFGALL )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      WRITE ( kunit_EFGALL )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      WRITE ( kunit_EFGALL )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      WRITE ( kunit_EFGALL )  ntype
      WRITE ( kunit_EFGALL )  ( atypei ( it ) , it = 1 , ntype )
      WRITE ( kunit_EFGALL )  ( natmi  ( it ) , it = 1 , ntype )
      WRITE ( kunit_EFGALL )  efg_t
    endif

  endif

  if ( lstatic ) CLOSE( kunit_EFGALL )

  return

END SUBROUTINE write_EFGALL


! *********************** SUBROUTINE write_EFALL ******************************
!
!>\brief
! write configuration (pos,vel) to CONTFF file
!
! ******************************************************************************
SUBROUTINE write_EFALL

  USE io,                       ONLY :  kunit_EFALL
  USE pim,                      ONLY :  lwfc
  USE coulomb,                  ONLY :  ef_t
  USE cell,                     ONLY :  kardir , periodicbc , dirkar
  USE control,                  ONLY :  lstatic
  USE config,                   ONLY :  system , natm , ntype , atype , itype , simu_cell, atypei, natmi

  implicit none

  ! local
  integer :: ia , it

  if ( ionode ) then

    if ( lstatic ) OPEN ( kunit_EFALL ,file = 'EFALL',STATUS = 'UNKNOWN')
      WRITE ( kunit_EFALL , * )  natm
      WRITE ( kunit_EFALL , * )  system
      WRITE ( kunit_EFALL , * )  simu_cell%A(1,1) , simu_cell%A(2,1) , simu_cell%A(3,1)
      WRITE ( kunit_EFALL , * )  simu_cell%A(1,2) , simu_cell%A(2,2) , simu_cell%A(3,2)
      WRITE ( kunit_EFALL , * )  simu_cell%A(1,3) , simu_cell%A(2,3) , simu_cell%A(3,3)
      WRITE ( kunit_EFALL , * )  ntype
      WRITE ( kunit_EFALL , * )  ( atypei ( it ) , it = 1 , ntype )
      WRITE ( kunit_EFALL , * )  ( natmi  ( it ) , it = 1 , ntype )
      WRITE ( kunit_EFALL ,'(a)') &
                                 '      ia type                   Ex                      Ey                      Ez'
      do ia = 1 , natm
        it = itype ( ia )
        if ( lwfc( it ) .ge. 0 ) then
          WRITE ( kunit_EFALL ,'(i8,2x,a3,3e24.16)') ia , atype ( ia ) , ef_t ( ia , 1 ) , ef_t ( ia , 2 ) , ef_t ( ia , 3 )
        endif
      enddo

  endif

  if ( lstatic ) CLOSE( kunit_EFALL )

  return

END SUBROUTINE write_EFALL

END MODULE field 
! ===== fmV =====

