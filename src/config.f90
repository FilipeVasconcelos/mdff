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

#include "symbol.h"

! ======= Hardware =======
!#define debug

! ======= Hardware =======

! *********************** MODULE CONF ******************************************
!> \brief 
!> This module should deal with everything related to the current configuration
!> positions, velocities, forces, density , box .... 
!> \author
!> FMV
! ******************************************************************************
MODULE config

  USE constants,                ONLY :  dp 
  USE io,                       ONLY :  ionode, stdout
  USE cell,                     ONLY :  celltype
  USE mpimdff,                  ONLY :  decomposition

  implicit none

  integer, PARAMETER                           :: ntypemax = 16      !< maximum number of types
  integer, PARAMETER                           :: vnlmax   = 2000    !< maximum number of types

  character(len=60), SAVE                      :: system             !< system name                                              

  integer                                      :: natm               !< number of atoms
  integer                                      :: ntype              !< number of types
  integer                                      :: npairs             !< number of types pairs
  integer, dimension(:),           allocatable :: itype              !< type of atome i array 
  logical, dimension(:),           allocatable :: ipolar             !< .eq. 1 if polar 
  integer, dimension(0:ntypemax)               :: natmi              !< number of atoms (per type)

  TYPE ( celltype )                            :: simu_cell          !< simulation cell
  real(kind=dp)                                :: tau_nonb ( 3 , 3 ) !< stress tensor ( lennard-jones , morse ... )
  real(kind=dp)                                :: tau_coul ( 3 , 3 ) !< stress tensor coulombic
  real(kind=dp)                                :: rho                !< density  

  real(kind=dp), dimension(:)    , allocatable :: rx  , ry  , rz     !< positions
  real(kind=dp), dimension(:)    , allocatable :: vx  , vy  , vz     !< velocities
  real(kind=dp), dimension(:)    , allocatable :: fx  , fy  , fz     !< forces
  real(kind=dp), dimension(:)    , allocatable :: fxs , fys , fzs    !< forces (previous step t-dt) beeman
  real(kind=dp), dimension(:)    , allocatable :: rxs , rys , rzs    !< previous positions for leap-frog integrator
  real(kind=dp), dimension(:)    , allocatable :: xs  , ys  , zs     !< last positions in verlet list
  real(kind=dp), dimension(:)    , allocatable :: rix , riy , riz    !< positions in the center of mass reference 


  real(kind=dp), dimension(:)      , allocatable :: massia             !< mass on ion 
  real(kind=dp), dimension(:)      , allocatable :: qia                !< charge on ion 
  real(kind=dp), dimension(:)      , allocatable :: quadia_nuc         !< quadrupolar moment on ion
  real(kind=dp), dimension(:,:)    , allocatable :: dipia              !< dipole on ion 
  real(kind=dp), dimension(:,:)    , allocatable :: dipia_wfc          !< induced dipole on ion from Wannier centers
  real(kind=dp), dimension(:,:,:)  , allocatable :: quadia             !< quadrupole on ion 
  real(kind=dp), dimension(:,:,:)  , allocatable :: poldipia           !< dipole polarisability on ion
  real(kind=dp), dimension(:,:,:,:), allocatable :: polquadia          !< quadrupole polarisability on ion
  real(kind=dp), dimension(:,:,:)  , allocatable :: invpoldipia           !< polarisation on ion

  real(kind=dp), dimension(:)    , allocatable   :: phi_coul_tot       !< coulombic potential 

  character(len=3), dimension(:) , allocatable   :: atype            !< atom type label  
  character(len=3), dimension(0:ntypemax)        :: atypei           !< type label (per type)
  character(len=3), dimension(:,:) , allocatable :: allowedmove      !< atom type label  

  TYPE(decomposition), SAVE :: atom_dec

  TYPE :: verlet_list 
    integer, dimension(:),           allocatable :: list, point        !< vnlist info
    real(kind=dp)                                :: cut  
    character(len=4)                             :: listname
  END TYPE

  TYPE( verlet_list ) :: verlet_vdw
  TYPE( verlet_list ) :: verlet_coul


  ! =====================================================
  !   type of positions coordinates 
  ! =====================================================
  character(len=60), SAVE :: coord_format_allowed(4)
  data coord_format_allowed / 'Direct' , 'D' , 'Cartesian' , 'C' /

CONTAINS

! *********************** SUBROUTINE config_init *******************************
!>\brief
!> set default values, read and check consistenstency of conifig parameters
!> \author
!> FMV
! ******************************************************************************
SUBROUTINE config_init 

  USE control,  ONLY :  calc , cutshortrange , cutlongrange

  implicit none

  if ( calc .ne. 'md' .and. calc .ne. 'dist' ) return 
  ! ===========================================================
  ! read initial configuration only for calc = md
  ! for opt, vib and efg configuations are read in other routines 
  ! ===========================================================
    CALL read_pos

#ifdef debug
  CALL print_config_sample(0,0)
#endif

  verlet_vdw%cut=cutshortrange
  verlet_coul%cut=cutlongrange

  npairs =  ntype * ( ntype + 1 ) / 2

  return
 
END SUBROUTINE config_init 


! *********************** SUBROUTINE config_print_info *************************
!
!>\brief
! print information to standard output about the starting configuration
!
! ******************************************************************************
SUBROUTINE config_print_info(kunit)

  implicit none

  ! global 
  integer , intent ( in ) :: kunit 

  ! local
  integer :: it , i

  if ( ionode ) then
    blankline(kunit)
    separator(kunit)
    blankline(kunit)
    WRITE ( kunit ,'(a)')            'CONFIG MODULE ... WELCOME'
    blankline(kunit)
    WRITE ( kunit ,'(a,a)')          'system                : ',system
    WRITE ( kunit ,'(a,i16)')        'natm                  = ',natm
    WRITE ( kunit ,'(a,i16)')        'ntype                 = ',ntype
    do it = 1 , ntype     
      WRITE ( kunit ,'(a,a,a,i16,f8.2,a1)') &
                          'n',atypei(it),'                  = ',natmi(it),DBLE(natmi(it))/DBLE(natm) * 100.0_dp,'%'
    enddo
    blankline(kunit)
    lseparator(kunit)
    WRITE ( kunit ,'(a)')            'direct     basis : '
    lseparator(kunit)
    WRITE ( kunit ,'(a,3f12.4)')     'a_vector              = ',simu_cell%A(1,1),simu_cell%A(2,1),simu_cell%A(3,1) 
    WRITE ( kunit ,'(a,3f12.4)')     'b_vector              = ',simu_cell%A(1,2),simu_cell%A(2,2),simu_cell%A(3,2) 
    WRITE ( kunit ,'(a,3f12.4)')     'c_vector              = ',simu_cell%A(1,3),simu_cell%A(2,3),simu_cell%A(3,3) 
    WRITE ( kunit ,'(a,3f12.4)')     'cell param.           = ',(simu_cell%ANORM(i),i=1,3)
    WRITE ( kunit ,'(a,3f12.4)')     'perpend. width        = ',simu_cell%WA,simu_cell%WB,simu_cell%WC
    WRITE ( kunit ,'(a,3f12.4)')     'angles                = ',simu_cell%ALPH,simu_cell%BET,simu_cell%GAMM
    WRITE ( kunit ,'(a,f12.4)')      'volume                = ',simu_cell%omega
    blankline(kunit)
    blankline(kunit)
    lseparator(kunit)
    WRITE ( kunit ,'(a)')            'reciprocal basis : '
    lseparator(kunit)
    WRITE ( kunit ,'(a,3f12.4)')     'a*_vector             = ',simu_cell%B(1,1),simu_cell%B(2,1),simu_cell%B(3,1) 
    WRITE ( kunit ,'(a,3f12.4)')     'b*_vector             = ',simu_cell%B(1,2),simu_cell%B(2,2),simu_cell%B(3,2) 
    WRITE ( kunit ,'(a,3f12.4)')     'c*_vector             = ',simu_cell%B(1,3),simu_cell%B(2,3),simu_cell%B(3,3) 
    blankline(kunit)
    WRITE ( kunit ,'(a,3f12.4)')     'cell param.           = ',(simu_cell%BNORM(i),i=1,3)
    WRITE ( kunit ,'(a,3f12.4)')     'perpend. width        = ',simu_cell%RWA,simu_cell%RWB,simu_cell%RWC
    WRITE ( kunit ,'(a,3f12.4)')     'angles                = ',simu_cell%RALPH,simu_cell%RBET,simu_cell%RGAMM
    WRITE ( kunit ,'(a,f12.4)')      'volume                = ',simu_cell%romega
    blankline(kunit)
    lseparator(kunit)
    blankline(kunit)
    
  endif 

  return
  
END SUBROUTINE config_print_info

! *********************** SUBROUTINE write_CONTFF ******************************
!
!>\brief
! write configuration (pos,vel) to CONTFF file
!
! ******************************************************************************
SUBROUTINE write_CONTFF

  USE io,                       ONLY :  kunit_CONTFF
  USE cell,                     ONLY :  kardir , periodicbc , dirkar

  implicit none

  ! local
  integer :: ia , it
  real(kind=dp), dimension (:) , allocatable :: xxx , yyy , zzz

  allocate ( xxx ( natm ) , yyy ( natm ) , zzz ( natm ) )

  xxx = rx
  yyy = ry
  zzz = rz
  ! ======================================
  !                  PBC
  ! ======================================
  CALL kardir     ( natm , xxx , yyy , zzz , simu_cell%B )
  CALL periodicbc ( natm , xxx , yyy , zzz )
  CALL dirkar     ( natm , xxx , yyy , zzz , simu_cell%A )
  
  if ( ionode ) then
  OPEN ( kunit_CONTFF ,file = 'CONTFF',STATUS = 'UNKNOWN')
      WRITE ( kunit_CONTFF,'(i8)') natm 
      WRITE ( kunit_CONTFF,'(a)') system
      WRITE ( kunit_CONTFF,'(3f20.12)') simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
      WRITE ( kunit_CONTFF,'(3f20.12)') simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
      WRITE ( kunit_CONTFF,'(3f20.12)') simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
      WRITE ( kunit_CONTFF,'(i4)') ntype 
      WRITE ( kunit_CONTFF,*) ( atypei(it) , it=1,ntype ) 
      WRITE ( kunit_CONTFF,*) ( natmi (it) , it=1,ntype ) 
      WRITE ( kunit_CONTFF,'(A)') 'Cartesian' 
      WRITE ( kunit_CONTFF,'(a,9e24.16)') ( atype ( ia ) , xxx ( ia ) , yyy ( ia ) , zzz ( ia ) , & 
                                                           vx  ( ia ) , vy  ( ia ) , vz  ( ia ) , &
                                                           fx  ( ia ) , fy  ( ia ) , fz ( ia )  , ia = 1 , natm )
  CLOSE (kunit_CONTFF)
  endif

  deallocate ( xxx , yyy , zzz ) 

  return

END SUBROUTINE write_CONTFF


! *********************** SUBROUTINE config_alloc ******************************
!> \brief
!> allocation of principal arrays of the calculation
! ******************************************************************************
SUBROUTINE config_alloc

  implicit none

  allocate( rx  ( natm ) , ry ( natm )  , rz ( natm ) )
  allocate( vx  ( natm ) , vy ( natm )  , vz ( natm ) )
  allocate( fx  ( natm ) , fy ( natm )  , fz ( natm ) )
  allocate( fxs ( natm ) , fys ( natm ) , fzs ( natm ) )
  allocate( rxs ( natm ) , rys ( natm ) , rzs ( natm ) )
  allocate( xs ( natm ) , ys ( natm ) , zs ( natm ) ) 
  allocate( atype ( natm ) )
  allocate( itype ( natm ) )
  allocate( allowedmove ( 3 , natm ) )
  allocate( verlet_vdw%list ( natm * vnlmax )  , verlet_vdw%point (  natm + 1 ) )
  allocate( verlet_coul%list ( natm * vnlmax ) , verlet_coul%point (  natm + 1 ) )
  allocate( qia ( natm ) )
  allocate( massia ( natm ) )
  allocate( quadia_nuc ( natm ) )
  allocate( dipia ( 3 , natm ) )
  allocate( dipia_wfc ( 3 , natm ) )
  allocate( quadia ( 3 , 3 , natm ) )
  allocate( poldipia ( 3 , 3  , natm ) )
  allocate( polquadia ( 3 , 3 , 3  , natm ) )
  allocate( invpoldipia ( 3 , 3 , natm ) )
  allocate( ipolar ( natm ) )
  allocate( phi_coul_tot ( natm ) ) !< only if we calculated coulombic interactions

  rx    = 0.0_dp
  ry    = 0.0_dp
  rz    = 0.0_dp
  vx    = 0.0_dp
  vy    = 0.0_dp
  vz    = 0.0_dp
  fx    = 0.0_dp
  fy    = 0.0_dp
  fz    = 0.0_dp
  fxs   = 0.0_dp
  fys   = 0.0_dp
  fzs   = 0.0_dp
  xs    = 0.0_dp
  ys    = 0.0_dp
  zs    = 0.0_dp
  rxs   = 0.0_dp
  rys   = 0.0_dp
  rzs   = 0.0_dp
  atype = ''
  verlet_vdw%list    = 0
  verlet_vdw%point   = 0
  verlet_coul%list   = 0
  verlet_coul%point  = 0
  qia        = 0.0_dp
  massia     = 1.0_dp
  quadia_nuc = 0.0_dp
  dipia      = 0.0_dp
  quadia     = 0.0_dp
  dipia_wfc = 0.0_dp
  poldipia     = 0.0_dp
  polquadia     = 0.0_dp
  ipolar    = .false. 
  phi_coul_tot = 0.0_dp

  return 
 
END SUBROUTINE config_alloc


! *********************** SUBROUTINE config_dealloc ****************************
!> \brief
!! Deallocate config quantities (see config_alloc)
! ******************************************************************************
SUBROUTINE config_dealloc

  USE control, ONLY : calc

  implicit none 
       
  ! tmp 
  if ( calc .eq.'rmc' .or. calc .eq.'stochio') return

  deallocate( rx  , ry  , rz )
  deallocate( vx  , vy  , vz )
  deallocate( fx  , fy  , fz )
  deallocate( fxs , fys , fzs )
  deallocate( rxs , rys , rzs )
  deallocate( xs , ys , zs )
  deallocate( atype )
  deallocate( itype )
  deallocate( allowedmove )
  deallocate( verlet_vdw%list , verlet_vdw%point )
  deallocate( verlet_coul%list , verlet_coul%point )
  deallocate( qia ) 
  deallocate( massia ) 
  deallocate( quadia_nuc ) 
  deallocate( dipia ) 
  deallocate( dipia_wfc ) 
  deallocate( quadia ) 
  deallocate( poldipia ) 
  deallocate( polquadia ) 
  deallocate( invpoldipia ) 
  deallocate( ipolar ) 
  deallocate( phi_coul_tot ) !< well only if we calculated coulombic interactions

  return 

END SUBROUTINE config_dealloc

! *********************** SUBROUTINE center_of_mass ****************************
!> \brief
!! for the moment still here but should be moved somewhere else
!! \param[in] ax, ay, az position vector
!! \param[out] com center of mass
!! \author
!! unknown
!! \note
!! adapted from quantum-espresso
!! \todo
!! should depend on mass 
! ******************************************************************************
SUBROUTINE center_of_mass ( ax , ay , az , com )

  implicit none

  ! global
  real(kind=dp) , intent ( in  ) :: ax ( natm ) , ay ( natm ) , az ( natm )
  real(kind=dp) , intent ( out ) :: com ( 0 : ntypemax , 3 )

  ! local
  integer :: ia , it 

  com = 0.0_dp

  do ia = 1 , natm
    it = itype ( ia )    
    com ( it , 1 ) = com ( it , 1 )  + ax ( ia ) ! * m 
    com ( it , 2 ) = com ( it , 2 )  + ay ( ia ) ! * m
    com ( it , 3 ) = com ( it , 3 )  + az ( ia ) ! * m 
    com ( 0  , 1 ) = com ( 0  , 1 )  + ax ( ia ) ! * m 
    com ( 0  , 2 ) = com ( 0  , 2 )  + ay ( ia ) ! * m
    com ( 0  , 3 ) = com ( 0  , 3 )  + az ( ia ) ! * m
  enddo

  do it = 0 , ntype
    com ( it , 1 )  = com ( it , 1 ) / DBLE ( natmi ( it ) )
    com ( it , 2 )  = com ( it , 2 ) / DBLE ( natmi ( it ) )
    com ( it , 3 )  = com ( it , 3 ) / DBLE ( natmi ( it ) )
  enddo

  return

END SUBROUTINE center_of_mass


! *********************** SUBROUTINE linear_momentum ***************************
!> \brief
!! Calculate the linear momentum (should be conserved along nve traj)
!! \note
!! not used so far
! ******************************************************************************
SUBROUTINE linear_momentum

  implicit none

  ! local
  integer :: ia
  real(kind=dp) :: Px, Py, Pz, normP

  do ia = 1 , natm
    Px = Px + vx ( ia )
    Py = Py + vy ( ia ) 
    Pz = Pz + vz ( ia )
  enddo

  normP = SQRT(Px*Px + Py*Py + Pz*Pz)

  return  

END SUBROUTINE linear_momentum

! *********************** SUBROUTINE angular_momentum **************************
!> \brief
!! Calculate the angular momentum (not conserved with pbc)
!! \note
!! not used so far
! ******************************************************************************
SUBROUTINE angular_momentum ( Lx , Ly , Lz , normL )

  implicit none

  ! local
  integer :: ia 
  real(kind=dp) :: Lx, Ly, Lz, normL

  Lx = 0.0_dp
  Ly = 0.0_dp
  Lz = 0.0_dp      
  do ia = 1 , natm
   Lx = Lx + ry ( ia ) * vz ( ia ) - rz ( ia ) * vy ( ia ) 
   Ly = Ly + rz ( ia ) * vx ( ia ) - rx ( ia ) * vz ( ia ) 
   Lz = Lz + rx ( ia ) * vy ( ia ) - ry ( ia ) * vx ( ia )
  enddo

  normL = SQRT( Lx * Lx + Ly * Ly + Lz * Lz)

  return

END SUBROUTINE angular_momentum


! *********************** SUBROUTINE ions_reference_positions ******************
!> \brief
!! Calculate the real position of atoms relative to the center of mass (cdm)
!! and store them in taui
!! cdmi: initial position of the center of mass (cdm) in cartesian coor.  
!! \note
!! not used so far ... probably not working as well 
! ******************************************************************************
SUBROUTINE ions_reference_positions

  implicit none
  real(kind=dp) :: com ( 0:ntypemax, 3 )
  integer  :: ia

  CALL center_of_mass ( rx , ry , rz , com )

  do ia = 1 , natm
    rix ( ia ) = rx ( ia ) - com ( 0 , 1 )
    riy ( ia ) = ry ( ia ) - com ( 0 , 1 )
    riz ( ia ) = rz ( ia ) - com ( 0 , 1 )
  enddo

  return 
 
END SUBROUTINE ions_reference_positions


! *********************** SUBROUTINE ions_displacement *************************
!> \brief
!! Calculate the sum of the quadratic displacements of the atoms in the ref.
!! of cdm respect to the initial positions.
!! \note
!! not used so far
!! tau_ref: starting position in center-of-mass ref. in real units
! ******************************************************************************
SUBROUTINE ions_displacement( dis, ax , ay , az )

  implicit none
 
  ! global  
  real(kind=dp), intent ( out ) :: dis(:)
  real(kind=dp), intent ( in  ) :: ax (:) , ay(:) , az(:)

  ! local
  real(kind=dp) :: rdist(3), r2, com(0:ntypemax,3)
  INTEGER  :: it, ia, isa

  ! =========================================================
  !  Compute the current value of cdm "Centro Di Massa"
  ! =========================================================
  CALL center_of_mass ( ax , ay , az , com )
 
  isa = 0
  do it = 1, ntype
    dis(it) = 0.0_dp
    r2      = 0.0_dp
    do ia = 1 , natmi(it)
      isa = isa + 1
      rdist ( 1 ) = rx (isa) - com ( 0 , 1 ) 
      rdist ( 2 ) = ry (isa) - com ( 0 , 2 )
      rdist ( 3 ) = rz (isa) - com ( 0 , 3 )
      r2 = r2 + ( rdist( 1 ) - rix(isa) )**2 + &
                ( rdist( 2 ) - riy(isa) )**2 + &
                ( rdist( 3 ) - riz(isa) )**2 
    enddo 
    dis(it) = dis(it) + r2 / DBLE(natmi(it))
  enddo
  
  return
 
END SUBROUTINE ions_displacement

! *********************** SUBROUTINE write_trajff_xyz **************************
!> \brief
!! write trajectory (pos, vel, for) to TRAJFF file
! ******************************************************************************
SUBROUTINE write_trajff_xyz

  USE control,                  ONLY :  itraj_format , trajff_data
  USE io,                       ONLY :  kunit_TRAJFF
  USE cell,                     ONLY :  periodicbc , kardir , dirkar

  implicit none

  ! local
  integer :: ia , it , i , j
  real(kind=dp), dimension (:) , allocatable :: xxx , yyy , zzz

  allocate ( xxx ( natm ) , yyy ( natm ) , zzz ( natm ) )

  xxx = rx
  yyy = ry
  zzz = rz

  ! ======================================
  !             PBC
  ! ======================================
  CALL kardir     ( natm , xxx , yyy , zzz , simu_cell%B )
  CALL periodicbc ( natm , xxx , yyy , zzz )
  CALL dirkar     ( natm , xxx , yyy , zzz , simu_cell%A )


  if ( ionode ) then
   if ( itraj_format .ne. 0) then
     WRITE ( kunit_TRAJFF , '(i6)' ) natm
     WRITE ( kunit_TRAJFF , '(a)' ) system
     do i = 1 , 3
       WRITE ( kunit_TRAJFF , '(3f20.12)' ) (simu_cell%A(i,j),j=1,3)
     enddo
     WRITE ( kunit_TRAJFF , '(i4)' ) ntype
     WRITE ( kunit_TRAJFF , * ) ( atypei ( it ) , it = 1 , ntype )
     WRITE ( kunit_TRAJFF , * ) ( natmi  ( it ) , it = 1 , ntype )
     WRITE ( kunit_TRAJFF,'(A)') 'Cartesian'
    if ( trajff_data .eq. 'rvf' ) then  
      do ia = 1 , natm
        WRITE ( kunit_TRAJFF ,'(a,9e20.12)' ) atype ( ia ) , xxx ( ia ) , yyy ( ia ) , zzz ( ia ) , vx( ia ) , vy( ia ) , vz( ia ) , fx( ia ) , fy( ia ) , fz( ia )
      enddo
    endif
    if ( trajff_data .eq. 'rvn' ) then
      do ia = 1 , natm
        WRITE ( kunit_TRAJFF ,'(a,9e20.12)' ) atype ( ia ) , xxx ( ia ) , yyy ( ia ) , zzz ( ia ) , vx( ia ) , vy( ia ) , vz( ia )
      enddo
    endif
    if ( trajff_data .eq. 'rnf' ) then
      do ia = 1 , natm
        WRITE ( kunit_TRAJFF ,'(a,9e20.12)' ) atype ( ia ) , xxx ( ia ) , yyy ( ia ) , zzz ( ia ) , fx( ia ) , fy( ia ) , fz( ia )
      enddo
    endif
    if ( trajff_data .eq. 'rnn' ) then
      do ia = 1 , natm
        WRITE ( kunit_TRAJFF ,'(a,9e20.12)' ) atype ( ia ) , xxx ( ia ) , yyy ( ia ) , zzz ( ia )
      enddo
    endif
   endif
   if ( itraj_format .eq. 0) then
     WRITE ( kunit_TRAJFF ) natm
     WRITE ( kunit_TRAJFF ) system
     do i = 1 , 3
        WRITE ( kunit_TRAJFF ) (simu_cell%A(i,j),j=1,3)
     enddo
     WRITE ( kunit_TRAJFF ) ntype
     WRITE ( kunit_TRAJFF ) ( atypei ( it ) , it = 1 , ntype )
     WRITE ( kunit_TRAJFF ) ( natmi  ( it ) , it = 1 , ntype )
     WRITE ( kunit_TRAJFF ) 'C'
    if ( trajff_data .eq. 'rvf' ) then  
      WRITE ( kunit_TRAJFF ) xxx , yyy , zzz , vx , vy , vz , fx , fy , fz
    endif
    if ( trajff_data .eq. 'rvn' ) then
      WRITE ( kunit_TRAJFF ) xxx , yyy , zzz , vx , vy , vz
    endif
    if ( trajff_data .eq. 'rnf' ) then
      WRITE ( kunit_TRAJFF ) xxx , yyy , zzz , fx , fy , fz
    endif
    if ( trajff_data .eq. 'rnn' ) then
      WRITE ( kunit_TRAJFF ) xxx , yyy , zzz 
    endif
   endif
  endif

  deallocate ( xxx  , yyy  , zzz )

  return

END SUBROUTINE write_trajff_xyz

SUBROUTINE read_traj_header ( kunit , iformat ) 

  implicit none

  integer, intent(in) :: kunit , iformat 
  integer            :: i , it
  character(len=60)  :: coord_format
  character(len=60)  :: coord_format1
  character(len=1)   :: coord_format2

  if ( iformat .ne. 0 ) then
    READ ( kunit , * ) natm
    READ ( kunit , * ) system
    READ ( kunit , * ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit , * ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit , * ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit , * ) ntype
    READ ( kunit , * ) ( atypei ( it ) , it = 1 , ntype )
    READ ( kunit , * ) ( natmi ( it )  , it = 1 , ntype )
    READ ( kunit , * ) coord_format1
    coord_format = coord_format1
  endif
  if ( iformat .eq. 0 ) then
    READ ( kunit ) natm
    READ ( kunit ) system
    READ ( kunit ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit ) ntype
    READ ( kunit ) ( atypei ( it ) , it = 1 , ntype )
    READ ( kunit ) ( natmi ( it )  , it = 1 , ntype )
    READ ( kunit ) coord_format2
    coord_format = coord_format2
  endif
  if ( ionode ) WRITE ( stdout ,'(A,20A3)' ) 'found type information on TRAJFF : ', atypei ( 1:ntype )



  ! ===================
  !  coord_format test
  ! ===================
  CALL check_allowed_tags ( size( coord_format_allowed ) , coord_format_allowed , coord_format , ' in  TRAJFF at line 9' , 'coord_format' )

  if ( ionode .and. &
       ( coord_format .eq. 'Direct' .or. coord_format .eq. 'D' ) ) &
       WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in TRAJFF'
  if ( ionode .and. &
       ( coord_format .eq. 'Cartesian' .or. coord_format .eq. 'C' ) ) &
       WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in cartesian coordinates in TRAJFF'

#ifdef debug
   write(*,*) simu_cell%A
#endif
  CLOSE ( kunit )

  return

END SUBROUTINE read_traj_header

SUBROUTINE read_traj ( kunit , iformat , csave ) 

  USE cell,             ONLY :  dirkar

  implicit none

  integer, intent(in) :: kunit , iformat
  character(len=3)    :: csave
  integer             :: ia , i , it
  character(len=60)   :: coord_format
  character(len=60)   :: coord_format1
  character(len=1)    :: coord_format2

  if ( iformat .ne. 0 ) then
    READ ( kunit , * ) natm
    READ ( kunit , * ) system
    READ ( kunit , * ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit , * ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit , * ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit , * ) ntype
    READ ( kunit , * ) ( atypei ( it ) , it = 1 , ntype )
    READ ( kunit , * ) ( natmi ( it )  , it = 1 , ntype )
    READ ( kunit , * ) coord_format1
    coord_format = coord_format1
    if ( csave .eq. 'rvf' ) then  
      do ia = 1 , natm
        READ ( kunit , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , vx( ia ) , vy( ia ) , vz( ia ) , fx( ia ) , fy( ia ) , fz( ia )
      enddo
    endif
    if ( csave .eq. 'rvn' ) then
      do ia = 1 , natm
        READ ( kunit , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , vx( ia ) , vy( ia ) , vz( ia )
      enddo
    endif
    if ( csave .eq. 'rnf' ) then
      do ia = 1 , natm
        READ ( kunit , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , fx( ia ) , fy( ia ) , fz( ia )
      enddo
    endif
    if ( csave .eq. 'rnn' ) then
      do ia = 1 , natm
        READ ( kunit , * ) atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia )
      enddo
    endif
  endif
  if ( iformat .eq. 0 ) then
    READ ( kunit ) natm
    READ ( kunit ) system
    READ ( kunit ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit ) ntype
    READ ( kunit ) ( atypei ( it ) , it = 1 , ntype )
    READ ( kunit ) ( natmi ( it )  , it = 1 , ntype )
    READ ( kunit ) coord_format2
    coord_format = coord_format2
    if ( csave .eq. 'rvf' ) then  
      READ ( kunit ) rx , ry , rz , vx , vy , vz , fx , fy , fz
    endif
    if ( csave .eq. 'rvn' ) then
      READ ( kunit ) rx , ry , rz , vx , vy , vz
    endif
    if ( csave .eq. 'rnf' ) then
      READ ( kunit ) rx , ry , rz , fx , fy , fz
    endif
    if ( csave .eq. 'rnn' ) then
      READ ( kunit ) rx , ry , rz 
    endif
  endif

  ! ===================
  !  coord_format test
  ! ===================
  CALL check_allowed_tags ( size ( coord_format_allowed ), coord_format_allowed, coord_format , ' in POSFF at line 9' , 'coord_format' )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  if ( coord_format .eq. 'Direct' .or. coord_format .eq.'D' ) CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE read_traj

END MODULE config
! ===== fmV =====
