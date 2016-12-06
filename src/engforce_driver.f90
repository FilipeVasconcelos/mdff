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
!#define debug
!#define debug_ES
!#define debug_ES_field_forces
!#define debug_ES_energy
!#define debug_ES_stress
!#define debug_ES_efg
!#define debug_ES_dir
!#define debug_scf_pola
!#define debug_wfc
!#define debug_morse
!#define debug_nmlj
!#define debug_nmlj_pbc
!#define debug_quadratic
!#define debug_para
!#define debug_mu
!#define debug_cg
!#define debug_extrapolate
!#define debug_print_dipff_scf
!#define debug_scf_kO_inner
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE engforce_driver 

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE kspace,                           ONLY :  kmesh 
  USE rspace,                           ONLY :  rmesh
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  USE tensors_rk,                       ONLY :  interaction_dd
  USE mpimdff

  implicit none

  integer :: cccc=0
  logical, PRIVATE :: symmetric_pot

  ! ============================================================  
  !                    force field type info 
  ! ============================================================  

  ! type dependent properties
  real(kind=dp)    :: mass     ( ntypemax )            !< masses ( not yet tested everywhere )
  real(kind=dp)    :: qch      ( ntypemax )            !< charges 
  real(kind=dp)    :: dip      ( ntypemax , 3 )        !< dipoles 
  real(kind=dp)    :: quad     ( ntypemax , 3 , 3 )    !< quadrupoles
  real(kind=dp)    :: poldip   ( ntypemax , 3 , 3 )    !< dipole     polarizability if ldip_polar( it ) = .true. 
  real(kind=dp)    :: poldip_iso ( ntypemax )          !< isotropic dipole polarizability if ldip_polar( it ) = .true.
  real(kind=dp)    :: polquad  ( ntypemax , 3 , 3 , 3 )!< quadrupole polarizability if lquad_polar( it ) = .true.
  real(kind=dp)    :: polquad_iso ( ntypemax )         !< isotropic quadrupole polarizability if ldip_polar( it ) = .true. 
  real(kind=dp)    :: quad_nuc ( ntypemax )            !< quadrupolar moment nucleus NMR
  logical          :: ldip_polar   ( ntypemax )                          !< induced moment from pola. Is this type of ion polarizable ?

  INTERFACE 

    SUBROUTINE multipole_ES ( ef , efg , mu , theta , task , damp_ind , &
                                            do_efield , do_efg , do_forces , do_stress , do_rec , do_dir , do_strucfact , use_ckrskr )
      real(selected_real_kind(15,300))     :: ef     ( : , : )
      real(selected_real_kind(15,300))     :: efg    ( : , : , : )
      real(selected_real_kind(15,300))     :: mu     ( : , : )
      real(selected_real_kind(15,300))     :: theta  ( : , : , : )
      logical  :: task   ( : )
      logical  :: damp_ind , do_efield , do_efg, do_forces, do_stress, do_rec , do_dir , do_strucfact , use_ckrskr

    END SUBROUTINE multipole_ES

  END INTERFACE        


CONTAINS

! *********************** SUBROUTINE engforce ***************************
!
!> \brief
!! this subroutine is the main driver to perform the potential energy, forces 
!! calculation 
!
! ******************************************************************************
SUBROUTINE engforce

  USE config,                   ONLY :  natm , ntype, system , simu_cell, atypei ,natmi, atype, itype , tau_coul, tau_nonb, tau
  USE control,                  ONLY :  lnmlj , lcoulomb , lbmhft , lbmhftd, lmorse , lharm , longrange, lnon_bonded, iefgall_format
  USE non_bonded,               ONLY :  engforce_nmlj_pbc, engforce_bmhftd_pbc
  USE coulomb,                  ONLY :  task_coul, multipole_ES, pair_thole, pair_thole_distance
  USE field,                    ONLY :  doefield, doefg
  USE io,                       ONLY :  kunit_DIPFF, kunit_EFGALL , kunit_EFALL

  implicit none

  ! local 
  real(kind=dp) , allocatable :: ef ( : , : ) , efg ( : , : , : )
  real(kind=dp) , allocatable :: mu ( : , : )
  real(kind=dp) , allocatable :: theta ( : , : , : )
  logical :: didpim


  ! test purpose only
  ! harmonic oscillator ( test purpose )
  !if ( lharm ) then
  !  CALL engforce_harm
  !endif

  tau =0.0_dp
  ! =================================
  !    n-m lennard-jones potential 
  ! =================================
  if ( lnmlj ) then
    CALL engforce_nmlj_pbc
    tau = tau + tau_nonb
  endif

  ! =================================
  !   bmft(d) potentials (d:damping) 
  ! =================================
  if ( lbmhftd .or. lbmhft )  then
    CALL engforce_bmhftd_pbc
    tau = tau + tau_nonb
  endif

  ! =================================
  !   coulombic potential 
  ! =================================
  if ( lcoulomb ) then

     allocate( ef(3,natm) , efg(3,3,natm) , mu(3,natm) , theta(3,3,natm))
     allocate ( pair_thole ( natm ) )
     allocate ( pair_thole_distance ( natm ) )
!     pair_thole = 0
!     pair_thole_distance = 0.0_dp
!     mu=0.0_dp
!     theta=0.0_dp

     ! this subroutine set the total dipole moment from :
     ! static          (if given in control file dip TODO: read from DIPFF )
     ! wannier centers (if given in POSFF )
     ! induced         (if polar are set in control file)
!     CALL get_dipole_moments ( mu , theta , didpim )
!     theta=0.0_dp

     ! ====================================
     ! get all the electrostatic quantities
     ! ====================================
     CALL multipole_ES ( ef , efg , mu , theta , task_coul , damp_ind=.true. , &
                         do_efield=doefield , do_efg=doefg , do_forces=.true. , &
                         do_stress=.true. , do_rec=.true. , do_dir=.true. , do_strucfact =.false. , use_ckrskr = didpim )
!     theta=0.0_dp
!     mu_t     = mu
!     theta_t  = theta
 !    ef_t     = ef
 !    efg_t    = efg

     deallocate( ef , efg , mu , theta )
     deallocate ( pair_thole )
     deallocate ( pair_thole_distance )

     tau = tau + tau_coul

  endif


  ! ===========================
  !   other potentials ...
  ! ===========================
  ! CALL engforce_<other>_pbc


  return

END SUBROUTINE engforce


END MODULE engforce_driver 
! ===== fmV =====

