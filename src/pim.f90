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
MODULE pim

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode

  implicit none

!  logical, SAVE     :: lwrite_dip_wfc    !< write dipoles from wannier centers to file
!  logical, SAVE     :: lwrite_dip        !< write dipoles 
!  logical, SAVE     :: lwrite_quad       !< write quadrupoles to QUADFF
!  logical, SAVE     :: lwrite_efg        !< write electric field gradient to EFGALL
!  logical, SAVE     :: lwrite_ef         !< write electric field s to EFALL
  logical, SAVE     :: ldip_wfc          !< calculate electrostatic contribution from dipolar momemt coming from wfc


  ! ============================================================  
  !                    force field type info 
  ! ============================================================  

  ! type dependent properties
  real(kind=dp)    :: poldip   ( ntypemax , 3 , 3 )    !< dipole     polarizability if ldip_polar( it ) = .true. 
  real(kind=dp)    :: poldip_iso ( ntypemax )          !< isotropic dipole polarizability if ldip_polar( it ) = .true.
  real(kind=dp)    :: polquad  ( ntypemax , 3 , 3 , 3 )!< quadrupole polarizability if lquad_polar( it ) = .true.
  real(kind=dp)    :: polquad_iso ( ntypemax )         !< isotropic quadrupole polarizability if ldip_polar( it ) = .true. 


  ! =====================================================
  !                 polarizability  
  ! =====================================================  
  real(kind=dp)    :: omegakO
  real(kind=dp)    :: conv_tol_ind                                   !< convergence tolerance of the scf induced dipole calculation
  integer          :: min_scf_pol_iter                               !< 
  integer          :: max_scf_pol_iter                               !< 
  integer          :: extrapolate_order                              !< 
  logical          :: ldip_polar   ( ntypemax )                          !< induced moment from pola. Is this type of ion polarizable ?
  logical          :: ldip_damping ( ntypemax , ntypemax , ntypemax) !< dipole damping 
  real(kind=dp)    :: pol_damp_b ( ntypemax, ntypemax,ntypemax )     !< dipole damping : parameter b [length]^-1
  real(kind=dp)    :: pol_damp_c ( ntypemax, ntypemax,ntypemax )     !< dipole damping : parameter c no units
  integer          :: pol_damp_k ( ntypemax, ntypemax,ntypemax )     !< dipole damping : Tang-Toennies function order
  logical          :: lquad_polar    ( ntypemax )                    !< induced quadrupole from pola


  character(len=4) :: algo_ext_dipole                  !< set the algorithm used to get induced moments from polarization 
  character(len=4) :: algo_ext_dipole_allowed(2)
  data                algo_ext_dipole_allowed       / 'poly','aspc'/

  character(len=11) :: algo_moment_from_pola            !< set the algorithm used to get induced moments from polarization 
  character(len=11) :: algo_moment_from_pola_allowed(6) !< set the algorithm used to get induced moments from polarization 
  data                 algo_moment_from_pola_allowed / 'scf' , 'scf_kO_v1' , 'scf_kO_v2' , 'scf_kO_v3' , 'scf_kO_v4_1' , 'scf_kO_v4_2'/  !! scf ( self consistent ) 
  
  integer          :: lwfc     ( ntypemax )            !< moment from wannier centers 
  real(kind=dp)    :: rcut_wfc                         !< radius cut-off for WFs searching

END MODULE pim
