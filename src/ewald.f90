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
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE ewald 

  USE constants,                        ONLY :  dp 
  USE mpimdff

  implicit none

CONTAINS

! *********************** SUBROUTINE multipole_ES ******************************
!> \brief
!! This subroutine calculates electric field, electric field gradient, 
!! potential energy, virial, electric potential and forces at ions in
!! a multipole expansion by Ewald summation
!
!> \param[in]  mu electric dipole at ions
!> \param[out] ef electric field
!
!> \todo
!! make it more condensed 
! ******************************************************************************
SUBROUTINE multipole_ES ( ef , efg , mu , theta , task , damp_ind , &
                          do_efield , do_efg , do_forces , do_stress , do_rec , do_dir , do_strucfact , use_ckrskr )

  USE ewald_dir,        ONLY :  multipole_ES_dir
  USE ewald_rec,        ONLY :  multipole_ES_rec
  USE control,          ONLY :  lsurf
  USE constants,        ONLY :  tpi , piroot, coul_unit, press_unit
  USE config,           ONLY :  natm, qia, rx ,ry ,rz, simu_cell, fx, fy, fz, tau, tau_coul, atype 
  USE coulomb,          ONLY :  alphaES
  USE thermodynamic,    ONLY :  u_coul, u_pol, pvirial_coul
  USE time,             ONLY :  fcoultimetot1 , fcoultimetot2
  USE tensors_rk,       ONLY :  interaction_dd
  USE dumb

  implicit none

  ! global 
  real(kind=dp)     :: ef     ( : , : )
  real(kind=dp)     :: efg    ( : , : , : )
  real(kind=dp)     :: mu     ( : , : )
  real(kind=dp)     :: theta  ( : , : , : )
  logical           :: task   ( : )
  logical           :: damp_ind , do_efield , do_efg, do_forces, do_stress, do_rec , do_dir , do_strucfact , use_ckrskr 

  ! local 
  integer         :: ia , ierr
  real(kind=dp)                                :: u_dir , u_rec , u_surf , u_self
  real(kind=dp)                                :: u_surf_qq , u_surf_qd , u_surf_dd
  real(kind=dp), dimension(:,:)  , allocatable :: ef_dir, ef_rec, ef_surf, ef_self
  real(kind=dp), dimension(:,:,:), allocatable :: efg_dir, efg_rec, efg_self
  real(kind=dp), dimension(:)    , allocatable :: fx_coul , fy_coul , fz_coul
  real(kind=dp), dimension(:)    , allocatable :: fx_dir , fy_dir , fz_dir
  real(kind=dp), dimension(:)    , allocatable :: fx_rec , fy_rec , fz_rec
  real(kind=dp), dimension(:)    , allocatable :: fx_surf , fy_surf , fz_surf
  real(kind=dp) :: tau_dir( 3 , 3 )
  real(kind=dp) :: tau_rec( 3 , 3 )
  real(kind=dp) :: qtot ( 3 ) , qsq , mutot ( 3 ) , musq , qmu_sum ( 3 ) , thetasq
  real(kind=dp) :: tpi_V, tpi_3V , fpi_3V , alpha2 , selfa , selfa2, selfa3
  real(kind=dp) :: ttt1, ttt2 
  real(kind=dp) :: u_self_1 , u_self_2, u_self_3

#ifdef debug_ES
        call print_config_sample(0,0)
#endif

  allocate( ef_dir  (3, natm)  , ef_rec(3 , natm)   , ef_surf(3,natm) ,ef_self(3,natm) )
  allocate( efg_dir (3,3,natm), efg_rec(3,3, natm), efg_self(3,3, natm) )
  allocate( fx_coul (natm)    , fy_coul (natm)   , fz_coul (natm) )
  allocate( fx_dir  (natm)    , fy_dir  (natm)   , fz_dir  (natm) )
  allocate( fx_rec  (natm)    , fy_rec  (natm)   , fz_rec  (natm) )
  allocate( fx_surf (natm)    , fy_surf (natm)   , fz_surf (natm) )
  ef_dir   = 0.0_dp;  ef_rec   = 0.0_dp;  ef_surf  = 0.0_dp; ef_self= 0.0_dp
  efg_dir  = 0.0_dp;  efg_rec  = 0.0_dp;  efg_self = 0.0_dp
  fx_dir   = 0.0_dp;  fy_dir   = 0.0_dp;  fz_dir   = 0.0_dp
  fx_rec   = 0.0_dp;  fy_rec   = 0.0_dp;  fz_rec   = 0.0_dp
  fx_surf  = 0.0_dp;  fy_surf  = 0.0_dp;  fz_surf  = 0.0_dp
  u_dir    = 0.0_dp;  u_rec    = 0.0_dp;  u_self   = 0.0_dp;  u_surf   = 0.0_dp
  tau_dir  = 0.0_dp;  tau_rec  = 0.0_dp; tau_coul =0.0_dp


  ! ==================================================
  !  total charge / moment / square 
  ! ==================================================
  ! note here that qtot is a dipole moment  !!!! 
  mutot = 0.0_dp
  musq  = 0.0_dp
  qtot  = 0.0_dp
  qsq   = 0.0_dp
  thetasq = 0.0_dp
  do ia = 1 , natm
    mutot ( 1 ) = mutot ( 1 ) + mu ( 1 , ia )
    mutot ( 2 ) = mutot ( 2 ) + mu ( 2 , ia )
    mutot ( 3 ) = mutot ( 3 ) + mu ( 3 , ia )
    musq = musq + ( mu ( 1 , ia ) * mu ( 1 , ia ) +  mu ( 2 , ia ) * mu ( 2 , ia ) +  mu ( 3 , ia ) * mu ( 3 , ia ) )
    qtot ( 1 ) = qtot ( 1 ) + qia ( ia ) * rx ( ia )
    qtot ( 2 ) = qtot ( 2 ) + qia ( ia ) * ry ( ia )
    qtot ( 3 ) = qtot ( 3 ) + qia ( ia ) * rz ( ia )
    qsq = qsq + qia ( ia ) * qia ( ia )
    thetasq= thetasq + (theta(1,1,ia)+theta(2,2,ia)+theta(3,3,ia))*qia(ia)
  enddo
  !WRITE(stdout,'(a,2f18.10)') 'thetasq',thetasq,theta(1,1,1)
  qmu_sum ( 1 ) = qtot ( 1 ) + mutot ( 1 )
  qmu_sum ( 2 ) = qtot ( 2 ) + mutot ( 2 )
  qmu_sum ( 3 ) = qtot ( 3 ) + mutot ( 3 )

  ! ===============
  !    constants
  ! ===============
  tpi_V  = tpi    / simu_cell%omega  ! 2pi / V
  tpi_3V = tpi_V  / 3.0_dp           ! 2pi / 3V 
  fpi_3V = tpi_3V * 2.0_dp           ! 4pi / 3V
  alpha2 = alphaES * alphaES
  selfa  = alphaES / piroot
  selfa2 = 2.0_dp * selfa * alpha2 / 3.0_dp
  selfa3 = 2.0_dp * selfa2         / 3.0_dp

  ! ==============================================
  !        direct space part
  ! ==============================================
  if ( do_dir ) then
#ifdef MPI
    ttt1 = MPI_WTIME(ierr)
#endif
      CALL multipole_ES_dir ( u_dir , ef_dir, efg_dir, fx_dir , fy_dir , fz_dir , tau_dir , mu , theta , task , damp_ind , & 
                              do_efield , do_efg , do_forces , do_stress )
#ifdef MPI
    ttt2 = MPI_WTIME(ierr)
    fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 )  
#endif
  endif


  ! ==============================================
  !        reciprocal space part
  ! ==============================================
  if ( do_rec ) then
#ifdef MPI
    ttt1 = MPI_WTIME(ierr)
#endif
    CALL multipole_ES_rec ( u_rec , ef_rec , efg_rec , fx_rec , fy_rec , fz_rec , tau_rec , mu , theta , task , & 
                            do_efield , do_efg , do_forces , do_stress , do_strucfact , use_ckrskr )
#ifdef MPI
    ttt2 = MPI_WTIME(ierr)
    fcoultimetot2 = fcoultimetot2 + ( ttt2 - ttt1 )  
#endif
  endif

  ! ====================================================== 
  !              Surface contribution ????? 
  ! ====================================================== 
  ! spherical symmetry
  ! electrostatic energy and virial
  ! qq
  u_surf_qq = qtot ( 1 ) * qtot ( 1 ) + qtot ( 2 ) * qtot ( 2 ) + qtot ( 3 ) * qtot ( 3 )
  u_surf_dd = mutot ( 1 ) * mutot ( 1 ) + mutot ( 2 ) * mutot ( 2 ) + mutot ( 3 ) * mutot ( 3 )
  u_surf_qd = 2.0_dp * ( qtot ( 1 ) * mutot ( 1 ) + qtot ( 2 ) * mutot ( 2 ) + qtot ( 3 ) * mutot ( 3 ) )
  u_surf    = u_surf_qq + u_surf_qd + u_surf_dd
  u_surf    = u_surf * tpi_3V

  ! potential, field , forces ( no contrib to efg ) 
  do ia = 1 , natm
    ef_surf ( 1 , ia ) = qmu_sum ( 1 )
    ef_surf ( 2 , ia ) = qmu_sum ( 2 )
    ef_surf ( 3 , ia ) = qmu_sum ( 3 )
    !fx_surf ( ia ) = qia ( ia ) * qmu_sum ( 1 )
    !fy_surf ( ia ) = qia ( ia ) * qmu_sum ( 2 )
    !fz_surf ( ia ) = qia ( ia ) * qmu_sum ( 3 )
  enddo
  !fx_surf  = - fx_surf  * fpi_3V
  !fy_surf  = - fy_surf  * fpi_3V
  !fz_surf  = - fz_surf  * fpi_3V
  ef_surf  = - ef_surf  * fpi_3V

  ! ====================================================== 
  !              Self contribution 
  ! ====================================================== 
  ! electrostic energy 
  u_self_1   =  - selfa  * qsq                ! q-q
  u_self_2   =  - selfa2 * musq               ! mu-mu
  u_self_3   =  - selfa3 * thetasq            ! theta-theta
  u_self     = u_self_1 + u_self_2 + u_self_3
  do ia = 1 , natm
    ef_self( 1 , ia ) = 2.0_dp * selfa2 * mu ( 1 , ia )
    ef_self( 2 , ia ) = 2.0_dp * selfa2 * mu ( 2 , ia )
    ef_self( 3 , ia ) = 2.0_dp * selfa2 * mu ( 3 , ia )
    efg_self ( 1 , 1 , ia ) =  - 2.0_dp * selfa2 * qia ( ia )
    efg_self ( 2 , 2 , ia ) =  - 2.0_dp * selfa2 * qia ( ia )
    efg_self ( 3 , 3 , ia ) =  - 2.0_dp * selfa2 * qia ( ia )
  enddo


  ! =====================================================
  !                  TOTAL and units
  !  TODO : electric field has not the output unit !!!!
  !         it keeps internal units because of induced_moment 
  !         subroutine. Make dipole, polarisation, electric 
  !         field more coherent !!!
  ! =====================================================

  if ( lsurf ) then
    u_coul   =      ( u_dir   + u_rec   + u_surf   + u_self  + u_pol  ) * coul_unit
    ef       =      ( ef_dir  + ef_rec  + ef_surf  + ef_self          ) 
    efg      =      ( efg_dir + efg_rec + efg_self                    ) * coul_unit
    tau_coul =      ( tau_dir + tau_rec                               ) * coul_unit / press_unit
    fx       = fx + ( fx_rec  + fx_dir  + fx_surf                     ) * coul_unit
    fy       = fy + ( fy_rec  + fy_dir  + fy_surf                     ) * coul_unit
    fz       = fz + ( fz_rec  + fz_dir  + fz_surf                     ) * coul_unit
    tau      = tau  + tau_coul
  else
    u_coul   =      ( u_dir   + u_rec   + u_self  + u_pol  ) * coul_unit
    ef       =      ( ef_dir  + ef_rec  + ef_self          ) 
    efg      =      ( efg_dir + efg_rec + efg_self         ) * coul_unit
    tau_coul =      ( tau_dir + tau_rec                    ) * coul_unit / press_unit
    fx       = fx + ( fx_rec  + fx_dir                     ) * coul_unit
    fy       = fy + ( fy_rec  + fy_dir                     ) * coul_unit
    fz       = fz + ( fz_rec  + fz_dir                     ) * coul_unit
    tau      = tau  + tau_coul
  endif

  
#ifdef debug_ES_field_forces
  WRITE ( stdout , '(a)' )     'Electric field at atoms :                           Forces at atoms :'
  do ia = 1 , natm
   WRITE ( stdout , '(i5,a3,a,3f18.10,10x,a,3f18.10)' ) &
   ia,atype(ia),' Efield   = ', ef ( 1 , ia )  , ef ( 2 , ia ) , ef ( 3 , ia ),'    f   = ', fx ( ia )  , fy ( ia ) , fz ( ia )
 enddo
 do ia = 1 , natm
   WRITE ( stdout , '(i5,a3,a,3f18.10,10x,a,3f18.10)' ) &
   ia,atype(ia),' ef_dir   = ', ef_dir ( 1, ia )  , ef_dir ( 2 , ia ) , ef_dir    ( 3  , ia ), '  f_dir   = ',fx_dir ( ia)  , fy_dir ( ia ) , fz_dir    ( ia  )
enddo
do ia = 1 , natm
  WRITE ( stdout , '(i5,a3,a,3f18.10,10x,a,3f18.10)' ) &
   ia,atype(ia),' ef_rec   = ', ef_rec ( 1, ia )  , ef_rec ( 2 , ia ) , ef_rec   ( 3 , ia ), '  f_rec   = ',fx_rec ( ia)  , fy_rec ( ia ) , fz_rec    ( ia  )
 enddo
 do ia = 1 , natm
   WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
   ia,atype(ia),' ef_surf  = ', ef_surf ( 1 , ia )  , ef_surf ( 2 , ia ) ,   ef_surf ( 3 , ia )
 enddo
 do ia = 1 , natm
   WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
   ia,atype(ia),' ef_self  = ', ef_self ( 1 , ia )  , ef_self ( 2 ,ia ) ,   ef_self ( 3 , ia )
 enddo
#endif

#ifdef debug_ES_energy
 WRITE ( stdout , '(6(a,f16.8))' ) ,' u_dir      = ', u_dir  * coul_unit , &
                                    ' u_rec      = ', u_rec  * coul_unit , &
                                    ' u_surf     = ', u_surf * coul_unit , & 
                                    ' u_self     = ', u_self * coul_unit , &
                                    ' u_pol      = ', u_pol  * coul_unit , &
                                    ' u_coul     = ', u_coul
  write(stdout , '(a)') 'self energies :'
  write(stdout , '(a,f12.6)') 'q-q         = ',u_self_1 
  write(stdout , '(a,f12.6)') 'µ-µ         = ',u_self_2 
  write(stdout , '(a,f12.6)') 'Θ-Θ         = ',u_self_3
#endif

#ifdef debug_ES_stress
  tau_dir  = tau_dir  / press_unit * coul_unit
  tau_rec  = tau_rec  / press_unit * coul_unit
  CALL print_tensor( tau_dir  ( : , : )     , 'TAU_DIR ' )
  CALL print_tensor( tau_rec  ( : , : )     , 'TAU_REC ' )
  CALL print_tensor( tau_coul ( : , : )     , 'TAU_COUL' )
#endif

#ifdef debug_ES_efg
  CALL print_tensor( efg_dir  ( : , : , 1 ) , 'EFG_DIRN' )
  CALL print_tensor( efg_rec  ( : , : , 1 ) , 'EFG_RECN' )
  CALL print_tensor( efg_self ( : , : , 1 ) , 'EFG_SELN' )
  CALL print_tensor( efg      ( : , : , 1 ) , 'EFG_TOTN' )
#endif



  deallocate( ef_dir  , ef_rec  , ef_surf ,ef_self)
  deallocate( efg_dir , efg_rec , efg_self)
  deallocate( fx_coul , fy_coul , fz_coul )
  deallocate( fx_dir  , fy_dir  , fz_dir  )
  deallocate( fx_rec  , fy_rec  , fz_rec  )
  deallocate( fx_surf , fy_surf , fz_surf )

  pvirial_coul = 1.0_dp / 3.0_dp * ( tau_coul(1,1) + tau_coul(2,2) + tau_coul(3,3) ) * press_unit

  return

END SUBROUTINE multipole_ES

END MODULE ewald 
! ===== fmV =====

