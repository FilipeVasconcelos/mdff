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
!#define debug_ES_rec
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE ewald_rec 

  USE constants,                        ONLY :  dp 
  USE mpimdff

  implicit none

CONTAINS

SUBROUTINE multipole_ES_rec ( u_rec , ef_rec, efg_rec , fx_rec , fy_rec , fz_rec , tau_rec , mu , theta , task , &
                              do_efield , do_efg , do_forces , do_stress , do_strucfact , use_ckrskr )

  USE constants,                ONLY :  imag, tpi
  USE config,                   ONLY :  natm, rx ,ry, rz, qia, simu_cell
  USE kspace,                   ONLY :  charge_density_k_q , charge_density_k_mu , struc_fact
  USE coulomb,                  ONLY :  alphaES, km_coul
!  USE time,                     ONLY :  fcoultimetot2_2
  USE dumb

  implicit none

  ! global
  real(kind=dp) :: u_rec
  real(kind=dp) :: ef_rec  (:,:)
  real(kind=dp) :: efg_rec (:,:,:)
  real(kind=dp) :: fx_rec  (:) , fy_rec (:) , fz_rec (:)
  real(kind=dp) :: tau_rec (:,:)
  real(kind=dp) :: mu      (:,:)
  real(kind=dp) :: theta   (:,:,:)
  logical       :: task(:), do_efield , do_efg , do_forces , do_stress , do_strucfact , use_ckrskr 

  ! local
  integer           :: ia , ik , ierr
  real(kind=dp)     :: qi
  real(kind=dp)     :: muix, muiy, muiz
  real(kind=dp)     :: thetaixx, thetaiyy, thetaizz
  real(kind=dp)     :: thetaixy, thetaixz, thetaiyz
  real(kind=dp)     :: kx   , ky   , kz , kk, Ak
  real(kind=dp)     :: rxi  , ryi  , rzi
  real(kind=dp)     :: fxij , fyij , fzij
  real(kind=dp)     :: str, k_dot_r ,  k_dot_mu , K_dot_Q , recarg, kcoe , rhonk_R , rhonk_I, recarg2
  real(kind=dp)     :: alpha2, tpi_V , fpi_V 
!  real(kind=dp)     :: rhonk_R_st_x , rhonk_R_st_y , rhonk_R_st_z 
!  real(kind=dp)     :: rhonk_I_st_x , rhonk_I_st_y , rhonk_I_st_z 
  real(kind=dp)     ,dimension (:), allocatable :: ckr , skr, tmp 
  real(kind=dp)     :: ckria , skria
  real(kind=dp)     :: tau_rec_tmp(3,3) 
  logical           :: ldip , lquad
  real(kind=dp)     :: ttt1,ttt2

  ! =================
  !  some constants 
  ! =================
  tpi_V  = tpi    / simu_cell%omega  ! 2pi / V
  fpi_V  = tpi_V  * 2.0_dp           ! 4pi / V
  alpha2 = alphaES * alphaES
  ldip = .false.
  if ( task(2) .or. task(3) ) ldip = .true.

!  if ( do_strucfact ) then
!  if ( .FALSE. ) then
!    CALL struc_fact ( km_coul )
!  endif
  
  !if ( .not. use_ckrskr ) allocate( ckr ( natm ) , skr (natm) ) 
!  if ( .TRUE. ) allocate( ckr ( natm ) , skr (natm) , tmp(natm)) 
  allocate( ckr ( natm ) , skr (natm) , tmp(natm)) 

  ! ==============================================
  !            reciprocal space part
  ! ==============================================
  kpoint : do ik = km_coul%kpt_dec%istart, km_coul%kpt_dec%iend
    
    if (km_coul%kptk(ik) .eq. 0.0_dp ) cycle

    kx     = km_coul%kptx(ik)
    ky     = km_coul%kpty(ik)
    kz     = km_coul%kptz(ik)
    kk     = km_coul%kptk(ik)
    Ak     = km_coul%Ak( ik )
    kcoe   = km_coul%kcoe(ik) 

!    if ( .not. use_ckrskr ) then
!    if ( .TRUE. ) then
#ifdef MPI
      ttt1 = MPI_WTIME(ierr)
#endif
      rhonk_R = 0.0_dp
      rhonk_I = 0.0_dp
      atom1 : do ia = 1, natm
        qi  = qia ( ia )
        rxi = rx  ( ia )
        ryi = ry  ( ia )
        rzi = rz  ( ia )
        k_dot_r  = ( kx * rxi + ky * ryi + kz * rzi )
        ckr(ia)  = COS(k_dot_r) 
        skr(ia)  = SIN(k_dot_r)
        rhonk_R    = rhonk_R + qi * ckr(ia) 
        rhonk_I    = rhonk_I + qi * skr(ia)  ! rhon_R + i rhon_I
        if ( .not. ldip ) cycle
        muix = mu ( 1 , ia )
        muiy = mu ( 2 , ia )
        muiz = mu ( 3 , ia )
        k_dot_mu = ( muix * kx + muiy * ky + muiz * kz )
        rhonk_R    = rhonk_R - k_dot_mu * skr(ia) 
        rhonk_I    = rhonk_I + k_dot_mu * ckr(ia) ! rhon_R + i rhon_I
  !      if ( .not. lquad ) cycle
  !      thetaixx = theta ( 1 , 1 , ia)
  !      thetaiyy = theta ( 2 , 2 , ia)
  !      thetaizz = theta ( 3 , 3 , ia)
  !      thetaixy = theta ( 1 , 2 , ia)
  !      thetaixz = theta ( 1 , 3 , ia)
  !      thetaiyz = theta ( 2 , 3 , ia)
  !      K_dot_Q =           thetaixx * kx * kx + thetaixy * kx * ky + thetaixz * kx * kz
  !      K_dot_Q = K_dot_Q + thetaixy * kx * ky + thetaiyy * ky * ky + thetaiyz * ky * kz
  !      K_dot_Q = K_dot_Q + thetaixz * kz * kx + thetaiyz * ky * kz + thetaizz * kz * kz
  !      K_dot_Q = K_dot_Q / 3.0_dp
  !      rhonk_R    = rhonk_R - K_dot_Q * ckr(ia)
  !      rhonk_I    = rhonk_I - K_dot_Q * skr(ia)  ! rhon_R + i rhon_I
      enddo atom1

#ifdef MPI
      ttt2 = MPI_WTIME(ierr)
      t12 = t12 + (ttt2-ttt1)
#endif
!    else 
!      rhonk_R = 0.0_dp
!      rhonk_I = 0.0_dp
!      km_coul%rhon_R = 0.0_dp
!      km_coul%rhon_I = 0.0_dp
!      if ( task(1) ) CALL charge_density_k_q  ( km_coul , ik )  
!      if ( ldip    ) CALL charge_density_k_mu ( km_coul , mu , ik )  
!      rhonk_R =  km_coul%rhon_R(ik)
!      rhonk_I =  km_coul%rhon_I(ik)
!    endif

#ifdef MPI
    ttt1 = MPI_WTIME(ierr)
#endif


    ! ===================================================
    ! potential energy 
    ! ===================================================
    str    = (rhonk_R*rhonk_R + rhonk_I*rhonk_I) * Ak  
    u_rec   = u_rec   + str


#ifdef debug    
    write(1000000,'(a,i,4e16.8)') 'qqqq',ik,rhonk_R,rhonk_I,str,Ak
#endif    

    atom2 : do ia = 1 , natm
!      if ( use_ckrskr ) then
!      if ( .FALSE. ) then
!        ckria = km_coul%ckr(ia,ik)  
!        skria = km_coul%skr(ia,ik)  
!      else
        ckria = ckr(ia)
        skria = skr(ia)
!      endif

      if ( do_efg .or. ( do_forces .and. ldip ) )  then
        recarg2 = Ak * ( rhonk_R * ckria + rhonk_I * skria )
      endif

      if ( do_efield .or. do_forces ) then
        recarg  = Ak * ( rhonk_I * ckria - rhonk_R * skria )
        fxij = kx * recarg
        fyij = ky * recarg
        fzij = kz * recarg
      endif

      if ( do_efield ) then
        ef_rec ( 1 , ia ) = ef_rec ( 1 , ia ) - fxij
        ef_rec ( 2 , ia ) = ef_rec ( 2 , ia ) - fyij
        ef_rec ( 3 , ia ) = ef_rec ( 3 , ia ) - fzij
      endif

      ! electric field gradient
      if ( do_efg ) then
        recarg2 = Ak * ( rhonk_R * ckria + rhonk_I * skria )
        efg_rec ( 1 , 1 , ia ) = efg_rec ( 1 , 1 , ia ) + kx * kx * recarg2
        efg_rec ( 2 , 2 , ia ) = efg_rec ( 2 , 2 , ia ) + ky * ky * recarg2
        efg_rec ( 3 , 3 , ia ) = efg_rec ( 3 , 3 , ia ) + kz * kz * recarg2
        efg_rec ( 1 , 2 , ia ) = efg_rec ( 1 , 2 , ia ) + kx * ky * recarg2
        efg_rec ( 1 , 3 , ia ) = efg_rec ( 1 , 3 , ia ) + kx * kz * recarg2
        efg_rec ( 2 , 3 , ia ) = efg_rec ( 2 , 3 , ia ) + ky * kz * recarg2
      endif

      if ( do_forces ) then
        qi  = qia ( ia )
        ! charges
        fx_rec ( ia ) = fx_rec ( ia ) - qi * fxij
        fy_rec ( ia ) = fy_rec ( ia ) - qi * fyij
        fz_rec ( ia ) = fz_rec ( ia ) - qi * fzij
        ! dipoles ( k_alpha * Ak * mu.k * recarg ) 
        if ( .not. ldip ) cycle
        muix = mu ( 1 , ia )
        muiy = mu ( 2 , ia )
        muiz = mu ( 3 , ia )
        k_dot_mu  =( muix * kx + muiy * ky + muiz * kz  ) * recarg2
        fx_rec ( ia ) = fx_rec ( ia ) + kx * k_dot_mu
        fy_rec ( ia ) = fy_rec ( ia ) + ky * k_dot_mu
        fz_rec ( ia ) = fz_rec ( ia ) + kz * k_dot_mu
   !     if ( .not. lquad ) cycle
   !     thetaixx = theta ( 1 , 1 , ia)
   !     thetaiyy = theta ( 2 , 2 , ia)
   !     thetaizz = theta ( 3 , 3 , ia)
   !     thetaixy = theta ( 1 , 2 , ia)
   !     thetaixz = theta ( 1 , 3 , ia)
   !     thetaiyz = theta ( 2 , 3 , ia)
   !     K_dot_Q =           thetaixx * kx * kx + thetaixy * kx * ky + thetaixz * kx * kz
   !     K_dot_Q = K_dot_Q + thetaixy * kx * ky + thetaiyy * ky * ky + thetaiyz * ky * kz
   !     K_dot_Q = K_dot_Q + thetaixz * kz * kx + thetaiyz * ky * kz + thetaizz * kz * kz
   !     K_dot_Q = K_dot_Q * recarg
   !     fx_rec ( ia ) = fx_rec ( ia ) + kx * K_dot_Q
   !     fy_rec ( ia ) = fy_rec ( ia ) + ky * K_dot_Q
   !     fz_rec ( ia ) = fz_rec ( ia ) + kz * K_dot_Q
      endif

    enddo atom2

    if ( do_stress ) then
     ! stress tensor symmetric !
     ! keep it out from the ia loop ! stupid bug ! 
     tau_rec(1,1) = tau_rec(1,1) + ( 1.0_dp - kcoe * kx * kx ) * str
     tau_rec(1,2) = tau_rec(1,2) -            kcoe * kx * ky   * str
     tau_rec(1,3) = tau_rec(1,3) -            kcoe * kx * kz   * str
     tau_rec(2,1) = tau_rec(2,1) -            kcoe * ky * kx   * str
     tau_rec(2,2) = tau_rec(2,2) + ( 1.0_dp - kcoe * ky * ky ) * str
     tau_rec(2,3) = tau_rec(2,3) -            kcoe * ky * kz   * str
     tau_rec(3,1) = tau_rec(3,1) -            kcoe * kz * kx   * str
     tau_rec(3,2) = tau_rec(3,2) -            kcoe * kz * ky   * str
     tau_rec(3,3) = tau_rec(3,3) + ( 1.0_dp - kcoe * kz * kz ) * str
! =================================================================================================================================
! Some differences inn stress tensor result for PIM were found ... seems to come from the reciprocal part
! the code show a further contribution ... but the implementation is quite different and is correct for all the other quantities
! checked so far with charges and dipoles. 
! CP2K
! ! The second one can be written in the following way
!          f0 = 2.0_dp * gauss
!          pv_tmp(1,1) = pv_tmp(1,1) + f0 * pw_grid%g(1,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(1,2) = pv_tmp(1,2) + f0 * pw_grid%g(1,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(1,3) = pv_tmp(1,3) + f0 * pw_grid%g(1,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(2,1) = pv_tmp(2,1) + f0 * pw_grid%g(2,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(2,2) = pv_tmp(2,2) + f0 * pw_grid%g(2,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(2,3) = pv_tmp(2,3) + f0 * pw_grid%g(2,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(3,1) = pv_tmp(3,1) + f0 * pw_grid%g(3,gpt) * REAL(summe_st(1,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(3,2) = pv_tmp(3,2) + f0 * pw_grid%g(3,gpt) * REAL(summe_st(2,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
!          pv_tmp(3,3) = pv_tmp(3,3) + f0 * pw_grid%g(3,gpt) * REAL(summe_st(3,gpt) * CONJG(summe_ef(gpt)),KIND=dp)
! =================================================================================================================================
     ! ( rhonk_R_stx + i rhonk_I_stx ) * ( rhonk_R - i rhonk_I )  = rhonk_R_stx * rhonk_R + rhonk_I_stx * rhonk_I 
!     f0 = 2.0_dp * Ak
!     tau_rec(1,1) = tau_rec(1,1) + f0 * kx * ( rhonk_R_st_x * rhonk_R + rhonk_I_st_x * rhonk_I ) 
!     tau_rec(1,2) = tau_rec(1,2) + f0 * kx * ( rhonk_R_st_y * rhonk_R + rhonk_I_st_y * rhonk_I )
!     tau_rec(1,3) = tau_rec(1,3) + f0 * kx * ( rhonk_R_st_z * rhonk_R + rhonk_I_st_z * rhonk_I )
!     tau_rec(2,1) = tau_rec(2,1) + f0 * ky * ( rhonk_R_st_x * rhonk_R + rhonk_I_st_x * rhonk_I )
!     tau_rec(2,2) = tau_rec(2,2) + f0 * ky * ( rhonk_R_st_y * rhonk_R + rhonk_I_st_y * rhonk_I )
!     tau_rec(2,3) = tau_rec(2,3) + f0 * ky * ( rhonk_R_st_z * rhonk_R + rhonk_I_st_z * rhonk_I )
!     tau_rec(3,1) = tau_rec(3,1) + f0 * kz * ( rhonk_R_st_x * rhonk_R + rhonk_I_st_x * rhonk_I )
!     tau_rec(3,2) = tau_rec(3,2) + f0 * kz * ( rhonk_R_st_y * rhonk_R + rhonk_I_st_y * rhonk_I )
!     tau_rec(3,3) = tau_rec(3,3) + f0 * kz * ( rhonk_R_st_z * rhonk_R + rhonk_I_st_z * rhonk_I )
   endif
#ifdef MPI
    ttt2 = MPI_WTIME(ierr)
#endif

  enddo kpoint

  ! "half" mesh
  u_rec   = u_rec   * 2.0_dp
  if ( do_efield ) ef_rec  = ef_rec  * 2.0_dp
  if ( do_efg    ) efg_rec = efg_rec * 2.0_dp
  if ( do_forces ) then
    fx_rec  = fx_rec  * 2.0_dp
    fy_rec  = fy_rec  * 2.0_dp
    fz_rec  = fz_rec  * 2.0_dp
  endif
  if ( do_stress) tau_rec = tau_rec * 2.0_dp

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_rec )
  if ( do_efield ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( ef_rec(1,:) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( ef_rec(2,:) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( ef_rec(3,:) , natm )
  endif
  if ( do_efg ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( 1 , 1 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( 2 , 2 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( 3 , 3 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( 1 , 2 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( 1 , 3 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_rec ( 2 , 3 , : ) , natm )
  endif
  if ( do_forces ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( fx_rec , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( fy_rec , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( fz_rec , natm )
  endif
  if ( do_stress ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( tau_rec ( 1 , : ) , 3 )
    CALL MPI_ALL_REDUCE_DOUBLE ( tau_rec ( 2 , : ) , 3 )
    CALL MPI_ALL_REDUCE_DOUBLE ( tau_rec ( 3 , : ) , 3 )
  endif

  ! Symmetrize the tensor
  if ( do_stress ) then
       tau_rec_tmp  = tau_rec
       tau_rec(1,2) = (tau_rec_tmp(1,2) + tau_rec_tmp(2,1))*0.5_dp
       tau_rec(1,3) = (tau_rec_tmp(1,3) + tau_rec_tmp(3,1))*0.5_dp
       tau_rec(2,3) = (tau_rec_tmp(2,3) + tau_rec_tmp(3,2))*0.5_dp
       tau_rec(2,1) = tau_rec(1,2)
       tau_rec(3,1) = tau_rec(1,3)
       tau_rec(3,2) = tau_rec(2,3)
  endif

  ! ======================================================
  ! remark on the unit :
  ! 1/(4*pi*epislon_0) = 1 => epsilon_0 = 1/4pi
  ! ======================================================
  u_rec   =   u_rec   * tpi_V
  if ( do_efield ) ef_rec  =   ef_rec  * fpi_V 
  if ( do_efg    ) efg_rec =   efg_rec * fpi_V
  if ( do_stress ) tau_rec =   tau_rec * tpi_V / simu_cell%omega
  if ( do_forces ) then
    fx_rec  =   fx_rec  * fpi_V
    fy_rec  =   fy_rec  * fpi_V
    fz_rec  =   fz_rec  * fpi_V
  endif

!  if ( .TRUE. ) deallocate( ckr , skr ) 
  deallocate( ckr , skr , tmp) 

  return

END SUBROUTINE multipole_ES_rec


END MODULE ewald_rec 
! ===== fmV =====

