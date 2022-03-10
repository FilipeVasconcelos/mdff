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
!#define debug_ES_dir
#define debug_thole
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE ewald_dir 

  USE constants,                        ONLY :  dp 
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  USE mpimdff

  implicit none

CONTAINS

SUBROUTINE multipole_ES_dir ( u_dir , ef_dir , efg_dir , fx_dir , fy_dir , fz_dir , tau_dir , mu , theta , & 
                              task , damp_ind , do_efield , do_efg , do_forces , do_stress )

  USE io,                       ONLY :  ioprint
  USE pim,                      ONLY :  ldip_damping, pol_damp_b, pol_damp_c, pol_damp_k
  USE control,                  ONLY :  lvnlist
  USE coulomb,                  ONLY :  alphaES, thole_functions, thole_param, thole_function_type, pair_thole, pair_thole_distance 
  USE config,                   ONLY :  natm, ntype,simu_cell, qia, rx ,ry ,rz ,itype , atom_dec, verlet_coul , atype, poldipia, ipolar
  USE constants,                ONLY :  piroot
  USE cell,                     ONLY :  kardir , dirkar
  USE tensors_rk,               ONLY :  tensor_rank0, tensor_rank1, tensor_rank2, tensor_rank3 , tensor_rank4 , tensor_rank5
  USE tt_damp,                  ONLY :  TT_damping_functions
  USE dumb
 
 
  implicit none

  ! global
  real(kind=dp) :: u_dir 
  real(kind=dp) :: ef_dir   ( : , : )
  real(kind=dp) :: efg_dir  ( : , : , : )
  real(kind=dp) :: fx_dir   ( : ) , fy_dir ( : ) , fz_dir ( : )
  real(kind=dp) :: tau_dir  ( : , : )
  real(kind=dp) :: mu       ( : , :  )
  real(kind=dp) :: theta    ( : , :  , : )
  logical       :: task ( : ) , damp_ind, do_efield , do_efg , do_forces , do_stress, store_interaction

  ! local 
  integer       :: ia , ja , ita, jta, j1 , jb ,je , i , j , k, l , m , it1,it2 , ierr
  real(kind=dp) :: qi, qj , qij , u_damp 
  real(kind=dp) :: mui(3)
  real(kind=dp) :: muj(3)
  real(kind=dp) :: thetai(3,3)
  real(kind=dp) :: thetaj(3,3)
  real(kind=dp) :: cutsq
  real(kind=dp) :: rxi  , ryi  , rzi
  real(kind=dp) :: rxj  , ryj  , rzj
  real(kind=dp) :: rij(3)
  real(kind=dp) :: sij(3)
  real(kind=dp) :: fij(3)
  real(kind=dp) :: d , d2 , d3  , d5 , d7 , d9
  real(kind=dp) :: dm1 , dm3 , dm5 , dm7 , dm9, dm11
  real(kind=dp) :: F0 , F1 , F2 , F3 , F4 , F5
  real(kind=dp) :: F1d , F2d 
  real(kind=dp) :: F1d2 , F2d2 
  real(kind=dp) :: alpha2 , alpha3 , alpha5 , alpha7 , alpha9 , expon 
  real(kind=dp), external :: errfc
  real(kind=dp), external :: errf
  real(kind=dp) :: fdamp , fdampdiff
  real(kind=dp) :: fdamp2 , fdampdiff2
  real(kind=dp) :: onesixth, twothird , sthole, uthole , vthole, vthole3, vthole4, ialpha, jalpha , F1thole, F2thole
  real(kind=dp) :: expthole, Athole , arthole, arthole2 , arthole3 , twopiroot, erfra , ra
  real(kind=dp) :: ttt1, ttt2 
  real(kind=dp) , dimension (:) , allocatable :: tmp 
  real(kind=dp) :: F1_dm3 , F1d_dm3 , F1d2_dm3 , F2_dm5 , F2d_dm5 , F2d2_dm5 , F3_dm7 , F4_dm9 , F5_dm11
  integer       :: cthole
  logical       :: ipol, jpol
  logical       :: ldamp 
  logical       :: charge_charge, charge_dipole, dipole_dipole, charge_quadrupole, dipole_quadrupole, quadrupole_quadrupole , dip_i  , dip_j
  
  TYPE ( tensor_rank0 ) :: T0
  TYPE ( tensor_rank1 ) :: T1
  TYPE ( tensor_rank2 ) :: T2
  TYPE ( tensor_rank3 ) :: T3
  TYPE ( tensor_rank4 ) :: T4
!  TYPE ( tensor_rank5 ) :: T5

  charge_charge         = task(1)
  charge_dipole         = task(2)
  dipole_dipole         = task(3)
  charge_quadrupole     = task(4)
  dipole_quadrupole     = task(5)
  quadrupole_quadrupole = task(6)

  !  few constants
  cutsq  = verlet_coul%cut * verlet_coul%cut !cutlongrange
  alpha2 = alphaES * alphaES
  alpha3 = alpha2  * alphaES
  alpha5 = alpha3  * alpha2
  alpha7 = alpha5  * alpha2
  alpha9 = alpha7  * alpha2

  onesixth = 1.0_dp / 6.0_dp
  twopiroot = 2.0_dp / piroot
  twothird  = onesixth * 4.0_dp  
 
  cthole = 0

#ifdef debug_ES_dir
    if ( ionode ) then
        write(stdout,'(a,e16.8)') 'debug multipole_ES_dir : cutsq',cutsq
        write(stdout,'(a,6l)')    'debug multipole_ES_dir : task',task
   endif
#endif

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  u_damp = 0.0_dp

  ion1 : do ia = atom_dec%istart , atom_dec%iend
    ! =====================================
    !  verlet list : ja index in point arrays
    ! =====================================
    if ( lvnlist ) then
      jb = verlet_coul%point( ia )
      je = verlet_coul%point( ia + 1 ) - 1
    else
      ! ====================================
      !         else all ja   
      ! ====================================
      jb = 1
      je = natm
    endif
    ita  = itype(ia)
    rxi = rx(ia)
    ryi = ry(ia)
    rzi = rz(ia)
    qi  = qia(ia)
    mui = mu ( : , ia )
    thetai = theta ( : , : , ia )
    ipol = ipolar(ia)
    ialpha = poldipia(1,1,ia)

    dip_i = any ( mui .ne. 0.0d0 ) 

    ion2 : do j1 = jb, je

      if ( lvnlist ) then
        ja = verlet_coul%list ( j1 )
      else
        ja = j1
      endif

      if ( ( lvnlist .and. ja .eq. ia ) .or. ( .not. lvnlist .and. ja .le. ia ) ) cycle

        fij = 0.0_dp
        jta  = itype(ja)
        qj   = qia(ja)
        muj = mu ( : , ja )
        thetaj = theta ( : , : , ja )
        dip_j = any ( muj .ne. 0.0d0 ) 
        qij  = qi * qj
        rxj  = rx(ja)
        ryj  = ry(ja)
        rzj  = rz(ja)
        rij(1) = rxi - rxj
        rij(2) = ryi - ryj
        rij(3) = rzi - rzj
        sij = rij - nint ( rij )
        rij=0.0_dp
        jpol = ipolar(ja)
        jalpha = poldipia(1,1,ja)
        
        do j=1, 3
          do k=1, 3
            rij(j) = rij(j) + sij(k) * simu_cell%A(j,k) 
          enddo
        enddo
        d2  = rij(1) * rij(1) + rij(2) * rij(2) + rij(3) * rij(3)
        if ( d2 .gt. cutsq ) cycle

          d    = SQRT ( d2 )
          d3   = d2 * d
          d5   = d3 * d2
          d7   = d5 * d2
          d9   = d7 * d2
          dm1  = 1.0_dp / d
          dm3  = dm1 / d2
          dm5  = dm3 / d2
          dm7  = dm5 / d2
          dm9  = dm7 / d2
          dm11 = dm9 / d2

          ! damping function 
          ldamp = .false.
          if ( ldip_damping(ita,ita,jta) .or. ldip_damping(jta,ita,jta) ) ldamp = .true.
          if ( .not. damp_ind ) ldamp = .false. 
          if ( ldamp ) then
            CALL TT_damping_functions(pol_damp_b(ita,ita,jta),pol_damp_c(ita,ita,jta),d,fdamp,fdampdiff,pol_damp_k(ita,ita,jta) )
            CALL TT_damping_functions(pol_damp_b(jta,ita,jta),pol_damp_c(jta,ita,jta),d,fdamp2,fdampdiff2,pol_damp_k(jta,ita,jta) )
          else
            fdamp = 1.0_dp
            fdamp2 = 1.0_dp
            fdampdiff = 0.0d0
            fdampdiff2 = 0.0d0
          endif

          expon = EXP ( - alpha2 * d2 )    / piroot
          F0    = errfc( alphaES * d )
          F1    = F0 +  2.0_dp * alphaES * d  * expon
          F2    = F1 +  4.0_dp * alpha3  * d3 * expon / 3.0_dp
          F3    = F2 +  8.0_dp * alpha5  * d5 * expon / 15.0_dp
          F4    = F3 + 16.0_dp * alpha7  * d7 * expon / 105.0_dp
          F5    = F4 + 32.0_dp * alpha9  * d9 * expon / 945.0_dp

          ! damping if no damping fdamp == 1 and fdampdiff == 0
          F1d   = - fdamp + 1.0d0 
          F2d   = F1d  + ( d / 3.0_dp ) * fdampdiff  ! recursive relation (10) in J. Chem. Phys. 133, 234101 (2010) 
          F1d2  = - fdamp2 + 1.0d0 
          F2d2  = F1d2 + ( d / 3.0_dp ) * fdampdiff2 ! recursive relation (10) in J. Chem. Phys. 133, 234101 (2010) 

          ! =========================================
          !   multipole interaction tensor rank = 0 
          ! =========================================
          T0%sca = dm1 * F0

          ! =========================================
          !   multipole interaction tensor rank = 1
          ! =========================================
          T1%a(:)  = - rij(:) * dm3
          if ( ldamp ) then
            T1%a_damp  = T1%a * F1d
            T1%a_damp2 = T1%a * F1d2
          endif
          T1%a = T1%a * F1
          
          ! =========================================
          !   multipole interaction tensor rank = 2
          ! =========================================
!          T2%ab = 0.0_dp
!          do j = 1 , 3
!            do k = 1 , 3
!              if ( j .gt. k ) cycle
!                T2%ab (j,k) = ( 3.0_dp * rij(j) * rij(k) * F2 ) * dm5
!              if ( j .eq. k ) T2%ab (j,j) = T2%ab (j,j) - F1 * dm3
!            enddo
!          enddo
!          T2%ab (2,1) = T2%ab (1,2)
!          T2%ab (3,1) = T2%ab (1,3)
!          T2%ab (3,2) = T2%ab (2,3)

          ! =========================================
          !   multipole interaction tensor rank = 2 
          !   nb of components = 6 => reduced = 5
          !   + damping 
          ! =========================================
          T2%ab = 0.0_dp
          T2%ab_damp = 0.0_dp
          T2%ab_damp2 = 0.0_dp
          F2_dm5   = 3.0_dp * F2   * dm5
          F2d_dm5  = 3.0_dp * F2d  * dm5
          F2d2_dm5 = 3.0_dp * F2d2 * dm5
          F1_dm3   = F1   * dm3
          F1d_dm3  = F1d  * dm3
          F1d2_dm3 = F1d2 * dm3
          do j = 1 , 3
            do k = 1 , 3
                                T2%ab (j,k) = rij(j) * rij(k) * F2_dm5
                if ( j .eq. k ) T2%ab (j,j) = T2%ab (j,j) - F1_dm3
                if ( ldamp ) then
                                  T2%ab_damp  (j,k) = rij(j) * rij(k) * F2d_dm5
                                  T2%ab_damp2 (j,k) = rij(j) * rij(k) * F2d2_dm5
                  if ( j .eq. k ) then
                    T2%ab_damp (j,j)  = T2%ab_damp (j,j)  - F1d_dm3
                    T2%ab_damp2 (j,j) = T2%ab_damp2 (j,j) - F1d2_dm3
                  endif
                endif
             enddo
          enddo


          ! =========================================
          !  "Thole" functions (linear)
          !  B.T. Thole, Chem. Phys. 59 p341 (1981)
          ! ========================================
          if ( thole_functions ) then

            ! if both ions have polarizability
            if ( ipol .and. jpol ) then
               ! =========================================
               !  "Thole" functions (linear)
               !  B.T. Thole, Chem. Phys. 59 p341 (1981)
               ! ========================================
              if ( thole_function_type .eq. 'linear' ) then
                ! linear = no ewald
                sthole = thole_param (ita,jta) * ( ialpha * jalpha ) ** onesixth
                ! check distance acoording to s = 1.662 (alpha_ialpha_j)^(1/6)
                ! rij < s
                if ( d .lt. sthole ) then
                  pair_thole(ia) = ja
                  pair_thole_distance (ia) = d
#ifdef debug_thole
                  io_print write(stdout,'(3a,2i5)') 'thole function for',atype(ia),atype(ja),ia,ja
                  io_print write(stdout,'(a,e12.6,a,e12.6)') 's =', sthole,'  r =',d
#endif
                  cthole = cthole + 1 
                  vthole = d / sthole
                  vthole3 = vthole * vthole * vthole
                  vthole4 = vthole3 * vthole
                  F2thole = vthole4
                  F1thole = ( 4.0_dp * vthole3 - 3.0_dp * vthole4 ) 
                  T2%ab_thole = 0.0_dp
                  do j = 1 , 3
                    do k = 1 , 3
                      if ( j .gt. k ) cycle
                        T2%ab_thole (j,k) = ( 3.0_dp * rij(j) * rij(k) * F2thole ) * dm5       ! v^4
                      if ( j .eq. k ) T2%ab_thole (j,j) = T2%ab_thole (j,j) - F1thole * dm3    ! (4v^3-3v^4)
                    enddo
                  enddo
                  T2%ab_thole (2,1) = T2%ab_thole (1,2)
                  T2%ab_thole (3,1) = T2%ab_thole (1,3)
                  T2%ab_thole (3,2) = T2%ab_thole (2,3)
                ! if rij >= s
                else
                  T2%ab_thole = T2%ab
                endif
               ! =========================================
               !  "Thole" functions (exp)
               !  B.T. Thole, Chem. Phys. 59 p341 (1981)
               ! ========================================
              else if ( thole_function_type .eq. 'expon1' ) then
                ! linear = ewald ?
                Athole =  ( ialpha * jalpha ) ** onesixth
                uthole = thole_param(ita,jta) * ( d / Athole ) ** 3.0_dp 
                expthole = EXP ( - uthole )  
                F2thole = F2 * ( 1.0_dp - ( 1.0_dp + uthole ) * expthole ) 
                F1thole = F1 * ( 1.0_dp - expthole ) 
!                F2thole = ( 1.0_dp - ( 1.0_dp + uthole ) * expthole ) 
!                F1thole = ( 1.0_dp - expthole ) 
                T2%ab_thole = 0.0_dp
                do j = 1 , 3
                  do k = 1 , 3
                    if ( j .gt. k ) cycle
                    T2%ab_thole (j,k) = ( 3.0_dp * rij(j) * rij(k) * F2thole ) * dm5     ! ( 1 - ( 1 + a * ( r / A )^3 )  * exp ( -a ( r / A) ^3 )    
                    if ( j .eq. k ) T2%ab_thole (j,j) = T2%ab_thole (j,j) - F1thole * dm3  !  1 - exp ( -a ( r / A) ^3 )
                  enddo
                enddo
                T2%ab_thole (2,1) = T2%ab_thole (1,2)
                T2%ab_thole (3,1) = T2%ab_thole (1,3)
                T2%ab_thole (3,2) = T2%ab_thole (2,3)
               ! =========================================
               !  "Thole" functions (exp)
               !  B.T. Thole, Chem. Phys. 59 p341 (1981)
               ! ========================================
              else if ( thole_function_type .eq. 'expon2' ) then
                ! expon2 = ewald 
                Athole =  ( ialpha * jalpha ) ** onesixth
                arthole = thole_param(ita,jta) * d * Athole
                arthole2 = arthole * arthole / 2.0_dp
                arthole3 = arthole2 * arthole / 3.0_dp
                expthole = exp ( -arthole ) 
                F2thole = F2 * ( 1.0_dp - ( 1.0_dp + arthole + arthole2 + arthole3 ) * expthole ) 
                F1thole = F1 * (  1.0_dp - ( 1.0_dp + arthole + arthole2 ) * expthole ) 
                !F2thole = ( 1.0_dp - ( 1.0_dp + arthole + arthole2 + arthole3 ) * expthole ) 
                !F1thole = (  1.0_dp - ( 1.0_dp + arthole + arthole2 ) * expthole ) 
                T2%ab_thole = 0.0_dp
                do j = 1 , 3
                  do k = 1 , 3
                    if ( j .gt. k ) cycle
                    T2%ab_thole (j,k) = ( 3.0_dp * rij(j) * rij(k) * F2thole ) * dm5     ! ( 1 - ( 1 + a * ( r / A )^3 )  * exp ( -a ( r / A) ^3 )    
                    if ( j .eq. k ) T2%ab_thole (j,j) = T2%ab_thole (j,j) - F1thole * dm3  !  1 - exp ( -a ( r / A) ^3 )
                  enddo
                enddo
                T2%ab_thole (2,1) = T2%ab_thole (1,2)
                T2%ab_thole (3,1) = T2%ab_thole (1,3)
                T2%ab_thole (3,2) = T2%ab_thole (2,3)
               ! =========================================
               !  "Thole" functions (gauss)
               !  B.T. Thole, Chem. Phys. 59 p341 (1981)
               ! ========================================
              else if ( thole_function_type .eq. 'gauss' ) then
                ! gauss = no ewald
                ra = d / thole_param(ita,jta)
                Athole = twopiroot * ra 
                expthole = Athole * EXP ( -ra*ra ) 
                erfra = errf ( ra ) 
                F2thole = F2 * ( erfra - expthole * ( 1.0_dp + twothird * ra ) ) 
                F1thole = F1 * ( erfra - expthole )
                T2%ab_thole = 0.0_dp
                do j = 1 , 3
                  do k = 1 , 3
                    if ( j .gt. k ) cycle
                    T2%ab_thole (j,k) = ( 3.0_dp * rij(j) * rij(k) * F2thole ) * dm5     ! ( 1 - ( 1 + a * ( r / A )^3 )  * exp ( -a ( r / A) ^3 )    
                    if ( j .eq. k ) T2%ab_thole (j,j) = T2%ab_thole (j,j) - F1thole * dm3  !  1 - exp ( -a ( r / A) ^3 )
                  enddo
                enddo
                T2%ab_thole (2,1) = T2%ab_thole (1,2)
                T2%ab_thole (3,1) = T2%ab_thole (1,3)
                T2%ab_thole (3,2) = T2%ab_thole (2,3)
              endif
            ! not both ions have polarizability
            else
              T2%ab_thole = T2%ab
            endif
          ! no thole functions
          else
            T2%ab_thole = T2%ab 
          endif

          ! =========================================
          !   multipole interaction tensor rank = 3  
          !   nb of components = 27 => reduced = 10
          ! =========================================
          F3_dm7 = F3 * dm7 * 15.0_dp
          T3%abc = 0.0_dp
          do i = 1 , 3
            do j = 1 , 3
              do k = 1 , 3
                T3%abc (i,j,k)               = - rij(i) * rij(j) * rij(k) * F3_dm7
                if ( i.eq.j ) T3%abc (i,j,k) = T3%abc (i,j,k) + rij(k) * F2_dm5
                if ( i.eq.k ) T3%abc (i,j,k) = T3%abc (i,j,k) + rij(j) * F2_dm5
                if ( j.eq.k ) T3%abc (i,j,k) = T3%abc (i,j,k) + rij(i) * F2_dm5
              enddo
            enddo
          enddo

          ! =========================================
          !   multipole interaction tensor rank = 4  
          !   nb of components = 81 => reduced = 15
          ! =========================================
          T4%abcd = 0.0_dp
          F4_dm9 = dm9 * F4 * 105.0_dp
          do i = 1 , 3
            do j = 1 , 3
              do k = 1 , 3
                do l = 1 , 3
                  T4%abcd  (i,j,k,l) = rij(i) * rij(j) * rij(k) * rij(l) * F4_dm9
                  if ( k.eq.l ) T4%abcd (i,j,k,l) = T4%abcd  (i,j,k,l) - rij(i)*rij(j) * F3_dm7
                  if ( j.eq.l ) T4%abcd (i,j,k,l) = T4%abcd  (i,j,k,l) - rij(i)*rij(k) * F3_dm7
                  if ( j.eq.k ) T4%abcd (i,j,k,l) = T4%abcd  (i,j,k,l) - rij(i)*rij(l) * F3_dm7
                  if ( i.eq.l ) T4%abcd (i,j,k,l) = T4%abcd  (i,j,k,l) - rij(j)*rij(k) * F3_dm7
                  if ( i.eq.k ) T4%abcd (i,j,k,l) = T4%abcd  (i,j,k,l) - rij(j)*rij(l) * F3_dm7
                  if ( i.eq.j ) T4%abcd (i,j,k,l) = T4%abcd  (i,j,k,l) - rij(k)*rij(l) * F3_dm7
                  if ( i .eq. j .and. k .eq. l ) T4%abcd  (i,j,k,l) = T4%abcd  (i,j,k,l) + F2_dm5
                  if ( i .eq. k .and. j .eq. l ) T4%abcd  (i,j,k,l) = T4%abcd  (i,j,k,l) + F2_dm5
                  if ( i .eq. l .and. j .eq. k ) T4%abcd  (i,j,k,l) = T4%abcd  (i,j,k,l) + F2_dm5
                enddo
              enddo
            enddo
          enddo
        
          ! =========================================
          !   multipole interaction tensor rank = 5  
          !   nb of components = 243 => reduced = ?
          ! =========================================
          !T5%abcde = 0.0_dp
          F5_dm11 = dm11 * F5 * 945.0_dp
          !do i = 1 , 3
          !  do j = 1 , 3
          !    do k = 1 , 3
          !      do l = 1 , 3
          !        do m = 1 , 3
          !        T5%abcde  (i,j,k,l,m) = rij(i) * rij(j) * rij(k) * rij(l) * rij(m) * F5_dm11
          !        if ( l.eq.m ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(i)*rij(j)*rij(k) * F4_dm9
          !        if ( k.eq.m ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(i)*rij(j)*rij(l) * F4_dm9
          !        if ( k.eq.l ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(i)*rij(j)*rij(m) * F4_dm9
          !        if ( j.eq.m ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(i)*rij(k)*rij(l) * F4_dm9
          !        if ( j.eq.l ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(i)*rij(k)*rij(m) * F4_dm9
          !        if ( j.eq.k ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(i)*rij(l)*rij(m) * F4_dm9
          !        if ( i.eq.m ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(j)*rij(k)*rij(l) * F4_dm9
          !        if ( i.eq.l ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(j)*rij(k)*rij(m) * F4_dm9
          !        if ( i.eq.k ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(j)*rij(l)*rij(m) * F4_dm9
          !        if ( i.eq.j ) T5%abcde (i,j,k,l,m) = T5%abcde (i,j,k,l,m) - rij(k)*rij(l)*rij(m) * F4_dm9
          !        if ( i .eq. j .and. k .eq. l ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(m) * F3_dm7
          !        if ( i .eq. j .and. l .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(k) * F3_dm7
          !        if ( i .eq. j .and. k .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(l) * F3_dm7
          !        if ( i .eq. k .and. j .eq. l ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(m) * F3_dm7
          !        if ( i .eq. k .and. l .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(j) * F3_dm7
          !        if ( i .eq. k .and. j .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(l) * F3_dm7
          !        if ( i .eq. l .and. k .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(j) * F3_dm7
          !        if ( i .eq. l .and. j .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(k) * F3_dm7
          !        if ( j .eq. k .and. l .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(i) * F3_dm7
          !        if ( j .eq. k .and. i .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(l) * F3_dm7
          !        if ( j .eq. k .and. i .eq. l ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(m) * F3_dm7
          !        if ( j .eq. l .and. k .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(i) * F3_dm7
          !        if ( j .eq. l .and. i .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(k) * F3_dm7
          !        if ( k .eq. l .and. j .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(i) * F3_dm7
          !        if ( k .eq. l .and. i .eq. m ) T5%abcde  (i,j,k,l,m) = T5%abcde (i,j,k,l,m) + rij(j) * F3_dm7
          !        enddo
          !      enddo
          !    enddo
          !  enddo
          !enddo



        ! note :
        ! les termes faisant intervenir les tenseurs d'interactions : T0,T2,T4...
        ! sont symétriques par changement de direction de l'interaction rij => rji
        ! alors que T1,T3,T5 ... ne le sont pas.
        ! en pratique : les sommes sur i ou j change de signe pour T1,T3,T5

        ! ===========================================================
        !                  charge-charge interaction
        ! ===========================================================
        qq : if ( charge_charge ) then
          ! energy
          u_dir = u_dir + qij * T0%sca

          ! electric field
          if ( do_efield ) then
            ef_dir ( : , ia )   = ef_dir ( : , ia ) - qj * T1%a(:) 
            ef_dir ( : , ja )   = ef_dir ( : , ja ) + qi * T1%a(:)
            if ( ldamp ) then
              ef_dir ( : , ia ) = ef_dir ( : , ia ) + qj * T1%a_damp(:) 
              ef_dir ( : , ja ) = ef_dir ( : , ja ) - qi * T1%a_damp2(:)
            endif
          endif
          
          ! electric field gradient
          if ( do_efg ) then
            efg_dir ( : , : , ia ) = efg_dir ( : , : , ia ) - qj * T2%ab ( : , : ) 
            efg_dir ( : , : , ja ) = efg_dir ( : , : , ja ) - qi * T2%ab ( : , : ) 
          endif

          ! forces
          if ( do_forces ) then
            fij(:)  = fij(:) + qij * T1%a(:)
          endif

        endif qq
        ! ===========================================================
        !                  dipole-dipole interaction
        ! ===========================================================
        dd : if ( dipole_dipole .and. ( dip_i .and. dip_j ) ) then

          ! remind "thole" interaction function : T2_thole = T2 (regular) for r>s (linear) 
          ! energy
          do j = 1, 3
            do k = 1, 3 
               u_dir = u_dir - mui(j) * T2%ab_thole(j,k) * muj(k) 
            enddo
          enddo

          ! remind thole interaction function : T2_thole = T2 (regular) for r>s  (linear)
          ! electric field
          if ( do_efield ) then
            do k = 1 , 3
              ef_dir ( : , ia ) = ef_dir ( : , ia ) + T2%ab_thole(:,k) * muj(k)
              ef_dir ( : , ja ) = ef_dir ( : , ja ) + T2%ab_thole(:,k) * mui(k) 
            enddo
          endif

          ! electric field gradient 
          if ( do_efg ) then
            do k = 1 , 3
              efg_dir ( : , : , ia ) = efg_dir ( : , : , ia ) + T3%abc (:,:,k) * muj(k)
              efg_dir ( : , : , ja ) = efg_dir ( : , : , ja ) - T3%abc (:,:,k) * mui(k)
            enddo
          endif

          ! forces
          if ( do_forces ) then
            do j = 1 , 3 
              do k = 1, 3
                fij(:) = fij(:) - mui(j) * T3%abc (j,:,k) * muj(k) 
              enddo
            enddo
          endif 

        endif dd

        ! ===========================================================
        !                  charge-dipole interaction
        ! ===========================================================

        qd : if ( charge_dipole .and. ( dip_i .or. dip_j ) ) then
          
          ! electrostatic energy
          do k = 1 , 3 
            u_dir = u_dir - qi *  T1%a(k) * muj(k)  
            u_dir = u_dir + qj *  T1%a(k) * mui(k) 
            if ( ldamp ) then
              u_dir = u_dir + qi *  T1%a_damp2 (k) * muj(k)  
              u_dir = u_dir - qj *  T1%a_damp  (k) * mui(k) 
            endif
          enddo


          ! forces
          ! remarque 1 : pour garder la construction : 
          !           f(ia) = f(ia) - fij 
          !           f(ja) = f(ja) + fij 
          ! en fin de boucle uniquement, 
          ! nous avons ici un changement de signe sur la somme sur j malgré
          ! le fait que le tenseur d'interaction soit paire (e.g T2 )
          ! remarque 2 : thole functions only apply to fields not forces 
          if ( do_forces ) then
            do j = 1 , 3
              do k = 1 , 3
                fij(k) = fij(k) - qi * T2%ab (k,j) * muj(j)
                fij(k) = fij(k) + qj * T2%ab (k,j) * mui(j)
              enddo
            enddo

          endif

          ! damping 
          if ( ldamp ) then
            ! forces
            if ( do_forces ) then
              !fij=0.0_dp
              do j = 1 , 3
                do k = 1 , 3
                  fij(k) = fij(k) - qj * T2%ab_damp  (j,k) * mui(j)
                  fij(k) = fij(k) + qi * T2%ab_damp2 (j,k) * muj(j)
                enddo
              enddo
            endif
          endif

        endif qd

        ! ===========================================================
        !                  charge-quadrupole interaction
        ! ===========================================================
        qquad : if ( charge_quadrupole ) then
          ! electrostatic energy
          do k = 1 , 3
            do j = 1 , 3
              u_dir = u_dir + qi *  T2%ab(k,j) * thetaj(k,j) / 3.0_dp
              u_dir = u_dir + qj *  T2%ab(k,j) * thetai(k,j) / 3.0_dp
            enddo
          enddo
        endif qquad
        
        ! ===========================================================
        !                  dipole-quadrupole interaction
        ! ===========================================================
        dquad : if ( dipole_quadrupole ) then
          ! electrostatic energy
          do l = 1 , 3
            do k = 1 , 3
              do j = 1 , 3
                u_dir = u_dir + mui(l) * T3%abc(l,k,j) * thetaj(k,j) / 3.0_dp
                u_dir = u_dir - muj(l) * T3%abc(l,k,j) * thetai(k,j) / 3.0_dp
              enddo
            enddo
          enddo
        endif dquad

        ! ===========================================================
        !                  quadrupole-quadrupole interaction
        ! ===========================================================
        quadquad : if ( quadrupole_quadrupole ) then

          ! energy
          do i = 1, 3
            do j = 1, 3
              do k = 1, 3
                do l = 1, 3
                  u_dir = u_dir + thetai(i,j) * T4%abcd(i,j,k,l) * thetaj(k,l) / 9.0_dp
                enddo
              enddo
            enddo
          enddo

          ! electric field
          if ( do_efield ) then
            do k = 1 , 3
              do j = 1 , 3
                ef_dir ( ia , :  ) = ef_dir ( ia , : ) - T3%abc(:,k,j) * thetaj(k,j)
                ef_dir ( ja , :  ) = ef_dir ( ja , : ) + T3%abc(:,k,j) * thetai(k,j)
              enddo
            enddo
          endif

          ! electric field gradient 
          if ( do_efg ) then
            do k = 1 , 3
              do j = 1 , 3
                efg_dir ( ia , : , : ) = efg_dir ( ia , : , :  ) + T4%abcd (:,:,k,j) * thetaj(k,j)
                efg_dir ( ja , : , : ) = efg_dir ( ja , : , :  ) - T4%abcd (:,:,k,j) * thetai(k,j)
              enddo
            enddo
          endif

        endif quadquad


        ! forces
        if ( do_forces ) then
          fx_dir ( ia ) = fx_dir ( ia ) - fij(1)
          fy_dir ( ia ) = fy_dir ( ia ) - fij(2)
          fz_dir ( ia ) = fz_dir ( ia ) - fij(3)
          fx_dir ( ja ) = fx_dir ( ja ) + fij(1)
          fy_dir ( ja ) = fy_dir ( ja ) + fij(2)
          fz_dir ( ja ) = fz_dir ( ja ) + fij(3)
        endif

        ! stress
        if ( do_stress ) then
          do j = 1, 3
            do k = 1, 3
              tau_dir(j,k) = tau_dir(j,k) - ( rij(j) * fij(k) + rij(k) * fij(j) )
            enddo
          enddo
        endif

    enddo ion2

  enddo ion1

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  allocate ( tmp ( natm ) ) 

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u_dir )
  if ( do_efield ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( 1, : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( 2, : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( ef_dir ( 3, : ) , natm )
  endif
  if ( do_efg ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( 1 , 1 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( 2 , 2 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( 3 , 3 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( 1 , 2 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( 1 , 3 , : ) , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( efg_dir ( 2 , 3 , : ) , natm )
  endif
  if ( do_forces ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( fx_dir , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( fy_dir , natm )
    CALL MPI_ALL_REDUCE_DOUBLE ( fz_dir , natm )
  endif
  if ( do_stress ) then
    CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 1 , : ) , 3 )
    CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 2 , : ) , 3 )
    CALL MPI_ALL_REDUCE_DOUBLE ( tau_dir ( 3 , : ) , 3 )
    tau_dir =   tau_dir / simu_cell%omega * 0.5_dp
  endif

  ! thole function overview 
  ! write it only when do_forces = .true. only 
  if ( do_forces ) then
    CALL MPI_ALL_REDUCE_INTEGER_SCALAR ( cthole )
    CALL MPI_ALL_REDUCE_INTEGER ( pair_thole , natm )
    CALL MPI_ALL_REDUCE_DOUBLE  ( pair_thole_distance , natm )
    if ( thole_functions .and. cthole .ne. 0 ) then
      io_printnode write( stdout , '(a,i5,a)') 'Thole damping used for ',cthole,' pairs'
      if ( ioprintnode ) then
        do ia = 1 ,natm
          ita= itype(ia)
          if ( pair_thole(ia) .ne. 0 ) then 
            jta= itype(pair_thole(ia))
            ialpha = poldipia(1,1,ia)
            jalpha = poldipia(1,1,pair_thole(ia))
            Athole = ( ialpha * jalpha ) ** onesixth 
            write(stdout, '(a,2i5,a5,a,a2,a,f8.5,a,f8.5,a,f8.5,a,f8.5,a)')  'pair : ',ia,pair_thole(ia), atype(ia) ,' -- ', atype(pair_thole(ia)),&
                  'distance = ',pair_thole_distance(ia),' catastrophe at distance = ',4.0_dp**onesixth*Athole,' thole param = ',thole_param(ita,jta) ,'(',thole_param(ita,jta)*Athole,')'
          endif
        enddo
      endif
    io_printnode WRITE( stdout , '(a)' ) '' 
    endif
  endif

  deallocate ( tmp ) 

  return

END SUBROUTINE multipole_ES_dir 

END MODULE ewald_dir 
! ===== fmV =====

