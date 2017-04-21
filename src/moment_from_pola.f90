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
MODULE moment_from_pola 

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  USE mpimdff

  implicit none

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

! *********************** SUBROUTINE get_dipole_moments *************************
!
!> \brief
!!  get dipole moments
!> \author
!! FMV
!
!> \date 
!! February 2014
!
! ******************************************************************************
SUBROUTINE get_dipole_moments ( mu , theta , didpim )

  USE coulomb,          ONLY :  mu_t, thole_functions, pair_thole , pair_thole_distance
  USE config,           ONLY :  natm , ntype , dipia , quadia , dipia_wfc, fx,fy,fz, atype 
  USE pim,              ONLY :  lwfc, ldip_wfc,algo_moment_from_pola
  USE field,            ONLY :  write_DIPFF
  USE io,               ONLY :  kunit_DIPFF

  implicit none

  ! global
  real(kind=dp) , intent ( out ) :: mu (:,:)
  real(kind=dp) , intent ( out ) :: theta (:,:,:)

  ! local
  integer :: it, ia ,ja 
  logical :: lwannier, lcata , didpim 
  real(kind=dp), dimension ( : )  , allocatable :: fx_save , fy_save, fz_save
  real(kind=dp), dimension ( : , :  ) , allocatable :: dipia_ind
  real(kind=dp), dimension ( : , : , :  ) , allocatable :: quadia_ind
  integer :: ierr 

  ! save total force fx , fy, fz as they are overwritted by moment_from_pola
  allocate ( fx_save(natm)  , fy_save (natm) , fz_save(natm)  )
  allocate ( dipia_ind ( 3, natm)  )
  allocate ( quadia_ind ( 3, 3, natm)  )
  fx_save = fx ; fy_save = fy ; fz_save = fz
  dipia_ind = 0.0d0

  ! ======================================
  !     induced moment from polarisation 
  ! ======================================
  if ( algo_moment_from_pola .eq. 'scf' ) then
    CALL moment_from_pola_scf    ( dipia_ind , quadia_ind , didpim )
  else if ( algo_moment_from_pola .eq. 'scf_kO_v1' ) then
    CALL moment_from_pola_scf_kO_v1 ( dipia_ind , quadia_ind , didpim )
  else if ( algo_moment_from_pola .eq. 'scf_kO_v2' ) then
    CALL moment_from_pola_scf_kO_v2 ( dipia_ind , quadia_ind , didpim )
  else if ( algo_moment_from_pola .eq. 'scf_kO_v3' ) then
    CALL moment_from_pola_scf_kO_v3 ( dipia_ind , quadia_ind , didpim )
  else if ( algo_moment_from_pola .eq. 'scf_kO_v4_1' .or. &
            algo_moment_from_pola .eq. 'scf_kO_v4_2'      ) then
    CALL moment_from_pola_scf_kO_v4 ( dipia_ind , quadia_ind , didpim )
  endif

  ! =========================================================
  !  Is there any Wannier center ? if yes lwannier = .TRUE.
  ! =========================================================
  lwannier = .false.
  do it = 1 , ntype
    if ( lwfc ( it ) .gt. 0 ) lwannier = .true.
  enddo

  if ( lwannier ) then
    ! ======================================
    !     induced moment from wannier center  
    ! ======================================
    CALL moment_from_WFc ( dipia_wfc )

    if ( ldip_wfc ) then
      if ( ionode ) write( stdout , '(A)' ) 'dipole contribution from wannier centers '
        ! ======================================
        !  total dipole mu :
        !   static +  induced  + "wannier"
        ! ======================================
        mu = dipia + dipia_ind + dipia_wfc
    else
        if ( ionode ) write( stdout , '(A)' ) 'no dipole contribution from wannier centers '
          mu = dipia + dipia_ind
    endif
  else
    mu    = dipia  + dipia_ind
    !theta = quadia + quadia_ind
    theta = 0.0_dp
  endif

  fx = fx_save
  fy = fy_save
  fz = fz_save

  deallocate ( fx_save, fy_save, fz_save ) 
  deallocate ( dipia_ind  )
  deallocate ( quadia_ind  )

  ! =====================================
  ! check for POLARIZATION CATASTROPH
  ! =====================================
  lcata = .false.
  cataloop : do ia = 1 , natm
    if ( any( mu(:, ia ).gt. 10_dp ) ) then
      lcata = .true. 
      if ( ionode ) then
        WRITE ( stdout , '(a)' )           '----------------------------------------------------------------------'
        WRITE ( stdout , '(a)' )           ' ERROR : POLARIZATION CATASTROPH'
        WRITE ( stdout , '(a)' )           ' ERROR : one (or more) dipole moment is extremely large ( > 10 [eA] )'
        WRITE ( stdout , '(a,i5,3e16.8)' ) ' ion   = ', ia , mu( : , ia )
        WRITE ( stdout , '(a)' )           '----------------------------------------------------------------------'
      endif
      exit cataloop
    endif      
  enddo cataloop
  if ( lcata ) then 
    OPEN  ( kunit_DIPFF , file = 'DIPFF',STATUS = 'UNKNOWN')
    mu_t = mu
    CALL    write_DIPFF 
    CLOSE ( kunit_DIPFF )
  endif
  if ( lcata .and. thole_functions ) then
    io_node print*,pair_thole
    do ja = 1 , natm
      if ( pair_thole(ja) .ne. 0 ) then
        if ( ionode ) then
          WRITE ( stdout, '(a)' )          'Ajust thole_param to avoid this issue'
          WRITE ( stdout, '(a,2i5,a5,a,a5,a,f12.8)') 'thole pair : ',ja,pair_thole(ja), atype(ja) ,'-', atype(pair_thole(ja)), 'distance = ',pair_thole_distance(ja)
          WRITE ( stdout, '(a,i5,3f12.8)') '  dipoles for ', ja , mu (:, ja ) 
          WRITE ( stdout, '(a,i5,3f12.8)') '  dipoles for ', pair_thole(ja) , mu (:,pair_thole(ja))
        endif
      endif
    enddo
  else if ( lcata ) then
    WRITE ( stdout , '(a)' )           'Use thole damping function : linear (recommended)'
  endif
  if ( lcata ) then
#ifdef MPI
    CALL MPI_FINALIZE(ierr)
#endif
    STOP
  endif

#ifdef debug_mu
        write( stdout,'(a)')  
#endif

  return

END SUBROUTINE get_dipole_moments


! *********************** SUBROUTINE moment_from_pola_scf **********************
!
!> \brief
!! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
!! The evaluation is done self-consistently starting from the the field due the
!! point charges only.
!
!> \param[out] mu_ind induced electric dipole from polarizabilities
!
!> \note
!! The stopping criteria is governed by conv_tol_ind
!
! ******************************************************************************
SUBROUTINE moment_from_pola_scf ( mu_ind , theta_ind , didpim ) 

  USE constants,        ONLY :  coul_unit
  USE coulomb,          ONLY :  alphaES, qch, dip, dipia_ind_t
  USE config,           ONLY :  natm , atype , fx , fy , fz , ntype , dipia , quadia , qia , ntypemax, poldipia , invpoldipia , itype 
  USE control,          ONLY :  longrange , calc
  USE thermodynamic,    ONLY :  u_pol, u_coul
  USE ewald,            ONLY :  multipole_ES
  USE pim,              ONLY :  ldip_polar, max_scf_pol_iter, min_scf_pol_iter, conv_tol_ind, extrapolate_order,algo_ext_dipole
  USE time,             ONLY :  time_moment_from_pola
#ifdef debug_print_dipff_scf  
  USE md,               ONLY :  itime
#endif

  implicit none

  ! global
  real(kind=dp) , intent (out) :: mu_ind ( : , : ) 
  real(kind=dp) , intent (out) :: theta_ind ( : , : , : ) 

  ! local
  integer :: ia , iscf , it , alpha, beta
  logical :: linduced , didpim
  real(kind=dp) :: tttt ,alphaES_save
  real(kind=dp) :: u_coul_stat , rmsd , u_coul_pol, u_coul_ind
  real(kind=dp) :: Efield( 3 , natm ) , Efield_stat ( 3 , natm ) , Efield_ind ( 3 , natm ), EfieldG_stat ( 3 , 3 , natm ) , EfieldG_ind ( 3 , 3 , natm )
  real(kind=dp) :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 
  logical       :: task_static (6), task_ind(6), ldip
  real(kind=dp) , dimension (:,:) , allocatable :: f_ind , mu_prev
#ifdef debug_print_dipff_scf
  integer :: kkkk
#endif 
#ifdef MPI
  dectime
#endif

#ifdef debug_print_dipff_scf
  kkkk=100000*itime+100000
  print*,kkkk,itime
#endif


  tttt=0.0_dp
#ifdef MPI
  statime
#endif
  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  didpim   = .false.
  do it = 1 , ntype
    if ( ldip_polar ( it ) ) linduced = .true.
  enddo
  if ( .not. linduced ) then
#ifdef debug
    write(stdout,'(a,e16.8)') 'quick return from moment_from_pola_scf',mu_ind(1,1)
#endif
    return
  endif
  didpim = linduced

  ! =============================================
  !  calculate static Efield ( charge + dipoles )
  ! =============================================
  u_pol = 0.0_dp
  u_coul = 0.0_dp
  Efield_stat = 0.0_dp
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  ldip=.false.
  task_static = .false.      
  task_static(1) = .true.      
  if ( any (dip .ne. 0.0d0 ) ) ldip = .true.
  if ( ldip ) then
    task_static = .true.        
  endif
  ! ==============================================
  !          main ewald subroutine 
  ! ==============================================
  ! do_strucfact = .true. : recalculate exp(ikr) for new r
  ! use_ckrskr  = .true. : use the store exp(ikr) 
  CALL  multipole_ES ( Efield_stat , EfieldG_stat , dipia , quadia , task_static , & 
                       damp_ind =.true. , do_efield=.true. , do_efg=.false. , & 
                       do_forces=.false. , do_stress = .false. , do_rec = .true. , do_dir = .true. , do_strucfact= .true. , use_ckrskr=.true. ) 
  u_coul_stat = u_coul 
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  iscf = 0
  rmsd = HUGE(0.0d0)
  ! =========================
  !  charges are set to zero 
  ! =========================
  alphaES_save = alphaES
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0_dp
  qia = 0.0_dp
  task_ind    = .false.
  task_ind(1) = .false. 
  task_ind(2) = .false. 
  task_ind(3) = .true. 
  io_printnode WRITE ( stdout ,'(a)') '' 
  u_pol = 0.0_dp

  lseparator_ioprintnode(stdout)
  io_printnode WRITE( stdout , '(a)' ) '                   running scf algo                          ' 
  lseparator_ioprintnode(stdout)
  ! =============================
  !           SCF LOOP
  ! =============================
  do while ( ( iscf < max_scf_pol_iter ) .and. ( rmsd .gt. conv_tol_ind )  .or. ( iscf < min_scf_pol_iter  ) )

    iscf = iscf + 1

    ! ==========================================================
    !  calculate mu_ind from Efield = Efield_stat + Efield_ind
    !  extrapolate from previous md steps :
    !  by ASPC (Always Stable Predictor Corrector) ( algo_ext_dipole .eq. 'aspc  ) by Kolafa J. Comp. Chem. v25 n3 (2003)
    ! ==========================================================
    if ( iscf.ne.1 .or. ( calc .ne. 'md' .and.  calc .ne. 'opt' ) ) then
      CALL induced_moment       ( Efield , mu_ind )  ! Efield in ; mu_ind and u_pol out

    else
      CALL extrapolate_dipole_aspc ( mu_ind , Efield , key=1 )  ! predictor
    endif
    ! ===============================================================
    !      u_pol = 1/ ( 2 alpha  ) * | mu_ind | ^ 2
    ! ---------------------------------------------------------------
    ! note on units :
    !    [ poldipia  ] = A ^ 3
    !    [ mu_ind ] = e A 
    !    [ u_pol  ] = e^2 / A ! electrostatic internal energy  
    ! ===============================================================
    u_pol = 0.0_dp
    do ia = 1 , natm
      do alpha = 1 , 3
        do beta = 1 , 3
          u_pol = u_pol + mu_ind ( alpha , ia ) * invpoldipia ( alpha , beta , ia ) * mu_ind ( beta , ia )
        enddo
      enddo
    enddo
    u_pol = u_pol * 0.5_dp



    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true., do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    u_coul_ind = u_coul 
    u_coul_pol = u_coul_stat-u_coul_ind

    Efield = Efield_stat + Efield_ind
    fx      = 0.0_dp
    fy      = 0.0_dp
    fz      = 0.0_dp

    ! ASPC corrector 
    if ( iscf.eq.1 .and. calc .eq. 'md' ) then
      CALL extrapolate_dipole_aspc ( mu_ind , Efield  , 2 ) ! corrector
    endif
    
    CALL get_rmsd_mu ( rmsd , Efield , mu_ind )

    ! output
    if ( calc .ne. 'opt' ) then
#ifdef debug
    io_printnode WRITE ( stdout ,'(a,i4,5(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' u_coul (qq)  = ', u_coul_stat, ' u_coul (dd)  = ', u_coul_ind,' u_coul_pol = ', u_coul_pol, ' rmsd       = ', rmsd
#endif
    io_printnode WRITE ( stdout ,'(a,i4,2(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' rmsd       = ', rmsd
    endif
 
#ifdef debug_print_dipff_scf
  do ia = 1 , natm
  !if ( mu_ind ( 1 , ia ) .eq. 0._dp ) cycle
  write(kkkk,'(3e16.8)') (mu_ind ( alpha , ia ),alpha=1,3)
  enddo
  kkkk = kkkk + 1 
#endif

  enddo ! end of SCF loop

  ! ===========================
  !  charge/force info is recovered
  ! ===========================
  alphaES = alphaES_save
  qch = qch_tmp
  qia = qia_tmp

  if ( ioprintnode .and.  calc .ne. 'opt' ) then
    blankline(stdout)
    WRITE ( stdout , '(a,i6,a)')            'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3,a,e10.3,a)') 'Electric field is converged at ',rmsd,' ( ',conv_tol_ind,' ) '
    blankline(stdout)
  endif
#ifdef debug_mu
  if ( ionode ) then
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm 
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
    enddo
    blankline(stdout)
  endif
#endif
  ! store induced dipole at t 
  dipia_ind_t(1,:,:) = mu_ind
  
#ifdef debug_extrapolate
    if ( ioprintnode ) then
      print*,'previous step',itime
      do ia=1,natm
        write(stdout,'(i5,<extrapolate_order+2>e16.8)') ia, mu_ind ( 1 , ia ), (dipia_ind_t ( t , 1 , ia ) , t=1, extrapolate_order+1 )
      enddo
    endif
#endif

#ifdef MPI
  stotime
  addtime(time_moment_from_pola)
#endif
  ! temporary
  theta_ind=0.0_dp

  return

END SUBROUTINE moment_from_pola_scf



! *********************** SUBROUTINE extrapolate_dipole_aspc *******************
! order of extrapolation is k+1 of Kolafa original derivation.
! as the zero order is taking juste 
! ******************************************************************************
SUBROUTINE extrapolate_dipole_aspc ( mu_ind , Efield , key )

  USE coulomb,          ONLY :  dipia_ind_t
  USE config,           ONLY :  atype, natm, invpoldipia
  USE pim,              ONLY :  extrapolate_order
  USE md,               ONLY :  itime

  implicit none
  ! global  
  real(kind=dp) , intent (out):: mu_ind (:,:)
  real(kind=dp) , intent (in) :: Efield (:,:)
  integer :: key
  ! local
  integer :: ia, k ,ext_ord, alpha, beta
  real(kind=dp)               :: mu_p ( 3 , natm )
  real(kind=dp)               :: mu_p_save (3 , natm )
  ! ASPC coefficient
  real(kind=dp) :: B_ASPC(6), w_ASPC

  ! switch to lower order of extrapolation 
  ! if time step is not large enough to store previous dipoles steps
  if ( itime .ge. extrapolate_order+1 ) then
    ext_ord = extrapolate_order
  else
    if (itime.eq.1.or.itime.eq.0)   then
      ext_ord = 0
    else
      ext_ord = itime - 1
    endif
  endif

!  write(*,'(a,2i,3f16.8)') '(1) key : ',key,ext_ord,mu_ind(1,:)

  SELECT CASE ( ext_ord )
  CASE DEFAULT
    io_node WRITE(stdout,'(a)') 'value of aspc extrapolation order not available should be 0<= extrapolate_order <= 4'
    STOP
  CASE(0)
    B_ASPC(1) =  1.0_dp
    W_ASPC    =  0.0_dp
  CASE(1)
    B_ASPC(1) =  2.0_dp
    B_ASPC(2) = -1.0_dp
    W_ASPC    =  2.0_dp/3.0_dp
  CASE(2)
    B_ASPC(1) =  2.5_dp
    B_ASPC(2) = -2.0_dp
    B_ASPC(3) =  0.5_dp
    W_ASPC    =  0.6_dp
  CASE(3)
    B_ASPC(1) =  2.8_dp
    B_ASPC(2) = -2.8_dp
    B_ASPC(3) =  1.2_dp
    B_ASPC(4) = -0.2_dp
    W_ASPC    =  4.0_dp/7.0_dp
  CASE(4)
    B_ASPC(1) =  3.0_dp
    B_ASPC(2) = -24.0_dp/7.0_dp
    B_ASPC(3) =  27.0_dp/14.0_dp
    B_ASPC(4) = -4.0_dp/7.0_dp
    B_ASPC(5) =  1.0_dp/14.0_dp
    W_ASPC    =  5.0_dp/9.0_dp
  CASE(5)
    B_ASPC(1) =  22.0_dp/7.0_dp
    B_ASPC(2) = -55.0_dp/14.0_dp
    B_ASPC(3) =  55.0_dp/21.0_dp
    B_ASPC(4) = -22.0_dp/21.0_dp
    B_ASPC(5) =  5.0_dp/21.0_dp
    B_ASPC(6) = -1.0_dp/42.0_dp
    W_ASPC    =  6.0_dp/11.0_dp
  END SELECT

  ! predictor
  if ( key .eq. 1 ) then
    mu_p = 0.0_dp
    do k=1 , ext_ord+1
      mu_p = mu_p + B_ASPC( k  ) * dipia_ind_t ( k , : , : )
    enddo
    mu_ind=mu_p
    do k = ext_ord + 1, 2, -1
      dipia_ind_t(k,:,:) = dipia_ind_t(k-1,:,:)
    enddo
!    write(*,'(a,2i,3f16.8)') '(2) key : ',key,ext_ord,mu_ind(1,:)
  endif

  ! corrector
  if ( key.eq.2) then
    mu_p_save = mu_ind
    CALL induced_moment ( Efield , mu_ind )
    mu_ind = w_ASPC * mu_ind + ( 1.0_dp - w_ASPC ) * mu_p_save
  endif

!#ifdef debug_mu
!      if ( ionode ) then
!        WRITE ( stdout , '(a)' )     'Induced dipoles at atoms from extrapolation: '
!        do ia = 1 , natm 
!          WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
!          ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
!        enddo
!        blankline(stdout)
!      endif
!#endif

  return

END SUBROUTINE extrapolate_dipole_aspc

! *********************** SUBROUTINE moment_from_WFc ***************************
!
!> \brief
!!  This routines evaluates the dipole moment induced by Wannier centers (Wfc).
!!  the listing of Wfc is done as for the verlet list
!
!> \description  
!!           wfc_point         : array of size natm+1 
!!                               gives the starting and finishing index of array
!!                               list for a given atom i
!!                               jbegin = point(i)  jend = point(i+1) - 1
!!           wfc_list          : index list of neighboring atoms
!!
!!  how to use it :
!!                          do ia = 1, natm
!!                            jbegin = wfc_point(i)
!!                            jend = wfc_point(i+1) - 1
!!                            do jvnl = jbegin , jend
!!                              ja = wfc_list ( jvnl ) 
!!                              then ia en ja are neighboors   
!!                            enddo
!!                          enddo
!!
!> \param[out] mu electric dipole 
!
!> \author
!! FMV
!
!> \date 
!! January 2013
!
! ******************************************************************************
SUBROUTINE moment_from_WFc ( mu )

  USE constants,                ONLY :  Debye_unit
  USE config,                   ONLY :  verlet_list , natm , ntype , itype , atype , simu_cell , rx , ry , rz 
  USE cell,                     ONLY :  kardir , dirkar 
  USE pim,                      ONLY :  lwfc, rcut_wfc, lwrite_dip_wfc
  USE io,                       ONLY :  kunit_DIPWFC

  implicit none

  ! global
  real(kind=dp), intent ( out )  :: mu ( 3 , natm ) 

  ! local
  integer :: it, jt , ia , ja , icount , k , jb , je , jwfc , ia_last
  integer, dimension ( : ) , allocatable :: cwfc
  TYPE( verlet_list ) :: verlet_wfc
  real(kind=dp) :: rijsq, cutsq
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: dmu
  logical          :: lwannier 

  ! ==========================================================================
  !  quick return 
  !  Is there any Wannier center ? if yes lwannier = .TRUE.
  !  redondant : cette condition est egalement evalué avant ( fast anyway )
  ! ==========================================================================
  lwannier = .false.
  do it = 1 , ntype
    if ( lwfc ( it ) .gt. 0 ) lwannier = .true.
  enddo
  if ( .not. lwannier ) then
    return
  endif

#ifdef debug_wfc
  WRITE ( stdout , '(a)' ) 'debug : in moment_from_WFs'
#endif

  cutsq = rcut_wfc * rcut_wfc
  mu = 0.0_dp

  allocate( verlet_wfc%list ( natm * 250 ) , verlet_wfc%point(  natm + 1 ) )
  allocate ( cwfc ( natm ) )
  cwfc = 0
  verlet_wfc%list  = 0
  verlet_wfc%point = 0

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  icount = 1
  do ia = 1 , natm 
    it = itype ( ia )
    ! loop only on non wannier-centres
    if ( lwfc ( it ) .lt. 0 )  cycle
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      k = 0    
      do ja = 1 , natm
        jt = itype (ja) 
        ! check distance to wannier-centres only
        if ( lwfc ( jt ) .ge. 0 ) cycle 
          rxij  = rxi - rx ( ja )
          ryij  = ryi - ry ( ja )
          rzij  = rzi - rz ( ja )
          sxij  = rxij - nint ( rxij )
          syij  = ryij - nint ( ryij )
          szij  = rzij - nint ( rzij )
          rxij  = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
          ryij  = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
          rzij  = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
          rijsq = rxij * rxij + ryij * ryij + rzij * rzij
          if ( rijsq .lt. rcut_wfc ) then
            cwfc ( ia ) = cwfc ( ia ) + 1
            icount = icount+1
            k = k + 1
            verlet_wfc%list( icount - 1 ) = ja
            ia_last = ia
#ifdef debug_wfc
          WRITE ( stdout ,'(a4,i4,a4,i4,2e17.8)') atype(ia),ia,atype(ja),ja,SQRT(rijsq),rcut_wfc
          WRITE ( stdout ,'(a,5e17.8)')          'distance ',rxij,ryij,rzij,SQRT(rijsq),rcut_wfc
#endif
          endif

      enddo
      verlet_wfc%point( ia ) = icount-k
  enddo
  verlet_wfc%point ( ia_last + 1 ) = icount


  ! ===========================================
  !  check if some wannier centers is missing
  ! ===========================================
  do ia = 1 , natm
    it = itype ( ia ) 
    if ( ( lwfc(it) .gt. 0 ) .and. ( cwfc ( ia ) .ne. lwfc ( it ) ) ) then
      WRITE ( * ,* ) 'ERROR in moment_from_WFc : wannier center is missing for atom ',ia , cwfc ( ia ) , lwfc ( it ) 
      STOP 
    endif
  enddo

  ! =============================================
  !  compute dipole moments from Wannier centers
  ! =============================================
  do ia = 1 , natm
    it = itype ( ia ) 
    jb = verlet_wfc%point ( ia )
    je = verlet_wfc%point ( ia + 1 ) - 1
    if ( lwfc ( it ) .lt. 0 ) cycle
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      jb = verlet_wfc%point ( ia )
      je = verlet_wfc%point ( ia + 1 ) - 1
      do jwfc = jb , je
        ja = verlet_wfc%list( jwfc )
        rxij  = rx ( ja ) - rxi
        ryij  = ry ( ja ) - ryi
        rzij  = rz ( ja ) - rzi
        sxij  = rxij - nint ( rxij )
        syij  = ryij - nint ( ryij )
        szij  = rzij - nint ( rzij )
        rxij  = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij  = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij  = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        mu ( 1 , ia ) = mu ( 1 , ia ) + rxij 
        mu ( 2 , ia ) = mu ( 2 , ia ) + ryij
        mu ( 3 , ia ) = mu ( 3 , ia ) + rzij
      enddo
  enddo

  ! charges of WFC
  mu = -2.0_dp * mu 

  if ( ionode .and. lwrite_dip_wfc ) then 
    OPEN ( UNIT = kunit_DIPWFC , FILE='DIPWFC' )     
    do ia= 1 , natm
      it = itype ( ia ) 
      if ( lwfc ( it ) .lt. 0 ) cycle
      WRITE ( kunit_DIPWFC , '(a,3e16.8)' ) atype( ia ) , mu ( 1 , ia ) , mu ( 2 , ia ) , mu ( 3 , ia )  
    enddo
    CLOSE ( kunit_DIPWFC ) 
  endif

  ! see constants.f90
  ! Debye_unit = 2.54174618782479816355_dp / bohr

  if ( ionode ) then
    WRITE ( stdout , '(a)' ) 'Dipoles from WFCs'
    WRITE ( stdout , '(a)' ) 'atom         |mu| [D]      |mu| [eA]'
    do ia = 1 , natm 
      it = itype ( ia )  
      if ( lwfc ( it ) .lt. 0 ) cycle
      dmu = mu ( 1 , ia ) * mu ( 1 , ia ) + mu ( 2 , ia ) * mu ( 2 , ia ) + mu ( 3 , ia ) *  mu ( 3 , ia )
      dmu = SQRT ( dmu ) 
      WRITE ( stdout , '(a,4x,2e16.8)' ) atype(ia), dmu * Debye_unit, dmu 
    enddo 
  endif

  deallocate ( verlet_wfc%list , verlet_wfc%point )
  deallocate ( cwfc ) 

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE moment_from_WFc


! *********************** SUBROUTINE induced_moment_inner ****************************
!> \brief
!! this subroutine calculates the induced moment from the total Electric field and the
!! polarizability tensor
!! Basically used in the SCF loop to get induced_moment
!! \f$ \mu_{i,\alpha} =  p_{i,\alpha,\beta} * E_{i,\beta} \f$
!> \param[in]  f_ind_ext 
!> \param[out] mu_ind induced electric dipole define at ion position 
!> \note
!! poldipia is the polarizability tensor
!! algorithm by kO ( Kirill Okhotnikov ) afternoon of 28/10/14 
! ******************************************************************************
SUBROUTINE induced_moment_inner ( f_ind_ext , mu_ind , theta_ind )

  USE pim,              ONLY :  ldip_polar, omegakO
  USE coulomb,          ONLY :  alphaES
  USE ewald,            ONLY :  multipole_ES
  USE config,           ONLY :  natm , itype , atypei, ntype , poldipia, invpoldipia, atype 

  implicit none

  ! global
  real(kind=dp) , intent ( in  ) :: f_ind_ext ( 3 , natm ) 
  real(kind=dp) , intent ( out ) :: mu_ind    ( 3 , natm ) 
  real(kind=dp) , intent ( out ) :: theta_ind    ( 3 , 3 , natm ) 

  ! local
  real(kind=dp), dimension(:,:), allocatable :: mu_prev , f_ind_total , zerovec , f_ind
  real(kind=dp) :: efg_dummy(3,3,natm)
  real(kind=dp) :: alphaES_save
  real(kind=dp) :: rmsd_inner, rmsd_ref 
  integer :: alpha , beta
  integer :: ia , it , ja  
  integer :: iscf_inner 
  logical :: task_ind(6)
  
  allocate ( mu_prev(3,natm) , f_ind_total ( 3 , natm ) , zerovec(3,natm) , f_ind (3,natm) )
  zerovec      = 0.0_dp
  f_ind        = 0.0_dp
  f_ind_total  = 0.0_dp
  alphaES_save = alphaES 
  rmsd_inner   = HUGE(0.0_dp)
  iscf_inner   = 1

  task_ind(1)  = .false.
  task_ind(2)  = .false.
  task_ind(3)  = .true.
  task_ind(4)  = .false.
  task_ind(5)  = .false.
  task_ind(6)  = .false.

  CALL get_rmsd_mu ( rmsd_ref , f_ind_ext , zerovec )

  f_ind_total = f_ind_ext
  ! ---------------------------------------------------------------
  ! \mu_{i,\alpha} =  alpha_{i,\alpha,\beta} * E_{i,\beta}
  ! ---------------------------------------------------------------
  ! note on units :
  ! everything are in internal units :
  !    [ f_ind_ext ] = e / A^2
  !    [ poldipia  ] = A ^ 3
  !    [ mu_ind ] = e A 
  ! ---------------------------------------------------------------
  mu_ind = 0.0_dp
  do ia = 1 , natm 
    it = itype(ia) 
    if ( .not. ldip_polar ( it ) ) cycle
    do alpha = 1 , 3 
      do beta = 1 , 3  
        mu_ind ( alpha , ia ) = mu_ind ( alpha , ia ) + poldipia ( alpha , beta  , ia ) * ( f_ind_total ( beta , ia ) )
      enddo
    enddo
  enddo
  mu_ind = omegakO * mu_ind  

!  do ia=1,natm
!    write(*,'(4f12.8)') mu_ind(1,ia), mu_ind(2,ia) , mu_ind(3,ia) ,omegakO
!  enddo
  !stop

  inner_loop : do while ( iscf_inner .le. 10 .and. rmsd_inner .gt. 0.1_dp * rmsd_ref )  
    mu_prev      = mu_ind 
   
    ! alpha is reduced for the short-range dipole-dipole interaction
    ! only real part is calculated in the inner loop
    alphaES = 0.001_dp
    ! note that do_strucfact , use_ckrskr has no effect as do_rec=.false.
    CALL  multipole_ES ( f_ind , efg_dummy, mu_ind , theta_ind , task_ind , &
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.false., do_dir=.true. , &
                         do_strucfact=.false. , use_ckrskr=.false. ) 

    f_ind_total = f_ind_ext + f_ind 
    CALL get_rmsd_mu ( rmsd_inner , f_ind_total , mu_ind )
    ! ---------------------------------------------------------------
    ! \mu_{i,\alpha} =  alpha_{i,\alpha,\beta} * E_{i,\beta}
    ! ---------------------------------------------------------------
    ! note on units :
    ! everything are in internal units :
    !    [ f_ind_ext ] = e / A^2
    !    [ poldipia  ] = A ^ 3
    !    [ mu_ind ] = e A 
    ! ---------------------------------------------------------------
    mu_ind = 0.0_dp
    do ia = 1 , natm 
      it = itype(ia) 
      if ( .not. ldip_polar ( it ) ) cycle
      do alpha = 1 , 3 
        do beta = 1 , 3  
          mu_ind ( alpha , ia ) = mu_ind ( alpha , ia ) + poldipia ( alpha , beta  , ia ) * ( f_ind_total ( beta , ia ) )
        enddo
      enddo
    enddo
    mu_ind = ( 1.0_dp - omegakO ) * mu_prev + omegakO * mu_ind  


#ifdef debug_scf_kO_inner
  io_printnode WRITE(stdout,'(a,i3,a,e16.8)')  'inner loop      step : ',iscf_inner,' rmsd_inner = ',rmsd_inner
#endif
  
    iscf_inner = iscf_inner + 1 

  enddo inner_loop 

  io_printnode WRITE(stdout,'(a,i3,a,e16.8)')  'inner loop converged in',iscf_inner-1,' steps    rmsd_inner = ',rmsd_inner
  ! restore alphaES
  alphaES = alphaES_save

  deallocate ( mu_prev , f_ind_total , zerovec , f_ind)

  return

END SUBROUTINE induced_moment_inner

! *********************** SUBROUTINE induced_moment ****************************
!> \brief
!! this subroutine calculates the induced moment from the total Electric field and the
!! polarizability tensor
!! Basically used in the SCF loop to get induced_moment
!! \f$ \mu_{i,\alpha} =  p_{i,\alpha,\beta} * E_{i,\beta} \f$
!> \param[in]  Efield electric field vector define at ion position 
!> \param[out] mu_ind induced electric dipole define at ion position 
!> \note
!! poldipia is the polarizability tensor
! ******************************************************************************
SUBROUTINE induced_moment ( Efield , mu_ind )

  USE pim,              ONLY :  ldip_polar
  USE config,           ONLY :  natm , itype , atypei, ntype , poldipia, atype 

  implicit none

  ! global
  real(kind=dp) , intent ( in  ) :: Efield ( 3 , natm ) 
  real(kind=dp) , intent ( out ) :: mu_ind ( 3 , natm ) 
  

  ! local 
  integer :: alpha , beta
  integer :: ia , it , ja  

  ! ---------------------------------------------------------------
  ! \mu_{i,\alpha} =  alpha_{i,\alpha,\beta} * E_{i,\beta}
  ! ---------------------------------------------------------------
  ! note on units :
  ! everything are in internal units :
  !    [ Efield ] = e / A^2
  !    [ poldipia  ] = A ^ 3
  !    [ mu_ind ] = e A 
  ! ---------------------------------------------------------------
  mu_ind = 0.0_dp
  do ia = 1 , natm 
    it = itype(ia) 
    if ( .not. ldip_polar ( it ) ) cycle
    do alpha = 1 , 3 
      do beta = 1 , 3  
        mu_ind ( alpha , ia ) = mu_ind ( alpha , ia ) + poldipia ( alpha , beta , ia ) * Efield ( beta , ia )  
      enddo
    enddo
  enddo

  return

END SUBROUTINE induced_moment

! *********************** SUBROUTINE get_rmsd_mu ****************************
!> \brief
! ***************************************************************************
SUBROUTINE get_rmsd_mu ( rmsd , g0 , mu )

  USE pim,              ONLY :  ldip_polar
  USE config,           ONLY :  natm , itype , poldipia

  implicit none
  ! global
  real(kind=dp) :: rmsd
  real(kind=dp) :: g0(3,natm)
  real(kind=dp) :: mu(3,natm)
  ! local
  integer :: alpha , ia ,  it
  integer :: npol

  rmsd = 0.0_dp
  npol=0
  do ia=1, natm
    it = itype( ia)
    if ( .not. ldip_polar( it ) ) cycle
    npol = npol + 1
    do alpha = 1 , 3
      if ( poldipia ( alpha , alpha  , ia ) .eq. 0.0_dp ) cycle
      rmsd  = rmsd  +  ( mu ( alpha , ia ) / poldipia ( alpha , alpha , ia ) - g0 ( alpha , ia ) ) ** 2.0_dp
    enddo
  enddo
  rmsd = SQRT ( rmsd /  REAL(npol,kind=dp) )

  return

END SUBROUTINE get_rmsd_mu

! *********************** SUBROUTINE moment_from_pola_scf_kO ********************
!
!> \brief
!! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
!! The evaluation is done self-consistently starting from the the field due the
!! point charges only.
!
!> \param[out] mu_ind induced electric dipole from polarizabilities
!
!> \note
!! The stopping criteria is governed by conv_tol_ind
!
! ******************************************************************************
SUBROUTINE moment_from_pola_scf_kO_v1 ( mu_ind , theta_ind , didpim )

  USE constants,        ONLY :  coul_unit
  USE config,           ONLY :  natm , atype , fx , fy , fz , ntype , dipia , quadia , qia , ntypemax, poldipia , invpoldipia , itype 
  USE control,          ONLY :  longrange , calc
  USE thermodynamic,    ONLY :  u_pol, u_coul
  USE time,             ONLY :  time_moment_from_pola

  USE coulomb,          ONLY :  alphaES, qch, dip, dipia_ind_t
  USE ewald,            ONLY :  multipole_ES
  USE pim,              ONLY :  ldip_polar, max_scf_pol_iter, min_scf_pol_iter, conv_tol_ind, extrapolate_order,algo_ext_dipole

#ifdef debug_print_dipff_scf  
  USE md,               ONLY :  itime
#endif

  implicit none

  ! global
  real(kind=dp) , intent (out) :: mu_ind ( : , : ) 
  real(kind=dp) , intent (out) :: theta_ind ( :, : , : ) 

  ! local
  integer :: ia , iscf , it , alpha, beta
  logical :: linduced , didpim
  real(kind=dp) :: alphaES_save
  real(kind=dp) :: u_coul_stat , rmsd , u_coul_pol, u_coul_ind
  real(kind=dp) :: Efield( 3 , natm ) , Efield_stat ( 3 , natm ) , Efield_ind ( 3 , natm ), EfieldG_stat ( 3 , 3 , natm ) , EfieldG_ind ( 3 , 3 , natm ) , Efield_ind_dir ( 3 , natm )
  real(kind=dp) :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 
  logical       :: task_static (3), task_ind(3), ldip
  real(kind=dp) , dimension (:,:) , allocatable :: f_ind , dmu_ind, zerovec
  real(kind=dp) , dimension (:,:,:) , allocatable :: dtheta_ind 
#ifdef debug_print_dipff_scf
  integer :: kkkk
#endif

#ifdef MPI
  dectime
#endif

#ifdef debug_print_dipff_scf
  kkkk=100000*itime+100000
#endif

#ifdef MPI
  statime
#endif
  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  didpim   = .false.
  do it = 1 , ntype
    if ( ldip_polar ( it ) ) linduced = .true.
  enddo
  if ( .not. linduced ) then
#ifdef debug
    write(stdout,'(a,e16.8)') 'quick return from moment_from_pola_scf',mu_ind(1,1)
#endif
    return
  endif
  didpim = linduced

  allocate ( f_ind ( 3 , natm ) , dmu_ind( 3 , natm ) , zerovec(3,natm) , dtheta_ind(3,3,natm) ) 
  f_ind = 0.0_dp
  dmu_ind = 0.0_dp
  zerovec = 0.0_dp

  ! ========================================================
  !  calculate static Efield ( charge + static dipoles )
  ! ========================================================
  u_pol = 0.0_dp
  u_coul = 0.0_dp
  Efield_stat = 0.0_dp
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  ! =============================================
  !  set coulombic task energy , forces (field) and virial
  ! =============================================
  ldip=.false.
  task_static = .false.      
  task_static(1) = .true.      
  if ( any (dip .ne. 0.0d0 ) ) ldip = .true.
  if ( ldip ) then
    task_static = .true.        
  endif
  ! ==============================================
  !          main ewald subroutine 
  ! ==============================================
  CALL  multipole_ES ( Efield_stat , EfieldG_stat , dipia , quadia , task_static , & 
                       damp_ind =.true. , do_efield=.true. , do_efg=.false. , & 
                       do_forces=.false. , do_stress = .false. , do_rec = .true. , do_dir = .true. , do_strucfact= .true. , use_ckrskr=.true. ) 
  u_coul_stat = u_coul 
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp
  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  iscf = 0
  rmsd = HUGE(0.0d0)
  ! =========================
  !  charges are set to zero 
  ! =========================
  alphaES_save = alphaES
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0_dp
  qia = 0.0_dp
  task_ind(1) = .false. 
  task_ind(2) = .false. 
  task_ind(3) = .true. 
  io_printnode WRITE ( stdout ,'(a)') '' 
  f_ind = Efield   
  Efield_ind_dir = 0.0_dp

  lseparator_ioprintnode(stdout)
  io_printnode WRITE( stdout , '(a)' ) '                running scf_kO_v1 algorithm                  '
  lseparator_ioprintnode(stdout)
  ! =============================
  !           SCF LOOP
  ! =============================
  do while ( ( iscf < max_scf_pol_iter ) .and. ( rmsd .gt. conv_tol_ind )  .or. ( iscf < min_scf_pol_iter  ) )

    iscf = iscf + 1

    ! ==========================================================
    !  calculate mu_ind from Efield = Efield_stat + Efield_ind
    ! ==========================================================
    if ( iscf.ne.1 .or. ( calc .ne. 'md' .and.  calc .ne. 'opt' ) ) then

          CALL induced_moment_inner ( f_ind , dmu_ind , dtheta_ind )         
          mu_ind = mu_ind + dmu_ind
          !print*,'here 1',dmu_ind(1,3)
    else

          !print*,'here 2',algo_ext_dipole
      if ( algo_ext_dipole .eq. 'poly' ) CALL extrapolate_dipole_poly ( mu_ind ) 
      if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield , key=1 )  ! predictor

#ifdef debug_mu
      if ( ionode ) then
        WRITE ( stdout , '(a)' )     'Induced dipoles at atoms from extrapolation: '
        do ia = 1 , natm 
          WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
          ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
        enddo
        blankline(stdout)
      endif
#endif

    endif
    ! ===============================================================
    !      u_pol = 1/ ( 2 alpha  ) * | mu_ind | ^ 2
    ! ---------------------------------------------------------------
    ! note on units :
    !    [ poldipia  ] = A ^ 3
    !    [ mu_ind ] = e A 
    !    [ u_pol  ] = e^2 / A ! electrostatic internal energy  
    ! ===============================================================
    u_pol = 0.0_dp
    do ia = 1 , natm
      do alpha = 1 , 3
        do beta = 1 , 3
          u_pol = u_pol + mu_ind ( alpha , ia ) * invpoldipia ( alpha , beta , ia ) * mu_ind ( beta , ia )
        enddo
      enddo
    enddo
    u_pol = u_pol * 0.5_dp


    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true., do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    u_coul_ind = u_coul 
    u_coul_pol = u_coul_stat-u_coul_ind

    Efield = Efield_stat + Efield_ind + Efield_ind_dir
    fx      = 0.0_dp
    fy      = 0.0_dp
    fz      = 0.0_dp

    ! ASPC corrector 
    if ( iscf.eq.1 .and. calc .eq. 'md' ) then
      if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield  , 2 ) ! corrector
    endif
    
    ! ===================
    !  stopping criteria
    ! ===================
    CALL get_rmsd_mu ( rmsd , Efield , mu_ind )
!kO
    do ia = 1 , natm
      do alpha = 1 , 3
        it = itype( ia)
        if ( .not. ldip_polar( it ) ) cycle
        f_ind(alpha,ia) = Efield(alpha,ia) - mu_ind(alpha,ia)/poldipia(alpha,alpha,ia)
      enddo
    enddo
!kO

    ! output
    if ( calc .ne. 'opt' ) then
#ifdef debug
    io_printnode WRITE ( stdout ,'(a,i4,5(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' u_coul (qq)  = ', u_coul_stat, ' u_coul (dd)  = ', u_coul_ind,' u_coul_pol = ', u_coul_pol, ' rmsd       = ', rmsd
#endif
    io_printnode WRITE ( stdout ,'(a,i4,2(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' rmsd       = ', rmsd
    endif

 
#ifdef debug_print_dipff_scf
  do ia = 1 , natm
  !if ( mu_ind ( 1 , ia ) .eq. 0._dp ) cycle
  write(kkkk,'(3e16.8)') (mu_ind ( alpha , ia ),alpha=1,3)
  enddo
  kkkk = kkkk + 1 
#endif

  enddo ! end of SCF loop

  ! ===========================
  !  charge/force info is recovered
  ! ===========================
  alphaES = alphaES_save
  qch = qch_tmp
  qia = qia_tmp

  if ( ioprintnode .and.  calc .ne. 'opt' ) then
    blankline(stdout)
    WRITE ( stdout , '(a,i6,a)')            'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3,a,e10.3,a)') 'Electric field is converged at ',rmsd,' ( ',conv_tol_ind,' ) '
    blankline(stdout)
  endif
#ifdef debug_mu
  if ( ionode ) then
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm 
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
    enddo
    blankline(stdout)
  endif
#endif
  ! store induced dipole at t 
  dipia_ind_t(1,:,:) = mu_ind
  
#ifdef debug_extrapolate
    if ( ioprintnode ) then
      print*,'previous step',itime
      do ia=1,natm
        write(stdout,'(i5,<extrapolate_order+2>e16.8)') ia, mu_ind ( 1 , ia ), (dipia_ind_t ( t , 1 , ia ) , t=1, extrapolate_order+1 )
      enddo
    endif
#endif


  deallocate ( f_ind , dmu_ind , zerovec ,dtheta_ind) 

#ifdef MPI
  stotime
  addtime(time_moment_from_pola)
#endif
  ! temporary
  theta_ind=0.0_dp

  return

END SUBROUTINE moment_from_pola_scf_kO_v1

! *********************** SUBROUTINE moment_from_pola_scf_kO ********************
!
!> \brief
!! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
!! The evaluation is done self-consistently starting from the the field due the
!! point charges only.
!
!> \param[out] mu_ind induced electric dipole from polarizabilities
!
!> \note
!! The stopping criteria is governed by conv_tol_ind
!
! ******************************************************************************
SUBROUTINE moment_from_pola_scf_kO_v2 ( mu_ind , theta_ind , didpim ) 

  USE constants,        ONLY :  coul_unit
  USE config,           ONLY :  natm , atype , fx , fy , fz , ntype , dipia , quadia , qia , ntypemax, poldipia , invpoldipia , itype 
  USE control,          ONLY :  longrange , calc
  USE thermodynamic,    ONLY :  u_pol, u_coul
  USE time,             ONLY :  time_moment_from_pola
  USE coulomb,          ONLY :  alphaES, qch, dip, dipia_ind_t
  USE ewald,            ONLY :  multipole_ES
  USE pim,              ONLY :  ldip_polar, max_scf_pol_iter, min_scf_pol_iter, conv_tol_ind, extrapolate_order,algo_ext_dipole

#ifdef debug_print_dipff_scf  
  USE md,               ONLY :  itime
#endif

  implicit none

  ! global
  real(kind=dp) , intent (out) :: mu_ind ( : , : ) 
  real(kind=dp) , intent (out) :: theta_ind ( : , : , : ) 

  ! local
  integer :: ia , iscf , it , alpha, beta
  logical :: linduced , didpim
  real(kind=dp) :: tttt , alphaES_save
  real(kind=dp) :: u_coul_stat , rmsd , u_coul_pol, u_coul_ind
  real(kind=dp) :: Efield( 3 , natm ) , Efield_stat ( 3 , natm ) , Efield_ind ( 3 , natm ), EfieldG_stat ( 3 , 3 , natm ) , EfieldG_ind ( 3 , 3 , natm ) , tmp ( 3 ,natm )
  real(kind=dp) :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 
  logical       :: task_static (3), task_ind(3), ldip 
  real(kind=dp) , dimension (:,:) , allocatable :: f_ind , dmu_ind, f_ind_prev , gmin , mu_prev, dmin, zerovec, g_ind
  real(kind=dp) , dimension (:,:,:) , allocatable :: dtheta_ind
#ifdef debug_print_dipff_scf
  integer :: kkkk
#endif
#ifdef MPI
  dectime
#endif

#ifdef debug_print_dipff_scf
  kkkk=100000*itime+100000
#endif

  tttt=0.0_dp
#ifdef MPI
  statime
#endif
  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  didpim   = .false.
  do it = 1 , ntype
    if ( ldip_polar ( it ) ) linduced = .true.
  enddo
  if ( .not. linduced ) then
#ifdef debug
    write(stdout,'(a,e16.8)') 'quick return from moment_from_pola_scf',mu_ind(1,1)
#endif
    return
  endif
  didpim = linduced

  allocate ( f_ind ( 3 , natm ) , dmu_ind( 3 , natm )  , f_ind_prev( 3 ,natm ) , gmin( 3 , natm ) , mu_prev( 3 , natm ) , dmin(3,natm) , zerovec(3,natm) , g_ind(3,natm ) ) 
  allocate ( dtheta_ind (3,3,natm) )
  f_ind = 0.0_dp
  dmu_ind = 0.0_dp
  f_ind_prev = 0.0_dp 
  gmin =0.0_dp
  mu_prev = 0.0_dp 
  dmin = 0.0_dp
  g_ind = 0.0_dp
  zerovec = 0.0_dp

  ! ========================================================
  !  calculate static Efield ( charge + static dipoles )
  ! ========================================================
  u_pol = 0.0_dp
  u_coul = 0.0_dp
  Efield_stat = 0.0_dp
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  ! =============================================
  !  set coulombic task energy , forces (field) and virial
  ! =============================================
  ldip=.false.
  task_static = .false.      
  task_static(1) = .true.      
  if ( any (dip .ne. 0.0d0 ) ) ldip = .true.
  if ( ldip ) then
    task_static = .true.        
  endif
  ! ==============================================
  !          main ewald subroutine 
  ! ==============================================
  CALL  multipole_ES ( Efield_stat , EfieldG_stat , dipia , quadia , task_static , & 
                       damp_ind =.true. , do_efield=.true. , do_efg=.false. , & 
                       do_forces=.false. , do_stress = .false. , do_rec = .true. , do_dir = .true. , do_strucfact= .true. , use_ckrskr=.true. ) 
  u_coul_stat = u_coul 
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp
  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  iscf = 0
  rmsd = HUGE(0.0d0)
  ! =========================
  !  charges are set to zero 
  ! =========================
  alphaES_save = alphaES
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0_dp
  qia = 0.0_dp
  task_ind(1) = .false. 
  task_ind(2) = .false. 
  task_ind(3) = .true. 
  io_printnode WRITE ( stdout ,'(a)') '' 
  f_ind = Efield   
  mu_prev = 0.0_dp

  lseparator_ioprintnode(stdout)
  io_printnode WRITE( stdout , '(a)' ) '                running scf_kO_v2 algorithm                  '
  lseparator_ioprintnode(stdout)
  ! =============================
  !           SCF LOOP
  ! =============================
  do while ( ( iscf < max_scf_pol_iter ) .and. ( rmsd .gt. conv_tol_ind )  .or. ( iscf < min_scf_pol_iter  ) )

    iscf = iscf + 1

    ! ==========================================================
    !  calculate mu_ind from Efield = Efield_stat + Efield_ind
    ! ==========================================================
    if ( iscf.ne.1 .or. ( calc .ne. 'md' .and.  calc .ne. 'opt' ) ) then

        if ( iscf .le. 2 ) then
          CALL induced_moment_inner ( f_ind , dmu_ind , dtheta_ind )         
          mu_ind = mu_ind + dmu_ind
        else
          CALL induced_moment_inner ( gmin , dmu_ind  , dtheta_ind )          
          mu_ind = mu_ind + dmu_ind
        endif

    else

      if ( algo_ext_dipole .eq. 'poly' ) CALL extrapolate_dipole_poly ( mu_ind ) 
      if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield , key=1 )  ! predictor

#ifdef debug_mu
      if ( ionode ) then
        WRITE ( stdout , '(a)' )     'Induced dipoles at atoms from extrapolation: '
        do ia = 1 , natm 
          WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
          ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
        enddo
        blankline(stdout)
      endif
#endif

    endif
    ! ===============================================================
    !      u_pol = 1/ ( 2 alpha  ) * | mu_ind | ^ 2
    ! ---------------------------------------------------------------
    ! note on units :
    !    [ poldipia  ] = A ^ 3
    !    [ mu_ind ] = e A 
    !    [ u_pol  ] = e^2 / A ! electrostatic internal energy  
    ! ===============================================================
    u_pol = 0.0_dp
    do ia = 1 , natm
      do alpha = 1 , 3
        do beta = 1 , 3
          u_pol = u_pol + mu_ind ( alpha , ia ) * invpoldipia ( alpha , beta , ia ) * mu_ind ( beta , ia )
        enddo
      enddo
    enddo
    u_pol = u_pol * 0.5_dp


    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true., do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    u_coul_ind = u_coul 
    u_coul_pol = u_coul_stat-u_coul_ind

    Efield = Efield_stat + Efield_ind
    fx      = 0.0_dp
    fy      = 0.0_dp
    fz      = 0.0_dp

    ! ASPC corrector 
    if ( iscf.eq.1 .and. calc .eq. 'md' ) then
      if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield  , 2 ) ! corrector
    endif
    
    ! ===================
    !  stopping criteria
    ! ===================
!kO
    do ia = 1 , natm
      do alpha = 1 , 3
        it = itype( ia)
        if ( .not. ldip_polar( it ) ) cycle
        f_ind(alpha,ia) = Efield(alpha,ia) - mu_ind(alpha,ia)/poldipia(alpha,alpha,ia)
      enddo
    enddo

    if ( iscf .ne. 1 ) then
      CALL find_min ( f_ind_prev , f_ind , (mu_ind - mu_prev) , gmin , dmin )
      do ia = 1 , natm
        do alpha = 1 , 3
          it = itype( ia)
          if ( .not. ldip_polar( it ) ) cycle
          g_ind(alpha,ia) = gmin(alpha,ia) + mu_ind(alpha,ia)/poldipia(alpha,alpha,ia)
        enddo
      enddo
      mu_ind = mu_prev + dmin
      CALL get_rmsd_mu ( rmsd , gmin , zerovec )
    endif
!kO

    ! output
    if ( calc .ne. 'opt' ) then
#ifdef debug
    io_printnode WRITE ( stdout ,'(a,i4,5(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' u_coul (qq)  = ', u_coul_stat, ' u_coul (dd)  = ', u_coul_ind,' u_coul_pol = ', u_coul_pol, ' rmsd       = ', rmsd
#endif
    io_printnode WRITE ( stdout ,'(a,i4,2(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' rmsd       = ', rmsd
    endif

#ifdef debug_print_dipff_scf
  do ia = 1 , natm
  if ( mu_ind ( 1 , ia ) .eq. 0._dp ) cycle
  write(kkkk,'(3e16.8)') (mu_ind ( alpha , ia ),alpha=1,3)
  enddo
  kkkk = kkkk + 1 
#endif

  ! save previous
  mu_prev = mu_ind
  f_ind_prev = gmin 

  enddo ! end of SCF loop

  ! ===========================
  !  charge/force info is recovered
  ! ===========================
  alphaES = alphaES_save
  qch = qch_tmp
  qia = qia_tmp

  if ( ioprintnode .and.  calc .ne. 'opt' ) then
    blankline(stdout)
    WRITE ( stdout , '(a,i6,a)')            'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3,a,e10.3,a)') 'Electric field is converged at ',rmsd,' ( ',conv_tol_ind,' ) '
    blankline(stdout)
  endif
#ifdef debug_mu
  if ( ionode ) then
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm 
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
    enddo
    blankline(stdout)
  endif
#endif
  ! store induced dipole at t 
  dipia_ind_t(1,:,:) = mu_ind
  
#ifdef debug_extrapolate
    if ( ioprintnode ) then
      print*,'previous step',itime
      do ia=1,natm
        write(stdout,'(i5,<extrapolate_order+2>e16.8)') ia, mu_ind ( 1 , ia ), (dipia_ind_t ( t , 1 , ia ) , t=1, extrapolate_order+1 )
      enddo
    endif
#endif


  deallocate ( f_ind , dmu_ind  , f_ind_prev , gmin , mu_prev , dmin , zerovec , g_ind) 

#ifdef MPI
  stotime
  addtime(time_moment_from_pola)
#endif


  ! temporary
  theta_ind=0.0_dp

  return

END SUBROUTINE moment_from_pola_scf_kO_v2

! *********************** SUBROUTINE moment_from_pola_scf_kO ********************
!
!> \brief
!! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
!! The evaluation is done self-consistently starting from the the field due the
!! point charges only.
!
!> \param[out] mu_ind induced electric dipole from polarizabilities
!
!> \note
!! The stopping criteria is governed by conv_tol_ind
!
! ******************************************************************************
SUBROUTINE moment_from_pola_scf_kO_v3 ( mu_ind , theta_ind , didpim ) 

  USE constants,        ONLY :  coul_unit
  USE config,           ONLY :  natm , atype , fx , fy , fz , ntype , dipia , quadia , qia , ntypemax, poldipia , invpoldipia , itype 
  USE control,          ONLY :  longrange , calc
  USE thermodynamic,    ONLY :  u_pol, u_coul
  USE time,             ONLY :  time_moment_from_pola
  USE coulomb,          ONLY :  alphaES, qch, dip, dipia_ind_t
  USE ewald,            ONLY :  multipole_ES
  USE pim,              ONLY :  ldip_polar, max_scf_pol_iter, min_scf_pol_iter, conv_tol_ind, extrapolate_order,algo_ext_dipole
#ifdef debug_print_dipff_scf  
  USE md,               ONLY :  itime
#endif

  implicit none

  ! global
  real(kind=dp) , intent (out) :: mu_ind ( : , : ) 
  real(kind=dp) , intent (out) :: theta_ind ( : , : , : ) 

  ! local
  integer :: ia , iscf , it , alpha, beta
  logical :: linduced , didpim
  real(kind=dp) :: tttt , alphaES_save
  real(kind=dp) :: u_coul_stat , rmsd , u_coul_pol, u_coul_ind
  real(kind=dp) :: Efield( 3 , natm ) , Efield_stat ( 3 , natm ) , Efield_ind ( 3 , natm ), EfieldG_stat ( 3 , 3 , natm ) , EfieldG_ind ( 3 , 3 , natm )
  real(kind=dp) :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 
  logical       :: task_static (3), task_ind(3), ldip , lfirst_inner
  real(kind=dp) , dimension (:,:) , allocatable :: f_ind , dmu_ind, f_ind_prev , gmin , mu_prev, dmin, zerovec, g_ind
#ifdef debug_print_dipff_scf
  integer :: kkkk
#endif
#ifdef MPI
  dectime
#endif

#ifdef debug_print_dipff_scf
  kkkk=100000*itime+100000
#endif

  tttt=0.0_dp
#ifdef MPI
  statime
#endif
  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  didpim   = .false.
  do it = 1 , ntype
    if ( ldip_polar ( it ) ) linduced = .true.
  enddo
  if ( .not. linduced ) then
#ifdef debug
    write(stdout,'(a,e16.8)') 'quick return from moment_from_pola_scf',mu_ind(1,1)
#endif
    return
  endif
  didpim = linduced

  allocate ( f_ind ( 3 , natm ) , dmu_ind( 3 , natm )  , f_ind_prev( 3 ,natm ) , gmin( 3 , natm ) , mu_prev( 3 , natm ) , dmin(3,natm) , zerovec(3,natm) , g_ind(3,natm ) ) 
  f_ind = 0.0_dp
  dmu_ind = 0.0_dp
  f_ind_prev = 0.0_dp 
  gmin =0.0_dp
  mu_prev = 0.0_dp 
  dmin = 0.0_dp
  g_ind = 0.0_dp
  zerovec = 0.0_dp

  ! ========================================================
  !  calculate static Efield ( charge + static dipoles )
  ! ========================================================
  u_pol = 0.0_dp
  u_coul = 0.0_dp
  Efield_stat = 0.0_dp
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  ! =============================================
  !  set coulombic task energy , forces (field) and virial
  ! =============================================
  ldip=.false.
  task_static = .false.      
  task_static(1) = .true.      
  if ( any (dip .ne. 0.0d0 ) ) ldip = .true.
  if ( ldip ) then
    task_static = .true.        
  endif
  ! ==============================================
  !          main ewald subroutine 
  ! ==============================================
  CALL  multipole_ES ( Efield_stat , EfieldG_stat , dipia , quadia , task_static , & 
                       damp_ind =.true. , do_efield=.true. , do_efg=.false. , & 
                       do_forces=.false. , do_stress = .false. , do_rec = .true. , do_dir = .true. , do_strucfact= .true. , use_ckrskr=.true. ) 
  u_coul_stat = u_coul 
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp
  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  iscf = 0
  rmsd = HUGE(0.0d0)
  ! =========================
  !  charges are set to zero 
  ! =========================
  alphaES_save = alphaES
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0_dp
  qia = 0.0_dp
  task_ind(1) = .false. 
  task_ind(2) = .false. 
  task_ind(3) = .true. 
  io_printnode WRITE ( stdout ,'(a)') '' 
  f_ind = Efield   
  mu_prev = 0.0_dp

  !  extrapolate from previous md steps :
  !   two methods : 
  !                  simple linear/polynomial extrapolation   ( algo_ext_dipole .eq. 'poly' ) : instable for extrapolate_order > 1  
  !                  ASPC (Always Stable Predictor Corrector) ( algo_ext_dipole .eq. 'aspc  ) by Kolafa J. Comp. Chem. v25 n3 (2003)
  if ( algo_ext_dipole .eq. 'poly' ) CALL extrapolate_dipole_poly ( mu_ind ) 
  if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield , key=1 )  ! predictor

  lfirst_inner = .true.
  if ( all ( mu_ind .eq. 0.0d0 ) ) lfirst_inner = .false.
  if ( extrapolate_order .lt. 1 )  lfirst_inner = .false. 
  
  if ( lfirst_inner ) then
    print*,' dipoles : l_firstinner=.true. '


  

    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , &
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true.,  &
                         do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    Efield = Efield_stat + Efield_ind



  else

    print*,'null dipoles or extrapolates : lfirst inner = .false. '
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true.,  &
                         do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    ! ASPC corrector 
    if ( calc .eq. 'md' .and. algo_ext_dipole .eq. 'aspc' ) then
      CALL extrapolate_dipole_aspc ( mu_ind , Efield_ind , key=2 ) ! corrector
    endif
    
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true.,  &
                         do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )
    
    Efield = Efield_stat + Efield_ind
    print*,'mu_ind',mu_ind(1,300)
  endif

  lseparator_ioprintnode(stdout)
  io_printnode WRITE( stdout , '(a)' ) '                running scf_kO_v3 algorithm                  '
  lseparator_ioprintnode(stdout)
  ! =============================
  !           SCF LOOP
  ! =============================
  outer_loop : do while ( ( iscf < max_scf_pol_iter ) .and. ( rmsd .gt. conv_tol_ind )  .or. ( iscf < min_scf_pol_iter  ) )

    iscf = iscf + 1

    CALL find_min ( f_ind_prev , f_ind , (mu_ind - mu_prev) , gmin , dmin )
    mu_ind = mu_ind + dmu_ind
    CALL get_rmsd_mu ( rmsd , gmin , mu_ind )

    ! ===============================================================
    !      u_pol = 1/ ( 2 alpha  ) * | mu_ind | ^ 2
    ! ---------------------------------------------------------------
    ! note on units :
    !    [ poldipia  ] = A ^ 3
    !    [ mu_ind ] = e A 
    !    [ u_pol  ] = e^2 / A ! electrostatic internal energy  
    ! ===============================================================
    u_pol = 0.0_dp
    do ia = 1 , natm
      do alpha = 1 , 3
        do beta = 1 , 3
          u_pol = u_pol + mu_ind ( alpha , ia ) * invpoldipia ( alpha , beta , ia ) * mu_ind ( beta , ia )
        enddo
      enddo
    enddo
    u_pol = u_pol * 0.5_dp


    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true., do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    u_coul_ind = u_coul 
    u_coul_pol = u_coul_stat-u_coul_ind

    Efield = Efield_stat + Efield_ind
    fx      = 0.0_dp
    fy      = 0.0_dp
    fz      = 0.0_dp

    do ia = 1 , natm
      do alpha = 1 , 3
        it = itype( ia)
        if ( .not. ldip_polar( it ) ) cycle
        f_ind(alpha,ia) = Efield(alpha,ia) - mu_ind(alpha,ia)/poldipia(alpha,alpha,ia)
      enddo
    enddo

!    if ( algo_moment_from_pola .eq. 'scf' ) then
!      CALL get_rmsd_mu ( rmsd , Efield , mu_ind )
!    else if ( algo_moment_from_pola .eq. 'scf_kO' ) then
      if ( iscf .ne. 1 ) then
        CALL find_min ( f_ind_prev , f_ind , (mu_ind - mu_prev) , gmin , dmin )
        do ia = 1 , natm
          do alpha = 1 , 3
            it = itype( ia)
            if ( .not. ldip_polar( it ) ) cycle
            g_ind(alpha,ia) = gmin(alpha,ia) + mu_ind(alpha,ia)/poldipia(alpha,alpha,ia)
          enddo
        enddo
        mu_ind = mu_prev + dmin
      CALL get_rmsd_mu ( rmsd , gmin , zerovec )
      endif
!    endif
!kO

!      CALL get_rmsd_mu ( rmsd , Efield , mu_ind )

    ! output
    if ( calc .ne. 'opt' ) then
#ifdef debug
    io_printnode WRITE ( stdout ,'(a,i4,5(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' u_coul (qq)  = ', u_coul_stat, ' u_coul (dd)  = ', u_coul_ind,' u_coul_pol = ', u_coul_pol, ' rmsd       = ', rmsd
#endif
    io_printnode WRITE ( stdout ,'(a,i4,2(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' rmsd       = ', rmsd
    endif
 
#ifdef debug_print_dipff_scf
  do ia = 1 , natm
  if ( mu_ind ( 1 , ia ) .eq. 0._dp ) cycle
  write(kkkk,'(3e16.8)') (mu_ind ( alpha , ia ),alpha=1,3)
  enddo
  kkkk = kkkk + 1 
#endif

  ! save previous
  mu_prev = mu_ind
  f_ind_prev = gmin 

  enddo outer_loop ! end of SCF loop

  ! ===========================
  !  field info is recovered
  ! ===========================
  alphaES = alphaES_save
  qch = qch_tmp
  qia = qia_tmp

  if ( ioprintnode .and.  calc .ne. 'opt' ) then
    blankline(stdout)
    WRITE ( stdout , '(a,i6,a)')            'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3,a,e10.3,a)') 'Electric field is converged at ',rmsd,' ( ',conv_tol_ind,' ) '
    blankline(stdout)
  endif
#ifdef debug_mu
  if ( ionode ) then
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm 
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
    enddo
    blankline(stdout)
  endif
#endif
  ! store induced dipole at t 
  dipia_ind_t(1,:,:) = mu_ind
  
#ifdef debug_extrapolate
    if ( ioprintnode ) then
      print*,'previous step',itime
      do ia=1,natm
        write(stdout,'(i5,<extrapolate_order+2>e16.8)') ia, mu_ind ( 1 , ia ), (dipia_ind_t ( t , 1 , ia ) , t=1, extrapolate_order+1 )
      enddo
    endif
#endif


  deallocate ( f_ind , dmu_ind  , f_ind_prev , gmin , mu_prev , dmin , zerovec , g_ind) 

#ifdef MPI
  stotime
  addtime(time_moment_from_pola)
#endif
  ! temporary
  theta_ind=0.0_dp

  return

END SUBROUTINE moment_from_pola_scf_kO_v3

! *********************** SUBROUTINE moment_from_pola_scf_v4 **********************
!
!> \brief
!! This routines evaluates the dipole moment induced by polarizabilities on atoms. 
!! The evaluation is done self-consistently starting from the the field due the
!! point charges only.
!
!> \param[out] mu_ind induced electric dipole from polarizabilities
!
!> \note
!! The stopping criteria is governed by conv_tol_ind
!
! ******************************************************************************
SUBROUTINE moment_from_pola_scf_kO_v4 ( mu_ind , theta_ind , didpim ) 

  USE constants,        ONLY :  coul_unit
  USE config,           ONLY :  natm , atype , fx , fy , fz , ntype , dipia , quadia , qia , ntypemax, poldipia , invpoldipia , itype 
  USE control,          ONLY :  longrange , calc
  USE thermodynamic,    ONLY :  u_pol, u_coul
  USE time,             ONLY :  time_moment_from_pola
  USE coulomb,          ONLY :  alphaES, qch, dip, dipia_ind_t
  USE ewald,            ONLY :  multipole_ES
  USE pim,              ONLY :  ldip_polar, max_scf_pol_iter, min_scf_pol_iter, conv_tol_ind, extrapolate_order,algo_ext_dipole,algo_moment_from_pola
#ifdef debug_print_dipff_scf  
  USE md,               ONLY :  itime
#endif

  implicit none

  ! global
  real(kind=dp) , intent (out) :: mu_ind ( : , : ) 
  real(kind=dp) , intent (out) :: theta_ind ( : , : , : ) 

  ! local
  integer :: ia , iscf , it , alpha, beta
  logical :: linduced , didpim
  real(kind=dp) :: tttt , alphaES_save
  real(kind=dp) :: u_coul_stat , rmsd , u_coul_pol, u_coul_ind
  real(kind=dp) :: Efield( 3 , natm ) , Efield_stat ( 3 , natm ) , Efield_ind ( 3 , natm ), EfieldG_stat ( 3 , 3 , natm ) , EfieldG_ind ( 3 , 3 , natm )
  real(kind=dp) :: qia_tmp ( natm )  , qch_tmp ( ntypemax ) 
  logical       :: task_static (3), task_ind(3), ldip
  real(kind=dp) , dimension (:,:) , allocatable :: F_ind , dmu
  real(kind=dp) , dimension (:,:,:) , allocatable :: F_h , mu_h  
#ifdef debug_print_dipff_scf
  integer :: kkkk
#endif 
#ifdef MPI
  dectime
#endif

#ifdef debug_print_dipff_scf
  kkkk=100000*itime+100000
  print*,kkkk,itime
#endif


  tttt=0.0_dp
#ifdef MPI
  statime
#endif
  ! =========================================================
  !  Is there any polarizability ? if yes linduced = .TRUE.
  ! =========================================================
  linduced = .false.
  didpim   = .false.
  do it = 1 , ntype
    if ( ldip_polar ( it ) ) linduced = .true.
  enddo
  if ( .not. linduced ) then
#ifdef debug
    write(stdout,'(a,e16.8)') 'quick return from moment_from_pola_scf',mu_ind(1,1)
#endif
    return
  endif
  didpim = linduced

  allocate ( F_h ( 3 , natm , 100 ) , mu_h ( 3 , natm , 100 ) , F_ind ( 3 , natm ) , dmu ( 3 , natm ) ) 
  F_ind = 0.0_dp
  F_h   = 0.0_dp
  mu_h  = 0.0_dp


  ! =============================================
  !  calculate static Efield ( charge + dipoles )
  ! =============================================
  u_pol = 0.0_dp
  u_coul = 0.0_dp
  Efield_stat = 0.0_dp
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  ! =============================================
  !  coulombic energy , forces (field) and virial
  ! =============================================
  ldip=.false.
  task_static = .false.      
  task_static(1) = .true.      
  if ( any (dip .ne. 0.0d0 ) ) ldip = .true.
  if ( ldip ) then
    task_static = .true.        
  endif
  ! ==============================================
  !          main ewald subroutine 
  ! ==============================================
  ! do_strucfact = .true. : recalculate exp(ikr) for new r
  ! use_ckrskr  = .true. : use the store exp(ikr) 
  CALL  multipole_ES ( Efield_stat , EfieldG_stat , dipia , quadia , task_static , & 
                       damp_ind =.true. , do_efield=.true. , do_efg=.false. , & 
                       do_forces=.false. , do_stress = .false. , do_rec = .true. , do_dir = .true. , do_strucfact= .true. , use_ckrskr=.true. ) 
  u_coul_stat = u_coul 
  fx      = 0.0_dp
  fy      = 0.0_dp
  fz      = 0.0_dp

  !F_h (:,:,1) = Efield_stat
  !mu_h(:,:,1) = 0.0_dp
  ! =============================================
  !  init total Efield to static only
  ! =============================================
  Efield = Efield_stat

  iscf = 0
  rmsd = HUGE(0.0d0)
  ! =========================
  !  charges are set to zero 
  ! =========================
  alphaES_save = alphaES
  qch_tmp = qch
  qia_tmp = qia
  qch = 0.0_dp
  qia = 0.0_dp
  task_ind(1) = .false. 
  task_ind(2) = .false. 
  task_ind(3) = .true. 
  io_printnode WRITE ( stdout ,'(a)') '' 
  u_pol = 0.0_dp

  lseparator_ioprintnode(stdout)
  io_printnode WRITE( stdout , '(a,a)' ) '                   running ',algo_moment_from_pola 
  lseparator_ioprintnode(stdout)
  ! =============================
  !           SCF LOOP
  ! =============================
  do while ( ( iscf < max_scf_pol_iter ) .and. ( rmsd .gt. conv_tol_ind )  .or. ( iscf < min_scf_pol_iter  ) )

    iscf = iscf + 1

    if ( iscf.ne.1 .or. ( calc .ne. 'md' .and.  calc .ne. 'opt' ) ) then

      CALL induced_moment       ( Efield , mu_ind )  ! Efield in ; mu_ind and u_pol out

    else
       if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield , key=1 )  ! predictor
       if ( all ( mu_ind .eq. 0.0_dp ) ) CALL induced_moment       ( Efield , mu_ind )
!       CALL induced_moment       ( Efield , mu_ind )
    endif

    ! ===============================================================
    !      u_pol = 1/ ( 2 alpha  ) * | mu_ind | ^ 2
    ! ---------------------------------------------------------------
    ! note on units :
    !    [ poldipia  ] = A ^ 3
    !    [ mu_ind ] = e A 
    !    [ u_pol  ] = e^2 / A ! electrostatic internal energy  
    ! ===============================================================
    u_pol = 0.0_dp
    do ia = 1 , natm
      do alpha = 1 , 3
        do beta = 1 , 3
          u_pol = u_pol + mu_ind ( alpha , ia ) * invpoldipia ( alpha , beta , ia ) * mu_ind ( beta , ia )
        enddo
      enddo
    enddo
    u_pol = u_pol * 0.5_dp

    ! ==========================================================
    !  calculate Efield_ind from mu_ind
    !  Efield_ind out , mu_ind in ==> charges and static dipoles = 0
    ! ==========================================================
    CALL  multipole_ES ( Efield_ind , EfieldG_ind , mu_ind , theta_ind , task_ind , & 
                         damp_ind = .false. , do_efield=.true. , do_efg = .false. , &
                         do_forces = .false. , do_stress =.false. , do_rec=.true., do_dir=.true. , do_strucfact=.false. , use_ckrskr = .true. )

    u_coul_ind = u_coul 
    u_coul_pol = u_coul_stat-u_coul_ind
    Efield = Efield_stat + Efield_ind 

    fx      = 0.0_dp
    fy      = 0.0_dp
    fz      = 0.0_dp

    do ia = 1 , natm
      do alpha = 1 , 3
        it = itype( ia)
        if ( .not. ldip_polar( it ) ) cycle
      
        if ( algo_moment_from_pola .eq. 'scf_kO_v4_1' ) then 
          F_ind(alpha,ia) = mu_ind(alpha,ia)/poldipia(alpha,alpha,ia) - Efield_ind(alpha,ia)  !!! find_min_ex
        else
          F_ind(alpha,ia) = mu_ind(alpha,ia)/poldipia(alpha,alpha,ia) - Efield(alpha,ia)       !!! find_min
        endif
      enddo
    enddo

    if ( algo_moment_from_pola .eq. 'scf_kO_v4_1' ) then
      CALL find_min_ex ( F_ind , mu_ind , F_h , mu_h , iscf , Efield_stat )
    else
      F_h(:,:,iscf)  = F_ind
      mu_h(:,:,iscf) = mu_ind 
      if ( iscf .ne. 1 ) CALL find_min ( F_h(:,:,iscf-1) , F_h(:,:,iscf) , (mu_h(:,:,iscf) - mu_h(:,:,iscf-1)) , F_ind , dmu  )
   
      mu_ind = mu_h(:,:,iscf-1) + dmu   !!! find_min
   endif

    do ia = 1 , natm
      do alpha = 1 , 3
        it = itype( ia)
        if ( .not. ldip_polar( it ) ) cycle
        Efield (alpha,ia) = - F_ind(alpha,ia) + mu_ind(alpha,ia)/poldipia(alpha,alpha,ia)    !!! for both 
      enddo
    enddo

    ! ASPC corrector 
    if ( iscf.eq.1 .and. calc .eq. 'md' ) then
      if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield  , 2 ) ! corrector
    endif


    CALL get_rmsd_mu ( rmsd , Efield , mu_ind )

    ! output
    if ( calc .ne. 'opt' ) then
#ifdef debug
    io_printnode WRITE ( stdout ,'(a,i4,5(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' u_coul (qq)  = ', u_coul_stat, ' u_coul (dd)  = ', u_coul_ind,' u_coul_pol = ', u_coul_pol, ' rmsd       = ', rmsd
#endif
    io_printnode WRITE ( stdout ,'(a,i4,2(a,e16.8))') &
    'scf = ',iscf,' u_pol = ',u_pol * coul_unit , ' rmsd       = ', rmsd
    endif
 
#ifdef debug_print_dipff_scf
  do ia = 1 , natm
  if ( mu_ind ( 1 , ia ) .eq. 0._dp ) cycle
  write(kkkk,'(3e16.8)') (mu_ind ( alpha , ia ),alpha=1,3)
  enddo
  kkkk = kkkk + 1 
#endif



  enddo ! end of SCF loop

  ! ===========================
  !  charge/force info is recovered
  ! ===========================
  alphaES = alphaES_save
  qch = qch_tmp
  qia = qia_tmp

  if ( ioprintnode .and.  calc .ne. 'opt' ) then
    blankline(stdout)
    WRITE ( stdout , '(a,i6,a)')            'scf calculation of the induced electric moment converged in ',iscf, ' iterations '
    WRITE ( stdout , '(a,e10.3,a,e10.3,a)') 'Electric field is converged at ',rmsd,' ( ',conv_tol_ind,' ) '
    blankline(stdout)
  endif
#ifdef debug_mu
  if ( ionode ) then
    WRITE ( stdout , '(a)' )     'Induced dipoles at atoms : '
    do ia = 1 , natm 
      WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
      ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
    enddo
    blankline(stdout)
  endif
#endif
  ! store induced dipole at t 
  dipia_ind_t(1,:,:) = mu_ind
  
#ifdef debug_extrapolate
    if ( ioprintnode ) then
      WRITE(stdout,'(a,i5)') 'previous step',itime
      do ia=1,natm
        write(stdout,'(i5,<extrapolate_order+2>e16.8)') ia, mu_ind ( 1 , ia ), (dipia_ind_t ( t , 1 , ia ) , t=1, extrapolate_order+1 )
      enddo
    endif
#endif

#ifdef MPI
  stotime
  addtime(time_moment_from_pola)
#endif

  deallocate ( F_h , mu_h , F_ind , dmu ) 
  ! temporary
  theta_ind=0.0_dp

  return

END SUBROUTINE moment_from_pola_scf_kO_v4


SUBROUTINE extrapolate_dipole_poly ( mu_ind )

  USE config,           ONLY :  natm
  USE md,               ONLY :  itime 
  USE coulomb,          ONLY :  dipia_ind_t
  USE pim,              ONLY :  extrapolate_order,ldip_polar

  implicit none

  ! global  
  real(kind=dp) , intent (out):: mu_ind ( 3 , natm ) 

  ! integer 
  real(kind=dp) :: err(3)
  integer :: xa (extrapolate_order) 
  integer :: ia , ie , t

  ! ====================================
  !   zero-th order : y(t+1) = y(t)
  ! ====================================
  if ( extrapolate_order .eq. 0 ) then
      mu_ind     ( : , : )     = dipia_ind_t ( 1 , : , : )
  endif

  ! =========================================
  !   first order : y(t+1) = 2 y(t) - y(t-1)
  ! =========================================
  if ( extrapolate_order .eq. 1 .and. itime .gt. extrapolate_order ) then
    mu_ind ( : , : ) = 2.0_dp * dipia_ind_t ( 1 , : , : ) - dipia_ind_t ( 2 , :, : ) 
    dipia_ind_t( 2 , : , : ) = dipia_ind_t ( 1 , : , : )
  else if ( extrapolate_order .eq. 1 .and. itime .le. extrapolate_order ) then
    print*,'zero order because itime <= extrapolate_order'
    ! first point : zeroth order
    mu_ind     ( : , : )     = dipia_ind_t ( 1 , : , : )
    dipia_ind_t( 2 , : , : ) = dipia_ind_t ( 1 , : , : )
    print*,'store mu(t-h)=mu(t)'
  endif


  ! =========================================
  !   >=2 order : call polint
  ! =========================================
  if (extrapolate_order .ge. 2 .and. itime .gt. extrapolate_order ) then 

    do t = 1 , extrapolate_order+1 
      xa(t) = 1 - t  
    enddo

#ifdef debug_extrapolate
!      write(stdout,'(a)') ' error in extrapolated dipole'
#endif
    do ia = 1, natm 
      CALL polint(xa , dipia_ind_t( : , 1 , ia ) , extrapolate_order+1 , mu_ind ( 1 , ia ) , err(1) )  
      CALL polint(xa , dipia_ind_t( : , 2 , ia ) , extrapolate_order+1 , mu_ind ( 2 , ia ) , err(2) )  
      CALL polint(xa , dipia_ind_t( : , 3 , ia ) , extrapolate_order+1 , mu_ind ( 3 , ia ) , err(3) )  
#ifdef debug_extrapolate
!      write(stdout,'(a,3e16.8)') ia, err(1), err(2) , err(3)
#endif
    enddo
  
    do ie = extrapolate_order+1, 1, -1
      dipia_ind_t(ie,:,:) = dipia_ind_t(ie-1,:,:)
    enddo

  else if ( extrapolate_order .ge. 2 .and. itime .le. extrapolate_order ) then

                                   
      if ( itime .eq. 1 ) mu_ind ( : , : ) =          dipia_ind_t ( 1 , : , : )                            ! y(t+1) = y(t) 
      if ( itime .eq. 2 ) mu_ind ( : , : ) = 2.0_dp * dipia_ind_t ( 1 , : , : ) - dipia_ind_t ( 2 , :, : ) ! y(t+1) = 2 y(t) - y(t-1) 
      if ( itime .gt. 2 ) then
        do t = 1 , itime+1 
          xa(t) = 1 - t  
        enddo
        do ia = 1, natm 
          CALL polint(xa , dipia_ind_t( : , 1 , ia ) , itime+1 , mu_ind ( 1 , ia ) , err(1) )  
          CALL polint(xa , dipia_ind_t( : , 2 , ia ) , itime+1 , mu_ind ( 2 , ia ) , err(2) )  
          CALL polint(xa , dipia_ind_t( : , 3 , ia ) , itime+1 , mu_ind ( 3 , ia ) , err(3) )  
        enddo
      endif

    do ie = extrapolate_order+1, 1, -1
      dipia_ind_t(ie,:,:) = dipia_ind_t(ie-1,:,:)
    enddo
  endif

  return 

END SUBROUTINE extrapolate_dipole_poly

! order of extrapolation is k+1 of Kolafa original derivation.
! as the zero order is taking juste 
SUBROUTINE extrapolate_dipole_aspc ( mu_ind , Efield , key ) 

  USE config,           ONLY :  atype, natm, invpoldipia
  USE pim,              ONLY :  extrapolate_order, ldip_polar
  USE coulomb,          ONLY :  dipia_ind_t
  USE md,               ONLY :  itime 

  implicit none
  ! global  
  real(kind=dp) , intent (out):: mu_ind (:,:) 
  real(kind=dp) , intent (in) :: Efield (:,:) 
  integer :: key 
  ! local
  integer :: ia, k ,ext_ord, alpha, beta
  real(kind=dp)               :: mu_p ( 3 , natm ) 
  real(kind=dp)               :: mu_p_save (3 , natm ) 
  ! ASPC coefficient
  real(kind=dp) :: B_ASPC(6), w_ASPC

  ! switch to lower order of extrapolation 
  ! if time step is not large enough to store previous dipoles steps
  if ( itime .ge. extrapolate_order+1 ) then
    ext_ord = extrapolate_order
  else
    if (itime.eq.1.or.itime.eq.0)   then
      ext_ord = 0 
    else
      ext_ord = itime - 1
    endif
  endif

!  write(*,'(a,2i,3f16.8)') '(1) key : ',key,ext_ord,mu_ind(1,:)

  SELECT CASE ( ext_ord ) 
  CASE DEFAULT
    io_node WRITE(stdout,'(a)') 'value of aspc extrapolation order not available should be 0<= extrapolate_order <= 4'
    STOP
  CASE(0)
    B_ASPC(1) =  1.0_dp 
    W_ASPC    =  0.0_dp
  CASE(1) 
    B_ASPC(1) =  2.0_dp
    B_ASPC(2) = -1.0_dp
    W_ASPC    =  2.0_dp/3.0_dp
  CASE(2)
    B_ASPC(1) =  2.5_dp
    B_ASPC(2) = -2.0_dp
    B_ASPC(3) =  0.5_dp
    W_ASPC    =  0.6_dp
  CASE(3)
    B_ASPC(1) =  2.8_dp
    B_ASPC(2) = -2.8_dp
    B_ASPC(3) =  1.2_dp
    B_ASPC(4) = -0.2_dp
    W_ASPC    =  4.0_dp/7.0_dp 
  CASE(4)
    B_ASPC(1) =  3.0_dp
    B_ASPC(2) = -24.0_dp/7.0_dp
    B_ASPC(3) =  27.0_dp/14.0_dp
    B_ASPC(4) = -4.0_dp/7.0_dp
    B_ASPC(5) =  1.0_dp/14.0_dp
    W_ASPC    =  5.0_dp/9.0_dp
  CASE(5)
    B_ASPC(1) =  22.0_dp/7.0_dp
    B_ASPC(2) = -55.0_dp/14.0_dp
    B_ASPC(3) =  55.0_dp/21.0_dp
    B_ASPC(4) = -22.0_dp/21.0_dp
    B_ASPC(5) =  5.0_dp/21.0_dp
    B_ASPC(6) = -1.0_dp/42.0_dp
    W_ASPC    =  6.0_dp/11.0_dp
  END SELECT
        
  ! predictor
  if ( key .eq. 1 ) then 
    mu_p = 0.0_dp
    do k=1 , ext_ord+1
      mu_p = mu_p + B_ASPC( k  ) * dipia_ind_t ( k , : , : )
    enddo
    mu_ind=mu_p
    do k = ext_ord + 1, 2, -1
      dipia_ind_t(k,:,:) = dipia_ind_t(k-1,:,:)
    enddo
!    write(*,'(a,2i,3f16.8)') '(2) key : ',key,ext_ord,mu_ind(1,:)
  endif

  ! corrector
  if ( key.eq.2) then 
    mu_p_save = mu_ind
    CALL induced_moment ( Efield , mu_ind ) 
    mu_ind = w_ASPC * mu_ind + ( 1.0_dp - w_ASPC ) * mu_p_save 
  endif

!#ifdef debug_mu
!      if ( ionode ) then
!        WRITE ( stdout , '(a)' )     'Induced dipoles at atoms from extrapolation: '
!        do ia = 1 , natm 
!          WRITE ( stdout , '(i5,a3,a,3f18.10)' ) &
!          ia,atype(ia),' mu_ind = ', mu_ind ( 1 , ia ) , mu_ind ( 2 , ia ) , mu_ind ( 3 , ia )
!        enddo
!        blankline(stdout)
!      endif
!#endif

  return

END SUBROUTINE extrapolate_dipole_aspc


SUBROUTINE polint(xa,ya,n,y,dy)

  implicit none
  INTEGER       :: n
  integer       :: xa(n)
  REAL(kind=dp) :: dy,y,ya(n)
  integer, PARAMETER :: NMAX=10
  INTEGER       :: i,m,ns
  REAL(kind=dp) ::  den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(1-xa(1))
  do i=1,n
    dift=abs(1-xa(i))
    if (dift.lt.dif) then
      ns=i
      dif=dift
    endif
  c = ya
  d = ya
  enddo
  !ns = n

  y=ya(ns)
  ns=ns-1
  do m=1,n-1
    do i=1,n-m
      ho=REAL(xa(i)-1,kind=dp)
      hp=REAL(xa(i+1)-1,kind=dp)
      w=c(i+1)-d(i)
      den=ho-hp
      if(den.eq.0.) stop 'failure in polint'
      den=w/den
      d(i)=hp*den
      c(i)=ho*den
    enddo
    if (2*ns.lt.n-m)then
      dy=c(ns+1)
    else
      dy=d(ns)
      ns=ns-1
    endif
    y=y+dy
  enddo

  return

END SUBROUTINE


SUBROUTINE find_min ( g0 , g1 , d1 , gmin , dmin )

  USE config,           ONLY :  natm , itype
  USE pim,              ONLY :  ldip_polar

  implicit none
  ! global 
  real(kind=dp) :: g0(3,natm)
  real(kind=dp) :: g1(3,natm)
  real(kind=dp) :: d1(3,natm)
  real(kind=dp) :: gmin(3,natm)
  real(kind=dp) :: dmin(3,natm)
  ! local 
  real(kind=dp) :: alpha
  real(kind=dp) :: g0sq, g0g1mg0
  real(kind=dp) :: g1mg0(3,natm)
  integer       :: ia, it, k 

  ! ===========================================
  !      alpha = - g0^2 / ( g0 ( g0 - g1 ) ) 
  ! ===========================================
  g1mg0   = g1 - g0
  g0sq    = 0.0_dp 
  g0g1mg0 = 0.0_dp
  do ia=1 , natm  
    it = itype( ia)
    if ( .not. ldip_polar( it ) ) cycle
    do k = 1 ,3  
      g0sq    = g0sq    + g0(k,ia)*g0   (k,ia)
      g0g1mg0 = g0g1mg0 + g0(k,ia)*g1mg0(k,ia)
    enddo
  enddo
  alpha = - g0sq / g0g1mg0 
  if ( abs(g0g1mg0) .lt. 1e-10) alpha = 1.0d0 

  gmin = g0 + alpha * g1mg0
  dmin =      alpha * d1

  return

END SUBROUTINE

SUBROUTINE find_min_ex ( F_ind , mu_ind , F_h , mu_h , iscf , Estat )

  USE config,           ONLY : natm

  implicit none

  ! global
  real(kind=dp) , intent (inout) :: F_ind (:,:) , mu_ind (:,:), Estat (:,:)
  real(kind=dp) , intent (inout) :: F_h  (:,:,:) , mu_h (:,:,:)
  integer       , intent ( in )  :: iscf

  ! local 
  integer :: i , j 
  real(kind=dp) , dimension ( :, : ) , allocatable :: A 
  real(kind=dp) , dimension ( : ) , allocatable    :: b, tmp1,tmp2,tmp3
  integer :: ipiv(iscf)
  INTEGER :: info
  REAL(KIND=DP), external :: ddot

  F_h (:,:,iscf) = F_ind
  mu_h(:,:,iscf) = mu_ind

  ! quick return 
  if ( iscf .eq. 1 ) return 
  allocate ( A ( iscf , iscf ) )      
  allocate ( b ( iscf ) )      
  allocate ( tmp1 ( natm ) , tmp2(natm) , tmp3 (natm ) ) 
  

  A = 0.0_dp
  b = 0.0_dp

  do j = 1 , iscf   
     do i = j , iscf
          A(i,j) =          DDOTF( F_h(:,:,i) , F_h(:,:,j) ) 
          A(j,i) = A(i,j)  
    enddo
  enddo

  do i = 1 , iscf
   b(i) =        DDOTF ( Estat(:,:) , F_h(:,:,i)  )
  enddo

  CALL DGESV( iscf , 1 , A , iscf, IPIV, b, iscf, INFO )

  if ( info .eq. 0 ) continue
  if ( info .gt. 0 ) print*,'DGESV trouble'

  tmp1=0.0_dp 
  tmp2=0.0_dp 
  tmp3=0.0_dp
  mu_ind = 0.0_dp 
  do i=1,iscf
    tmp1(:) = tmp1(:)  + b(i) * F_h(1,:,i)
    tmp2(:) = tmp2(:)  + b(i) * F_h(2,:,i)
    tmp3(:) = tmp3(:)  + b(i) * F_h(3,:,i)
    mu_ind(1,:)     = mu_ind(1,:) + b(i) * mu_h(1,:,i)
    mu_ind(2,:)     = mu_ind(2,:) + b(i) * mu_h(2,:,i)
    mu_ind(3,:)     = mu_ind(3,:) + b(i) * mu_h(3,:,i)
  enddo

  ! G = P - Eq
  F_ind(1,:) = - Estat(1,:) + tmp1 
  F_ind(2,:) = - Estat(2,:) + tmp2
  F_ind(3,:) = - Estat(3,:) + tmp3
 
  deallocate ( A )      
  deallocate ( b )      
  deallocate ( tmp1 , tmp2 , tmp3 ) 

  return 

END SUBROUTINE find_min_ex

REAL(KIND=DP) FUNCTION DDOTF(v1,v2)

  USE config,           ONLY :  natm, itype
  USE pim,              ONLY :  ldip_polar

  implicit none
  real(kind=dp) :: v1(3,natm)
  real(kind=dp) :: v2(3,natm)
  integer       :: ia, it, k
  REAL(KIND=DP), external :: ddot

  DDOTF=0.0_dp 
  do ia=1 , natm
    it = itype( ia)
    if ( .not. ldip_polar( it ) ) cycle
    do k = 1 ,3
      DDOTF = DDOTF + v1(k,ia)*v2(k,ia) 
    enddo
  enddo
  

END FUNCTION


END MODULE moment_from_pola 
! ===== fmV =====
