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
#define debug_ES
#define debug_ES_field_forces
#define debug_ES_energy
#define debug_ES_stress
#define debug_ES_efg
#define debug_ES_dir
! ======= Hardware =======

! *********************** MODULE field *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE coulomb 

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE kspace,                           ONLY :  kmesh 
  USE rspace,                           ONLY :  rmesh
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  USE tensors_rk,                       ONLY :  interaction_dd
  USE mpimdff

  implicit none

  integer :: cccc=0

  logical, SAVE     :: lautoES           !< auto-determination of Ewald parameter from epsw ( accuracy)


  ! ewald sum related 
  real(kind=dp)    :: epsw                             !< accuracy of the ewald sum 
  real(kind=dp)    :: alphaES                          !< Ewald sum parameter 
  integer          :: kES(3)                           !< kmax of ewald sum in reciprocal space
  TYPE ( kmesh )   :: km_coul                          !< kpoint mesh ( see kspace.f90 )
  logical          :: task_coul(6)                     !< q-q, q-d, d-d q-Q d-Q and Q-Q tasks
  real(kind=dp)    :: qch      ( ntypemax )            !< charges 
  real(kind=dp)    :: dip      ( ntypemax , 3 )        !< dipoles 
  real(kind=dp)    :: quad     ( ntypemax , 3 , 3 )    !< quadrupoles
  real(kind=dp)    :: quad_nuc ( ntypemax )            !< quadrupolar moment nucleus NMR

  ! direct sum
  integer          :: ncelldirect                      !< number of cells  in the direct summation
  TYPE ( rmesh )   :: rm_coul                          !< real space mesh ( see rspace.f90 )
  logical          :: doefield , doefg

 ! THOLE FUNCTION RELATED
  integer      , dimension (:) , allocatable :: pair_thole
  real(kind=dp), dimension (:) , allocatable :: pair_thole_distance
  logical          :: thole_functions                                !< use thole functions correction on dipole-dipole interaction
  real(kind=dp)    :: thole_param (ntypemax,ntypemax)                !< use thole functions correction on dipole-dipole interaction
  character(len=6) :: thole_function_type
  character(len=6) :: thole_function_type_allowed(4)
  data                thole_function_type_allowed / 'linear','expon1','expon2', 'gauss'/


  real(kind=dp), dimension ( : , : )     , allocatable :: mu_t         !< total dipoles at ions
  real(kind=dp), dimension ( : , : , : ) , allocatable :: theta_t      !< total quadupole at ions
  real(kind=dp), dimension ( : , : )     , allocatable :: ef_t         !< electric field vector
  real(kind=dp), dimension ( : , : , : ) , allocatable :: efg_t        !< electric field gradient tensor
  real(kind=dp), dimension ( : , : , : ) , allocatable :: dipia_ind_t  !< induced dipole at ions 

CONTAINS

! *********************** SUBROUTINE coulomb_init ********************************
!> \brief
!! coulombic force field initialisation
! ******************************************************************************
SUBROUTINE coulomb_init

  USE control,  ONLY :  longrange

  implicit none

  ! local
  character(len=132)    :: filename
  integer               :: ioerr

  namelist /coulombtag/   lautoES , &
                          epsw    , &
                          alphaES , &  
                          kES     , & 
                          doefield, &
                          doefg   

  ! ================================
  ! defaults values for field tags 
  ! ================================
  CALL coulomb_default_tag

  ! ================================
  ! read field tags
  ! ================================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
    READ ( stdin , coulombtag, iostat=ioerr)
    if ( ioerr .lt. 0 )  then
      io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : coulombtag section is absent'
      STOP
    elseif ( ioerr .gt. 0 )  then
      io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : coulombtag wrong tag',ioerr
      STOP
    endif
  CLOSE ( stdin )

  ! ================================
  ! check field tags values
  ! ================================
  CALL coulomb_check_tag

  ! ===============================================
  !  this routines generates the ewald parameters
  ! ===============================================
  if ( longrange .eq. 'ewald' ) CALL ewald_param

  ! =====================================  
  ! initialize constant parameters and arrays
  ! =====================================  
  CALL initialize_coulomb


  ! ================================
  !  print coulomb info
  ! ================================
  CALL coulomb_print_info(stdout)

  return

END SUBROUTINE coulomb_init


! *********************** SUBROUTINE coulomb_default_tag *************************
!> \brief
! ******************************************************************************
SUBROUTINE coulomb_default_tag

  implicit none

  ! =================
  !  default values
  ! =================

  ! direct convergence
  ncelldirect   =  2
  ! ewald convergence
  kES(1)        = 10
  kES(2)        = 10
  kES(3)        = 10
  alphaES       = 1.0_dp
  epsw          = 1e-6
  lautoES       = .false.

  ! field
  doefield      = .false. ! calculate electric field ( it is internally swicth on for induced polarization calculation )
  doefg         = .false. ! electric field gradient

  task_coul = .false.

  return

END SUBROUTINE coulomb_default_tag


! *********************** SUBROUTINE coulomb_check_tag ***************************
!> \brief
!! check coulomb tag values
! ******************************************************************************
SUBROUTINE coulomb_check_tag

  USE pim,                      ONLY :  ldip_polar, lquad_polar
  USE control,                  ONLY :  lbmhftd , lbmhft , lcoulomb
  USE config,                   ONLY :  ntype , natm , dipia
  USE tt_damp,                  ONLY :  maximum_of_TT_expansion , get_TT_damp

  implicit none

  ! local
  integer :: i, it, it2
  logical :: ldamp , ldip, lqua, lqch

  lqch = .false.
  ldip = .false.
  lqua = .false.

  ! static moments     
  if ( any(qch .ne.0.0_dp) ) lqch =.true. !charge
  if ( any(dip .ne.0.0_dp) ) ldip =.true. !dipoles     
  if ( any(quad.ne.0.0_dp) ) lqua =.true. !quadrupoles     

  ! dipole polarizabilities
  do it = 1 , ntype
    if ( ldip_polar(it) )   ldip = .true.
  enddo
  ! quadrupole polarizabilities
  do it = 1 , ntype
    if ( lquad_polar(it) )  lqua = .true.
  enddo

  ! set tasks
  if ( lqch ) task_coul(1) = .true.   ! q-q
  if ( ldip ) then
    if ( lqch ) task_coul(2) = .true. ! q-d
    task_coul(3) = .true.             ! d-d
  endif
  if ( lqua ) then
    if ( lqch ) task_coul(4) = .true. ! q-Q
    if ( ldip ) task_coul(5) = .true. ! d-Q
    task_coul(6) = .true.             ! Q-Q
  endif

  return 

END SUBROUTINE coulomb_check_tag

! *********************** SUBROUTINE coulomb_print_info **************************
!> \brief
!! print force field information to standard output
! ******************************************************************************
SUBROUTINE coulomb_print_info ( kunit )

  USE constants,        ONLY :  pi, pisq
  USE control,          ONLY :  longrange, cutlongrange
  USE config,           ONLY :  simu_cell, natmi, ntype

  implicit none

  ! global
  integer :: kunit

  ! local
  integer :: i, it, it1, it2
  real(kind=dp)      :: rcut2 , kmax2 , alpha2 , ereal , ereci(3) , ereci2(3) , qtot , qtot2 


  qtot   = 0.0_dp
  qtot2  = 0.0_dp
  do it = 1 , ntype
    qtot  = qtot  +   qch(it) * natmi ( it )
    qtot2 = qtot2 + ( qch(it) * natmi ( it ) ) * ( qch(it) * natmi ( it ) )
  enddo

  if ( ionode ) then
    separator(kunit)
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'COULOMB MODULE ... WELCOME'
    blankline(kunit)
    lseparator(kunit)
    WRITE ( kunit ,'(a)')           'coulombic interaction : '
    WRITE ( kunit ,'(a,l)')         'task : charge-charge        ', task_coul(1)
    WRITE ( kunit ,'(a,l)')         'task : charge-dipole        ', task_coul(2)
    WRITE ( kunit ,'(a,l)')         'task : dipole-dipole        ', task_coul(3)
    WRITE ( kunit ,'(a,l)')         'task : charge-quadrupole    ', task_coul(4)
    WRITE ( kunit ,'(a,l)')         'task : dipole-quadrupole    ', task_coul(5)
    WRITE ( kunit ,'(a,l)')         'task : quadrupole-quadrupole', task_coul(6)

    lseparator(kunit)
    blankline(kunit)
    if ( task_coul(1) .and. .not. task_coul(3)) then
      WRITE ( kunit ,'(a)')             '        qi qj   '
      WRITE ( kunit ,'(a)')             ' Vij = -------  '
      WRITE ( kunit ,'(a)')             '         rij    '
    endif
    if ( task_coul(1) .and. task_coul(3)) then
      WRITE ( kunit ,'(a)')             '                      -->   -->               /                  -->  -->     -->   -->   \'
      WRITE ( kunit ,'(a)')             '        qi qj      qi muj . rij          1    | -->   -->      ( mui .rij ) ( rij . muj )  |'
      WRITE ( kunit ,'(a)')             ' Vij = ------- +  ---------------- +  ------- | mui . muj - 3 ---------------------------  |'
      WRITE ( kunit ,'(a)')             '         rij          rij^3            rij^3  \                            rij^2          /'
      WRITE ( kunit ,'(a)')
      WRITE ( kunit ,'(a)')
      lseparator(kunit)


    endif
    if ( .not. task_coul(1) .and. task_coul(3)) then
      WRITE ( kunit ,'(a)')             '                /                  -->  -->     -->   -->   \'
      WRITE ( kunit ,'(a)')             '           1    | -->   -->      ( mui .rij ) ( rij . muj )  |'
      WRITE ( kunit ,'(a)')             ' Vij =  ------- | mui . muj - 3 ---------------------------  |'
      WRITE ( kunit ,'(a)')             '         rij^3  \                            rij^2          /'
    endif
    blankline(kunit)
    ! =================================
    !         Direct Summation
    ! =================================
    if ( longrange .eq. 'direct' )  then
      WRITE ( kunit ,'(a)')           'direct summation'
      WRITE ( kunit ,'(a)')           'cubic cutoff in real space'
      WRITE ( kunit ,'(a,i10)')       '-ncelldirect ... ncelldirect     = ',ncelldirect
      WRITE ( kunit ,'(a,i10)')       'total number of cells            = ',( 2 * ncelldirect + 1 ) ** 3
      WRITE ( kunit ,'(a,f10.5)')     'radial cutoff                    = ',cutlongrange
    endif
    ! =================================
    !         Ewald Summation
    ! =================================
    if ( longrange .eq. 'ewald' )  then
      if ( lautoES ) then
        WRITE ( kunit ,'(a)')         'ewald summation parameters ( automatic see field.f90 tmp construction )'
        WRITE ( kunit ,'(a,f10.5)')   'alpha                            = ',alphaES
        WRITE ( kunit ,'(a,f10.5)')   'cut-off (real)                   = ',cutlongrange
        WRITE ( kunit ,'(a,3i10)')    'kmax                             = ',(kES(i),i=1,3)
        blankline(kunit)
        WRITE ( kunit ,'(a,e12.5)')   'relative error (user defined)    : ',epsw
      else
       alpha2 = alphaES * alphaES
        do i=1,3
          kmax2 = pi * kES(i) / simu_cell%ANORM(i) / alphaES
          kmax2 = kmax2 * kmax2
          ereci(i)  = EXP ( - kmax2 ) / kmax2 * ( SQRT ( REAL(kES(i),kind=dp) ) / alphaES / simu_cell%ANORM(i) / simu_cell%ANORM(i) )
          ereci2(i) = qtot2 * alphaES / pisq / ( SQRT ( REAL(kES(i),kind=dp) * REAL(kES(i),kind=dp) * REAL( kES(i),kind=dp ) ) ) &
          * EXP ( - ( pi * REAL(kES(i),kind=dp) / alphaES / simu_cell%ANORM(i) ) ** 2 )
        enddo
        rcut2  = cutlongrange * cutlongrange
        ereal  = EXP ( - alpha2 * rcut2 ) / alpha2  / rcut2 * SQRT ( cutlongrange / 2.0_dp / simu_cell%omega )
        WRITE ( kunit ,'(a)')         'ewald summation parameters ( from input file )'
        WRITE ( kunit ,'(a,f10.5)')   'alpha                            = ',alphaES
        WRITE ( kunit ,'(a,f10.5)')   'cut-off (short range)            = ',cutlongrange
        WRITE ( kunit ,'(a,3i10)')    'kmax                             = ',(kES(i),i=1,3)
        blankline(kunit)
        WRITE ( kunit ,'(a,e12.5)')   'relative error in real space       with alphaES from input : ',ereal
        WRITE ( kunit ,'(a,3e12.5)')  'relative error in reciprocal space with alphaES from input : ',(ereci(i),i=1,3)
        WRITE ( kunit ,'(a,3e12.5)')  'relative error in reciprocal space with alphaES from input : ',(ereci2(i),i=1,3)
        blankline(kunit)
        blankline(kunit)
      endif
    endif
  endif


  return

END SUBROUTINE coulomb_print_info 

! *********************** SUBROUTINE gen_ewald_param ***************************
!> \brief
!! automatic determination of ewald parameter from epsw 
!> \note
!! there is several methods ( we follow dl_poly )
! ******************************************************************************
SUBROUTINE ewald_param
   
  USE constants,                ONLY :  pi , pisq
  USE config,                   ONLY :  simu_cell 
  USE control,                  ONLY :  cutlongrange

  implicit none

  !local
  real(kind=dp) :: aaa , aaa2 , rcut 
  real(kind=dp) :: alpha , eps , tol , tol1
  integer       :: nc ( 3 )
  integer       :: kmax_1 , kmax_2 , kmax_3

  rcut =  0.5_dp * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC)
  !CALL estimate_alpha( aaa , epsw ,rcut )
  !CALL accur_ES_frenkel_smit( epsw , aaa2 , rcut , nc )
  ! dl_poly like
  if ( min (cutlongrange,rcut) .ne. cutlongrange ) &
  WRITE ( stdout , '(a)') 'WARNING : cutlongrange will be changed according to simu_cell%W*'
  cutlongrange = min (cutlongrange,rcut)
  eps=min(abs(epsw),0.5_dp)
  tol=sqrt(abs(log(eps*cutlongrange)))
  alpha=sqrt(abs(log(eps*cutlongrange*tol)))/cutlongrange
  tol1=sqrt(-log(eps*cutlongrange*(2.0_dp*tol*alpha)**2))
  kmax_1=nint(0.25_dp+simu_cell%ANORM(1)*alpha*tol1/pi)
  kmax_2=nint(0.25_dp+simu_cell%ANORM(2)*alpha*tol1/pi)
  kmax_3=nint(0.25_dp+simu_cell%ANORM(3)*alpha*tol1/pi)
  if ( lautoES ) then
    ! dl_poly like
    io_node write ( stdout , '(a)' ) 'automatic Ewald Sum '  
    alphaES = alpha
    kES(1)=kmax_1
    kES(2)=kmax_2
    kES(3)=kmax_3
    ! frenkel_smit like
    !write ( * , * ) 'auto ES frenkel_smit like'  
    !alphaES = aaa2
    !kES=nc
    ! qe like
    !alphaES = aaa
    !kES=nc
  endif
  

  return

END SUBROUTINE ewald_param

! *********************** SUBROUTINE initialize_coulomb ************************
!> \brief
!! this subroutine initialize the common quantities for charged particules.
!! real space and reciprocal space summation
!> \note
!! It is used by EFG and Coulombic subroutines 
!> \warning
!! EFG routines and + Coulombic forces routines has be used together 
! ******************************************************************************
SUBROUTINE initialize_coulomb

  USE pim,      ONLY  : 
  USE config,   ONLY  : natm , natmi , ntype , qia , itype
  USE control,  ONLY  : longrange
  USE rspace,   ONLY  : direct_sum_init
  USE kspace,   ONLY  : kpoint_sum_init , kpoint_sum_init_BZ

  implicit none

  ! local
  integer :: nk , ncmax , j , k

  allocate ( ef_t ( 3 , natm ) )
  allocate ( efg_t ( 3 , 3 , natm ) )
!  allocate ( dipia_ind_t ( extrapolate_order+1, 3 , natm ) )
  allocate ( mu_t ( 3 , natm ) )
  allocate ( theta_t ( 3 , 3 , natm ) )
  ef_t        = 0.0_dp
  efg_t       = 0.0_dp
  dipia_ind_t = 0.0_dp
  mu_t        = 0.0_dp

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
   rm_coul%meshlabel='rm_coul'
   ncmax = ( 2 * ncelldirect + 1 ) ** 3
   rm_coul%ncmax=ncmax
   rm_coul%ncell=ncelldirect
   allocate( rm_coul%boxxyz( 3 , ncmax ) , rm_coul%lcell( ncmax ) , rm_coul%rr ( ncmax ) )
   CALL direct_sum_init ( rm_coul )
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then
    km_coul%meshlabel='km_coul'
    km_coul%kmax(1) = kES(1)
    km_coul%kmax(2) = kES(2)
    km_coul%kmax(3) = kES(3)

    ! full kpt
!    nk = ( 2 * km_coul%kmax(1) + 1 ) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )   

!    half 
!    nk  = ( km_coul%kmax(1) + 1 )  * ( km_coul%kmax(2) + 1 ) * ( km_coul%kmax(3) + 1 )
    !nk = nk - 1

!   with symmetry
    nk = km_coul%kmax(3) + km_coul%kmax(2) * ( 2 * km_coul%kmax(3) + 1 ) + km_coul%kmax(1) * ( 2 * km_coul%kmax(2) + 1 ) * ( 2 * km_coul%kmax(3) + 1 )
    km_coul%nk = nk
    allocate ( km_coul%kptk( nk ) , km_coul%kptx(nk), km_coul%kpty(nk), km_coul%kptz(nk) )
    allocate ( km_coul%Ak      ( nk ) )
    allocate ( km_coul%kcoe    ( nk ) )
!    allocate ( km_coul%ckr    ( natm , nk ) )
!    allocate ( km_coul%skr    ( natm , nk ) )
!    allocate ( km_coul%rhon_R    ( nk ) )
!    allocate ( km_coul%rhon_I    ( nk ) )

!    allocate ( km_coul%rhon    ( nk ) )
!    allocate ( km_coul%expikr  ( natm , nk ) )
!    allocate ( km_coul%expikm  ( natm , nk ) )
    CALL kpoint_sum_init ( km_coul , alphaES )
  endif


  return

END SUBROUTINE initialize_coulomb

! *********************** SUBROUTINE finalize_coulomb **************************
!> \brief
!! Deallocate main quanties used during coulombic calculation
! ******************************************************************************
SUBROUTINE finalize_coulomb

  USE control,  ONLY :  longrange , lcoulomb , calc

  implicit none

  if ( .not. lcoulomb .or. calc.ne.'md') return

  deallocate ( ef_t )
  deallocate ( efg_t )
  deallocate ( dipia_ind_t )
  deallocate ( mu_t )
  deallocate ( theta_t )

  ! ============
  !  direct sum
  ! ============
  if ( longrange .eq. 'direct' ) then
    deallocate( rm_coul%boxxyz , rm_coul%lcell , rm_coul%rr )
  endif

  ! ============
  !  ewald sum
  ! ============
  if ( longrange .eq. 'ewald')  then
    deallocate( km_coul%kptk , km_coul%kptx, km_coul%kpty,km_coul%kptz )
    deallocate ( km_coul%Ak    )
    deallocate ( km_coul%kcoe  )
 !   deallocate ( km_coul%rhon_R    )
 !   deallocate ( km_coul%rhon_I    )
 !   deallocate ( km_coul%ckr    )
 !   deallocate ( km_coul%skr    )
  endif


END SUBROUTINE finalize_coulomb

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

  USE control,          ONLY :  lsurf
  USE constants,        ONLY :  tpi , piroot, coul_unit, press_unit
  USE config,           ONLY :  natm, qia, rx ,ry ,rz, simu_cell, fx, fy, fz, tau, tau_coul, atype 
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


SUBROUTINE multipole_ES_dir ( u_dir , ef_dir , efg_dir , fx_dir , fy_dir , fz_dir , tau_dir , mu , theta , & 
                              task , damp_ind , do_efield , do_efg , do_forces , do_stress )

  USE pim,                      ONLY :  ldip_damping, pol_damp_b, pol_damp_c, pol_damp_k
  USE control,                  ONLY :  lvnlist
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
                  io_print write(stdout,'(3a,2i)') 'thole function for',atype(ia),atype(ja),ia,ja
                  io_print write(stdout,'(a,f,a,f)') 's =', sthole,'  r =',d
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

SUBROUTINE multipole_ES_rec ( u_rec , ef_rec, efg_rec , fx_rec , fy_rec , fz_rec , tau_rec , mu , theta , task , &
                              do_efield , do_efg , do_forces , do_stress , do_strucfact , use_ckrskr )

  USE constants,                ONLY :  imag, tpi
  USE config,                   ONLY :  natm, rx ,ry, rz, qia, simu_cell
  USE kspace,                   ONLY :  charge_density_k_q , charge_density_k_mu , struc_fact
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


END MODULE coulomb 
! ===== fmV =====

