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
  real(kind=dp)     :: utail               !< long-range correction (energy) of short-range interaction 
  real(kind=dp)     :: ptail               !< long-range correction (virial) of short-range interaction 
  character(len=60) :: ctrunc                                                     !< truncation of nmlj
  character(len=60) :: ctrunc_allowed(3)                                          !< truncation of nmlj 
  data                 ctrunc_allowed / 'notrunc', 'linear' , 'quadratic' /       !< see initialize_param_nmlj
  integer           :: trunc                                                      !< integer definition of truncation 

  logical, SAVE     :: lKA               !< use Kob-Andersen model for BMLJ                        
  logical, SAVE     :: lautoES           !< auto-determination of Ewald parameter from epsw ( accuracy)
  logical, SAVE     :: lwrite_dip_wfc    !< write dipoles from wannier centers to file
  logical, SAVE     :: lwrite_dip        !< write dipoles 
  logical, SAVE     :: lwrite_quad       !< write quadrupoles to QUADFF
  logical, SAVE     :: lwrite_efg        !< write electric field gradient to EFGALL
  logical, SAVE     :: lwrite_ef         !< write electric field s to EFALL
  logical, SAVE     :: ldip_wfc          !< calculate electrostatic contribution from dipolar momemt coming from wfc
  logical, SAVE     :: lquiet            !< internal stuff 
  logical, SAVE     :: symmetric_pot     !< symmetric potential ( default .true. but who knows ?)
  integer, dimension (:) , allocatable :: pair_thole
  real(kind=dp), dimension (:) , allocatable :: pair_thole_distance


  ! ============================================================  
  !                         Lennard - Jones
  ! ============================================================  
  !
  !              eps    /    / sigma*\ q         / sigma*\ p  \
  !     V  =   ------- |  p | ------- |   -  q  | ------- |    |      sigma* = 2^(1/6)*sigma
  !             q - p   \    \   r   /           \   r   /    /
  !
  ! main parameters
  real(kind=dp) :: qlj     ( ntypemax , ntypemax )
  real(kind=dp) :: plj     ( ntypemax , ntypemax )
  real(kind=dp) :: epslj   ( ntypemax , ntypemax )
  real(kind=dp) :: sigmalj ( ntypemax , ntypemax )
  real(kind=dp) :: rcutsq  ( ntypemax , ntypemax )  
  real(kind=dp) :: sigsq   ( ntypemax , ntypemax ) 
  real(kind=dp) :: epsp    ( ntypemax , ntypemax )
  real(kind=dp) :: fc      ( ntypemax , ntypemax )
  real(kind=dp) :: uc      ( ntypemax , ntypemax )
  real(kind=dp) :: uc1     ( ntypemax , ntypemax )
  real(kind=dp) :: uc2     ( ntypemax , ntypemax )
  real(kind=dp) :: testtab ( ntypemax , ntypemax )

  ! ============================================================  
  !                            Morse
  ! ============================================================  
  !
  !     V  =  
  !
  !
  ! main parameters
  real(kind=dp) :: rhomor  ( ntypemax , ntypemax )
  real(kind=dp) :: epsmor  ( ntypemax , ntypemax )
  real(kind=dp) :: sigmamor( ntypemax , ntypemax )
  real(kind=dp) :: rs      ( ntypemax , ntypemax )
  real(kind=dp) :: fm      ( ntypemax , ntypemax )

  ! ============================================================  
  !                            BMHFTD
  ! ============================================================  
  !
  !     V  =   A*exp(-B*r) - f_6*(r) C      f_8(r)*D      f_order(r)=1-exp(-BD * r) * \sum_{k=0}^order (BD * r)^k / k! 
  !                         -----------  - ----------
  !                            r ^ 6         r ^ 8
  !
  ! main parameters
  real(kind=dp) :: Abmhftd  ( ntypemax , ntypemax )
  real(kind=dp) :: Bbmhftd  ( ntypemax , ntypemax )
  real(kind=dp) :: Cbmhftd  ( ntypemax , ntypemax )
  real(kind=dp) :: Dbmhftd  ( ntypemax , ntypemax )
  real(kind=dp) :: BDbmhftd ( ntypemax , ntypemax )


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

  ! THOLE FUNCTION RELATED
  logical          :: thole_functions                                !< use thole functions correction on dipole-dipole interaction
  real(kind=dp)    :: thole_param (ntypemax,ntypemax)                !< use thole functions correction on dipole-dipole interaction
  character(len=6) :: thole_function_type
  character(len=6) :: thole_function_type_allowed(4)
  data                thole_function_type_allowed / 'linear','expon1','expon2', 'gauss'/


  character(len=4) :: algo_ext_dipole                  !< set the algorithm used to get induced moments from polarization 
  character(len=4) :: algo_ext_dipole_allowed(2)
  data                algo_ext_dipole_allowed       / 'poly','aspc'/

  character(len=11) :: algo_moment_from_pola            !< set the algorithm used to get induced moments from polarization 
  character(len=11) :: algo_moment_from_pola_allowed(6) !< set the algorithm used to get induced moments from polarization 
  data                 algo_moment_from_pola_allowed / 'scf' , 'scf_kO_v1' , 'scf_kO_v2' , 'scf_kO_v3' , 'scf_kO_v4_1' , 'scf_kO_v4_2'/  !! scf ( self consistent ) 
  
  integer          :: lwfc     ( ntypemax )            !< moment from wannier centers 
  real(kind=dp)    :: rcut_wfc                         !< radius cut-off for WFs searching

  ! ewald sum related 
  real(kind=dp)    :: epsw                             !< accuracy of the ewald sum 
  real(kind=dp)    :: alphaES                          !< Ewald sum parameter 
  integer          :: kES(3)                           !< kmax of ewald sum in reciprocal space
  TYPE ( kmesh )   :: km_coul                          !< kpoint mesh ( see kspace.f90 )
  logical          :: task_coul(6)                     !< q-q, q-d, d-d q-Q d-Q and Q-Q tasks

  ! direct sum
  integer          :: ncelldirect                      !< number of cells  in the direct summation
  TYPE ( rmesh )   :: rm_coul                          !< real space mesh ( see rspace.f90 )
  logical          :: doefield , doefg

  real(kind=dp), dimension ( : , : )     , allocatable :: mu_t         !< total dipoles at ions
  real(kind=dp), dimension ( : , : , : ) , allocatable :: theta_t      !< total quadupole at ions
  real(kind=dp), dimension ( : , : )     , allocatable :: ef_t         !< electric field vector
  real(kind=dp), dimension ( : , : , : ) , allocatable :: efg_t        !< electric field gradient tensor
  real(kind=dp), dimension ( : , : , : ) , allocatable :: dipia_ind_t  !< induced dipole at ions 

CONTAINS

! *********************** SUBROUTINE field_default_tag *************************
!> \brief
!! set default values to field tags
! ******************************************************************************
SUBROUTINE field_default_tag

  implicit none

  ! =================
  !  default values
  ! =================

  ctrunc        = 'notrunc' 
  symmetric_pot = .true.

  ! LJ
  lKA           = .false.
  epslj         = 0.0_dp
  sigmalj       = 0.0_dp
  qlj           = 12.0_dp
  plj           = 6.0_dp

  ! morse
  epsmor        = 0.0_dp 
  sigmamor      = 0.0_dp 
  rhomor        = 0.0_dp 

  ! bmhftd
  Abmhftd  = 0.0_dp 
  Bbmhftd  = 0.0_dp
  Cbmhftd  = 0.0_dp
  Dbmhftd  = 0.0_dp
  BDbmhftd = 0.0_dp

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
  qch           = 0.0_dp  ! charge
  quad_nuc      = 0.0_dp  ! quadrupolar moment
  dip           = 0.0_dp  ! dipolar moment
  doefield      = .false. ! calculate electric field ( it is internally swicth on for induced polarization calculation )
  doefg         = .false. ! electric field gradient
  lwrite_dip    = .false.            
  lwrite_quad   = .false.            
  lwrite_ef     = .false.            
  lwrite_efg    = .false.            

  ! polarization
  ldip_polar        = .false. 
  poldip        = 0.0_dp
  polquad       = 0.0_dp
  pol_damp_b    = 0.0_dp
  pol_damp_c    = 0.0_dp
  pol_damp_k    = 0
  ldip_damping  = .false.
  conv_tol_ind  = 1e-6
  min_scf_pol_iter = 3
  max_scf_pol_iter = 100
  extrapolate_order = 0 
  algo_ext_dipole = 'aspc'
  algo_moment_from_pola = 'scf'
  thole_functions       = .false.
  thole_function_type   = 'linear'
  thole_param           = 1.662_dp

  omegakO               = 0.7_dp 

  ! wannier centers related
  lwfc           = 0
  lwrite_dip_wfc = .false.            
  ldip_wfc       = .true.            
  rcut_wfc       = 0.5_dp

  mass           = 1.0_dp

  task_coul = .false.

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
  integer :: i, it, it2
  logical :: ldamp , ldip, lqua, lqch

  ! ========
  !  ctrunc
  ! ========
  CALL check_allowed_tags ( size ( ctrunc_allowed ), ctrunc_allowed, ctrunc, 'in fieldtag','ctrunc' ) 

  if ( ctrunc .eq. 'notrunc'   ) trunc = 0
  if ( ctrunc .eq. 'linear'    ) trunc = 1
  if ( ctrunc .eq. 'quadratic' ) trunc = 2

  ! =================================================
  !  KOB-ANDERSEN MODEL --- PhysRevE 51-4626 (1995) 
  ! =================================================
  if ( lKA .and. ntype .ne. 2 ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR fieldtag lKA should be used with 2 differents types'
    STOP 
  endif
  if ( lKA ) then
    sigmalj ( 1 , 1 ) = 1.0_dp
    sigmalj ( 2 , 2 ) = 0.88_dp
    sigmalj ( 1 , 2 ) = 0.8_dp
    sigmalj ( 2 , 1 ) = 0.8_dp
    epslj   ( 1 , 1 ) = 1.0_dp
    epslj   ( 2 , 2 ) = 0.5_dp
    epslj   ( 1 , 2 ) = 1.5_dp
    epslj   ( 2 , 1 ) = 1.5_dp
  endif

  ! =================================================
  ! symetrization of input potentials !
  ! =================================================
  if ( symmetric_pot ) then
    do it = 1 , ntype
      ! nmlj
      epslj   (:,it)   = epslj   (it,:) 
      sigmalj (:,it)   = sigmalj (it,:)
      plj     (:,it)   = plj     (it,:)
      qlj     (:,it)   = qlj     (it,:)
      ! bhmftd
      Abmhftd(:,it)  = Abmhftd(it,:)
      Bbmhftd(:,it)  = Bbmhftd(it,:)
      Cbmhftd(:,it)  = Cbmhftd(it,:)
      Dbmhftd(:,it)  = Dbmhftd(it,:)
      BDbmhftd(:,it) = BDbmhftd(it,:)
      !thole param
      thole_param(:,it) = thole_param(it,:)  
      do it2 = 1 ,ntype
        ! dip_damping
        ldip_damping(it2,:,it) = ldip_damping(it2,it,:)
        pol_damp_b  (it2,:,it) = pol_damp_b  (it2,it,:)
        pol_damp_c  (it2,:,it) = pol_damp_c  (it2,it,:)
        pol_damp_k  (it2,:,it) = pol_damp_k  (it2,it,:)
      enddo
    enddo
  endif

  if ( any( pol_damp_k .gt. maximum_of_TT_expansion ) ) then
    io_node  WRITE ( stdout , '(a,i5)' ) 'ERROR Tang-Toennieng expansion order too large (i.e pol_damp_k) max = ',maximum_of_TT_expansion 
    STOP
  endif
  ldamp=.false.
  if ( any ( ldip_damping )  )  ldamp = .true.

  if ( ldamp .or. lbmhftd ) CALL get_TT_damp

  ! ===============================
  !       coulombic tasks
  ! ===============================
  if ( lcoulomb ) then

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

  endif


  ! =======================
  !  algo_moment_from_pola 
  ! =======================
  CALL check_allowed_tags( size( algo_moment_from_pola_allowed), algo_moment_from_pola_allowed, algo_moment_from_pola, 'in fieldtag','algo_moment_from_pola' ) 

  ! =======================
  !  thole_function_type 
  ! =======================
  CALL check_allowed_tags( size(thole_function_type_allowed), thole_function_type_allowed, thole_function_type, 'in fieldtag','thole_function_type' ) 
   
  if ( thole_function_type .eq. 'expon1' ) thole_param = 0.572_dp
  if ( thole_function_type .eq. 'expon2' ) thole_param = 1.9088_dp
  if ( thole_function_type .eq. 'gauss' ) thole_param = 0.957_dp


  ! =================================================
  ! isotropic values for polarizabilities in input 
  ! set the tensors
  ! =================================================
  do it = 1, ntype
    if ( poldip_iso(it) .ne. 0.0_dp ) then
      poldip(it,:,:) = 0.0_dp
      poldip(it,1,1) = poldip_iso(it)
      poldip(it,2,2) = poldip_iso(it)
      poldip(it,3,3) = poldip_iso(it)
    endif
    if ( polquad_iso(it) .ne. 0.0_dp ) then
      polquad(it,:,:,:) = 0.0_dp
      polquad(it,1,1,1) = polquad_iso(it)
      polquad(it,2,2,2) = polquad_iso(it)
      polquad(it,3,3,3) = polquad_iso(it)
    endif
    !write(stdout,'(a,i,4e16.8)')'debug',it,polquad(it,1,1,1), polquad(it,2,2,2),polquad(it,3,3,3) ,polquad_iso(it)
  enddo




  return 

END SUBROUTINE field_check_tag

! *********************** SUBROUTINE field_init ********************************
!> \brief
!! force field initialisation
! ******************************************************************************
SUBROUTINE field_init

  USE control,                  ONLY :  calc , lnmlj , lcoulomb , lmorse , longrange, non_bonded

  implicit none

  ! local
  character(len=132)    :: filename
  integer               :: ioerr

  namelist /fieldtag/    lKA           , &       
                         ctrunc        , &
                         symmetric_pot , &
                         ncelldirect   , &
                         kES           , &
                         alphaES       , &
                         qlj           , & 
                         plj           , & 
                         sigmalj       , &
                         epslj         , &
                         sigmamor      , &
                         epsmor        , &
                         rhomor        , &
                         mass          , &
                         doefield      , &
                         doefg         , &
                         qch           , &
                         quad_nuc      , &
                         dip           , &
                         quad          , &
                         poldip        , &  
                         poldip_iso    , &  
                         polquad       , &  
                         polquad_iso   , &  
                         pol_damp_b    , &  
                         pol_damp_c    , &  
                         pol_damp_k    , &  
                         extrapolate_order , &
                         conv_tol_ind  , &  
                         min_scf_pol_iter, &
                         max_scf_pol_iter, &
                         algo_moment_from_pola, &
                         algo_ext_dipole , &
                         thole_functions, &
                         thole_function_type, &
                         thole_param   , &
                         omegakO       , &
                         epsw          , &  
                         lautoES       , &  
                         lwfc          , &            
                         lwrite_dip_wfc, &            
                         lwrite_dip    , &            
                         lwrite_quad   , &            
                         lwrite_ef     , &            
                         lwrite_efg    , &            
                         ldip_wfc      , &            
                         rcut_wfc      , &            
                         ldip_polar    , &
                         ldip_damping  , &
                         lquad_polar   , &
                         Abmhftd       , &  
                         Bbmhftd       , &
                         Cbmhftd       , &
                         Dbmhftd       , &
                         BDbmhftd 
  
                 

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
  if ( longrange .eq. 'ewald' .and. lcoulomb ) CALL ewald_param
 
  ! =====================================  
  !  if efg print field info and return 
  ! =====================================  
  if ( calc .eq. 'efg' .or. calc .eq. 'rmc' ) then 
    CALL field_print_info(stdout,quiet=.true.)
    return
  endif

  ! ================================
  ! initialize constant parameters
  ! ================================
  if ( non_bonded )    then
    CALL initialize_param_non_bonded
  endif

  if ( lcoulomb ) then
    CALL initialize_coulomb
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

  USE config,           ONLY :  natm , ntype , atype , atypei , natmi , simu_cell , rhoN, rho , massia 
  USE control,          ONLY :  calc , cutshortrange , lnmlj , lmorse , lbmhft , lbmhftd , lcoulomb , longrange , lreducedN , cutlongrange
  USE constants,        ONLY :  pi , pisq, rho_unit

  implicit none

  !local
  logical , optional :: quiet
  integer            :: kunit, it , it1 , it2 , i , j , ia 
  real(kind=dp)      :: rcut2 , kmax2 , alpha2 , ereal , ereci(3) , ereci2(3) , qtot , qtot2, total_mass
  logical            :: linduced, ldamp
  real(kind=dp)      :: Athole, dist_cata , dist_corr, mu_sum(3) , theta_sum(3,3)

  if ( ( present ( quiet ) .and. quiet ) .and. .not. lquiet ) then 
    lquiet = .true.
  else if ( ( present ( quiet ) .and. quiet ) .and. lquiet ) then
    return
  endif

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
  linduced = .false.
  do it = 1 , ntype
    if ( ldip_polar(it) )  linduced = .true.
  enddo
  ldamp = .false.
  if ( any ( ldip_damping )  )  ldamp = .true.


  if ( ionode ) then
    total_mass = 0.0_dp
    do ia = 1 , natm 
      total_mass = total_mass + massia(ia)
    enddo
    separator(kunit)    
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'FIELD MODULE ... WELCOME'
    blankline(kunit)
    lseparator(kunit) 
    WRITE ( kunit ,'(a,f12.4,a)')       'total mass            = ', total_mass ,' a.m '
    WRITE ( kunit ,'(a,2f12.4,a)')      'density               = ', rhoN , rho * rho_unit ,' g/cm^3 '
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
    if ( linduced ) then 
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'polarizabilities on atoms'
      if ( ldamp ) &
      WRITE ( kunit ,'(a)')             'electric field damping applied to polarizable atoms' 

      lseparator(kunit) 
      do it1 = 1 , ntype
        if ( ldip_polar( it1 ) ) then
          WRITE ( kunit ,'(a,a2,a,f12.4)')'polarizability on type ', atypei(it1),' : ' 
          WRITE ( kunit ,'(3f12.4)')      ( poldip ( it1 , 1 , j ) , j = 1 , 3 ) 
          WRITE ( kunit ,'(3f12.4)')      ( poldip ( it1 , 2 , j ) , j = 1 , 3 ) 
          WRITE ( kunit ,'(3f12.4)')      ( poldip ( it1 , 3 , j ) , j = 1 , 3 ) 
          blankline(kunit)
          !if ( ldamp ) then
            WRITE ( kunit ,'(a)') 'damping functions : '
            WRITE ( kunit ,'(a)') '                    b           c              k'
            do it2 = 1 ,ntype 
              !if ( ldip_damping(it1,it1,it2 ) )  &
              WRITE ( kunit ,'(a,a,a,a,2f12.4,i2,a,2f12.4,i2,a)') atypei(it1),' - ',atypei(it2), ' : ' ,& 
                                                    pol_damp_b(it1,it1,it2),pol_damp_c(it1,it1,it2),pol_damp_k(it1,it1,it2),&
                                              ' ( ',pol_damp_b(it1,it2,it1),pol_damp_c(it1,it2,it1),pol_damp_k(it1,it2,it1),' ) '
            enddo
          !endif


        else
          WRITE ( kunit ,'(a,a2)')      'no polarizability on type ', atypei(it1)
        endif
        lseparator(kunit) 
        blankline(kunit)
      enddo
    endif
    ! =================================
    !      LONG RANGE INTERACTIONS 
    ! =================================
    if ( lcoulomb )    then 
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
        if ( thole_functions ) then
        lseparator(kunit)
        WRITE ( kunit ,'(a)')             ' Thole function for dipole-dipole interaction ' 
        WRITE ( kunit ,'(a,a)')           ' type  : ',thole_function_type
        WRITE ( kunit ,'(a)')             ' thole parameter        : a '
        WRITE ( kunit ,'(a)')             ' catastrophe distance   : dc =   ( 4 alpha_i alpha_j )^(1/*6)'
        WRITE ( kunit ,'(a)')             ' thole distance cut-off : dt = a (   alpha_i alpha_j )^(1/*6)' 
        lseparator(kunit)
        WRITE ( kunit ,'(a)')             'pair               a              dc              dt'
        do it1 = 1 , ntype
          do it2 = it1 , ntype
            Athole    = ( poldip (it1,1,1) * poldip (it2,1,1) ) ** (1.0_dp / 6.0_dp )
            dist_cata = ( 4.0_dp ) ** (1.0_dp / 6.0_dp ) * Athole
            dist_corr = Athole * thole_param ( it1 , it2 )
            WRITE ( kunit ,120)   atypei(it1),'-',atypei(it2),'    ',thole_param ( it1 , it2 ), dist_cata , dist_corr 
          enddo
        enddo
        lseparator(kunit)

        
        endif
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
    ! =================================
    !     SHORT RANGE INTERACTIONS 
    ! =================================
    !       LENNARD-JONES
    ! =================================
    if ( lnmlj )       then     
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'lennard-jones           '
      lseparator(kunit) 
      blankline(kunit)
      WRITE ( kunit ,'(a)')             '       eps    /    / sigma* \ q       / sigma*  \ p  |'
      WRITE ( kunit ,'(a)')             ' V = ------- |  p | ------- |    - q | -------- |    |'   
      WRITE ( kunit ,'(a)')             '      q - p   \    \   r    /         \    r    /    /'
      blankline(kunit)
      WRITE ( kunit ,'(a,f10.5)')       'cutoff      = ',cutshortrange
      WRITE ( kunit ,'(a,a)')           'truncation  = ',ctrunc
      if ( .not. lreducedN ) &
      WRITE ( kunit ,'(a,2f20.9)')      'long range correction (energy)   : ',utail
      WRITE ( kunit ,'(a,2f20.9)')      'long range correction (pressure) : ',ptail
      blankline(kunit)
      blankline(kunit)
      do it1 = 1 , ntype
        do it2 = it1 , ntype 
          lseparator(kunit) 
          WRITE ( kunit ,'(a,a,a,a)')   atypei(it1),'-',atypei(it2),' interactions:'    
          lseparator(kunit) 
          if ( it1 .eq. it2 ) then
            WRITE ( kunit ,100)         'sigma                                = ',sigmalj ( it1 , it2 )
            WRITE ( kunit ,100)         'eps                                  = ',epslj   ( it1 , it2 )
            WRITE ( kunit ,100)         'q                                    = ',qlj     ( it1 , it2 )
            WRITE ( kunit ,100)         'p                                    = ',plj     ( it1 , it2 )
          else
            WRITE ( kunit ,110)         'sigma                                = ',sigmalj ( it1 , it2 ) , '( ',sigmalj ( it2 , it1 ), ' )' 
            WRITE ( kunit ,110)         'eps                                  = ',epslj   ( it1 , it2 ) , '( ',epslj   ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'q                                    = ',qlj     ( it1 , it2 ) , '( ',qlj     ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'p                                    = ',plj     ( it1 , it2 ) , '( ',plj     ( it2 , it1 ), ' )'
          endif
          if ( trunc .eq. 1 ) then
            WRITE ( kunit ,'(a,f10.5)') 'shift correction (linear)            : ',uc(it1,it2)
          else if ( trunc .eq. 2 ) then
            WRITE ( kunit ,'(a,2f10.5)')'shift correction (quadratic)         : ',uc(it1,it2),uc2(it1,it2)
          endif
        enddo
      enddo
    endif
    ! =================================
    !      BUCKINGHAM - MORSE 
    ! =================================
    if ( lmorse )       then
      blankline(kunit)
      lseparator(kunit) 
      WRITE ( kunit ,'(a)')             'Morse  potential'
      lseparator(kunit) 
      blankline(kunit)
      WRITE ( kunit ,'(a)')             ' V = eps * exp ( -2 rho ( r - sigma ) ) - 2 eps * exp ( -rho ( r - sigma) ) '
      do it1 = 1 , ntype
        do it2 = it1 , ntype
          lseparator(kunit) 
          WRITE ( kunit ,'(a,a,a,a)')   atypei(it1),'-',atypei(it2),' interactions:'
          lseparator(kunit) 
          if ( it1 .eq. it2 ) then
            WRITE ( kunit ,100)         'sigma                                = ',sigmamor ( it1 , it2 )
            WRITE ( kunit ,100)         'eps                                  = ',epsmor   ( it1 , it2 )
            WRITE ( kunit ,100)         'rho                                  = ',rhomor   ( it1 , it2 )
          else
            WRITE ( kunit ,110)         'sigma                                = ',sigmamor ( it1 , it2 ) , '( ',sigmamor ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'eps                                  = ',epsmor   ( it1 , it2 ) , '( ',epsmor   ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'rho                                  = ',rhomor   ( it1 , it2 ) , '( ',rhomor   ( it2 , it1 ), ' )'
          endif
        enddo
      enddo
    endif
    if ( lbmhftd .or. lbmhft ) then
      blankline(kunit)
      lseparator(kunit)
      WRITE ( kunit ,'(a)')             'BMHFTD potential'
      WRITE ( kunit ,'(a)')             'Born-Huggins-Meyer-Fumi-Tossi + Damping'
      lseparator(kunit)
      blankline(kunit)
      WRITE ( kunit ,'(a)')             '                                  C              D  '
      WRITE ( kunit ,'(a)')             ' V = A  exp ( - B  r ) - f_6(r) ----- - f_8(r) -----'
      WRITE ( kunit ,'(a)')             '                                 r^6            r^8 '
      blankline(kunit)
      WRITE ( kunit ,'(a)')             '                             n                      '
      WRITE ( kunit ,'(a)')             '                            ----   (BD * r)^k       '
      WRITE ( kunit ,'(a)')             ' f_n(r) = 1 - exp(-BD r ) * \     ------------      '
      WRITE ( kunit ,'(a)')             '                            /___       k!           '
      WRITE ( kunit ,'(a)')             '                             k=0                    '
      blankline(kunit)
      do it1 = 1 , ntype
        do it2 = it1 , ntype
          lseparator(kunit)
          WRITE ( kunit ,'(a,a,a,a)')   atypei(it1),'-',atypei(it2),' interactions:'
          lseparator(kunit)
          if ( it1 .eq. it2 ) then
            WRITE ( kunit ,100)         'A                                = ',Abmhftd ( it1 , it2 )
            WRITE ( kunit ,100)         'B                                = ',Bbmhftd ( it1 , it2 )
            WRITE ( kunit ,100)         'C                                = ',Cbmhftd ( it1 , it2 )
            WRITE ( kunit ,100)         'D                                = ',Dbmhftd ( it1 , it2 )
            if ( lbmhftd ) &
            WRITE ( kunit ,100)         'BD                               = ',BDbmhftd ( it1 , it2 )
          else
            WRITE ( kunit ,110)         'A                                = ',Abmhftd ( it1 , it2 ) , '( ', Abmhftd ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'B                                = ',Bbmhftd ( it1 , it2 ) , '( ', Bbmhftd ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'C                                = ',Cbmhftd ( it1 , it2 ) , '( ', Cbmhftd ( it2 , it1 ), ' )'
            WRITE ( kunit ,110)         'D                                = ',Dbmhftd ( it1 , it2 ) , '( ', Dbmhftd ( it2 , it1 ), ' )'
            if ( lbmhftd ) &
            WRITE ( kunit ,110)         'BD                               = ',BDbmhftd ( it1 , it2 ), '( ', BDbmhftd ( it2 , it1 ), ' )'
          endif
        enddo
      enddo
    endif


    blankline(kunit)
    blankline(kunit)
    separator(kunit)    
    blankline(kunit)
  endif

  return
100 FORMAT(a,e16.8) 
110 FORMAT(a,e16.8,a,e16.8,a) 
120 FORMAT(a,a,a,a,3e16.8) 

END SUBROUTINE field_print_info


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

! *********************** SUBROUTINE initialize_param_non_bonded ***************
!> \brief
!! initialisation of principal parameters for lennard-jones and morse potentials
! ******************************************************************************
SUBROUTINE initialize_param_non_bonded

  USE constants,                ONLY :  tpi , press_unit
  USE config,                   ONLY :  ntypemax , natm , natmi , rhoN , atype , itype  , ntype , simu_cell
  USE control,                  ONLY :  skindiff , cutshortrange , calc

  implicit none

  ! local
  integer :: it,jt
  real(kind=dp) :: one13, one16, two16, rskinmax
  real(kind=dp) :: rcut3 ( ntypemax , ntypemax )
  real(kind=dp) :: rskin ( ntypemax , ntypemax ) 
  real(kind=dp) :: rskinsq ( ntypemax , ntypemax ) 
  real(kind=dp) :: ut ( ntypemax , ntypemax ) 
  real(kind=dp) :: pt ( ntypemax , ntypemax ) 

  real(kind=dp) :: rcut ( ntype , ntype ) 
  real(kind=dp) :: ppqq ( ntype , ntype )
  real(kind=dp) :: pp ( ntype , ntype )
  real(kind=dp) :: qq ( ntype , ntype )
  real(kind=dp) :: pp3 ( ntype , ntype ) 
  real(kind=dp) :: qq3 ( ntype , ntype ) 
  real(kind=dp) :: sr2 ( ntype , ntype ) 
  real(kind=dp) :: sr ( ntype , ntype ) 
  real(kind=dp) :: srp ( ntype , ntype ) 
  real(kind=dp) :: srq ( ntype , ntype ) 

  rskinmax = 0.0_dp
  utail    = 0.0_dp
  ptail    = 0.0_dp
  rcut3    = 0.0_dp
  rcut     = 0.0_dp
  rskin    = 0.0_dp
  rskinsq  = 0.0_dp
  ut       = 0.0_dp
  ppqq     = 0.0_dp
  pp       = 0.0_dp
  qq       = 0.0_dp
  sr2      = 0.0_dp
  sr       = 0.0_dp
  srp      = 0.0_dp
  srq      = 0.0_dp
  
  do it = 1 , ntype
    do jt = 1 , ntype
      pp ( it , jt ) = plj ( it , jt ) 
      qq ( it , jt ) = qlj ( it , jt ) 
    enddo
  enddo
  pp3 = pp - 3.0_dp
  qq3 = qq - 3.0_dp

  one13 = (1.0_dp / 3.0_dp)
  one16 = (1.0_dp / 6.0_dp)
  two16 = 2.0_dp **  one16

  ! ==================================================================================================
  ! TAIL ( checked september 2011) :
  ! be careful two minus are vanishing
  !
  !     2 pi rc^3 espilon NA NB    /     p         / sigma* \ q           q         / sigma*  \ p  \
  !    -------------------------- |  ----------   | ------- |    --   ----------   | -------- |    |
  !         ( q - p )   V          \ ( q - 3 )     \   rc   /          ( p - 3 )    \  rc     /    /
  !
  !  virial the same with qp in the first term
  !
  ! rskinmax ??????
  ! ==================================================================================================
  
  
  ! ==================================================================================================
  !  simple truncation  
  !
  !  ctrunc = 'linear'
  !  trunc = 1 
  !
  !
  !           eps    /    / sigma* \ q         / sigma* \ p  \
  !  V  =   ------- |  p | ------- |    -  q  | --------|    |   -  c1   with  sigma* = 2^(1/6)*sigma
  !          q - p   \    \   r    /           \    r   /    /
  !
  ! 
  !          eps      /     /  sigma* \ q        /  sigma* \ p  \
  !  c1 = ---------  |   p |  -------- |    - q |  -------- |   |      with rc = cutshortrange 
  !         q - p     \     \   rc    /          \   rc    /    / 
  ! 
  ! ==================================================================================================
  

  ! ==================================================================================================
  !  truncation presented in J. Chem. Phys. 135 (2011) , Sengupta, Vasconcelos, Affouard, Sastry
  ! 
  !  ctrunc = 'quadratic'
  !  trunc  = 2    
  ! 
  !
  !           eps    /    / sigma* \ q         / sigma* \ p  \
  !  V  =   ------- |  p | ------- |    -  q  | --------|    |   +  c1 r^2  -  c2    with  sigma* = 2^(1/6)*sigma
  !          q - p   \    \   r   /            \    r   /    /
  !
  !
  ! Sage (January 2013) 
  !     
  !              eps p  q           /  / sigma* \ q       /  sigma* \ p \
  !    c1 =  --------------------- |  | --------|    -   | --------- |  |
  !          2  ( q - p ) * rc^2    \  \   rc   /         \   rc    /   /
  ! 
  !    and
  !     
  !           epsilon     /           / sigma* \ q               / sigma* \ p  \
  !    c2 =  ----------- |  (2p+qp)  | --------|     - (2q+pq)  | --------|   |
  !          2 ( q - p )  \           \   rc   /                 \   rc   /   /       
  !
  ! ==================================================================================================
  do it = 1 , ntype 
    do jt = 1 , ntype

          ! intermediate terms 
          rcut   ( it , jt ) = cutshortrange                               ! rc
          rcutsq ( it , jt ) = rcut ( it , jt ) * rcut   ( it , jt )       ! rc^2 
          rcut3  ( it , jt ) = rcut ( it , jt ) * rcutsq ( it , jt )       ! rc^3
          ppqq   ( it , jt ) = pp   ( it , jt ) * qq     ( it , jt )       ! p x q
          rskin  ( it , jt ) = rcut ( it , jt ) + skindiff

         if ( rskin ( it , jt ) .gt. rskinmax ) rskinmax = rskin ( it , jt )
           rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
           ! eps / q - p   
           epsp ( it , jt )= epslj( it , jt )/( qq ( it , jt )-pp ( it , jt ) )
           ! sigma*^2 
           sigsq ( it , jt )= two16 * two16 * sigmalj ( it , jt ) * sigmalj ( it , jt )
           ! sigma*^2 / rc ^2
           sr2 ( it , jt ) = sigsq ( it , jt )/ rcutsq ( it , jt )
           ! sigma* / rc 
           sr( it ,jt )    = SQRT ( sr2( it , jt )  )
           ! (sigma* / rc ) ^ p
           srp( it , jt )  = sr ( it , jt ) ** pp ( it , jt )
           ! (sigma* / rc ) ^ q 
           srq( it , jt )  = sr ( it , jt ) ** qq ( it , jt )
           ! trunc = 1
           uc ( it , jt )  = epsp ( it , jt ) * ( pp ( it , jt ) * srq( it , jt ) - qq ( it , jt ) * srp( it , jt ) )
           ! trunc = 2
           ! c1
           uc1 ( it , jt ) = epsp ( it , jt ) *  ppqq ( it , jt ) / ( 2.0_dp * rcutsq ( it , jt ) ) * ( srq( it , jt ) - srp( it , jt ) ) 
           ! c2
           uc2 ( it , jt ) = 0.5_dp * epsp( it , jt )  * (  &
                         ( 2.0_dp * pp ( it , jt ) + ppqq ( it , jt )  ) * srq( it , jt ) - &
                         ( 2.0_dp * qq ( it , jt ) + ppqq ( it , jt )  ) * srp( it , jt ) )
           ! for the virial
           fc ( it , jt ) =  ppqq ( it , jt ) * epsp ( it , jt ) /  sigsq ( it , jt ) 
           ! morse
           fm ( it , jt ) = - 2.0_dp * epsmor ( it , jt ) * rhomor ( it , jt ) * EXP ( rhomor ( it , jt ) * sigmamor ( it , jt ) )  
           rs ( it , jt ) = EXP ( rhomor ( it , jt ) * sigmamor( it , jt ) ) 
           ! tail energy
           ut ( it , jt ) = epsp ( it , jt ) * ( pp ( it , jt ) * srq ( it , jt ) / qq3 ( it , jt ) - &
                                                 qq ( it , jt ) * srp ( it , jt ) / pp3 ( it , jt )  )       
           pt ( it , jt ) = ppqq ( it , jt ) * epsp ( it , jt ) * ( srq ( it , jt ) / qq3 ( it , jt ) - &
                                                                    srp ( it , jt ) / pp3 ( it , jt )  )       

           ut ( it , jt ) = ut ( it , jt ) * rcut3 ( it , jt ) * tpi 
           pt ( it , jt ) = pt ( it , jt ) * rcut3 ( it , jt ) * tpi

           if ( ( natmi ( it ) .ne. 0 ) .and. ( natmi ( jt ) .ne. 0 ) ) &
           utail = utail + ut ( it , jt ) * natmi ( it ) * natmi ( jt ) / simu_cell%omega
           if ( ( natmi ( it ) .ne. 0 ) .and. ( natmi ( jt ) .ne. 0 ) ) &
           ptail = ptail + pt ( it , jt ) * natmi ( it ) * natmi ( jt ) / simu_cell%omega


!#ifdef debug
!  WRITE ( stdout , '(2i6,7e16.6)' ) it , jt , uc ( it , jt )  , epsp ( it , jt ) , pp ( it , jt ) , qq ( it , jt ) , srq( it , jt ) , srp( it , jt ) ,rcutsq ( it , jt ) 
!#endif
    enddo
  enddo
           ptail = ptail /press_unit/simu_cell%omega/3.0d0 ! virial to pressure
           io_node write(stdout,'(a,e16.6)') 'long range correction (init) energy   ',utail
           io_node write(stdout,'(a,e16.6)') 'long range correction (init) pressure ',ptail
 
#ifdef debug_quadratic
  do it = 1 , ntype 
    WRITE ( stdout , '(a,<ntype>e16.6)' ) 'uc1 ',(uc1(it,jt),jt=1,ntype)
    WRITE ( stdout , '(a,<ntype>e16.6)' ) 'uc2 ',(uc2(it,jt),jt=1,ntype)
  enddo
#endif

  return

END SUBROUTINE initialize_param_non_bonded


! *********************** SUBROUTINE engforce_driver ***************************
!
!> \brief
!! this subroutine is the main driver to perform the potential energy, forces 
!! calculation 
!
! ******************************************************************************
SUBROUTINE engforce_driver 

  USE config,                   ONLY :  natm , ntype, system , simu_cell, atypei ,natmi, atype, itype 
  USE control,                  ONLY :  lnmlj , lcoulomb , lbmhft , lbmhftd, lmorse , lharm , longrange, non_bonded, iefgall_format
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

  ! =================================
  !    n-m lennard-jones potential 
  ! =================================
  if ( lnmlj )                     CALL engforce_nmlj_pbc

  ! =================================
  !   bmft(d) potentials (d:damping) 
  ! =================================
  if ( lbmhftd .or. lbmhft )       CALL engforce_bmhftd_pbc

  ! =================================
  !   coulombic potential 
  ! =================================
  if ( lcoulomb ) then

     allocate( ef(3,natm) , efg(3,3,natm) , mu(3,natm) , theta(3,3,natm))     
     allocate ( pair_thole ( natm ) ) 
     allocate ( pair_thole_distance ( natm ) ) 
     pair_thole = 0
     pair_thole_distance = 0.0_dp
     mu=0.0_dp
     theta=0.0_dp

     ! this subroutine set the total dipole moment from :
     ! static          (if given in control file dip TODO: read from DIPFF )
     ! wannier centers (if given in POSFF )
     ! induced         (if polar are set in control file)
     CALL get_dipole_moments ( mu , theta , didpim )
     theta=0.0_dp
 
     ! ====================================
     ! get all the electrostatic quantities
     ! ====================================
     CALL multipole_ES ( ef , efg , mu , theta , task_coul , damp_ind=.true. , &
                         do_efield=doefield , do_efg=doefg , do_forces=.true. , &
                         do_stress=.true. , do_rec=.true. , do_dir=.true. , do_strucfact =.false. , use_ckrskr = didpim )
     theta=0.0_dp
     mu_t     = mu
     theta_t  = theta
     ef_t     = ef
     efg_t    = efg

     deallocate( ef , efg , mu , theta )      
     deallocate ( pair_thole )  
     deallocate ( pair_thole_distance )  
   
  endif
  ! ===========================
  !   other potentials ...
  ! ===========================
  ! CALL engforce_<other>_pbc


  return

END SUBROUTINE engforce_driver

! *********************** SUBROUTINE engforce_nmlj_pbc *************************
!
!> \brief
!! total potential energy forces for each atoms for a nmlj potential with 
!! periodic boundaries conditions, with or without vnlist.
!
!> \author
!! F.Affouard / FMV
!
!> \note
!! adapted from F. Affouard code. Parallelized in december 2008
!
! ******************************************************************************
SUBROUTINE engforce_nmlj_pbc 

  USE constants,                ONLY :  press_unit
  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz, tau_nonb ,  &
                                        atype , itype , verlet_vdw , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist 
  USE thermodynamic,            ONLY :  u_lj , vir_lj , write_thermo
  USE time,                     ONLY :  forcetimetot 
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer          :: ia , ja , it , jt , j1 , je , jb , ierr 
  integer          :: p1 , p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij 
  real(kind=dp) :: sxij , syij , szij 
  real(kind=dp) :: sr2 , rijsq , srp , srq
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: ptwo ( ntype , ntype )
  real(kind=dp) :: qtwo ( ntype , ntype )
  real(kind=dp) :: ttt1 , ttt2 
  real(kind=dp) :: u , vir 

#ifdef debug_nmlj
  WRITE ( stdout , '(a,2i6)' ) 'debug : atom decomposition ',atom_dec%istart,atom_dec%iend
  do ia=1,natm
    WRITE ( stdout , '(a,i6,a,a)' )  'debug : atype ',ia,'',atype(ia)
    WRITE ( stdout , '(a,i6,a,i4)' ) 'debug : itype ',ia,'',itype(ia)
  enddo
#endif
#ifdef MPI
  ttt1 = MPI_WTIME(ierr) ! timing info
#endif
  cccc = cccc + 1

  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp
  tau_nonb = 0.0d0
 
  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5_dp
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5_dp
    enddo
  enddo

!  if ( lvnlist ) CALL vnlistcheck ( verlet_vdw )

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    ! =====================================
    !  verlet list : ja index in point arrays
    ! =====================================
    if ( lvnlist ) then
      jb = verlet_vdw%point( ia )
      je = verlet_vdw%point( ia + 1 ) - 1
    else
    ! ====================================
    !         else all ja   
    ! ====================================
      jb = 1 
      je = natm
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = verlet_vdw%list ( j1 )
      else 
        ja = j1
      endif
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
      !if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and.  ( ( ia .gt. ja .and. ( MOD ( ia + ja , 2 ) .eq. 0 ) ) .or. &
      !                                                                ( ia .lt. ja .and. ( MOD ( ia + ja , 2 ) .ne. 0 ) ) ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          ! potential energy ( truncated )  
          if ( trunc .eq. 0 ) then
            u = u  + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp )
          endif
          if ( trunc .eq. 1 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) - uc(p1,p2)
          endif
          if ( trunc .eq. 2 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
          endif
          wij = fc(p1,p2) * (srq-srp) * sr2
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          !virial 
          vir = vir + wij * rijsq
          ! forces 
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
          ! stress tensor
          tau_nonb(1,1) = tau_nonb(1,1) + (rxij * fxij + rxij * fxij ) * 0.5_dp
          tau_nonb(1,2) = tau_nonb(1,2) + (rxij * fyij + ryij * fxij ) * 0.5_dp
          tau_nonb(1,3) = tau_nonb(1,3) + (rxij * fzij + rzij * fxij ) * 0.5_dp 
          tau_nonb(2,1) = tau_nonb(2,1) + (ryij * fxij + rxij * fyij ) * 0.5_dp
          tau_nonb(2,2) = tau_nonb(2,2) + (ryij * fyij + ryij * fyij ) * 0.5_dp
          tau_nonb(2,3) = tau_nonb(2,3) + (ryij * fzij + rzij * fyij ) * 0.5_dp
          tau_nonb(3,1) = tau_nonb(3,1) + (rzij * fxij + rxij * fzij ) * 0.5_dp
          tau_nonb(3,2) = tau_nonb(3,2) + (rzij * fyij + ryij * fzij ) * 0.5_dp
          tau_nonb(3,3) = tau_nonb(3,3) + (rzij * fzij + rzij * fzij ) * 0.5_dp
        endif
      endif
    enddo
  enddo
  tau_nonb = tau_nonb / simu_cell%omega / press_unit
  vir = vir/3.0_dp

#ifdef MPI
  ttt2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot + ( ttt2 - ttt1 )
#endif

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 
  
  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm ) 
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm ) 

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 1, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 2, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 3, : ) , 3  )

  u_lj = u 
  vir_lj = vir
  
  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE engforce_nmlj_pbc

! *********************** SUBROUTINE engforce_nmlj_nopbc ***********************
!
!> \brief
!! total potential energy forces for each atoms for a nmlj potential with 
!! *NO* periodic boundaries conditions, with or without vnlist (lvnlist=.TRUE.OR.FALSE.)
!
! ******************************************************************************
SUBROUTINE engforce_nmlj_nopbc 

  USE constants,                ONLY :  press_unit
  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz , itype , ntype , tau_nonb , simu_cell , atom_dec , verlet_vdw
  USE control,                  ONLY :  lvnlist
  USE thermodynamic,            ONLY :  u_lj , vir_lj

  implicit none

  ! local
  integer :: ia , ja , it, jt, j1, je, jb !, ierr
  integer :: p1, p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij , sr2 , rijsq , srp , srq
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: ptwo ( ntype , ntype )
  real(kind=dp) :: qtwo ( ntype , ntype )
  real(kind=dp) :: u, vir

  u = 0.0_dp
  vir = 0.0_dp
  fx = 0.0_dp
  fy = 0.0_dp
  fz = 0.0_dp
  tau_nonb = 0.0_dp

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5_dp
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5_dp
    enddo
  enddo

!  if ( lvnlist ) CALL vnlistcheck ( verlet_vdw )
  do ia = atom_dec%istart, atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = verlet_vdw%point ( ia )
      je = verlet_vdw%point ( ia + 1 ) - 1
    else
      jb = ia 
      je = atom_dec%iend
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = verlet_vdw%list(j1)
      else
        ja = j1
      endif
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          if ( trunc .eq. 1 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq - qlj(p1,p2) * srp ) - uc(p1,p2)
          endif
          if ( trunc .eq. 2 ) then
            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
          endif
          wij = fc(p1,p2) * (srq-srp) * sr2
          vir = vir + wij * rijsq
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
          ! stress tensor
          tau_nonb(1,1) = tau_nonb(1,1) + rxij * fxij
          tau_nonb(1,2) = tau_nonb(1,2) + rxij * fyij
          tau_nonb(1,3) = tau_nonb(1,3) + rxij * fzij
          tau_nonb(2,1) = tau_nonb(2,1) + ryij * fxij
          tau_nonb(2,2) = tau_nonb(2,2) + ryij * fyij
          tau_nonb(2,3) = tau_nonb(2,3) + ryij * fzij
          tau_nonb(3,1) = tau_nonb(3,1) + rzij * fxij
          tau_nonb(3,2) = tau_nonb(3,2) + rzij * fyij
          tau_nonb(3,3) = tau_nonb(3,3) + rzij * fzij
        endif
      endif
    enddo
  enddo
  tau_nonb = tau_nonb / simu_cell%omega / press_unit
  vir = vir / 3.0_dp

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u ) 
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir ) 
  
  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb ( 1 , : )  , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb ( 2 , : )  , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb ( 3 , : )  , 3  )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_nmlj_nopbc

! *********************** SUBROUTINE engforce_nmlj_pbc_noshift *****************
!
!> \brief
!! same as engforce_nmlj_pbc but with no shift in the potential 
!
!> \note
!! if ok it should be merged !
!
!> \todo
!! test it !
!
! ******************************************************************************
SUBROUTINE engforce_nmlj_pbc_noshift 

  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz , itype , verlet_vdw , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist 
  USE thermodynamic,            ONLY :  u_lj , vir_lj
  USE time,                     ONLY :  forcetimetot

  implicit none

  ! local
  integer       :: ia , ja , it , jt , j1 , je , jb , ierr
  integer       :: p1 , p2
  real(kind=dp) :: rxi , ryi, rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: sr2 , rijsq , srp , srq
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: ptwo ( ntype , ntype )
  real(kind=dp) :: qtwo ( ntype , ntype )
  real(kind=dp) :: ttt1 , ttt2
  real(kind=dp) :: u , vir
 
#ifdef MPI 
  ttt1 = MPI_WTIME(ierr) ! timing info
#endif 
  
  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp
  

  do jt = 1 , ntype
    do it = 1 , ntype
      ptwo ( it , jt ) = plj ( it , jt ) * 0.5_dp
      qtwo ( it , jt ) = qlj ( it , jt ) * 0.5_dp
    enddo
  enddo

!  if ( lvnlist ) CALL vnlistcheck ( verlet_vdw )
  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = verlet_vdw%point ( ia )
      je = verlet_vdw%point ( ia + 1 ) - 1
    else
      jb = ia 
      je = atom_dec%iend 
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = verlet_vdw%list ( j1 )
      else
        ja = j1
      endif 
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia )  ) then
        rxij  = rxi - rx ( ja )
        ryij  = ryi - ry ( ja )
        rzij  = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1    = itype ( ia )
        p2    = itype ( ja )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
          sr2 = sigsq(p1,p2) / rijsq
          srp = sr2 ** (ptwo(p1,p2))
          srq = sr2 ** (qtwo(p1,p2))
          u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq - qlj(p1,p2) * srp ) 
          wij = fc(p1,p2) * (srq-srp) * sr2
          vir = vir + wij * rijsq
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
        endif
      endif
    enddo
  enddo
  vir = vir/3.0_dp

#ifdef MPI 
  ttt2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot + ( ttt2 - ttt1 )
#endif

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  u_lj = u
  vir_lj = vir

  return

END SUBROUTINE engforce_nmlj_pbc_noshift

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

  USE config,   ONLY  : natm , natmi , ntype , qia , itype
  USE control,  ONLY  : longrange
  USE rspace,   ONLY  : direct_sum_init
  USE kspace,   ONLY  : kpoint_sum_init , kpoint_sum_init_BZ

  implicit none

  ! local
  integer :: nk , ncmax , j , k

  allocate ( ef_t ( 3 , natm ) )
  allocate ( efg_t ( 3 , 3 , natm ) )
  allocate ( dipia_ind_t ( extrapolate_order+1, 3 , natm ) )
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

  USE config,           ONLY : natm , itype , atypei, ntype , poldipia, invpoldipia, atype 

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
  logical :: task_ind(3)
  
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

  USE config,           ONLY : natm , itype , atypei, ntype , poldipia, atype 

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
  USE config,           ONLY :  natm, qia, rx ,ry ,rz, simu_cell, fx, fy, fz, tau_coul, atype 
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
  else
    u_coul   =      ( u_dir   + u_rec   + u_self  + u_pol  ) * coul_unit
    ef       =      ( ef_dir  + ef_rec  + ef_self          ) 
    efg      =      ( efg_dir + efg_rec + efg_self         ) * coul_unit
    tau_coul =      ( tau_dir + tau_rec                    ) * coul_unit / press_unit
    fx       = fx + ( fx_rec  + fx_dir                     ) * coul_unit
    fy       = fy + ( fy_rec  + fy_dir                     ) * coul_unit
    fz       = fz + ( fz_rec  + fz_dir                     ) * coul_unit
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

  USE control,                  ONLY :  lvnlist
  USE config,                   ONLY :  natm, ntype,simu_cell, qia, rx ,ry ,rz ,itype , atom_dec, verlet_coul , atype, poldipia, ipolar
  USE constants,                ONLY :  piroot
  USE cell,                     ONLY :  kardir , dirkar
  USE tensors_rk,               ONLY :  tensor_rank0, tensor_rank1, tensor_rank2, tensor_rank3 , tensor_rank4 , tensor_rank5
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
        if ( .not. lquad ) cycle
        thetaixx = theta ( 1 , 1 , ia)
        thetaiyy = theta ( 2 , 2 , ia)
        thetaizz = theta ( 3 , 3 , ia)
        thetaixy = theta ( 1 , 2 , ia)
        thetaixz = theta ( 1 , 3 , ia)
        thetaiyz = theta ( 2 , 3 , ia)
        K_dot_Q =           thetaixx * kx * kx + thetaixy * kx * ky + thetaixz * kx * kz
        K_dot_Q = K_dot_Q + thetaixy * kx * ky + thetaiyy * ky * ky + thetaiyz * ky * kz
        K_dot_Q = K_dot_Q + thetaixz * kz * kx + thetaiyz * ky * kz + thetaizz * kz * kz
        K_dot_Q = K_dot_Q / 3.0_dp
        rhonk_R    = rhonk_R - K_dot_Q * ckr(ia)
        rhonk_I    = rhonk_I - K_dot_Q * skr(ia)  ! rhon_R + i rhon_I
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
        if ( .not. lquad ) cycle
        thetaixx = theta ( 1 , 1 , ia)
        thetaiyy = theta ( 2 , 2 , ia)
        thetaizz = theta ( 3 , 3 , ia)
        thetaixy = theta ( 1 , 2 , ia)
        thetaixz = theta ( 1 , 3 , ia)
        thetaiyz = theta ( 2 , 3 , ia)
        K_dot_Q =           thetaixx * kx * kx + thetaixy * kx * ky + thetaixz * kx * kz
        K_dot_Q = K_dot_Q + thetaixy * kx * ky + thetaiyy * ky * ky + thetaiyz * ky * kz
        K_dot_Q = K_dot_Q + thetaixz * kz * kx + thetaiyz * ky * kz + thetaizz * kz * kz
        K_dot_Q = K_dot_Q * recarg
        fx_rec ( ia ) = fx_rec ( ia ) + kx * K_dot_Q
        fy_rec ( ia ) = fy_rec ( ia ) + ky * K_dot_Q
        fz_rec ( ia ) = fz_rec ( ia ) + kz * K_dot_Q
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


SUBROUTINE engforce_bmhftd_pbc

  USE constants,                ONLY :  press_unit
  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz, tau_nonb ,  &
                                        atype , itype , verlet_vdw , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist , lbmhftd
  USE thermodynamic,            ONLY :  u_bmhft , vir_bmhft , write_thermo
  USE time,                     ONLY :  forcetimetot 
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer       :: ia , ja , j1 , je , jb !, ierr
  integer       :: p1 , p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: rijsq , d , erh 
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: f6 , f8
  real(kind=dp) :: u , vir 
  real(kind=dp) :: ir2, ir6 , ir7 , ir8 , ir9 , ir6d ,ir8d 
  real(kind=dp) :: fdiff6, fdiff8

  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp
  tau_nonb = 0.0d0

!  if ( lvnlist ) CALL vnlistcheck ( verlet_vdw )
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    ! =====================================
    !  verlet list : ja index in point arrays
    ! =====================================
    if ( lvnlist ) then
      jb = verlet_vdw%point( ia )
      je = verlet_vdw%point( ia + 1 ) - 1
    else
      ! ====================================
      !         else all ja   
      ! ====================================
      jb = 1
      je = natm
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = verlet_vdw%list ( j1 )
      else
        ja = j1
      endif
      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and. ja .gt. ia ) ) then
!      if ( ( lvnlist .and. ja .ne. ia ) .or. ( .not. lvnlist .and.  ( ( ia .gt. ja .and. ( MOD ( ia + ja , 2 ) .eq. 0 ) ) .or. &
!                                                                      ( ia .lt. ja .and. ( MOD ( ia + ja , 2 ) .ne. 0 ) ) ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if ( Abmhftd(p1,p2) .eq. 0.0_dp ) cycle
        if ( rijsq .lt. rcutsq(p1,p2) ) then

          ir2 = 1.0_dp / rijsq
          d = SQRT(rijsq)
          erh = Abmhftd(p1,p2) * EXP ( - Bbmhftd(p1,p2) * d )
          ir6 = ir2 * ir2 * ir2 
          ir8 = ir6 * ir2
          ir6 = ir6 * Cbmhftd(p1,p2) 
          ir8 = ir8 * Dbmhftd(p1,p2) 
          if ( lbmhftd ) then
            CALL TT_damping_functions ( BDbmhftd(p1,p2), 1.0_dp , d , f6 , fdiff6, order=6 ) 
            CALL TT_damping_functions ( BDbmhftd(p1,p2), 1.0_dp , d , f8 , fdiff8, order=8 ) 
            ir6d = ir6 * f6
            ir8d = ir8 * f8
          else
            ir6d = ir6
            ir8d = ir8
          endif
          ir7 = 6.0_dp * ir6d / d 
          ir9 = 8.0_dp * ir8d / d
          u = u + erh - ir6d - ir8d 
          wij  =  Bbmhftd(p1,p2) * erh - ir7 - ir9 + ir6 * fdiff6 + ir8 * fdiff8 
          wij  = wij / d 
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          vir  = vir + ( rxij * fxij + ryij * fyij + rzij * fzij )
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
          ! stress tensor
          tau_nonb(1,1) = tau_nonb(1,1) + ( rxij * fxij + rxij * fxij )
          tau_nonb(1,2) = tau_nonb(1,2) + ( rxij * fyij + ryij * fxij ) 
          tau_nonb(1,3) = tau_nonb(1,3) + ( rxij * fzij + rzij * fxij )
          tau_nonb(2,1) = tau_nonb(2,1) + ( ryij * fxij + rxij * fyij ) 
          tau_nonb(2,2) = tau_nonb(2,2) + ( ryij * fyij + ryij * fyij ) 
          tau_nonb(2,3) = tau_nonb(2,3) + ( ryij * fzij + rzij * fyij ) 
          tau_nonb(3,1) = tau_nonb(3,1) + ( rzij * fxij + rxij * fzij )
          tau_nonb(3,2) = tau_nonb(3,2) + ( rzij * fyij + ryij * fzij )
          tau_nonb(3,3) = tau_nonb(3,3) + ( rzij * fzij + rzij * fzij ) 
        endif
     endif
   enddo 
  enddo
  vir = vir/3.0_dp
  tau_nonb = tau_nonb / simu_cell%omega / press_unit * 0.5_dp

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 1, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 2, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 3, : ) , 3  )

  ! ======================================
  !         direct to cartesian 
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  u_bmhft = u
  vir_bmhft = vir

  return

END SUBROUTINE engforce_bmhftd_pbc

! *********************** SUBROUTINE engforce_morse_pbc ************************
!
!> \brief
!! total potential energy forces for each atoms for a morse potential with 
!! periodic boundaries conditions, with or without vnlist
!! (lvnlist=.TRUE.OR.FALSE.)
!
!> \warning
!! not fully tested
!
! ******************************************************************************
SUBROUTINE engforce_morse_pbc 

  USE constants,                ONLY :  press_unit
  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz, atype , itype , verlet_vdw, & 
                                        ntype , simu_cell , atom_dec , tau_nonb
  USE control,                  ONLY :  lvnlist  
  USE thermodynamic,            ONLY :  u_morse , vir_morse
  USE time,                     ONLY :  forcetimetot
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer          :: ia , ja , j1 , je , jb , ierr
  integer          :: p1 , p2
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: rijsq , d , erh , erh2 
  real(kind=dp) :: wij , fxij , fyij , fzij
  real(kind=dp) :: forcetime1 , forcetime2 
  real(kind=dp) :: u , vir , u2
  real(kind=dp) :: expon 

#ifdef debug_morse
  if ( ionode ) then
  do ia=1, natm
    WRITE ( stdout , '(a,a)' )  'debug : atype ',atype(ia)
    WRITE ( stdout , '(a,i4)' ) 'debug : itype ',itype(ia)
  enddo
  endif
#endif

#ifdef MPI
  forcetime1 = MPI_WTIME(ierr) ! timing info
#endif

  u   = 0.0_dp
  vir = 0.0_dp
  fx  = 0.0_dp
  fy  = 0.0_dp
  fz  = 0.0_dp

!  if ( lvnlist ) CALL vnlistcheck ( verlet_vdw )
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    if ( lvnlist ) then
      jb = verlet_vdw%point( ia )
      je = verlet_vdw%point( ia + 1 ) - 1
    else
      jb = 1 
      je = natm !atom_dec%iend
    endif
    do j1 = jb, je
      if ( lvnlist ) then
        ja = verlet_vdw%list ( j1 )
      else
        ja = j1
      endif
      if ( ja .ne. ia ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq= rxij * rxij + ryij * ryij + rzij * rzij
        p1   = itype ( ja )
        p2   = itype ( ia )
        if ( rijsq .lt. rcutsq(p1,p2) ) then
           d   = SQRT ( rijsq ) 
           erh   = EXP ( - d * rhomor(p1,p2) ) 
           expon = erh * rs(p1,p2) - 1.0_dp
           expon = expon * expon - 1.0_dp
           u     = u  + epsmor(p1,p2) * expon 
           expon = EXP( rhomor(p1,p2) * ( 1.0_dp - d / sigmamor(p1,p2) )  ) 
           u2    = u2 + ( expon * ( expon - 2.0_dp ) ) * epsmor(p1,p2)
        !  if ( ia .eq. 1 )  print*,d,erh,erh2,expon,u
!          if ( trunc .eq. 1 ) then
!            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) - uc(p1,p2)
!          endif
!          if ( trunc .eq. 2 ) then
!            u =  u + epsp(p1,p2) * ( plj(p1,p2) * srq -qlj(p1,p2) * srp ) + uc1(p1,p2) * rijsq - uc2(p1,p2)
!          endif
          wij  = fm(p1,p2) * ( erh - rs(p1,p2) * erh2 ) 
          vir  = vir + wij * rijsq
          fxij = wij * rxij
          fyij = wij * ryij
          fzij = wij * rzij
          fx ( ia ) = fx ( ia ) + fxij
          fy ( ia ) = fy ( ia ) + fyij
          fz ( ia ) = fz ( ia ) + fzij
          fx ( ja ) = fx ( ja ) - fxij
          fy ( ja ) = fy ( ja ) - fyij
          fz ( ja ) = fz ( ja ) - fzij
          ! stress tensor
          tau_nonb(1,1) = tau_nonb(1,1) + rxij * fxij
          tau_nonb(1,2) = tau_nonb(1,2) + rxij * fyij
          tau_nonb(1,3) = tau_nonb(1,3) + rxij * fzij
          tau_nonb(2,1) = tau_nonb(2,1) + ryij * fxij
          tau_nonb(2,2) = tau_nonb(2,2) + ryij * fyij
          tau_nonb(2,3) = tau_nonb(2,3) + ryij * fzij
          tau_nonb(3,1) = tau_nonb(3,1) + rzij * fxij
          tau_nonb(3,2) = tau_nonb(3,2) + rzij * fyij
          tau_nonb(3,3) = tau_nonb(3,3) + rzij * fzij

        endif
      endif
    enddo
  enddo
  vir = vir/3.0_dp
  tau_nonb = tau_nonb / simu_cell%omega / press_unit

#ifdef MPI
  forcetime2 = MPI_WTIME(ierr) ! timing info
  forcetimetot = forcetimetot+(forcetime2-forcetime1)
#endif

  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( u   )
  CALL MPI_ALL_REDUCE_DOUBLE_SCALAR ( vir )

  CALL MPI_ALL_REDUCE_DOUBLE ( fx , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fy , natm )
  CALL MPI_ALL_REDUCE_DOUBLE ( fz , natm )

  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 1, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 2, : ) , 3  )
  CALL MPI_ALL_REDUCE_DOUBLE ( tau_nonb( 3, : ) , 3  )

  u_morse = u
  vir_morse = vir

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE engforce_morse_pbc

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
  USE config,           ONLY :  natm , atype , fx , fy , fz , ntype , dipia , quadia , qia , ntypemax, poldipia , invpoldipia , itype 
  USE control,          ONLY :  longrange , calc
  USE thermodynamic,    ONLY :  u_pol, u_coul
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
  logical       :: task_static (3), task_ind(3), ldip
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
    !   two methods : 
    !                  simple linear/polynomial extrapolation   ( algo_ext_dipole .eq. 'poly' ) : instable for extrapolate_order > 1  
    !                  ASPC (Always Stable Predictor Corrector) ( algo_ext_dipole .eq. 'aspc  ) by Kolafa J. Comp. Chem. v25 n3 (2003)
    ! ==========================================================
    if ( iscf.ne.1 .or. ( calc .ne. 'md' .and.  calc .ne. 'opt' ) ) then

      CALL induced_moment       ( Efield , mu_ind )  ! Efield in ; mu_ind and u_pol out

    else

      if ( algo_ext_dipole .eq. 'poly' ) CALL extrapolate_dipole_poly ( mu_ind ) 
      if ( algo_ext_dipole .eq. 'aspc' ) CALL extrapolate_dipole_aspc ( mu_ind , Efield , key=1 )  ! predictor


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

  USE config,           ONLY : natm , ntype , dipia , quadia , dipia_wfc, fx,fy,fz, atype 
  USE io,               ONLY : kunit_DIPFF

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

! *********************** SUBROUTINE TT_damping_functions **********************
!> \brief
!! Tang-Toennies damping function
!> \author
!! FMV
!> \date 
!! February 2014
! ******************************************************************************
SUBROUTINE TT_damping_functions(b,c,r,f,fd,order)

  USE tt_damp,          ONLY : E_TT ! Tang-Toennies coefficients define once up to order 8 

  implicit none

  ! global 
  integer :: order
  real(kind=dp) :: b , c , r, f , fd  ! f damping function , fd first derivative
 

  ! local 
  integer :: k 
  real(kind=dp) :: expbdr , br
  
  br = b * r 
  expbdr = EXP(-br) * c

  f = E_TT(order)
  do k=order-1,1,-1
    f = f * br + E_TT(k)
  enddo 
  f = f * br + E_TT(0)
  f = 1.0_dp - f * expbdr

  ! derivative of f checked on sage 18/02/14 worksheet BMHFTD
  fd = E_TT(order) * ( br ) ** ( order ) * expbdr * b 

  return

END SUBROUTINE TT_damping_functions

SUBROUTINE extrapolate_dipole_poly ( mu_ind )

  USE config,           ONLY :  natm
  USE md,               ONLY :  itime 

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


! *********************** SUBROUTINE write_DIPFF ******************************
!
!>\brief
! write dipoles at ions to DIPFF file
!
! ******************************************************************************
SUBROUTINE write_DIPFF 

  USE io,                       ONLY :  kunit_DIPFF
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

SUBROUTINE find_min ( g0 , g1 , d1 , gmin , dmin )

  USE config,           ONLY :  natm , itype

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

SUBROUTINE get_rmsd_mu ( rmsd , g0 , mu ) 

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

  USE config,           ONLY : natm, itype

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


! *********************** SUBROUTINE write_EFGALL ******************************
!
!>\brief
! write Electric Field Gradient at ions to EFGALL file
!
! ******************************************************************************
SUBROUTINE write_EFGALL

  USE io,                       ONLY :  kunit_EFGALL
  USE cell,                     ONLY :  kardir , periodicbc , dirkar
  USE control,                  ONLY :  lstatic, iefgall_format
  USE config,                   ONLY :  system , natm , ntype , atype , itype , simu_cell, atypei, natmi

  implicit none

  ! local
  integer :: ia , it

  if ( ionode .and. doefg ) then

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
  USE cell,                     ONLY :  kardir , periodicbc , dirkar
  USE control,                  ONLY :  lstatic
  USE config,                   ONLY :  system , natm , ntype , atype , itype , simu_cell, atypei, natmi

  implicit none

  ! local
  integer :: ia , it

  if ( ionode .and. doefg ) then

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

