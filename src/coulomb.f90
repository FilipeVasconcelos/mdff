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

  logical, PRIVATE :: symmetric_pot
  logical, SAVE     :: lautoES           !< auto-determination of Ewald parameter from epsw ( accuracy)
  logical           :: doefield , doefg

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

  namelist /coulombtag/   lautoES      , &
                          epsw         , &
                          alphaES      , &  
                          kES          , &
                          symmetric_pot, & 
                          doefield     , &
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
  thole_functions       = .false.
  thole_function_type   = 'linear'
  thole_param           = 1.662_dp

  task_coul = .false.

  return

END SUBROUTINE coulomb_default_tag


! *********************** SUBROUTINE coulomb_check_tag ***************************
!> \brief
!! check coulomb tag values
! ******************************************************************************
SUBROUTINE coulomb_check_tag

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

  if ( symmetric_pot ) then
    !thole param
    thole_param(:,it) = thole_param(it,:)  
  endif

  ! =======================
  !  thole_function_type 
  ! =======================
  CALL check_allowed_tags( size(thole_function_type_allowed), thole_function_type_allowed, thole_function_type, 'in fieldtag','thole_function_type' ) 
   
  if ( thole_function_type .eq. 'expon1' ) thole_param = 0.572_dp
  if ( thole_function_type .eq. 'expon2' ) thole_param = 1.9088_dp
  if ( thole_function_type .eq. 'gauss' ) thole_param = 0.957_dp

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

  USE config,   ONLY  : natm , natmi , ntype , qia , itype
  USE control,  ONLY  : longrange
  USE rspace,   ONLY  : direct_sum_init
  USE kspace,   ONLY  : kpoint_sum_init , kpoint_sum_init_BZ

  implicit none

  ! local
  integer :: nk , ncmax , j , k

  allocate ( ef_t ( 3 , natm ) )
  allocate ( efg_t ( 3 , 3 , natm ) )
  allocate ( dipia_ind_t ( 5 , 3 , natm ) ) ! 5 => max extrapolate_order+1
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

   allocate ( pair_thole ( natm ) )
   allocate ( pair_thole_distance ( natm ) )

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

   deallocate ( pair_thole  )
   deallocate ( pair_thole_distance )

   return

END SUBROUTINE finalize_coulomb

END MODULE coulomb 
! ===== fmV =====

