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

! *********************** MODULE non_bonded *****************************************
!> \brief
!! Module related to force-field calculation
! ******************************************************************************
MODULE non_bonded 

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  !USE kspace,                           ONLY :  kmesh 
  !USE rspace,                           ONLY :  rmesh
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  !USE tensors_rk,                       ONLY :  interaction_dd
  USE mpimdff

  implicit none

  integer :: cccc=0
  logical, PRIVATE :: symmetric_pot
  real(kind=dp)     :: utail               !< long-range correction (energy) of short-range interaction 
  real(kind=dp)     :: ptail               !< long-range correction (virial) of short-range interaction 
  logical, SAVE     :: lKA               !< use Kob-Andersen model for BMLJ                        

  character(len=60) :: ctrunc                                                     !< truncation of nmlj
  character(len=60) :: ctrunc_allowed(3)                                          !< truncation of nmlj 
  data                 ctrunc_allowed / 'notrunc', 'linear' , 'quadratic' /       !< see initialize_param_nmlj
  integer           :: trunc                                                      !< integer definition of truncation 

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


CONTAINS

! *********************** SUBROUTINE non_bonded_init ********************************
!> \brief
!! non bonded force field initialisation
! ******************************************************************************
SUBROUTINE non_bonded_init

!  USE control,                  ONLY :  calc , lnmlj , lcoulomb , lmorse , longrange, lnon_bonded

  implicit none

  ! local
  character(len=132)    :: filename
  integer               :: ioerr

  namelist /non_bondedtag/               &
                         lKA           , &       
                         ctrunc        , &
                         symmetric_pot , &
                         qlj           , & 
                         plj           , & 
                         sigmalj       , &
                         epslj         , &
                         sigmamor      , &
                         epsmor        , &
                         rhomor        , &
                         Abmhftd       , &
                         Bbmhftd       , &
                         Cbmhftd       , &
                         Dbmhftd       , &
                         BDbmhftd


  ! ================================
  ! defaults values for non bonded tags 
  ! ================================
  CALL non_bonded_default_tag

  ! ================================
  ! read non_bonded tags
  ! ================================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
    READ ( stdin , non_bondedtag, iostat=ioerr)
    if ( ioerr .lt. 0 )  then
      io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : non_bondedtag section is absent'
      STOP
    elseif ( ioerr .gt. 0 )  then
      io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : non_bondedtag wrong tag',ioerr
      STOP
    endif
  CLOSE ( stdin )

  ! ================================
  ! check non_bonded tags values
  ! ================================
  CALL non_bonded_check_tag

  ! ================================
  ! initialize constant parameters
  ! ================================
  CALL initialize_param_non_bonded

  ! ================================
  !  print non_bonded info
  ! ================================
  CALL non_bonded_print_info(stdout)

END SUBROUTINE non_bonded_init


! *********************** SUBROUTINE non_bonded_default_tag *************************
!> \brief
!! set default values to non_bonded tags
! ******************************************************************************
SUBROUTINE non_bonded_default_tag

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

  return

END SUBROUTINE non_bonded_default_tag


! *********************** SUBROUTINE non_bonded_check_tag ***************************
!> \brief
!! check non_bonded tag values
! ******************************************************************************
SUBROUTINE non_bonded_check_tag

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
  CALL check_allowed_tags ( size ( ctrunc_allowed ), ctrunc_allowed, ctrunc, 'in non_bondedtag','ctrunc' ) 

  if ( ctrunc .eq. 'notrunc'   ) trunc = 0
  if ( ctrunc .eq. 'linear'    ) trunc = 1
  if ( ctrunc .eq. 'quadratic' ) trunc = 2

  ! =================================================
  !  KOB-ANDERSEN MODEL --- PhysRevE 51-4626 (1995) 
  ! =================================================
  if ( lKA .and. ntype .ne. 2 ) then
    io_node WRITE ( stdout , '(a)' ) 'ERROR non_bondedtag lKA should be used with 2 differents types'
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
    enddo
  endif

  if ( lbmhftd ) CALL get_TT_damp

  return 

END SUBROUTINE non_bonded_check_tag

! *********************** SUBROUTINE non_bonded_print_info ***************************
!> \brief
! ******************************************************************************
SUBROUTINE non_bonded_print_info(kunit)

  USE control,  ONLY :  lnmlj, lmorse, lbmhft, lbmhftd, cutshortrange, lreducedN
  USE config,   ONLY :  ntype, atypei

  implicit none

  ! global
  integer            :: kunit

  !local 
  integer            :: it1 , it2 


  if ( ionode ) then
    separator(kunit)    
    blankline(kunit)
    WRITE ( kunit ,'(a)')               'NON BONDED MODULE ... WELCOME'
    blankline(kunit)
    lseparator(kunit) 
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

END SUBROUTINE non_bonded_print_info

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


! *********************** SUBROUTINE engforce_bmhftd_pbc *****************

SUBROUTINE engforce_bmhftd_pbc

  USE constants,                ONLY :  press_unit
  USE config,                   ONLY :  natm , rx , ry , rz , fx , fy , fz, tau_nonb ,  &
                                        atype , itype , verlet_vdw , ntype , simu_cell , atom_dec
  USE control,                  ONLY :  lvnlist , lbmhftd
  USE thermodynamic,            ONLY :  u_bmhft , vir_bmhft , write_thermo
  USE time,                     ONLY :  forcetimetot 
  USE cell,                     ONLY :  kardir , dirkar
  USE tt_damp,                  ONLY :  TT_damping_functions 

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
            CALL TT_damping_functions ( BDbmhftd(p1,p2), 1.0_dp , d , f6 , fdiff6, 6 ) 
            CALL TT_damping_functions ( BDbmhftd(p1,p2), 1.0_dp , d , f8 , fdiff8, 8 ) 
            ir6d = ir6 * f6
            ir8d = ir8 * f8
          else
            !print*,'here in non_bonded'
            fdiff6 = 0.0_dp
            fdiff8 = 0.0_dp
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

END MODULE non_bonded 
! ===== fmV =====

