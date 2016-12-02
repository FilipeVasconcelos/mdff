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
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
! ===== fmV =====

! ======= Hardware =======
#include "symbol.h"
!#define debug
!#define debug_invert
!#define debug_invert2
!#define debug_chisq
! ======= Hardware =======

! *********************** MODULE rmc ******************************************
!> \brief
!! Module related to reverve monte carlo simuations 
! ******************************************************************************
MODULE rmc 

  USE io,               ONLY :  ionode, stdout, stdin, stderr
  USE constants,        ONLY :  dp
  USE mpimdff

  implicit none

  integer       :: acceptance_max      !< maximum number of acceptance configurations
  integer       :: nprint              !< print period
  integer       :: efg_start           !< calculate efg from step efg_start
  real(kind=dp) :: cut_rmc             !< radial cut-off 
  real(kind=dp) :: resolution_gr       !< resolution in g(r) distribution 
  real(kind=dp) :: max_displacement    !< maximum dispalcement
  real(kind=dp) :: temp_rmc(5)         !< rmc temperature
  logical       :: restart_rmc         !< restart from current POSFF
  logical       :: lrmc_gr             !<  
  logical       :: lrmc_var            !<  
  logical       :: lrmc_u              !<  
  logical       :: lrmc_vzz            !<  
  logical       :: lrmc_eta            !<  
  real(kind=dp) :: resolution_u_efg    !< resolution in U      (efg tensor)
  real(kind=dp) :: resolution_vzz_efg  !< resolution in vzz    (efg tensor)
  real(kind=dp) :: resolution_eta_efg  !< resolution in eta    (efg tensor)
  real(kind=dp) :: u_efg_min           !< minimum value of U   (efg tensor)
  real(kind=dp) :: vzz_efg_min         !< minimum value of vzz (efg tensor)
  real(kind=dp) , dimension ( :, : , :) , allocatable :: dab
  real(kind=dp) , dimension ( :, : , :) , allocatable :: dab_updated

  !character(len=60) :: rmc_allowed(4)
  !data rmc_allowed  / 'rmc'mc-invert-efg' /

CONTAINS

! *********************** SUBROUTINE rmc_init *******************************
!
!> \brief
!! initialize reverse monte carlo parameters
!
! ******************************************************************************
SUBROUTINE rmc_init

  USE config,                   ONLY :  simu_cell
  USE control,                  ONLY :  calc

  implicit none

  ! local
  integer            :: ioerr 
  character(len=132) :: filename

  namelist /rmctag/   nprint              , &
                      lrmc_gr             , &
                      lrmc_var            , &
                      lrmc_vzz            , &
                      lrmc_eta            , &
                      lrmc_u              , &
                      efg_start           , &
                      max_displacement    , &
                      cut_rmc             , &
                      temp_rmc            , &
                      acceptance_max      , &
                      restart_rmc         , &
                      u_efg_min           , &
                      vzz_efg_min         , &
                      resolution_u_efg    , &
                      resolution_vzz_efg  , &
                      resolution_eta_efg  , &
                      resolution_gr  

  CALL rmc_default_tag
  
  ! =================
  !  read grtag tags
  ! =================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , rmctag , iostat=ioerr )
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : rmctag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : rmctag wrong tag',ioerr
   STOP
  endif
  CLOSE ( stdin )

  CALL rmc_check_tag

  CALL rmc_print_info(stdout)

  return 
 
END SUBROUTINE rmc_init

! *********************** SUBROUTINE rmc_default_tag ****************************
!
!> \brief
!! set default values to rmc tag
!
! ******************************************************************************
SUBROUTINE rmc_default_tag

  USE config,           ONLY : simu_cell

  implicit none

  ! ===============
  !  default value
  ! ===============
  nprint           = 1
  restart_rmc      = .false.
  max_displacement = 0.1_dp
  efg_start        = 0
  lrmc_gr          = .false.          
  lrmc_var         = .false.          
  lrmc_u           = .false.          
  lrmc_vzz         = .false.         
  lrmc_eta         = .false.          

  return

END SUBROUTINE rmc_default_tag

SUBROUTINE rmc_check_tag

  implicit none

  return

END SUBROUTINE


! *********************** SUBROUTINE rmc_print_info *****************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE rmc_print_info(kunit)

  USE control,                  ONLY :  calc

  implicit none
 
  ! local
  integer :: kunit

   if ( ionode ) then
     blankline(kunit)
     WRITE ( kunit ,'(a)')                 'RMC MODULE ... WELCOME'
     WRITE ( kunit ,'(a,i5)')              'rmc calculation'
     if ( lrmc_var ) WRITE ( kunit ,'(a)') 'rmc-invert algorithm (add references!!)'
     WRITE ( kunit ,'(a,i5)')              ' acceptance max : ',acceptance_max
     blankline(kunit)
   endif 
  return

END SUBROUTINE rmc_print_info

! *********************** SUBROUTINE rmc_print_info *****************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE rmc_main

  USE constants,                ONLY : pi , dzero
  USE io,                       ONLY : kunit_RMCFF, kunit_GRTFF, &
                                       kunit_POSFF , kunit_RMCLOG, kunit_TRAJFF , kunit_DTIBUFF , kunit_DTETAFF , kunit_DTVZZFF
  USE config,                   ONLY : natm, ntype , simu_cell, natmi, atypei, atype, allowedmove, rx, ry , rz, &
                                       config_alloc, coord_format_allowed, write_CONTFF, system , &
                                       vx, vy, vz, fx, fy, fz, rhoN , config_print_info , itype , atom_dec , write_trajff_xyz
  USE radial_distrib,           ONLY : nbins, npairs , read_grtff, gr_main , gr_alloc , gr , resg, cutgr
  USE control,                  ONLY : lcoulomb
  USE cell,                     ONLY : lattice, dirkar, kardir
  USE efg,                      ONLY : PANU , PANeta, PANvzz , resu , reseta , resvzz , umin , vzzmin , &
                                       read_dtibuff , read_dtvzzff , read_dtetaff , dibUtot , dibetatot , dibvzztot , &
                                       mu , efg_alloc , efg_mesh_alloc , multipole_efg_es , efg_ia , nmr_convention
  USE field,                    ONLY : field_init , km_coul 
  USE time,                     ONLY : rmcgrtimetot_comm

  implicit none

  integer :: kacc , kall , kprint
  integer :: i , ia, it, isq , igr , mp , bin , ui
  integer :: random_particule , nmax , ierr

  real(kind=dp),     dimension (     : , : ) , allocatable :: grr_exp
  real(kind=dp),     dimension (     : , : ) , allocatable :: grr_calc
  real(kind=dp),     dimension ( : , : , : ) , allocatable :: dibU_exp
  real(kind=dp),     dimension ( : , : , : ) , allocatable :: dibU_calc
  real(kind=dp),     dimension (     : , : ) , allocatable :: dibvzz_exp
  real(kind=dp),     dimension (     : , : ) , allocatable :: dibvzz_calc
  real(kind=dp),     dimension (     : , : ) , allocatable :: dibeta_exp
  real(kind=dp),     dimension (     : , : ) , allocatable :: dibeta_calc
  real(kind=dp)     :: x , y , z , rr
  ! TODO : create a type structure for chisq
  character(len=3)  :: label_chisq(6) 
  data label_chisq / 'rmc' , 'var' , 'efg' , 'vzz' , 'eta' ,'tot' /
  real(kind=dp)     :: chisq(6)       ! array of chisq : chisq_rmc, chisq_var, chisq_uefg, chisq_eta, chisq_vzz ,total
  real(kind=dp)     :: chisq_old(6)   ! array of chisq same order as above the last one is the total
  real(kind=dp)     :: delta_chisq(6) ! array of delta same order as above the last one is the total
  logical           :: accept(6)
  real(kind=dp)     :: metro, randmetro, randpart
  real(kind=dp)     :: exec_loops,timing_loop1,timing_loop2,timing_all_start
  character(len=60) :: cpos
  character(len=20) :: FMT 
  

#ifdef MPI
  timing_all_start = MPI_WTIME(ierr)
#endif
  kall = 1 
  kacc = 1
  kprint = 1

  ! main quantities
  chisq        = 0.0_dp

  if ( .not. restart_rmc ) then 
    ! ==============================================================
    ! starting from scratch some configuration information
    ! are readed from RMCFF file : number of atoms, types, cell ...
    ! ==============================================================
    OPEN(UNIT=kunit_RMCFF,FILE='RMCFF')
      READ( kunit_RMCFF, * ) natm
      READ( kunit_RMCFF, * ) system 
      READ( kunit_RMCFF ,* ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
      READ( kunit_RMCFF ,* ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
      READ( kunit_RMCFF ,* ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
      READ( kunit_RMCFF, * ) ntype
      READ( kunit_RMCFF, * ) (atypei(it),it=1,ntype)
      READ( kunit_RMCFF, * ) (natmi(it),it=1,ntype)
    CLOSE(UNIT=kunit_RMCFF) 
    CALL lattice ( simu_cell ) 
    rhoN = REAL ( natm , kind = dp ) / simu_cell%omega
    npairs = ntype * ( ntype + 1 ) / 2
    CALL config_alloc
    CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms' )
    CALL typeinfo_init
  else
    ! ===============================================================
    !  if we restart from previous job we read POSFF file
    ! ===============================================================
    OPEN ( kunit_POSFF , file = 'POSFF' )
    READ ( kunit_POSFF ,* ) natm
    READ ( kunit_POSFF ,* ) system
    READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
    READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
    READ ( kunit_POSFF ,* ) simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
    READ ( kunit_POSFF ,* ) ntype 
    READ ( kunit_POSFF ,* ) ( atypei ( it ) , it = 1 , ntype )
    IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on POSFF : ', atypei ( 1:ntype )
    READ( kunit_POSFF ,*) ( natmi ( it ) , it = 1 , ntype )
    READ( kunit_POSFF ,*) cpos

    ! =============
    !  check cpos
    ! =============
    CALL check_allowed_tags( size( coord_format_allowed ), coord_format_allowed, cpos, 'in POSFF at line 9','' )

    if ( cpos .eq. 'Direct' .or. cpos .eq. 'D' ) then 
      io_node WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in direct coordinates in POSFF'
    else if ( cpos .eq. 'Cartesian' .or. cpos .eq. 'C' ) then 
      io_node WRITE ( stdout      ,'(A,20A3)' ) 'atomic positions in cartesian coordinates in POSFF'
    endif

    CALL lattice ( simu_cell )
    rhoN = REAL ( natm , kind = dp ) / simu_cell%omega
    CALL config_alloc

    READ  ( kunit_POSFF , * ) ( atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , &
                                               vx ( ia ) , vy ( ia ) , vz ( ia ) , &
                                               fx ( ia ) , fy ( ia ) , fz ( ia ) , ia = 1 , natm )

    if ( cpos .eq. 'Direct' .or. cpos .eq. 'D' ) then
      CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
    else if ( cpos .eq. 'Cartesian' .or. cpos .eq. 'C' ) then
    endif
    CLOSE ( kunit_POSFF )
    CALL typeinfo_init
    CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms' )
    npairs = ntype * ( ntype + 1 ) / 2
  endif

  ! =============
  !  init rmc_gr
  ! =============
  if ( lrmc_gr ) then
    ! ================================================
    !  - Experimental pdfs readed from GRTFF.exp.
    !      - At the moment we are stuck to follow the grid 
    !        given by GRTFF.exp
    !  - Experimental EFG's from DTIBUFF.exp 
    !   or from DTVZZFF.exp and DTETAFF.exp. 
    ! ================================================
    nbins = int ( cut_rmc / resolution_gr ) + 1
    resg  = resolution_gr
    cutgr = cut_rmc
    nmax = maxval(natmi(1:))
    allocate ( grr_exp  ( 0 : npairs , 0 : nbins-1 ) )
    allocate ( grr_calc ( 0 : npairs , 0 : nbins-1 ) )
    CALL read_GRTFF( grr_exp , 'GRTFF.exp')
#ifdef debug
  do igr=0,nbins-1
#ifdef GFORTRAN
    write(FMT,* ) npairs+2
    write(*,'('// ADJUSTL(FMT) //'e20.8)') ( grr_exp( mp , igr ), mp=0,npairs ) 
#else
    write(*,'(<npairs+2>e20.8)') ( grr_exp( mp , igr ), mp=0,npairs ) 
#endif
  enddo
#endif
  endif
  ! ==============
  !  init rmc_var 
  ! ==============
  if ( lrmc_var ) then
    allocate ( dab ( nmax, ntype , natm ) ) 
    allocate ( dab_updated ( nmax, ntype , natm ) ) 
    dab         = 0.0_dp
    dab_updated = 0.0_dp
  endif
  ! =============
  !  init rmc_u 
  ! =============
  if ( lrmc_u ) then
    resu   = resolution_U_efg
    umin   = u_efg_min
    vzzmin = vzz_efg_min
    PANU   = int ((2.0_dp*ABS (u_efg_min))/resolution_U_efg)
    allocate ( dibU_exp    ( 6 , 0:ntype , 0:PANU ) ) 
    allocate ( dibU_calc   ( 6 , 0:ntype , 0:PANU ) ) 
    CALL read_DTIBUFF( dibU_exp   , 'DTIBUFF.exp')
    if ( ionode ) write ( stdout , * ) 'reading U efg component from DTIBUFF.exp'
  endif
  ! ==============
  !  init rmc_vzz 
  ! ==============
  if ( lrmc_vzz ) then
    resvzz = resolution_vzz_efg
    PANvzz = int ((2.0_dp*ABS (vzz_efg_min ))/ resolution_vzz_efg )
    allocate ( dibvzz_exp  (     0:ntype , 0:PANvzz ) ) 
    allocate ( dibvzz_calc (     0:ntype , 0:PANvzz ) ) 
    CALL read_DTVZZFF( dibvzz_exp , 'DTVZZFF.exp')
  endif
  ! ==============
  !  init rmc_eta
  ! ==============
  if ( lrmc_eta ) then
    reseta = resolution_eta_efg
    PANeta = int ( 1.0_dp / resolution_eta_efg )
    allocate ( dibeta_exp  (     0:ntype , 0:PANeta ) ) 
    allocate ( dibeta_calc (     0:ntype , 0:PANeta ) ) 
    CALL read_DTETAFF( dibeta_exp , 'DTETAFF.exp')
  endif
 
  OPEN (unit = kunit_TRAJFF ,file = 'TRAJFF')
  OPEN (unit = kunit_RMCLOG ,file = 'RMCLOG')

  ! ================================================
  !              "calc" structure
  ! ================================================
  if ( .not. restart_rmc ) then 
    ! generate random structure
    if ( ionode ) write( stdout , '(a)') "random positions in the box : (direct coordinates)"
    do ia = 1 , natm 
      CALL RANDOM_NUMBER(x)
      CALL RANDOM_NUMBER(y)
      CALL RANDOM_NUMBER(z)
      rx(ia) = x 
      ry(ia) = y
      rz(ia) = z
    enddo      
    CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
  endif 
  CALL config_print_info( stdout )

  if ( lrmc_gr ) then
    ! alloc main histogram for g(r) calculation
    CALL gr_alloc
    CALL rmc_gr ( grr_calc )
!#ifdef debug
  do igr=0,nbins-1
    rr  = ( REAL( igr , kind = dp ) + 0.5_dp ) * resolution_gr
    write(*,'(3e20.8)') rr, grr_calc ( 0 , igr ) , grr_exp ( 0 , igr )
  enddo
!#endif
  endif

  ! =============================================
  !  chi² EFG : U_i Vzz eta
  ! =============================================
  if ( ( lrmc_u .or. lrmc_vzz .or. lrmc_eta ) .and. kacc .ge. efg_start ) then 
    CALL field_init
    CALL typeinfo_init
    CALL efg_alloc
    CALL efg_mesh_alloc
    if ( lcoulomb ) CALL do_split ( km_coul%nk , myrank , numprocs , km_coul%kpt_dec ,'k-pts' )
    CALL rmc_efg ( dibU_calc , dibeta_calc , dibvzz_calc )
    if ( lrmc_u )   CALL eval_chisq_Uefg ( dibU_exp   , dibU_calc   , chisq(3) )
    if ( lrmc_vzz ) CALL eval_chisq_vzz  ( dibvzz_exp , dibvzz_calc , chisq(4) )
    if ( lrmc_eta ) CALL eval_chisq_eta  ( dibeta_exp , dibeta_calc , chisq(5) )
  endif

  ! ================================================
  !  chi first/input structre
  ! ================================================
  if ( lrmc_gr  ) CALL eval_chisq_rmc ( grr_exp , grr_calc , chisq(1) )
  if ( lrmc_var ) CALL eval_chisq_var_distance ( chisq(2) )
  chisq(6)     = sum( chisq(1:5) )    
  chisq_old    = chisq
  delta_chisq = 0.0_dp

#ifdef debug_chisq
    write(stdout,'(a)') '-------------------------------------------------------------------------------------'
    write(stdout,'(a)') 'label     current          old             delta            metropolis         "temp"'
    write(stdout,'(a)') '-------------------------------------------------------------------------------------'
    do isq=1,5
      write(stdout,'(a,3e16.8,L16,e16.8)') label_chisq(isq), chisq(isq) , chisq_old(isq) , delta_chisq(isq), accept(isq) , temp_rmc(isq)
    enddo
    write(stdout,'(a,3e16.8,L16,e16.8)') label_chisq(6), chisq(6) , chisq_old(6) , delta_chisq(6) , accept(isq)
    write(stdout,'(a)') ''
#endif


  ! ===================================================================
  !   main RMC loop
  ! ===================================================================
  kall = 1 
  kacc = 1
  kprint = 1
#ifdef MPI
  timing_loop1 = MPI_WTIME(ierr)
  exec_loops = timing_loop1 - timing_all_start
#endif
  if ( ionode ) then
    lseparator(stdout)
    write( stdout , '(a,a10,8(a16))' )       'accepted  chi : ',label_chisq(6),(label_chisq(isq),isq=1,5),' delta_tot ',' rate '
    lseparator(stdout)
    write( stdout       ,200) kacc,' : ', chisq(6), (chisq(isq) , isq=1,5 ) , delta_chisq(6), REAL(kacc,kind=dp)/REAL(kall,kind=dp),' cpu : ',exec_loops
    write( kunit_RMCLOG ,200) kacc,' : ', chisq(6), (chisq(isq) , isq=1,5 ) , delta_chisq(6), REAL(kacc,kind=dp)/REAL(kall,kind=dp),' cpu : ',exec_loops
  endif 

  rmcloop : do while ( kacc .lt. acceptance_max ) 
    
    kall = kall + 1
    ! =============================================
    !  init all chi^2 quantities 
    ! =============================================
    chisq = 0.0_dp

    ! ==========================================
    !          move random particule 
    ! ==========================================
    CALL RANDOM_NUMBER(HARVEST = randpart )
    random_particule = int ( randpart * natm ) + 1
    CALL rmc_move(  max_displacement , random_particule , x , y , z ) 

    ! ===========================================
    !  g(r) calcualtion   
    ! ===========================================
    if ( lrmc_gr ) CALL rmc_gr ( grr_calc )

#ifdef debug
  do igr=0,nbins-1
    rr  = ( REAL( igr , kind = dp ) + 0.5_dp ) * resolution_gr
    write(*,'(3e20.8)') rr, grr_calc ( 0 , igr ) , grr_exp ( 0 , igr )
  enddo
#endif

    ! =============================================
    !  chi² g(r)
    ! =============================================
    CALL eval_chisq_rmc ( grr_exp , grr_calc , chisq(1) )

    ! =============================================
    !  chi² EFG : U_i Vzz eta
    ! =============================================
    if ( ( lrmc_u .or. lrmc_vzz .or. lrmc_eta ) .and. kacc .ge. efg_start ) then 
     ! ===========================================
     !  EFG calculation   
     ! ===========================================
     CALL rmc_efg ( dibU_calc , dibeta_calc , dibvzz_calc )
     if ( lrmc_u )   CALL eval_chisq_Uefg ( dibU_exp   , dibU_calc   , chisq(3) )
     if ( lrmc_vzz ) CALL eval_chisq_vzz  ( dibvzz_exp , dibvzz_calc , chisq(4) )
     if ( lrmc_eta ) CALL eval_chisq_eta  ( dibeta_exp , dibeta_calc , chisq(5) )
    endif
 
    ! ===========================================
    ! chi_var (i.e RMC-INVERT) g(r)
    ! ===========================================
    if ( lrmc_var )  CALL eval_chisq_var_distance (  chisq(2) )
    !if ( lrmc_var )  CALL eval_chisq_var_distance_simple_update ( random_particule , chisq(2) )
    chisq(6)     = sum( chisq(1:5) )    ! total
    delta_chisq = chisq - chisq_old

    ! ===============================
    ! metropolis-hastings algorithm
    ! ===============================
    accept=.FALSE.
    do isq = 1 , 6
      if ( delta_chisq(isq) .lt. 0.0_dp ) then
        accept(isq) = .TRUE.
      else
        metro = EXP ( - delta_chisq(isq) * 0.5_dp / temp_rmc(isq) )
        CALL RANDOM_NUMBER(HARVEST = randmetro )
        if ( metro .ge. randmetro ) then
          accept(isq) = .TRUE.
        else
          accept(isq) = .FALSE.
        endif
      endif
    enddo

#ifdef debug_chisq
    write(stdout,'(a)') '-------------------------------------------------------------------------------------'
    write(stdout,'(a)') 'label     current          old             delta            metropolis         "temp"'
    write(stdout,'(a)') '-------------------------------------------------------------------------------------'
    do isq=1,5
      write(stdout,'(a,3e16.8,L16,e16.8)') label_chisq(isq), chisq(isq) , chisq_old(isq) , delta_chisq(isq), accept(isq) , temp_rmc(isq)
    enddo
    write(stdout,'(a,3e16.8,L16,e16.8)') label_chisq(6), chisq(6) , chisq_old(6) , delta_chisq(6) , accept(isq)
    write(stdout,'(a)') ''
#endif

    !!!! ========== !!!!
    !!!!   succes   !!!!
    !!!! ========== !!!!
    acc : if ( accept(6) ) then

      kacc = kacc + 1
      chisq_old = chisq
      dab = dab_updated

      ! =========
      !   output
      ! =========
      if ( ionode .and. mod( kacc , nprint) .eq. 0.0_dp ) then

        ! print header each 50 accpeted configuration
        kprint = kprint + 1 
        if ( ionode .and. mod( kprint , 50) .eq. 0.0_dp ) then
          lseparator(stdout)
          write( stdout , '(a,a10,8(a16))' )       'accepted  chi : ',label_chisq(6),(label_chisq(isq),isq=1,5),' delta_tot ',' rate '
          lseparator(stdout)
          lseparator(kunit_RMCLOG)
          write( kunit_RMCLOG , '(a,a10,8(a16))' ) 'accepted  chi : ',label_chisq(6),(label_chisq(isq),isq=1,5),' delta_tot ',' rate '
          lseparator( kunit_RMCLOG )
        endif
 
        ! print to standard output and log file
#ifdef MPI
        timing_loop2 = MPI_WTIME(ierr)
        exec_loops = timing_loop2 - timing_loop1
        timing_loop1 = MPI_WTIME(ierr)
#endif
        write( stdout       ,200) kacc,' : ', chisq(6), (chisq(isq) , isq=1,5 ) , delta_chisq(6), REAL(kacc,kind=dp)/REAL(kall,kind=dp),' cpu : ',exec_loops
        write( kunit_RMCLOG ,200) kacc,' : ', chisq(6), (chisq(isq) , isq=1,5 ) , delta_chisq(6), REAL(kacc,kind=dp)/REAL(kall,kind=dp),' cpu : ',exec_loops

        ! print distribution files
        ! GRTFF

!        if ( lrmc_gr ) then
          CALL rmc_gr ( grr_calc )
          OPEN(UNIT=kunit_GRTFF,FILE='GRTFF.calc')
          do igr=0, nbins-1
            rr  = ( REAL( igr , kind = dp ) + 0.5_dp ) * resolution_gr
#ifdef GFORTRAN
            WRITE ( FMT , * ) npairs+2
            WRITE ( kunit_GRTFF ,'('// ADJUSTL(FMT) //'e20.10)') rr , ( grr_calc ( mp , igr ) , mp = 0 , npairs )
#else
            WRITE ( kunit_GRTFF ,'(<npairs+2>e20.10)') rr , ( grr_calc ( mp , igr ) , mp = 0 , npairs )
#endif
          enddo
          CLOSE(UNIT=kunit_GRTFF)
!        endif

        ! EFG distribution related 
        if (  ( lrmc_u .or. lrmc_vzz .or. lrmc_eta ) .and. kacc .ge. efg_start ) then
          ! DTIBUFF
          if ( lrmc_u ) then
            OPEN(UNIT=kunit_DTIBUFF,FILE='DTIBUFF.calc')
            do bin=0, PANU
#ifdef GFORTRAN
              WRITE ( FMT , * ) 6*ntype+7
              WRITE (kunit_DTIBUFF ,'('// ADJUSTL(FMT) //'f15.8)') u_efg_min  + REAL ( bin * resolution_u_efg , kind = dp ) , ( ( dibU_calc(ui,it,bin) , ui=1,6) , it=0,ntype )
#else
              WRITE (kunit_DTIBUFF ,'(<6*ntype+7>f15.8)') u_efg_min  + REAL ( bin * resolution_u_efg , kind = dp ) , ( ( dibU_calc(ui,it,bin) , ui=1,6) , it=0,ntype )
#endif
            enddo
            CLOSE(UNIT=kunit_DTIBUFF)
          endif
          ! DTVZZFF
          if ( lrmc_vzz ) then
            OPEN(UNIT=kunit_DTVZZFF,FILE='DTVZZFF.calc')
            do bin=0, PANvzz
#ifdef GFORTRAN
              WRITE ( FMT , * ) 6*ntype+7
              WRITE (kunit_DTVZZFF ,'('// ADJUSTL(FMT) //'f15.8)') vzz_efg_min  + REAL ( bin * resolution_vzz_efg , kind = dp ) , ( dibvzz_calc(it,bin) , it=0,ntype )
#else
              WRITE (kunit_DTVZZFF ,'(<6*ntype+7>f15.8)') vzz_efg_min  + REAL ( bin * resolution_vzz_efg , kind = dp ) , ( dibvzz_calc(it,bin) , it=0,ntype )
#endif
            enddo
            CLOSE(UNIT=kunit_DTVZZFF)
          endif
          ! DTETAFF
          if ( lrmc_eta ) then
            OPEN(UNIT=kunit_DTETAFF,FILE='DTETAFF.calc')
#ifdef GFORTRAN
            WRITE ( FMT , * ) ntype+2
            WRITE (kunit_DTETAFF,'('// ADJUSTL(FMT) //'f15.8)')  dzero, ( dzero , it = 0 , ntype )
#else
            WRITE (kunit_DTETAFF,'(<ntype+2>f15.8)')  dzero, ( dzero , it = 0 , ntype )
#endif
            do bin=0, PANeta-1
#ifdef GFORTRAN
            WRITE ( FMT , * ) 6*ntype+7 
            WRITE (kunit_DTETAFF ,'('// ADJUSTL(FMT) //'f15.8)') REAL ( bin * resolution_eta_efg , kind = dp ) , ( dibeta_calc(it,bin), it=0,ntype )
#else
            WRITE (kunit_DTETAFF ,'(<6*ntype+7>f15.8)') REAL ( bin * resolution_eta_efg , kind = dp ) , ( dibeta_calc(it,bin), it=0,ntype )
#endif
            enddo
            CLOSE(UNIT=kunit_DTETAFF) 
          endif
        endif 
        CALL write_CONTFF
        if ( kacc .gt. (acceptance_max-10000) ) CALL write_trajff_xyz 
      endif

    else
    !!!! ============================================= !!!!
    !!!!  try again we get back to the previous config !!!!
    !!!! ============================================= !!!!
#ifdef debug
      write( stdout ,970 ) "debug : before rejecting",random_particule,rx(random_particule),ry(random_particule),rz(random_particule)
#endif
      CALL kardir ( natm , rx , ry , rz , simu_cell%B )
      rx(random_particule) = rx(random_particule) - x 
      ry(random_particule) = ry(random_particule) - y 
      rz(random_particule) = rz(random_particule) - z
      CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
#ifdef debug
      write( stdout , 970 ) "debug : after rejecting",random_particule,rx(random_particule),ry(random_particule),rz(random_particule)
#endif

     endif acc

  enddo rmcloop
  ! =====================
  !  end of main loop
  ! =====================

  CLOSE (unit = kunit_TRAJFF )
  CLOSE (unit = kunit_RMCLOG )

970 FORMAT(a,i6,3f12.6)
200 FORMAT(i9,a3,8e16.8,a7,f6.2)
310 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,' min =  ',F10.5,' max = ',F10.5)    


END SUBROUTINE rmc_main

! *********************** SUBROUTINE rmc_print_info *****************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE eval_chisq_rmc ( grr_exp , grr_calc , chisq_rmc )

  USE radial_distrib,    ONLY :  nbins, npairs

  implicit none

  integer :: mp , bin
  real(kind=dp) :: chisq_rmc
  real(kind=dp) , dimension ( 0 :npairs , 0: nbins-1) :: grr_exp, grr_calc
 
  chisq_rmc = 0.0_dp
  do bin=0, nbins-1
    do mp=0, npairs
      chisq_rmc = chisq_rmc + ( grr_calc( mp , bin ) - grr_exp( mp , bin ) )**2.0_dp
    enddo
  enddo
  
  return

END SUBROUTINE eval_chisq_rmc

! *********************** SUBROUTINE rmc_print_info ****************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE eval_chisq_Uefg ( dibU_exp , dibU_calc , chisq_uefg )

  USE efg,    ONLY :  PANU
  USE config,    ONLY :  ntype

  implicit none

  integer :: mp , bin, ui
  real(kind=dp) :: chisq_uefg
  real(kind=dp) , dimension ( 6 , 0:ntype , 0:PANU ) :: dibU_exp, dibU_calc

  chisq_uefg = 0.0_dp
  do bin=0, PANU
    do mp=0, ntype
      do ui =1 , 5 
        chisq_uefg = chisq_uefg + ( dibU_calc( ui , mp , bin ) - dibU_exp( ui , mp , bin ) )**2.0_dp
      enddo
    enddo
  enddo

  return

END SUBROUTINE eval_chisq_Uefg

! *********************** SUBROUTINE rmc_print_info ****************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE eval_chisq_vzz ( dibvzz_exp , dibvzz_calc , chisq_vzz )

  USE efg,    ONLY :  PANvzz
  USE config,    ONLY :  ntype

  implicit none

  integer :: mp , bin
  real(kind=dp) :: chisq_vzz
  real(kind=dp) , dimension ( 0:ntype , 0:PANvzz ) :: dibvzz_exp, dibvzz_calc

  chisq_vzz = 0.0_dp
  do bin=0, PANvzz
    do mp=0, ntype
      chisq_vzz = chisq_vzz + ( dibvzz_calc( mp , bin ) - dibvzz_exp( mp , bin ) )**2.0_dp
    enddo
  enddo

  return

END SUBROUTINE eval_chisq_vzz

! *********************** SUBROUTINE rmc_print_info ****************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE eval_chisq_eta ( dibeta_exp , dibeta_calc , chisq_eta )

  USE efg,    ONLY :  PANeta
  USE config,    ONLY :  ntype

  implicit none

  integer :: mp , bin
  real(kind=dp) :: chisq_eta
  real(kind=dp) , dimension ( 0:ntype , 0:PANeta ) :: dibeta_exp, dibeta_calc

  chisq_eta = 0.0_dp
  do bin=0, PANeta
    do mp=0, ntype
      chisq_eta = chisq_eta + ( dibeta_calc( mp , bin ) - dibeta_exp( mp , bin ) )**2.0_dp
    enddo
  enddo

  return

END SUBROUTINE eval_chisq_eta



! *********************** SUBROUTINE eval_chisq_var_distance *******************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE eval_chisq_var_distance ( chisq_var )

  USE control,          ONLY :  myrank
  USE config,           ONLY :  natm, ntype , itype , natmi , rx , ry , rz, simu_cell, rx, ry, rz 
  USE cell,             ONLY :  kardir, dirkar
  USE time,             ONLY :  chisqvartimetot , chisqvartime2 , chisqvartime3 , chisqvartime4 , chisqvartime5 

  implicit none

  ! global
  real(kind=dp) ,intent (out) :: chisq_var

  ! local
  integer       :: nmax
  integer       :: ia, ja , it , n , it1 ,it2, ierr 
  real(kind=dp) :: rijsq, d 
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij

  real(kind=dp) , dimension ( :,:,: ) , allocatable :: dab_sorted
  real(kind=dp) , dimension ( :,:,: ) , allocatable :: dab_sorted_shifted
  real(kind=dp) , dimension ( :,:,: ) , allocatable :: dab_aver
  integer       , dimension ( :,: )   , allocatable :: nm 
  integer       , dimension ( :,: )   , allocatable :: np 
  real(kind=dp)     :: ttt1 , ttt2 , ttt3 ,ttt4, ttt5
!  real(kind=dp)     :: ttt2_1
  ! merging working array
  real(kind=dp) , dimension (:)       , allocatable :: tab
  integer       , dimension (:)       , allocatable :: lab
  integer       , dimension (:)       , allocatable :: ltab

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif

  ! nm(it,ia) gives the number of neigbours of type 
  ! it for an atom ia 
  ! np(it,ia) gives nmax - nm(it,ia)
  ! this will be useful to reshape arrays after sorting.
  nmax = maxval(natmi(1:))
  allocate ( nm ( ntype, natm ) )      
  allocate ( np ( ntype, natm ) )      

  do ia = 1 , natm
    do it = 1 , ntype
      if ( it .eq. itype(ia) ) then
        nm(it,ia) = natmi ( itype(ia) ) - 1
        np(it,ia) = nmax - nm(it,ia)
      else
        nm(it,ia) = natmi ( it ) 
        np(it,ia) = nmax - nm(it,ia)
      endif
    enddo
  enddo

#ifdef debug_invert
   write( stdout, '(a)' ) 'debug : nm,np'
   do ia = 1 , natm
    do it = 1 , ntype
      write( stdout, '(a,2i)') 'debug : ', nm(it,1),np(it,1)
      write( stdout, '(a,2i)') 'debug : ', nm(it,1500),np(it,1500)
    enddo
   enddo
#endif

  allocate ( dab_sorted ( nmax, ntype , natm ) ) 
  allocate ( dab_sorted_shifted ( nmax, ntype , natm ) ) 
  allocate ( dab_aver ( nmax , ntype, ntype ) ) 
  allocate ( tab ( ( nmax + 1) / 2 ) ) 
  allocate ( lab ( nmax ) ) 
  allocate ( ltab( ( nmax + 1) / 2 ) ) 
  dab        = 0.0_dp
  dab_sorted = 0.0_dp
  dab_sorted_shifted = 0.0_dp
  dab_aver   = 0.0_dp
  tab=0.0_dp
  lab=0
  ltab=0

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = 1 , natm 
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    do it = 1 , ntype
      n = 0
      do ja = 1 , natm  
        if ( ja .ne. ia ) then
          if ( itype(ja) .eq. it ) then
             n = n + 1
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
             d = sqrt(rijsq)
             dab( n , it , ia ) = d
          endif
        endif
      enddo
    enddo
  enddo

!  do it = 1, ntype
!    do n=1,nmax
!      CALL MPI_ALL_REDUCE_DOUBLE ( dab( n , it , : ) , natm ) 
!    enddo
!  enddo
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  chisqvartime2 = chisqvartime2 + ( ttt2 - ttt1 )
#endif
#ifdef debug_invert
  print*,'debug : avant sorting dab'      
  do ia=1, natm
    do it=1, ntype
      write ( stdout , '(a)' )    '                ia       it      nm(it,ia)'
      write ( stdout , '(a,3i)' ) '          ',ia,it,nm (it,ia)
      do n=1, nmax
        write( stdout, '(a,i,f)') 'debug : ',myrank,dab( n , it , ia )
      enddo
    enddo
  enddo
#endif 

  ! ============================================
  ! merge sort of dab array => dab_sorted
  ! dab is conserved because it will be just 
  ! updated at random_particule  
  ! ============================================
  dab_sorted = dab
  dab_aver(:,:,:) = 0.0_dp
  do ia = 1 , natm
    do it = 1 , ntype
      call merge_sort ( dab_sorted (:,it,ia) , nmax , tab , lab , ltab )
      !call QsortC( dab_sorted (:,it,ia) )
      dab_sorted_shifted(1:nm(it,ia),it,ia) = dab_sorted (np(it,ia)+1:nmax,it,ia)
      do it2 = 1, ntype
         if (itype(ia).eq.it2 ) then
           do n=1,nmax
             if (dab_sorted_shifted( n , it , ia ) .ne. 0.0_dp) then
               dab_aver(n,it,it2) = dab_aver(n,it,it2) + dab_sorted_shifted (n,it,ia)
             endif
           enddo
         endif
      enddo
    enddo
  enddo

!   do it = 1, ntype
!    do n=1,nmax
!      CALL MPI_ALL_REDUCE_DOUBLE ( dab_sorted_shifted ( n , it , : ) , natm )
!    enddo
!  enddo
!#ifdef debug_invert
#ifdef MPI
  ttt3 = MPI_WTIME(ierr)
  chisqvartime3 = chisqvartime3 + ( ttt3 - ttt2 )
#endif
!#endif

  ! we do not need this anymore
  deallocate ( tab  ) 
  deallocate ( lab  ) 
  deallocate ( ltab ) 

#ifdef debug_invert
  write (stdout , '(a)' ) 'debug : apres sorting'      
  do ia=1, natm
    do it=1, ntype
      write ( stdout , '(a)' )    '          ia,it,nm (it,ia)'
      write ( stdout , '(a,3i)' ) 'debug : ',ia,it,nm (it,ia)
      do n=1, nm(it,ia) 
       write ( stdout , '(a,2f)' ) 'debug : ' , dab_sorted( n , it , ia ) , dab_sorted_shifted ( n , it , ia ) 
      enddo
    enddo
  enddo
#endif

  dab_sorted = 0.0_dp
  dab_sorted = dab_sorted_shifted
  deallocate ( dab_sorted_shifted ) 

!  dab_aver(:,:,:) = 0.0_dp
!  do it2=1, ntype
!    do ia=1, natm
!      do it1=1, ntype 
!        if (itype(ia).eq.it2 ) then
!          do n=1,nmax
!            if (dab_sorted( n , it1 , ia ) .ne. 0.0_dp) then
!              dab_aver(n,it1,it2) = dab_aver(n,it1,it2) + dab_sorted(n,it1,ia)
!            endif
!          enddo
!        endif
!      enddo
!    enddo
!  enddo

! removed to be included directly is the chi loops
!  do it1=1,ntype
!    do it2=1,ntype
!      do n=1,nmax
!        dab_aver(n,it2,it1) = dab_aver(n,it2,it1) / REAL(natmi(it1), kind=dp)
!      enddo
!    enddo
!  enddo 

#ifdef debug_invert
  write(stdout, '(a)') 'debug : apres aver (not normalized => in the loop where chi is calcualated) '      
  do it2=1, ntype 
    do it1=1, ntype
      write ( stdout , '(a)' )    '          it2,it1'
      write ( stdout , '(a,2i)' ) 'debug : ',it2,it1
      do n=1,nmax
        if ( dab_aver(n,it2,it1).ne. 0.0_dp ) then
          write(stdout, '(a,f)' ) 'debug : ',dab_aver( n , it2 , it1 ) 
        endif
      enddo
    enddo
  enddo

  print*,'---------------------------------------------------------------------------------'
  print*,'it = 1'
  print*,'ia = 1'
  print*,dab_sorted(:,1,1)
  print*,'ia = 2'
  print*,dab_sorted(:,1,2)
  print*,'ia = 3'
  print*,dab_sorted(:,1,3)
  print*,'ia = 4'
  print*,dab_sorted(:,1,4)
  print*,'ia = 5'
  print*,dab_sorted(:,1,5)
  print*,'ia = 6'
  print*,dab_sorted(:,1,6)
  print*,'average it = 1-1'
  print*,dab_aver(:,1,1) / REAL(natmi(1), kind=dp)
  print*,'---------------------------------------------------------------------------------'

  do ia=1,natm
    it1 = itype(ia)
    do it2=1,ntype
      do n=1,nmax
        if ( dab_sorted(n,it2,ia) .ne. 0.0_dp) then
          write(stdout, '(a,2f)' ) 'debug :', dab_sorted(n,it2,ia),dab_aver(n,it2,it1)
        endif
      enddo
    enddo
  enddo
#endif
#ifdef MPI
  ttt4 = MPI_WTIME(ierr)
  chisqvartime4 = chisqvartime4 + ( ttt4 - ttt3 )
#endif
!#endif

  chisq_var = 0.0_dp
 do ia=1,natm
    it1 = itype(ia)
    do it2=1,ntype
      do n=1,nmax
        if (dab_aver(n,it2,it1).ne.0.0_dp) then 
          chisq_var = chisq_var + ( dab_sorted (n,it2,ia) - (dab_aver(n,it2,it1))/ REAL(natmi(it1), kind=dp) ) ** 2.0_dp / (dab_aver(n,it2,it1)/REAL(natmi(it1), kind=dp))**2.0_dp 
        endif
      enddo
    enddo
  enddo

  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  deallocate ( nm )
  deallocate ( dab_aver  ) 
  deallocate ( dab_sorted  ) 

#ifdef MPI
  ttt5 = MPI_WTIME(ierr)
!#ifdef debug_invert
  chisqvartime5 = chisqvartime5 + ( ttt5 - ttt4 )
!#endif
  chisqvartimetot = chisqvartimetot + ( ttt5 - ttt1 ) 
#endif

  return

END SUBROUTINE eval_chisq_var_distance

! *********************** SUBROUTINE eval_chisq_var_distance_simple_update *****
!
!> \brief
!
! ******************************************************************************
SUBROUTINE eval_chisq_var_distance_simple_update ( random_particule , chisq_var )
 
  USE config,           ONLY :  natm, ntype , itype , natmi , rx , ry , rz, simu_cell, rx, ry, rz 
  USE cell,             ONLY :  kardir, dirkar
  USE time,             ONLY :  chisqvartimetotu, chisqvartime2u , chisqvartime3u , chisqvartime4u

  implicit none

  ! global
  integer       ,intent (in)  :: random_particule
  real(kind=dp) ,intent (out) :: chisq_var

  ! local
  integer            :: ia, ja , it , n , it1 ,it2, ierr
  integer            :: nmax
  real(kind=dp)      :: rijsq, d 
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: ttt1!,ttt2 
!  real(kind=dp) :: ttt3,ttt4
  real(kind=dp) :: ttt5!,ttt6

  real(kind=dp) , dimension ( :,:,: ) , allocatable :: dab_sorted
  real(kind=dp) , dimension ( :,:,: ) , allocatable :: dab_sorted_shifted
  real(kind=dp) , dimension ( :,:,: ) , allocatable :: dab_aver
  integer       , dimension ( :,: )   , allocatable :: nm 
  integer       , dimension ( :,: )   , allocatable :: np 
  ! merging working array
!  real(kind=dp) , dimension ( (nmax + 1) / 2)       :: tab
!  integer , dimension ( nmax )                      :: lab
!  integer , dimension ( (nmax + 1) / 2)             :: ltab
  real(kind=dp) , dimension (:)       , allocatable  :: tab
  integer       , dimension (:)       , allocatable  :: lab
  integer       , dimension (:)       , allocatable  :: ltab
  real(kind=dp) , external :: QuickSort

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif

  ! nm(it,ia) gives the number of neigbours of type 
  ! it for an atom ia 
  ! np(it,ia) gives nmax - nm(it,ia)
  ! this will be useful to reshape arrays after sorting.
  nmax = maxval(natmi(1:))
  allocate ( nm ( ntype, natm ) )      
  allocate ( np ( ntype, natm ) )      

  do ia = 1 , natm
    do it = 1 , ntype
      if ( it .eq. itype(ia) ) then
        nm(it,ia) = natmi ( itype(ia) ) - 1
        np(it,ia) = nmax - nm(it,ia)
      else
        nm(it,ia) = natmi ( it ) 
        np(it,ia) = nmax - nm(it,ia)
      endif
    enddo
  enddo

#ifdef debug_invert
   write( stdout, '(a)' ) 'debug : nm,np'
   do ia = 1 , natm
    do it = 1 , ntype
      write( stdout, '(a,2i)') 'debug : ', nm(it,1),np(it,1)
      write( stdout, '(a,2i)') 'debug : ', nm(it,1500),np(it,1500)
    enddo
   enddo
#endif

 
  allocate ( dab_sorted ( nmax, ntype , natm ) ) 
  allocate ( dab_sorted_shifted ( nmax, ntype , natm ) ) 
  allocate ( dab_aver ( nmax , ntype, ntype ) ) 
  allocate ( tab ( ( nmax + 1) / 2 ) ) 
  allocate ( lab ( nmax ) ) 
  allocate ( ltab( ( nmax + 1) / 2 ) ) 
  tab=0.0_dp
  lab=0
  ltab=0
  dab_updated = 0.0_dp
  dab_sorted = 0.0_dp
  dab_sorted_shifted = 0.0_dp

  dab_updated = dab 

  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  ! this works !!!
  do ia = 1 , natm
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    do it = 1 , ntype
      n = 0
      do ja = 1 , natm
        if ( ja .ne. ia ) then
          if ( itype(ja) .eq. it ) then
            n = n + 1
            ! only update ia , ja == random_particule
            if ( ia .eq. random_particule .or. ja .eq. random_particule ) then 
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
              d = sqrt(rijsq)
              dab_updated ( n , it , ia ) = d
            endif
          endif
        endif
      enddo
    enddo
  enddo

  ! merging 
!  do it = 1, ntype
!    do n=1,nmax
!      CALL MPI_ALL_REDUCE_DOUBLE ( dab_updated( n , it , : ) , natm )
!    enddo
!  enddo


#ifdef debug_invert
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  chisqvartime2u = chisqvartime2u + ( ttt2 - ttt1 ) 
#endif
  write(stdout ,'(a)' ) 'debug : avant sorting dab'      
  do ia=1, natm
    do it=1, ntype
      print*,ia,it,nm (it,ia)
      do n=1, nmax
        print*,dab( n , it , ia ),dab_updated(n,it,ia)
      enddo
    enddo
  enddo
#endif 

  !merge sort seems faster
  dab_sorted = dab_updated
  dab_sorted_shifted = 0.0_dp
  dab_aver(:,:,:) = 0.0_dp
  do ia = 1 , natm
    do it = 1 , ntype
!      call merge_sort ( dab_sorted (:,it,ia) , nmax , tab , lab , ltab )
      dab_sorted_shifted(1:nm(it,ia),it,ia) = dab_sorted (np(it,ia)+1:nmax,it,ia)
      do it2 = 1, ntype
         if (itype(ia).eq.it2 ) then
           do n=1,nmax 
             if (dab_sorted_shifted( n , it , ia ) .ne. 0.0_dp) then           
               dab_aver(n,it,it2) = dab_aver(n,it,it2) + dab_sorted_shifted (n,it,ia)  
             endif
           enddo
         endif
      enddo
    enddo
  enddo
 
  deallocate ( dab_sorted )
  deallocate ( tab  ) 
  deallocate ( lab  ) 
  deallocate ( ltab ) 


#ifdef debug_invert
#ifdef MPI
  ttt3 = MPI_WTIME(ierr)
  chisqvartime3u = chisqvartime3u + ( ttt3 - ttt2 ) 
#endif
  print*,'apres sorting'
  do ia=1, natm
    do it=1, ntype
      print*,ia,it,nm (it,ia)
      do n=1, nm(it,ia)
        print*,n,dab_sorted( n , it , ia ) , dab_sorted_shifted ( n , it , ia )
      enddo
    enddo
  enddo
#endif
! -------------------------------------------------------------------
  ! ========================================================
  ! merveilleuse fonction de fortran : eoshift 
  ! sort put zeros at the beginning of the array
  ! we shift to the left until the first index value 
  ! is not null ( commented out )
  !
  ! update : the reshaping version is clearly more 
  ! efficient for large problems
  ! ========================================================
!  do ia=1, natm
!    do it=1, ntype
!      do while ( dab_sorted (1,it,ia) .eq. 0.0_dp  )
!        dab_sorted (:,it,ia) = EOSHIFT(dab_sorted (:,it,ia), np(it,ia), BOUNDARY=0.0_dp, DIM=1)
!         dab_sorted_shifted(1:nm(it,ia),it,ia) = dab_sorted (np(it,ia)+1:nmax,it,ia)
!         print*,'sorted',dab_sorted(:,it,ia)
!         print*,'shifted',dab_sorted_shifted(:,it,ia)
!      enddo 
!    enddo
!  enddo
! -------------------------------------------------------------------
!  dab_sorted = 0.0_dp
!  dab_sorted = dab_sorted_shifted
!  deallocate ( dab_sorted_shifted ) 
!  ttt6 = MPI_WTIME(ierr)
!  chisqvartime5 = chisqvartime5 + ( ttt6 - ttt2 ) 

!  dab_aver(:,:,:) = 0.0_dp
!  do it2=1, ntype
!    do ia=1, natm
!      do it1=1, ntype 
!        if (itype(ia).eq.it2 ) then
!          do n=1,nmax
!            if (dab_sorted( n , it1 , ia ) .ne. 0.0_dp) then
!              dab_aver(n,it1,it2) = dab_aver(n,it1,it2) + dab_sorted (n,it1,ia)
!            endif
!          enddo
!        endif
!      enddo
!    enddo
!  enddo

#ifdef debug_invert
  print*,'apres aver (not normalized => in the loop where chi is calcualated) '      
  do it2=1, ntype 
    do it1=1, ntype
      print*,it2,it1
      do n=1,nmax
        if ( dab_aver(n,it2,it1).ne. 0.0_dp ) then
          print*,dab_aver( n , it2 , it1 ) 
        endif
      enddo
    enddo
  enddo

  print*,'---------------------------------------------------------------------------------'
  print*,'it = 1'
  print*,'ia = 1'
  print*,dab_sorted_shifted (:,1,1)
  print*,'ia = 2'
  print*,dab_sorted_shifted (:,1,2)
  print*,'ia = 3'
  print*,dab_sorted_shifted (:,1,3)
  print*,'ia = 4'
  print*,dab_sorted_shifted (:,1,4)
  print*,'ia = 5'
  print*,dab_sorted_shifted (:,1,5)
  print*,'ia = 6'
  print*,dab_sorted_shifted (:,1,6)
  print*,'average it = 1-1'
  print*,dab_aver(:,1,1)/REAL(natmi(1), kind=dp)
  print*,'---------------------------------------------------------------------------------'

  do ia=1,natm
    it1 = itype(ia)
    do it2=1,ntype
      do n=1,nmax
        if ( dab_sorted_shifted (n,it2,ia) .ne. 0.0_dp) then
          print*,ia,it2,dab_sorted_shifted (n,it2,ia),dab_aver(n,it2,it1)
        endif
      enddo
    enddo
  enddo
#endif

  chisq_var = 0.0_dp
  do ia=1,natm
    it1 = itype(ia)
    do it2=1,ntype
      do n=1,nmax
        if (  dab_sorted_shifted (n,it2,ia) .ne. 0.0_dp) then
          chisq_var = chisq_var + ( dab_sorted_shifted (n,it2,ia) - ( dab_aver(n,it2,it1) / REAL(natmi(it1), kind=dp) )  ) ** 2.0_dp / ( dab_aver(n,it2,it1) / REAL(natmi(it1), kind=dp) ) **2.0_dp 
        endif
      enddo
    enddo
  enddo
#ifdef debug_invert
#ifdef MPI
  ttt4 = MPI_WTIME(ierr)
  chisqvartime4u = chisqvartime4u + (ttt4 - ttt3) 
#endif
#endif

  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  deallocate ( nm )
  deallocate ( np )
  deallocate ( dab_aver  ) 
  deallocate ( dab_sorted_shifted  ) 

#ifdef MPI
  ttt5 = MPI_WTIME(ierr)
  chisqvartimetotu = chisqvartimetotu + ( ttt5 - ttt1 ) 
#endif

  return

END SUBROUTINE eval_chisq_var_distance_simple_update

! *********************** SUBROUTINE rmc_move **********************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE rmc_move ( max_displacement , random_particule , x1 , y1 , z1 )

  USE config,           ONLY :  natm , simu_cell , rx , ry ,rz
  USE cell,             ONLY :  kardir, dirkar

  implicit none

  integer :: random_particule
  real(kind=dp) :: x ,  y  , z 
  real(kind=dp) :: x1 , y1 , z1 
  real(kind=dp) :: max_displacement 

#ifdef debug    
    write(*,970) "debug : before moving",random_particule,rx(random_particule),ry(random_particule),rz(random_particule)
#endif

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  CALL RANDOM_NUMBER(HARVEST = x)
  CALL RANDOM_NUMBER(HARVEST = y)
  CALL RANDOM_NUMBER(HARVEST = z)
  x1 = ( max_displacement / simu_cell%ANORM(1) ) * ( 2.0_dp * x - 1.0_dp ) 
  y1 = ( max_displacement / simu_cell%ANORM(2) ) * ( 2.0_dp * y - 1.0_dp ) 
  z1 = ( max_displacement / simu_cell%ANORM(3) ) * ( 2.0_dp * z - 1.0_dp ) 
  rx(random_particule) = rx(random_particule) + x1
  ry(random_particule) = ry(random_particule) + y1
  rz(random_particule) = rz(random_particule) + z1

  if ( rx(random_particule) .gt. 1.0_dp   .or. & 
         rx(random_particule) .lt. 0.0_dp ) then 
     rx ( random_particule ) = rx ( random_particule ) - NINT ( rx ( random_particule ) )
  endif
  if ( ry(random_particule) .gt. 1.0_dp   .or. & 
         ry(random_particule) .lt. 0.0_dp ) then 
     ry ( random_particule ) = ry ( random_particule ) - NINT ( ry ( random_particule ) )
  endif
  if ( rz(random_particule) .gt. 1.0_dp   .or. & 
         rz(random_particule) .lt. 0.0_dp ) then 
     rz ( random_particule ) = rz ( random_particule ) - NINT ( rz ( random_particule ) )
  endif

  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

#ifdef debug    
    write(*,970) "debug : after moving",random_particule,rx(random_particule),ry(random_particule),rz(random_particule)
970  FORMAT(a,i6,3f12.6)
#endif

  return

END SUBROUTINE

! *********************** SUBROUTINE rmc_move **********************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE rmc_gr ( grr_calc ) 

  USE constants,                ONLY :  pi
  USE radial_distrib,           ONLY :  nbins , npairs , gr_main , gr , resg , cutgr
  USE config,                   ONLY :  natm, natmi , ntype , simu_cell
  USE time,                     ONLY :  rmcgrtimetot_comm

  implicit none
 
  ! global
  real(kind=dp) ::  grr_calc ( 0 : npairs , 0 : nbins-1 )

  ! local 
  integer :: igr , ierr, it1 , it2 , mp 
  real(kind=dp) :: ttt1,ttt2
  real(kind=dp) :: rr , vol


  grr_calc = 0.0_dp

  ! ===========================================
  !  g(r) calcualtion and array merging  
  ! ===========================================
  CALL gr_main 
#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif
  CALL MPI_ALL_REDUCE_INTEGER ( gr(:,0,0), nbins )
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      CALL MPI_ALL_REDUCE_INTEGER ( gr(:,it1,it2), nbins )
    enddo
  enddo
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  rmcgrtimetot_comm = rmcgrtimetot_comm + ( ttt2 - ttt1 )
#endif
  
  ! ============================================
  !        g(r) normalization
  ! ============================================
  do igr = 0 , nbins-1
    rr  = ( REAL( igr , kind = dp ) + 0.5_dp ) * resolution_gr
    vol = 4.0_dp * pi * ( resolution_gr * rr * rr + ( resolution_gr**3 ) / 12.0_dp)
    vol = vol / simu_cell%omega
    ! all - all pairs 
    grr_calc ( 0 , igr ) = REAL ( gr( igr , 0 , 0 ) , kind = dp ) / ( vol * natm * natm )
    ! type pairs
    mp = 1
    do it1 = 1 , ntype
      do it2 = it1 , ntype
        if ( mp .lt. 0 .and. mp .gt. npairs ) then
          WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
          STOP
        endif
        grr_calc ( mp , igr ) = REAL ( gr ( igr ,  it1 , it2 ) , kind = dp ) / REAL ( vol * natmi ( it1 ) * natmi ( it2 ) , kind = dp )
        mp = mp + 1
      enddo
    enddo
  enddo
  gr = 0

  return

END SUBROUTINE

! *********************** SUBROUTINE rmc_move **********************************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE rmc_efg ( dibU_calc , dibeta_calc , dibvzz_calc ) 

  USE config,                   ONLY :  natm, ntype , itype , natmi
  USE field,                    ONLY :  alphaES , lwfc, km_coul
  USE efg,                      ONLY :  PANU , PANeta, PANvzz , resu , reseta , resvzz , umin , vzzmin , &
                                        read_dtibuff , read_dtvzzff , read_dtetaff , dibUtot , dibetatot , dibvzztot , &
                                        mu , efg_alloc , efg_mesh_alloc , multipole_efg_es , efg_ia , nmr_convention

  implicit none

  ! global
  real(kind=dp),     dimension ( 6 , 0:ntype , 0:PANU )   :: dibU_calc
  real(kind=dp),     dimension (     0:ntype , 0:PANvzz ) :: dibvzz_calc
  real(kind=dp),     dimension (     0:ntype , 0:PANeta ) :: dibeta_calc

  ! local
  integer :: bin
  integer :: ia, it, ui, ku, keta, kvzz
  real(kind=dp) :: uk,vzzk,etak
  real(kind=dp),     dimension (     : , : ) , allocatable :: U
  real(kind=dp),     dimension (     : , : ) , allocatable :: nmr
  real(kind=dp)     :: sq3 , sq32 

  ! variables DSYEV
  integer :: ifail
  integer, parameter                             :: lwork = 6
  real(kind=dp)                                  :: w(3)
  real(kind=dp)                                  :: work(3 * lwork)


  allocate ( U        ( natm , 5 ) )                       ! Czjzek U component
  allocate ( nmr      ( natm , 4 ) )

  sq3 = SQRT ( 3.0_dp )
  sq3 = 1.0_dp / sq3
  sq32 = sq3 * 0.5_dp

  CALL multipole_efg_ES ( km_coul , alphaES , mu )

  do ia = 1, natm

    ! =================
    !  diagonalisation
    ! =================
    CALL DSYEV ( 'N' , 'U' , 3 , efg_ia(ia,:,:) , 3 , w , work , 3 * lwork , ifail )
    if ( ifail .ne. 0 ) then
      io_node WRITE ( stderr , * ) 'ERROR: DSYEV, STOP in rmc MODULE'
      STOP
    endif
    CALL nmr_convention( w , nmr ( ia , : )  , ia )

    U ( ia , 1 ) = efg_ia ( ia, 3 , 3 ) * 0.5_dp
    U ( ia , 2 ) = efg_ia ( ia, 1 , 3 ) * sq3
    U ( ia , 3 ) = efg_ia ( ia, 2 , 3 ) * sq3
    U ( ia , 4 ) = efg_ia ( ia , 1 , 2 ) * sq3
    U ( ia , 5 ) = ( efg_ia ( ia, 1 , 1 ) - efg_ia ( ia, 2 , 2 ) ) * sq32
    ! ======================== 
    !  quadrupolar parameters
    ! ======================== 
    vzzk = (nmr(ia,3)-vzz_efg_min) / resolution_vzz_efg
    etak = nmr(ia,4) / resolution_eta_efg
    kvzz = int(vzzk) + 1
    keta = int(etak)
    ! ====================== 
    !  test out of bound
    ! ====================== 
    if ( kvzz .lt. 0 .or. kvzz .gt. PANvzz ) then
      cycle
      !if ( ionode ) &
      !WRITE ( stderr , '(a,a,i,f)' ) 'ERROR: out of bound distribvzz in rmc
      !MODULE ',atype(ia), kvzz, nmr(ia,3)
      !STOP
    endif
    if ( keta .lt. 0 .or. keta .gt. PANeta ) then
      cycle
      io_node WRITE ( stderr , '(a)' ) 'ERROR: out of bound distribeta in rmc MODULE '
      io_node WRITE ( stderr ,210) ia, keta, nmr(ia,4), nmr(ia,4), &
                                             nmr(ia,1), nmr(ia,2), &
                                             nmr(ia,3)
      STOP
    endif
    ! distribution per type 
    dibvzztot(0,kvzz) = dibvzztot(0,kvzz) + 1
    dibetatot(0,keta) = dibetatot(0,keta) + 1
    ! type specific 
    do it=1,ntype
      if ( itype(ia) .eq. it ) then
        dibvzztot(it,kvzz) = dibvzztot(it,kvzz) + 1
        dibetatot(it,keta) = dibetatot(it,keta) + 1
      endif
    enddo

    if ( lrmc_u ) then
      ! ========
      !    Ui
      ! ========
      do ui=1,5
        uk = (U(ia,ui) - u_efg_min )/resolution_u_efg
        ku = int(uk) + 1
        ! ====================== 
        !  test out of bound
        ! ====================== 
        if (ku.lt.0.or.ku.gt.PANU) then
          cycle
          !io_node WRITE ( stderr , * ) 'ERROR: out of bound dibU1'
          !io_node WRITE ( stderr ,310) ia,ku,U(ia,ui),umin,ABS (umin)
        endif
        ! all types
        dibUtot(ui,0,ku) = dibUtot(ui,0,ku) + 1
        ! type specific
        do it=1,ntype
          if (itype(ia).eq.it) then
            dibUtot(ui,it,ku) = dibUtot(ui,it,ku) + 1
          endif
        enddo
        ! average Uk,k>1
        if ( ui .ne. 1 ) then
          dibUtot(6,0,ku) = dibUtot(6,0,ku) + 1
          do it=1,2
            if (itype(ia).eq.it) then
              dibUtot(6,it,ku) = dibUtot(6,it,ku) + 1
            endif
          enddo
        endif
      enddo 
    endif
  enddo 
  if ( lrmc_u ) then
    do bin=0,PANU
      do ui=1,5
        dibU_calc(ui,0,bin)  = REAL( dibUtot(ui,0,bin),kind=dp) / ( resolution_u_efg * REAL( natm ,kind = dp ) )
        do it=1,ntype
          dibU_calc(ui,it,bin) = REAL( dibUtot(1,0,bin),kind=dp) / ( resolution_u_efg * REAL( natmi(it) ,kind = dp ) )
        enddo
      enddo
    enddo
  endif
  if ( lrmc_eta ) then
    do bin = 0 , PANeta-1
      dibeta_calc(0,bin) = REAL( dibetatot( 0 , bin ),kind=dp )  / ( resolution_eta_efg * REAL( natm , kind=dp ) )
      do it=1,ntype
        dibeta_calc(it,bin) = REAL( dibetatot( it , bin ),kind=dp )  / ( resolution_eta_efg * REAL( natmi(it) , kind=dp ) )
      enddo
    enddo
  endif
  if ( lrmc_vzz ) then
    do bin = 0 , PANvzz
      dibvzz_calc(0,bin) = REAL( dibvzztot( 0 , bin ),kind=dp )  / ( resolution_vzz_efg * REAL( natm , kind=dp ) )
      do it=1,ntype
        dibvzz_calc(it,bin) = REAL( dibvzztot( it , bin ),kind=dp )  / ( resolution_vzz_efg * REAL( natmi(it) , kind=dp ) )
      enddo
    enddo
  endif
  
  dibUtot = 0
  dibvzztot = 0
  dibetatot = 0

  return

210 FORMAT('atom = ',I6,' k = ',I12, 'value = ',F14.8,&
           ' min =  0.0_dp    max =  1.0_dp',' vaa{a = xx,yy,zz}, eta = ',4F14.8)

END SUBROUTINE


END MODULE rmc
! ===== fmV =====
