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
!#define debug_symmetry
! ======= Hardware =======

! *********************** MODULE vib *******************************************
!
!> \brief
!! phonons, band , dos related module
!
!> \todo
!! check if any shape of cell would work on this module
!
! ******************************************************************************
MODULE vib

  USE constants,                ONLY : dp
  USE mpimdff

  implicit none

  integer          :: nconf          !< number of configurations in ISCFF to be analysed
  integer          :: ngconf         !< number of configurations that would be generated 
  integer          :: ncell          !< number of primitiv cell in each direction
  integer          :: imod           !< id calc.eq.'vib+gmod' 
  integer          :: PANdos         !< ( internal ) number of bins
  integer          :: nkphon         !< number of kpoint between ks and kf used when calc = 'vib+band'
                                     !< or (nkphon+1) * (nkphon+1) * (nkphon+1) generated in IBZKPTFF when calc='vib+dos'
  character(len=3) :: path           !< path name
  real(kind=dp)    :: ks(3)          !< start kpoint
  real(kind=dp)    :: kf(3)          !< final kpoint    
  real(kind=dp)    :: resdos         !< resolution in density of states     
  real(kind=dp)    :: omegamax       !< maximum value in dos
  logical          :: lwrite_vectff  !< write the vector field
  real(kind=dp)    :: tempmod        !< temperature at which we populate the mode imod


  character(len=60) :: vib_allowed(5) ! vib allowed flags                     
  ! allowed vib = 'vib' , 'vib+fvib' , 'vib+gmod' , 'vib+band' , 'vib+dos'
  data vib_allowed / 'vib' , 'vib+fvib' , 'vib+gmod' , 'vib+band' , 'vib+dos' /



CONTAINS

! *********************** SUBROUTINE vib_init **********************************
!
!> \brief
!! initialize vib calculation
!
! ******************************************************************************
SUBROUTINE vib_init

  USE io,                  ONLY :  stdin , stdout , stderr , ionode

  implicit none

  integer               :: ioerr
  character(len=132)    :: filename

  namelist /vibtag/ lwrite_vectff , &
                    nconf         , & 
                    ngconf        , & 
                    ncell         , & 
                    imod          , & 
                    tempmod       , & 
                    nkphon        , &
                    resdos        , &
                    omegamax      , &
                    path          , &
                    ks            , &
                    kf   
            
  ! ======================
  !  vibtag default values
  ! ======================
  CALL vib_default_tag
  ! ======================
  !  read vibtag namelist
  ! ======================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , vibtag,iostat=ioerr)
  if ( ioerr .lt. 0 )  then
    io_node WRITE ( stderr , '(a)') 'ERROR reading input_file : vibtag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node WRITE ( stderr , '(a,i8)') 'ERROR reading input_file : vibtag wrong tag',ioerr
    STOP
  endif
  CLOSE  ( stdin )

  ! ======================
  !  check vib tags
  ! ======================
  CALL vib_check_tag
  
  PANdos = INT ( omegamax / resdos , kind=dp ) + 1

  ! ======================
  !  vibtag print info
  ! ======================
  CALL vib_print_info(stdout)


  return
  
END SUBROUTINE vib_init

! *********************** SUBROUTINE vib_default_tag ***************************
!
!> \brief
!! set default values to vib tag
!
! ******************************************************************************
SUBROUTINE vib_default_tag

  implicit none

  ! =================
  !  default values
  ! =================
  lwrite_vectff = .false.
  nconf        = 0
  ngconf       = 0
  ncell        = 0
  imod         = 4
  resdos       = 1.0_dp             
  omegamax     = 100.0_dp

  return

END SUBROUTINE vib_default_tag


! *********************** SUBROUTINE vib_check_tag *****************************
!
!> \brief
!! check vib tag values
!
! ******************************************************************************
SUBROUTINE vib_check_tag

  USE io,                  ONLY :  stderr , ionode

  implicit none

  if ( ncell .eq. 0 ) then
    io_node WRITE ( stderr , * ) 'ERROR in vib module, ncell was not set'
    STOP
  endif

  return

END SUBROUTINE vib_check_tag

! *********************** SUBROUTINE vib_print_info ****************************
!
!> \brief
!! print information about vib calculation
!
! ******************************************************************************
SUBROUTINE vib_print_info(kunit)

  USE control,                  ONLY :  calc 
  USE config,                   ONLY :  natm , ntype, rho
  USE io,                       ONLY :  ionode

  implicit none

  ! local
  integer :: kunit

  if ( ionode ) then
                      blankline(kunit)
                      WRITE ( kunit ,'(a)')       'NORMAL MODE ANALYSIS             '
                      WRITE ( kunit ,'(a)')       'Hessian calc. and generation of config.' 
                      WRITE ( kunit ,'(a)')       'for a given mode'
                      blankline(kunit)
     if (calc.eq.'vib + fvib')  then  
                      WRITE ( kunit ,'(a)')       'Fvibcalc:'
                      WRITE ( kunit ,'(a)')       'This program reads in the vibrational frequencies '
                      WRITE ( kunit ,'(a)')       'and calculates the frequency dependent part of    ' 
                      WRITE ( kunit ,'(a)')       'the vibrational free energy fvib.                 ' 
                      WRITE ( kunit ,'(a)')       'It also histograms the fvib according to the energy'
                      WRITE ( kunit ,'(a)')       'of the minima. ener are energies that are READ  in,'
                      WRITE ( kunit ,'(a)')       'eigval are eigen values, fvib is calculated in the' 
                      WRITE ( kunit ,'(a)')       'program, and so is fvibhist.'
                      WRITE ( kunit ,'(a)')       'Original author S. Sastry JNCASR and probably F. Affouard'
     endif
                      WRITE ( kunit ,'(a)')       'periodic boundary conditions '
                      blankline(kunit)
                      WRITE ( kunit ,'(a)')       'configuration (at equilibrium) file : ISCFF'
                      WRITE ( kunit ,'(a,i5)')    'number of configurations in file    = ',nconf
                      blankline(kunit)
                      WRITE ( kunit ,'(a)')       'START HESSIAN CALCULATION            '
                      WRITE ( kunit ,'(a)')       'save eingenvalues                     : EIGFF'
    if ( lwrite_vectff )  & 
                      WRITE ( kunit ,'(a)')       'save eingenvectors                    : VECTFF'
                      blankline(kunit)
  endif



  return

END SUBROUTINE vib_print_info

! *********************** SUBROUTINE vib_main **********************************
!
!> \brief
!! main program of the vib calculation.
!! this subroutine reads configuration from ISCFF (usually optimized structures)
!
! ******************************************************************************
SUBROUTINE vib_main 

  USE config,                   ONLY :  system , natm , natmi , rx , ry , rz ,  & 
                                        atype , atypei , itype , simu_cell , & 
                                        rho , ntype, config_alloc , coord_format_allowed , atom_dec , read_traj_header , read_traj
  USE control,                  ONLY :  calc , iscff_format , iscff_data
  USE io,                       ONLY :  ionode , stdout , stderr , kunit_ISCFF , kunit_EIGFF , kunit_VECTFF , & 
                                        kunit_DOSFF , kunit_MODFF, kunit_DOSKFF , kunit_IBZKPTFF
  USE thermodynamic,            ONLY :  u_tot , pressure_tot , calc_thermo
  USE cell,                     ONLY :  lattice , dirkar
  USE field,                    ONLY :  field_init , engforce_driver
  USE time,                     ONLY :  vibtimetot , diaghessiantimetot
  USE kspace,                   ONLY :  kmesh , kpoint_sum_init_BZ

  implicit none

  ! local
  integer                                     :: i , j , ik , ierr , ikd , ja , im , ic , ka , nk 
  integer, dimension(:), allocatable          :: dostab
  integer, dimension(:), allocatable          :: dostabtot
  integer, dimension(:,:), allocatable        :: dostabktot
  real(kind=dp), dimension(:,:), allocatable  :: hess
  real(kind=dp), dimension(:), allocatable    :: work, deig, ipiv !,deig_reorder(3 * n)
  real(kind=dp),dimension (:,:), allocatable  :: eigenk
  real(kind=dp)                               :: ak
  TYPE(kmesh)                                 :: km 
!  real(kind=dp)                               :: pressure0, pot0 removed 28/05/13
  character(len=1)                            :: jobz, uplo
  character(len=20)                           :: FMT
  integer*4                                   :: info, lwork
  real(kind=dp)                               :: ttt1 , ttt2 
  real(kind=dp)                               :: ttt1p , ttt2p 
  real(kind=dp)                               :: ttt1d , ttt2d

#ifdef MPI
  ttt1 = MPI_WTIME( ierr )
#endif
 
  OPEN ( UNIT = kunit_EIGFF , FILE = 'EIGFF'  )
  OPEN ( UNIT = kunit_VECTFF, FILE = 'VECTFF' )
  OPEN ( UNIT = kunit_DOSFF , FILE = 'DOSFF'  )
  OPEN ( UNIT = kunit_DOSKFF, FILE = 'DOSKFF'  )

  if (calc.ne.'vib+band'.and.calc.ne.'vib+dos') then 
  io_node WRITE ( kunit_VECTFF ,'(a)') & 
  '       #rx               ry             dx              dy              rx              rz              dx             dz              ry              rz              dy             dz           eigenvalue'   
  endif

  ! ====================================================
  !  read main config parameters from ISCFF and reopen
  ! ====================================================
  if ( iscff_format .ne. 0 ) OPEN ( UNIT = kunit_ISCFF , FILE = 'ISCFF'  )
  if ( iscff_format .eq. 0 ) OPEN ( UNIT = kunit_ISCFF , FILE = 'ISCFF',form='unformatted'  )
  CALL read_traj_header( kunit_ISCFF , iscff_format ) 
  if ( iscff_format .ne. 0 ) OPEN ( UNIT = kunit_ISCFF , FILE = 'ISCFF'  )
  if ( iscff_format .eq. 0 ) OPEN ( UNIT = kunit_ISCFF , FILE = 'ISCFF',form='unformatted'  )
  

  CALL lattice ( simu_cell )
  rho = REAL ( natm , kind=dp ) / simu_cell%omega 
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms')
  CALL field_init
  CALL typeinfo_init

  ! ======================================================
  ! allocate main quantities hessian, eigenvalues arrays
  ! ======================================================
  allocate( hess ( 3 * natm , 3 * natm )          )
  allocate( work ( 9 * natm )                     )
  allocate( deig ( 3 * natm ) , ipiv ( 3 * natm ) )
  allocate( dostabtot ( 0 : PANdos + 1 )          )
  allocate( dostab    ( 0 : PANdos + 1 )          )
  allocate( dostabktot( 0 : PANdos + 1,0:3*ntype )      )
  dostabktot = 0
  dostabtot = 0

  CALL print_general_info ( stdout )

  io_node blankline(stdout)
  io_node WRITE ( stdout ,'(a,i9,a,i9)')  'hessian dimension ', 3 * natm , ' x ' , 3 * natm

  ! ===========================================
  !  LOOP OVER CONFIGURATIONS (AT EQUILIBRIUM)
  ! ===========================================

  conf : do ic = 1, nconf
#ifdef MPI
    ttt1p = MPI_WTIME( ierr )
#endif

    if ( ionode ) WRITE ( stdout ,'(a,i5)') 'read config',ic

    CALL read_traj ( kunit_ISCFF , iscff_format , iscff_data )

    CALL lattice ( simu_cell )
    rho = natm / simu_cell%omega

    CALL typeinfo_init

    ! ===============================================
    ! do we need forces or energy at this point ??? 
    ! ===============================================
    !CALL engforce_driver 

    !CALL calc_thermo
    !pot0      = u_tot
    !pressure0 = pressure_tot 

    ! ========================
    !  get the hessian matrix
    ! ========================
    CALL hessian(hess)

    ! ===========================
    !  real space dos: 3N modes
    ! ===========================
    if (calc.ne.'vib+band'.and.calc.ne.'vib+dos') then

      ! ===========================
      !  diagonalisation of hess
      ! ===========================
      io_node WRITE( stdout ,'(a)') 'start diagonalization'
      lwork = 9 * natm
      jobz = 'V'
      uplo = 'U'
#ifdef MPI
      ttt1d = MPI_WTIME( ierr )  
#endif
      CALL DSYEV( jobz , uplo , 3 * natm , hess ,3 * natm , deig , work , lwork , info )
      if ( info .ne. 0 ) then
        io_node WRITE ( stderr ,'(a,i5)') 'ERROR in vib_main : improper termination. of DSYEV info = ',info
        STOP 
      else
#ifdef MPI
        ttt2d = MPI_WTIME( ierr )  
        diaghessiantimetot = diaghessiantimetot + ( ttt2d - ttt1d )
#endif
        io_node WRITE ( stderr ,'(a)')    'hessian diagonalisation ok'
      endif

      ! ===========================================
      !  WRITE frequencies^2 (hessian eigenvalues) 
      ! ===========================================
      do im = 1,3 * natm
        io_node WRITE ( kunit_EIGFF ,'(f24.10)') deig(im)
      enddo
      ! =========================================
      !  DOS: density of states  of the 3N modes
      ! =========================================
      dostab = 0
      do i = 4,3 * natm
        ak = (SQRT (deig(i)))/resdos
        ka = INT (ak) + 1
        if (ka.gt.PANdos + 1.or.ka.lt.0) then
           WRITE ( stderr , * ) 'ERROR out of bound in dostab'
           WRITE ( stderr , * ) i,ka,PANdos + 1,ak,deig(i)
           STOP
        endif
        dostab(ka) = dostab(ka) + 1
        dostabtot(ka) = dostabtot(ka) + 1
      enddo
      dostab(0) = 3
      dostabtot(0) = dostabtot(0) + 3

      do i = 0 , PANdos + 1
        io_node WRITE ( kunit_DOSFF ,'(3f16.8)') &
        REAL( i , kind=dp ) * resdos, REAL(dostab(i),kind=dp)    / ( 3.0_dp * REAL(natm,kind=dp) * resdos ) , &
                                      REAL(dostabtot(i),kind=dp) / ( 3.0_dp * REAL(natm,kind=dp) * resdos * REAL( ic ,kind = dp) )
      enddo
      io_node blankline(kunit_DOSFF )
      io_node blankline(kunit_DOSFF )
      io_node WRITE ( stdout ,'(a)') 'dos ok' 


      ! ====================
      !  WRITE vector field
      ! ====================
      if ( lwrite_vectff ) then
        do im = 1,3 * natm
          do ja = 1, natm
            if ( ionode ) &
            WRITE  ( kunit_VECTFF ,'(13f16.8)') rx(ja),ry(ja),hess(ja,im),hess(ja + natm,im), &
                                                ry(ja),rz(ja),hess(ja + natm,im),hess(ja + 2 * natm,im), &
                                                rx(ja),rz(ja),hess(ja,im),hess(ja + 2 * natm,im), &
                                                deig(im)
          enddo
          io_node blankline(kunit_VECTFF)
          io_node blankline(kunit_VECTFF)
        enddo
      endif

    endif 
    ! ==============================================
    !  beware hess is changing when diagonalized!!! 
    ! ==============================================
    !  using reciprocal space we get the complete : 
    ! ==============================================
    !  band structure    
    ! ==============================================
    if ( calc .eq.'vib+band' ) CALL band(hess)
    ! ==============================================
    ! total dos for more than 3N k - vector
    ! ==============================================
    if ( calc .eq.'vib+dos'  ) then
   
      ! ============================================ 
      !  generate k-point mesh in reciprocal lattice 
      !  write kpoints on IBZKPTFF file 
      ! ============================================ 
!      CALL gene_IBZKPTFF 
      nk = ( ( nkphon + 1 ) *( nkphon + 1 ) * ( nkphon + 1 ) )
      ! kmax could be read in control.F
      km%kmax(1) = nkphon
      km%kmax(2) = nkphon
      km%kmax(3) = nkphon
      km%nk = nk
      km%meshlabel='km-dos'
      allocate ( km%kptx ( nk ) ,  km%kpty ( nk ),  km%kptz ( nk ) , km%kptk(nk) )
      CALL do_split ( km%nk , myrank , numprocs , km%kpt_dec , 'kpts ')
      CALL kpoint_sum_init_BZ ( km , 1.0_dp )

      ! why the hell did I keet it 
      ! ============================================
      !  read number of kpoints in IBZKPTFF
      ! ============================================
!      OPEN ( UNIT = kunit_IBZKPTFF ,FILE = 'IBZKPTFF')
!      READ ( kunit_IBZKPTFF , * ) cccc
!      READ ( kunit_IBZKPTFF , * ) nk
!      READ ( kunit_IBZKPTFF , * ) cccc
!      CLOSE( kunit_IBZKPTFF     )

      ! ============================================
      !  allocate eigenvalues 3 modes per kpoints
      ! ============================================
      allocate( eigenk ( nk , 3 * ntype  ) )
      eigenk = 0.0_dp
       
      ! ============================================
      !   diagonalisation of the 3*nk states
      ! ============================================
      CALL doskpt( hess , eigenk , km )
      ! ============================================
      !  DOS: density of states of the 3*nk modes
      ! ============================================
      if ( ionode ) write(  stdout , '(a)' ) 'dos of the 3*nk modes'
      k_loop : do ik = 1, nk
        ! ==========================================
        !    loop on modes
        ! ==========================================
        do ikd = 1 , 3 * ntype 
          if (eigenk(ik,ikd).lt.0.0_dp ) then
            cycle 
          endif
          ak = (SQRT(eigenk(ik,ikd)))/resdos
          ka = INT (ak) + 1
          !write (stdout, '(i9,2f12.6,2i9)') ,ik,ak,eigenk(ik,ikd),ka,PANdos
          if (ka.gt.PANdos + 1.or.ka.lt.0) then
             WRITE ( stderr , * ) 'ERROR out of bound in eigenk'
             WRITE ( stderr , '(3i12,2f12.5)' ) ik,ka,PANdos + 1,ak,eigenk(ik,ikd)
             STOP
          endif
          ! ========================================
          !   bin the dos for
          !   index 0 : all modes
          !   index ikd : separately
          ! ========================================
          dostabktot(ka,0)   = dostabktot(ka,0) + 1
          dostabktot(ka,ikd) = dostabktot(ka,ikd) + 1
        enddo

      enddo k_loop 

      deallocate(eigenk) 
      deallocate(km%kptx,km%kpty,km%kptz,km%kptk) 
        
    endif ! vib+dos

#ifdef MPI
    ttt2p = MPI_WTIME( ierr )
    io_node WRITE ( stdout , 110 ) 'config : ',ic,' VIB  ', ttt2p - ttt1p  
#endif

  enddo conf 

  ! ===========================================
  !    write density of states to DOSKFF file
  ! ===========================================
  do i = 0,PANdos + 1
#ifdef GFORTRAN
    WRITE(FMT,*) 6*ntype+3
    io_node WRITE ( kunit_DOSKFF ,'('// ADJUSTL(FMT) //'f16.8)') &
#else
    io_node WRITE ( kunit_DOSKFF ,'(<6*ntype+3>f16.8)') &
#endif
    REAL ( i, kind = dp ) * resdos , ( REAL ( dostabktot(i,j) , kind = dp ) / ( 3.0 * nk * resdos * nconf ), j = 0 , 3*ntype )
  enddo
  io_node blankline(kunit_DOSKFF )
  io_node blankline(kunit_DOSKFF )
  io_node WRITE ( stdout ,'(a)')       'dos (k) ok'

  CLOSE ( kunit_DOSKFF ) 
  CLOSE ( kunit_DOSFF ) 
  CLOSE ( kunit_VECTFF )
  CLOSE ( kunit_EIGFF )
  CLOSE ( kunit_ISCFF )

  ! ===================================
  !  generate config from a given mode
  ! ===================================
  if (calc.eq.'vib+gmod') then
    if ( ionode ) then
      blankline(stdout)
      blankline(stdout)
      WRITE ( stdout ,'(a,i4,a,f8.3)')  'generate ', ngconf ,' configurations at temperature  = ', tempmod
      WRITE ( stdout ,'(a)')            'save configurations in file          : MODFF'
      blankline(stdout)
      WRITE ( stdout ,'(a)')            'Method:'
      WRITE ( stdout ,'(a)')            'gaussian distributions of normal coordinates'
      WRITE ( stdout ,'(a)')            'equation 11 of J.Phys.Chem.B 2005, 109, 7245'
      blankline(stdout)
    endif

    !generation of modes imod
    OPEN (unit = kunit_MODFF ,FILE = 'MODFF')
    CALL generate_modes(deig,hess,kunit_MODFF)
    CLOSE ( kunit_MODFF )
  endif 
  ! ====================================
  !  start vibrational free energy calc
  ! ====================================
  if (calc.eq.'vib+fvib') then
    if ( ionode ) then
      blankline(stdout)
      WRITE ( stdout ,'(a)')       'read thermo (Fvibcalc)              : ISTHFF'
      WRITE ( stdout ,'(a)')       'WARNING !!!! DEV. STATE !!!! WARNING !!! NEED TO BE TESTED!!!'
    endif
    CALL fvibcalc
    STOP
  endif

  deallocate(dostabktot)
  deallocate(dostab)
  deallocate(dostabtot)
  deallocate(deig)
  deallocate(work)
  deallocate(hess) 

#ifdef MPI
  ttt2 = MPI_WTIME( ierr )
  vibtimetot = vibtimetot + ( ttt2 - ttt1 ) 
#endif

  return

110   FORMAT(2X,A8,I8,A20,' :  cpu time',F9.2)

END SUBROUTINE vib_main

! *********************** SUBROUTINE hessian ***********************************
!
!
! ******************************************************************************
SUBROUTINE hessian ( hess )

  USE config,           ONLY :  natm , rx , ry , rz , itype , ntype , simu_cell 
  USE io,               ONLY :  ionode , stdout , stderr
  USE field,            ONLY :  rcutsq , sigsq , epsp , fc , uc , plj , qlj
  USE cell,             ONLY :  kardir , dirkar
  USE time,             ONLY :  hessiantimetot

  implicit none

  ! global
  real(kind=dp) :: hess(3 * natm,3 * natm)

  ! local
  integer       :: it , jt , p1, p2, ia, ja , ierr
  real(kind=dp) :: rxi , ryi, rzi 
  real(kind=dp) :: rxij , ryij , rzij , rijsq
  real(kind=dp) :: sxij, syij, szij
  real(kind=dp) :: sr2 , sr , srp2 , srp4 , srq2 , srq4 , c1 , c2        
  real(kind=dp) :: txx , txy , tyy , tyz , tzz , txz
  real(kind=dp) :: np2 (ntype , ntype ) 
  real(kind=dp) :: nq2 (ntype , ntype ) 
  real(kind=dp) :: np4 (ntype , ntype ) 
  real(kind=dp) :: nq4 (ntype , ntype ) 
  real(kind=dp) :: ttt1 , ttt2 
#ifdef debug_symmetry
  real(kind=dp) , allocatable  :: tmpdebug( : , : )       
  real(kind=dp) :: fnnb , snnb , tnnb
  integer       :: si , sj , ita , jta 
#endif

#ifdef debug_symmetry
  allocate ( tmpdebug ( 3 * ntype , 3 * ntype ) ) 
#endif
#ifdef MPI
  ttt1 = MPI_WTIME(ierr) ! timing info
#endif

  do it = 1 , ntype
    do jt = 1 , ntype
      np2 ( it , jt ) = plj ( it , jt ) 
      nq2 ( it , jt ) = qlj ( it , jt ) 
      np4 ( it , jt ) = plj ( it , jt ) 
      nq4 ( it , jt ) = qlj ( it , jt ) 
      np2 ( it , jt ) = np2 ( it , jt ) + 2
      nq2 ( it , jt ) = nq2 ( it , jt ) + 2
      np4 ( it , jt ) = np4 ( it , jt ) + 4
      nq4 ( it , jt ) = nq4 ( it , jt ) + 4
    enddo
  enddo

  hess = 0.0_dp

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = 1 , natm 

    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )

    do ja = ia + 1 , natm

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
        sr2 = sigsq(p1,p2)/rijsq
        sr = SQRT (sr2)
        srp2 = sr ** np2(p1,p2)
        srp4 = sr ** np4(p1,p2)
        srq2 = sr ** nq2(p1,p2)
        srq4 = sr ** nq4(p1,p2)

        ! ===========
        !  Hessian 
        ! ===========
        c1 = ( fc ( p1 , p2 ) / sigsq ( p1 , p2 ) ) * &
        ( ( plj ( p1 , p2 ) + 2.0_dp ) * srp4 - ( qlj ( p1 , p2 ) + 2.0_dp ) * srq4 )
        c2 =  fc(p1,p2) * (srq2 - srp2)

        hess( ia , ja )                       = rxij * rxij * c1 + c2  ! x,x
        hess( ia , natm + ja )                = rxij * ryij * c1       ! x,y
        hess( ia , 2 * natm + ja )            = rxij * rzij * c1       ! x,z  

        hess( natm + ia , ja )                = ryij * rxij * c1       ! x,y
        hess( natm + ia , natm + ja )         = ryij * ryij * c1 + c2  ! y,y  
        hess( natm + ia , 2 * natm + ja )     = ryij * rzij * c1       ! y,z

        hess( 2 * natm + ia , ja )            = rzij * rxij * c1       ! x,z
        hess( 2 * natm + ia , natm + ja )     = rzij * ryij * c1       ! y,z
        hess( 2 * natm + ia , 2 * natm + ja ) = rzij * rzij * c1 + c2  ! z,z 

        ! ====================== 
        !  hessian is symmetric 
        ! ====================== 

        hess( ja , ia )                       = hess( ia , ja ) 
        hess( natm + ja , ia )                = hess( ia , natm + ja ) 
        hess( 2 * natm + ja , ia )            = hess( ia , 2 * natm + ja ) 

        hess( ja , natm + ia )                = hess( natm + ia , ja ) 
        hess( natm + ja , natm + ia )         = hess( natm + ia , natm + ja ) 
        hess( 2 * natm + ja , natm + ia )     = hess( natm + ia , 2 * natm + ja ) 

        hess( ja , 2 * natm + ia )            = hess( 2 * natm + ia , ja ) 
        hess( natm + ja , 2 * natm + ia )     = hess( 2 * natm + ia , natm + ja) 
        hess( 2 * natm + ja , 2 * natm + ia ) = hess( 2 * natm + ia , 2 * natm + ja ) 

      endif
    enddo
  enddo

!=================================================================================================================
! elements diagonaux:
! "appelle abusivement self - force. il traduit simplement le fait que
! lorsqu'un atome est deplace, il ressent une force exercee par l'ensemble des autres atom du cristal.   
! La forme de ce terme se clarifie quand on s'aperercoit que, deplacer un atome dans une direction
! donnee, est equivalent a deplacer l'ensemble du cristal a l'exception de cet atome dans la direction opposee
!
! Notes de Cours de Physique des materiaux
! P. GHOSEZ, J - Y. RATY
! Universite de Liege
! http://www.phythema.ulg.ac.be/Teaching/Cours/Cours - Phys_Mat/Notes - Phys_Mat.pdf
!=================================================================================================================

  do ia = 1, natm
    txx = 0.0_dp
    txy = 0.0_dp
    txz = 0.0_dp
    tyy = 0.0_dp
    tyz = 0.0_dp
    tzz = 0.0_dp
    do ja = 1, natm
      if ( ( ia .eq. ja ) .and. ( hess ( ia , ja ) .ne. 0.0_dp ) ) then
        io_node WRITE ( stderr , * )' ERROR : in hessian subroutine force constant is null', ia , ja , hess(ia,ja)
      endif
      if ( ia .eq. ja ) cycle ! zero anyway 
      txx = txx + hess ( ia , ja )
      txy = txy + hess ( ia , natm + ja )
      txz = txz + hess ( ia , 2 * natm + ja )
      tyy = tyy + hess ( natm + ia , natm + ja )
      tyz = tyz + hess ( natm + ia , 2 * natm + ja )
      tzz = tzz + hess ( 2 * natm + ia , 2 * natm + ja )
       
    enddo
    ! diagonal = - sum i,j
    hess ( ia , ia )                       = - txx
    hess ( ia , natm + ia)                 = - txy
    hess ( ia , 2 * natm + ia )            = - txz
    hess ( natm + ia , natm + ia )         = - tyy
    hess ( natm + ia , 2 * natm + ia )     = - tyz
    hess ( 2 * natm + ia , 2 * natm + ia ) = - tzz
   
    !symmetry
    hess ( natm + ia , ia ) = hess ( ia , natm + ia )
    hess ( 2 * natm + ia , ia ) = hess ( ia , 2 * natm + ia )
    hess ( 2 * natm + ia , natm + ia) = hess ( natm + ia , 2 * natm + ia )
        

    if ( hess ( ia , ia ) .eq. 0.0_dp .and. ionode )                       &
    WRITE ( stdout , * ) 'in forcehes zero at ', ia
    if ( hess ( ia + natm , ia + natm ) .eq. 0.0_dp .and. ionode )         &
    WRITE ( stdout , * ) 'in forcehes zero at ', natm + ia
    if ( hess ( ia + 2 * natm , ia + 2 * natm ) .eq. 0.0_dp .and. ionode ) & 
    WRITE ( stdout , * ) 'in forcehes zero at ', 2 * natm + ia
  enddo


#ifdef debug_symmetry
 ! =========================================================================== 
 !  this debug unit, check the symmetrical form of the force-constant matrix 
 !  in case of a fcc structrure c.f Table 3.1 of "Thermal vibration in
 !  Crystallography" by B.T.M. Willis and A.W. Pryor Cambridge Univ. Press 1975
 ! =========================================================================== 

  if ( ntype .eq. 1 ) then
    lseparator(stdout) 
    write ( stdout , '(a)' ) ' monoatomic fcc'
    ! the 3 first squared fcc distances
    ! first nearest neighbour fcc
    fnnb = simu_cell%ANORM(1)*SQRT(2.0_dp)/2.0_dp/ncell
    fnnb = fnnb*fnnb
    ! second nearest neighbour fcc
    snnb = simu_cell%ANORM(1)/ncell
    snnb = snnb*snnb
    ! third nearest neighbour fcc
    tnnb = SQRT(6.0_dp)*simu_cell%ANORM(1)/ncell/2.0_dp
    tnnb = tnnb*tnnb
  endif
  
  if ( ntype .eq. 2 ) then
    lseparator(stdout) 
    write ( stdout , '(a)' ) 'NaCl fcc'
    ! the 3 first squared NaCl distances
    ! first nearest neighbour NaCl 
    fnnb = simu_cell%ANORM(1)/2.0_dp/ncell
    fnnb = fnnb*fnnb
    ! second nearest neighbour NaCl
    snnb = SQRT(2.0_dp)*simu_cell%ANORM(1)/2.0_dp/ncell
    snnb = snnb*snnb
    ! third nearest neighbour NaCl
    tnnb = SQRT(3.0_dp)*simu_cell%ANORM(1)/ncell/2.0_dp
    tnnb = tnnb*tnnb
  endif
  

  do ia = 1 , natm

    ita = itype ( ia ) 
    si = 3 * ita - 3
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )

    do ja = 1 , natm

      jta  = itype ( ja ) 
      sj   = 3 * jta - 3
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

      if ( rijsq .lt. rcutsq(ita,jta) ) then

        tmpdebug = 0.0_dp
        tmpdebug ( 1 + si , 1 + sj ) = hess ( ja            , ia            )
        tmpdebug ( 2 + si , 2 + sj ) = hess ( ja + natm     , ia + natm     )
        tmpdebug ( 3 + si , 3 + sj ) = hess ( ja + 2 * natm , ia + 2 * natm )
        tmpdebug ( 1 + si , 2 + sj ) = hess ( ja            , ia + natm     )
        tmpdebug ( 1 + si , 3 + sj ) = hess ( ja            , ia + 2 * natm )
        tmpdebug ( 2 + si , 3 + sj ) = hess ( ja + natm     , ia + 2 * natm )
        tmpdebug ( 2 + si , 1 + sj ) = hess ( ia + natm     , ja            )
        tmpdebug ( 3 + si , 1 + sj ) = hess ( ia + 2 * natm , ja            )
        tmpdebug ( 3 + si , 2 + sj ) = hess ( ia + 2 * natm , ja + natm     )

        ! =============================================================================================================
        !                                           check the monoatomic fcc
        ! =============================================================================================================
        if ( ntype .eq. 1 ) then
          ! ============================
          ! first nearest neighbour fcc
          ! ============================
          if ( abs ( rijsq-fnnb ) .lt. 0.1d0 ) then 
            lseparator(stdout) 
            write ( stdout , '(a,i6,a,i6)' ) 'ia  = ', ia ,' ja  = ',ja
            write ( stdout , '(a,i6,a,i6)' ) 'ita = ', ita,' jta = ',jta
            write ( stdout , '(a,2f12.5)' ) 'ia , ja are first neighbour',SQRT( rijsq ),SQRT(fnnb )
            write ( stdout , '(a)' ) 'by symmetry the force constant matrix is :'
            if ( rxij .eq. 0.0 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| b   0   0  |' 
              write ( stdout , '(a)' ) '| 0   a   g  |' 
              write ( stdout , '(a)' ) '| 0   g   a  |' 
              io_node blankline(stdout)
            endif
            if ( ryij .eq. 0.0 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   0   g  |' 
              write ( stdout , '(a)' ) '| 0   b   0  |' 
              write ( stdout , '(a)' ) '| g   0   a  |' 
              io_node blankline(stdout)
            endif
            if ( rzij .eq. 0.0 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   g   0  |' 
              write ( stdout , '(a)' ) '| g   a   0  |' 
              write ( stdout , '(a)' ) '| 0   0   b  |' 
              io_node blankline(stdout)
            endif
            CALL print_tensor ( tmpdebug , 'HESSFCC ' )
          ! ============================
          ! second nearest neighbour fcc
          ! ============================
          else if ( abs ( rijsq-snnb ) .lt. 0.1d0 ) then 
            lseparator(stdout) 
            write ( stdout , '(a,i6,a,i6)' ) 'ia  = ', ia ,' ja  = ',ja
            write ( stdout , '(a,i6,a,i6)' ) 'ita = ', ita,' jta = ',jta
            write ( stdout , '(a,2f12.5)'  ) 'ia , ja are second neighbour',SQRT( rijsq ),SQRT(snnb)
            write ( stdout , '(a)' ) 'by symmetry the force constant matrix is :'
            if ( rxij .ne. 0.0 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   0   0  |' 
              write ( stdout , '(a)' ) '| 0   b   0  |' 
              write ( stdout , '(a)' ) '| 0   0   b  |' 
              io_node blankline(stdout)
            endif
            if ( ryij .ne. 0.0 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| b   0   0  |' 
              write ( stdout , '(a)' ) '| 0   a   0  |' 
              write ( stdout , '(a)' ) '| 0   0   b  |' 
              io_node blankline(stdout)
            endif
            if ( rzij .ne. 0.0 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| b   0   0  |' 
              write ( stdout , '(a)' ) '| 0   b   0  |' 
              write ( stdout , '(a)' ) '| 0   0   a  |' 
              io_node blankline(stdout)
            endif
            CALL print_tensor ( tmpdebug , 'HESSFCC ' )
          ! ============================
          ! third nearest neighbour fcc
          ! ============================
          else if ( abs ( rijsq-tnnb ) .lt. 0.1d0 ) then 
            lseparator(stdout) 
            write ( stdout , '(a,i6,a,i6)' ) 'ia  = ', ia ,' ja  = ',ja
            write ( stdout , '(a,i6,a,i6)' ) 'ita = ', ita,' jta = ',jta
            write ( stdout , '(a,2f12.5)' ) 'ia , ja are third neighbour',SQRT( rijsq ),SQRT(tnnb)
            write ( stdout , '(a)' ) 'by symmetry the force constant matrix is :'
            if ( abs ( ryij )  .eq. abs ( rzij ) ) then
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   b   b  |' 
              write ( stdout , '(a)' ) '| b   g   d  |' 
              write ( stdout , '(a)' ) '| b   d   g  |' 
              io_node blankline(stdout)
            endif
            if ( abs ( rxij ) .eq. abs ( rzij ) ) then
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| g   b   d  |' 
              write ( stdout , '(a)' ) '| b   a   b  |' 
              write ( stdout , '(a)' ) '| d   b   g  |' 
              io_node blankline(stdout)
            endif
            if ( abs ( rxij ) .eq. abs ( ryij ) ) then
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| g   d   b  |' 
              write ( stdout , '(a)' ) '| d   g   b  |' 
              write ( stdout , '(a)' ) '| b   b   a  |' 
              io_node blankline(stdout)
            endif
            CALL print_tensor ( tmpdebug , 'HESSFCC ' )
          endif
        endif

        ! =============================================================================================================
        !                    NaCl rock-salt structure ( i.e 2 atoms in the primitiv cell )
        ! =============================================================================================================
        if ( ntype .eq. 2 ) then
          ! ============================
          ! first nearest neighbour fcc
          ! ============================
          if ( abs ( rijsq-fnnb ) .lt. 0.1d0 ) then 
            lseparator(stdout) 
            write ( stdout , '(a,i6,a,i6)' ) 'ia  = ', ia ,' ja  = ',ja
            write ( stdout , '(a,i6,a,i6)' ) 'ita = ', ita,' jta = ',jta
            write ( stdout , '(a,i6,a,i6)' ) 'si  = ', si ,' sj  = ',sj
            write ( stdout , '(a,2f12.5)' ) 'ia , ja are first neighbour',SQRT( rijsq ),SQRT(fnnb )
            write ( stdout , '(a)' ) 'by symmetry the force constant matrix is :'
            if ( abs( rxij ) .gt. 1e-6 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   0   0  |' 
              write ( stdout , '(a)' ) '| 0   b   0  |' 
              write ( stdout , '(a)' ) '| 0   0   b  |' 
              io_node blankline(stdout)
            endif
            if ( dabs( ryij ) .gt. 1e-6 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| b   0   0  |' 
              write ( stdout , '(a)' ) '| 0   a   0  |' 
              write ( stdout , '(a)' ) '| 0   0   b  |' 
              io_node blankline(stdout)
            endif
            if ( dabs( rzij ) .gt. 1e-6 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| b   0   0  |' 
              write ( stdout , '(a)' ) '| 0   b   0  |' 
              write ( stdout , '(a)' ) '| 0   0   a  |' 
              io_node blankline(stdout)
            endif
            CALL print_tensor_nxn ( tmpdebug , 'HESSFCC ' , 3*ntype)
          ! ============================
          ! second nearest neighbour fcc
          ! ============================
          else if ( abs ( rijsq-snnb ) .lt. 0.1d0 ) then 
            lseparator(stdout) 
            write ( stdout , '(a,i6,a,i6)' ) 'ia  = ', ia ,' ja  = ',ja
            write ( stdout , '(a,i6,a,i6)' ) 'ita = ', ita,' jta = ',jta
            write ( stdout , '(a,i6,a,i6)' ) 'si  = ', si ,' sj  = ',sj
            write ( stdout , '(a,2f12.5)'  ) 'ia , ja are second neighbour',SQRT( rijsq ),SQRT(snnb)
            write ( stdout , '(a)' ) 'by symmetry the force constant matrix is :'
            if ( abs( rxij ) .lt. 1e-6 ) then 
              blankline(stdout)
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| b   0   0  |' 
              write ( stdout , '(a)' ) '| 0   a   g  |' 
              write ( stdout , '(a)' ) '| 0   g   a  |' 
              io_node blankline(stdout)
            endif
            if ( abs( ryij ) .lt. 1e-6 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   0   g  |' 
              write ( stdout , '(a)' ) '| 0   b   0  |' 
              write ( stdout , '(a)' ) '| g   0   a  |' 
              io_node blankline(stdout)
            endif
            if ( abs( rzij ) .lt. 1e-6 ) then 
              io_node blankline(stdout)
              write ( stdout , '(a)' ) '| a   g   0  |' 
              write ( stdout , '(a)' ) '| g   a   0  |' 
              write ( stdout , '(a)' ) '| 0   0   b  |' 
              io_node blankline(stdout)
            endif
            CALL print_tensor_nxn ( tmpdebug , 'HESSFCC ' , 3*ntype)
          ! ============================
          ! third nearest neighbour fcc
          ! ============================
          else if ( abs ( rijsq-tnnb ) .lt. 0.1d0 ) then 
            lseparator(stdout) 
            write ( stdout , '(a,i6,a,i6)' ) 'ia  = ', ia ,' ja  = ',ja
            write ( stdout , '(a,i6,a,i6)' ) 'ita = ', ita,' jta = ',jta
            write ( stdout , '(a,i6,a,i6)' ) 'si  = ', si ,' sj  = ',sj
            write ( stdout , '(a,2f12.5)' ) 'ia , ja are third neighbour',SQRT( rijsq ),SQRT(tnnb)
            write ( stdout , '(a)' ) 'by symmetry the force constant matrix is :'
            io_node blankline(stdout)
            write ( stdout , '(a)' ) '| a   b   b  |' 
            write ( stdout , '(a)' ) '| b   a   b  |' 
            write ( stdout , '(a)' ) '| b   b   a  |' 
            io_node blankline(stdout)
            CALL print_tensor_nxn ( tmpdebug , 'HESSFCC ' , 3*ntype)
          endif
        endif
      endif
    enddo
  enddo

  deallocate ( tmpdebug )
  stop
#endif

  ! ======================================
  !         direct to kartesian 
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

#ifdef MPI
  ttt2 = MPI_WTIME(ierr) ! timing info
  hessiantimetot = hessiantimetot + ( ttt2 - ttt1 )
#endif

  return

END SUBROUTINE hessian

! *********************** SUBROUTINE fvibcalc **********************************
!
!> \brief
!! This program reads in the vibrational frequencies and calculates 
!! the frequency dependent part of the vibrational free energy fvib. 
!! It also histograms the fvib according to the energy of the minima.
!! ener are energies that are READ  in, eigval are eigen values, fvib 
!! is calculated in the program, and so is fvibhist.
!
! ******************************************************************************
SUBROUTINE fvibcalc

  USE config,   ONLY :  natm
  USE io,  ONLY :  ionode , stdout , kunit_ISTHFF , kunit_EIGFF , kunit_VIBFF

  implicit none

!lGocal
  integer            :: i, j
  integer,parameter  :: nbin = 100, ncmax = 100000
  real(kind=dp)      :: ener(ncmax), eigval(3 * natm), fvib(ncmax), fvibhist(nbin,3)
  character(len=100) ::   header
  integer            :: ncs , nskip ,ifv, ieg, ind
  real(kind=dp)      :: emin, emax
  real(kind=dp)      :: enerind, pres, de, dum, fvibav, pofet,atmp 

  ! trash  
  integer            :: iiii
  real(kind=dp)      :: xxxx 

!# config                 eIS                grad                Pres  Iter  Neng
!u_initial     Press_initial
!       1     - 3.947162344773      0.000000000037      6.708240000000    20    70 - 2.603606441227   13.556241500514
!       2     - 3.947162344773      0.000000000049      6.708240000000    21    74 - 2.328527224939   15.013373656527

  ncs = 1
  emin = 999999.0_dp
  emax = - 999999.0_dp
  OPEN (unit = kunit_ISTHFF ,file = 'ISTHFF',status = 'old')
  READ ( kunit_ISTHFF , * ) header
  READ ( kunit_ISTHFF , * ) header
  DO i = 1,nconf
    READ ( kunit_ISTHFF , * ) iiii, enerind, xxxx, pres, iiii, iiii, xxxx
    if ( enerind .gt. emax) emax = enerind
    if ( enerind .lt. emin) emin = enerind
  ENDDO
  CLOSE ( kunit_ISTHFF )

  de = ABS (emax - emin)/ REAL ( nbin , kind = dp )

!  PRINT  * ,' *  - ISthermo file = '
!  PRINT  * ,' *  - vib      file = ',F_IN2
!  PRINT  * ,' - fvibn     file = ',F_IN3
  !print  * ,'nconf = ',nconf
!  print  * ,'emin  = ',emin
!  print  * ,'emax  = ',emax
!  print  * ,'de    = ',de

  OPEN (unit = kunit_ISTHFF ,file = 'ISTHFF',status = 'old')
  READ ( kunit_ISTHFF , * ) header
  READ ( kunit_ISTHFF , * ) header
  OPEN (unit = kunit_EIGFF ,file = 'EIGFF',status = 'old')
  OPEN (unit = kunit_VIBFF ,file = 'VIBFF',status = 'unknown')

! =========================================
!  specify minimum energy and the bin size 
! =========================================

  do i = 1, nconf
    fvib(i) = 0.0
  enddo

  do i = 1, nbin
    fvibhist(i,1) = emin + (i - .5) * de  
    fvibhist(i,2) = 0.0
    fvibhist(i,3) = 0.0
  enddo

  dum = 10.0_dp
  ifv = 0

  do i = 1, nconf
    READ ( kunit_ISTHFF , * ) iiii, enerind, xxxx, pres, iiii, iiii, xxxx
    if ( mod ( i , ncs ) .eq. 0 )then
      ifv = ifv + 1
      ener(ifv) = enerind 

      do ieg = 1,3 * natm 
      READ ( kunit_EIGFF , * ) eigval(ieg)
      enddo

      nskip = 3
      do j = 4, 3 * natm 
        if ( eigval ( j ) .lt. 0.0_dp ) then
          nskip = nskip + 1 
        else 
          fvib ( ifv ) = fvib ( ifv ) + log ( eigval ( j ) )
        endif
      enddo 

      if ( nskip .gt. 3 .and. ionode ) WRITE ( stdout , * ) ifv , nskip 

      fvib ( ifv ) = fvib ( ifv ) * ( 3.0_dp * natm )/ ( ( 2.0_dp * natm ) * ( 3.0_dp * natm - nskip * 1.0_dp ) )
    endif 
  enddo

  CLOSE ( kunit_EIGFF )
  CLOSE ( kunit_ISTHFF )
  

  do i = 1, nconf
    atmp = (ener(i) - emin)/de
    ind = INT (atmp)  
    ind = ind + 1
    fvibhist(ind,2) = fvibhist(ind,2) + 1.0_dp
    fvibhist(ind,3) = fvibhist(ind,3) + fvib(i)
  enddo 

  do i = 1, nbin 
    if ( fvibhist ( i , 2 ) .ne. 0.0_dp ) then
      io_node WRITE ( kunit_VIBFF ,'(3(2x,e14.6))') &
      fvibhist(i,1),fvibhist(i,3)/fvibhist(i,2),fvibhist(i,2)/(de * nconf)
    endif
  enddo


  fvibav = 0.0_dp
  pofet = 0.0_dp

  do i = 1, nbin 

    if ( fvibhist ( i , 2 ) .ne. 0.0_dp ) then
      fvibav = fvibav + fvibhist ( i , 3 )
      pofet  = pofet  + fvibhist ( i , 2 )
    endif 
  enddo

  io_node WRITE ( stdout , * ) fvibav/pofet, pofet

  return

END SUBROUTINE fvibcalc

! *********************** SUBROUTINE generate_modes ****************************
!
!> \brief
!! generate configuration from a given mode
!
! ******************************************************************************
SUBROUTINE generate_modes ( deig , hess , kunit )

  USE constants,                ONLY :  dzero
  USE config,                   ONLY :  system , natm , ntype , rx , ry , rz , atype , atypei , natmi , simu_cell
  USE io,                  ONLY :  ionode , stdout

  implicit none

  ! global
  real(kind=dp) :: hess(3 * natm,3 * natm)
  real(kind=dp) :: deig(3 * natm)
  integer :: kunit
  ! local
  integer :: ia , igconf , it  , i , j
  real(kind=dp), dimension(:), allocatable :: rrx , rry , rrz
  real(kind=dp), dimension(:), allocatable :: qx , qxd , xsq
  real(kind=dp) :: dexpo, omeg, G1

  allocate( rrx ( natm ) , rry ( natm ) , rrz ( natm ) )
  allocate( qx ( 3 * natm ) , qxd ( 3 * natm ) , xsq ( 3 * natm ) )  ! normal coordinates

  ! ===============================
  !  inverse of eigenvectors matrix
  ! ===============================
  !CALL DGETRF( 3 * natm, 3 * natm, hess, 3 * natm, ipiv, info )
  !CALL DGETRI( 3 * natm, hess, 3 * natm, ipiv, work, lwork, info )

  ! ===========================================
  !  transformations to get normal coordinates
  ! ===========================================
!  qx = 0.0_dp
!  do i = 1,3 * natm
!    do ja = 1 , natm
!      qx ( imod ) = qx ( imod ) + &
!      hess ( ja , imod ) * rx ( ja ) + &
!      hess ( ja + natm , imod ) * ry ( ja ) + &
!      hess ( ja + 2 * natm , imod) * rz ( ja )
!    enddo
!  enddo

  ! ===========================================
  !  following MgO mauri paper 
  !  xsq(i): width of gaussian distribution 
  !  on q(i) coordinates (normal coordinates) 
  !  of mode i
  ! ===========================================

  do imod = 4,3 * natm ! loop over modes
    omeg  = SQRT ( deig ( imod ) )
    dexpo = EXP  ( - omeg / tempmod )
    xsq ( imod ) = 1.0_dp + dexpo
    xsq ( imod ) = xsq ( imod ) / ( 2.0_dp * omeg * ( 1.0_dp - dexpo ) )
    WRITE ( 6000 , * ) imod, xsq ( imod )
  enddo

  mode : do imod = 4 , 3 * natm
  ! =============================================
  !  generate ngconf configurations for each imod
  ! =============================================
  gconf : do igconf = 1 , ngconf

    io_node WRITE ( stdout ,'(a,i6,a,i6)') 'conf = ',igconf,' of mode ',imod

        ! ========================================================================
        ! gaussian distribution: 
        !   mean value: equilibrium position rx,ry,rz - > qx, qy, qz 
        !   width x^2: equation (11) of   J. Phys. Chem. B 2005, 109, 7245 - 7250
        ! ========================================================================
        CALL boxmuller_basic(G1, 0.0_dp, xsq(imod)) 
        qxd(imod) = G1

      ! ==========================
      ! inverse transformations 
      ! to get cartesian coordinates
      ! ==========================
      rrx = 0.0_dp
      rry = 0.0_dp
      rrz = 0.0_dp
      do ia = 1 , natm
          rrx(ia) = rrx(ia) + hess(ia           ,imod) * qxd(imod)
          rry(ia) = rry(ia) + hess(ia + natm    ,imod) * qxd(imod)
          rrz(ia) = rrz(ia) + hess(ia + 2 * natm,imod) * qxd(imod)
      enddo

      ! ===========================================================
      !            write configuration to file 
      ! ===========================================================
      if ( ionode ) then
        WRITE ( kunit , '(i6)' ) natm
        WRITE ( kunit , '(a)' ) system
        do i = 1 , 3
          WRITE ( kunit , '(3f20.12)' ) (simu_cell%A(i,j),j=1,3)
        enddo
        WRITE ( kunit , '(i4)' ) ntype
        WRITE ( kunit , * ) ( atypei ( it ) , it = 1 , ntype )
        WRITE ( kunit , * ) ( natmi  ( it ) , it = 1 , ntype )
 
        do ia = 1 , natm
          WRITE ( kunit , 200 ) atype(ia), rx(ia) + rrx(ia),ry(ia) + rry(ia),rz(ia) + rrz(ia) , &
                                           dzero , dzero ,dzero , &
                                           dzero , dzero ,dzero 
        enddo
      endif

    enddo gconf

  enddo mode

  deallocate(rrx,rry,rrz)
  deallocate(qx,qxd,xsq)

  return 

 200 FORMAT(A2,9E20.12)

END SUBROUTINE generate_modes



! *********************** SUBROUTINE band **************************************
!
!> \brief
!! calculates dispersion curve (band) in a given direction
!
! ******************************************************************************
SUBROUTINE band ( hess )

  USE config,           ONLY :  natm , rx , ry , rz , itype , simu_cell , ntype
  USE control,          ONLY :  calc
  USE constants,        ONLY :  tpi , pi
  USE io,          ONLY :  ionode , stdout , stderr , kunit_DOSKFF
  USE field,            ONLY :  rcutsq , sigsq , epsp , fc , uc 
  USE kspace,           ONLY :  kpath , get_kpath
  USE cell,             ONLY :  dirkar , kardir
  USE time,             ONLY :  bandtimetot

  implicit none

  ! global  
  real(kind=dp)    :: hess(3 * natm,3 * natm)
  ! local
  integer          :: ik , ia , ja , p1 , p2 , ierr
  integer*4        :: info ,  lwork
  real(kind=dp)    :: rxi , ryi , rzi 
  real(kind=dp)    :: rxij , ryij , rzij , rijsq
  real(kind=dp)    :: sxij, syij, szij
  real(kind=dp)    :: sphase
  real(kind=dp)    :: kx , ky , kz , kr
  real(kind=dp), allocatable  :: tmphess ( : , : ) , hessij ( : , : )
  real(kind=dp)    :: ww ( 3 ) , work ( 9 )
  real(kind=dp)    :: ttt1 , ttt2 
  TYPE (kpath)     :: kp
  character(len=1) :: jobz , uplo

#ifdef MPI
  ttt1 = MPI_WTIME(ierr) ! timing info
#endif

  ! =====================================================
  ! WARNING : 
  ! we assume ntype to be 
  ! the number of atoms in the primitiv cell 
  ! in case of several inequivalent position 
  ! of the same type, one could define a further 
  ! identical type .
  ! example :
  ! let say there is 3 atoms in the primitiv cell
  ! 2 differents types
  ! in POSFF :
  ! --------------------------------------------------
  ! 3
  ! cell  
  ! 3 
  ! A A B 
  ! 1 1 1
  ! A 0.0 0.0 0.0    0.0 0.0 0.0   0.0 0.0 0.0
  ! A 0.2 0.2 0.2    0.0 0.0 0.0   0.0 0.0 0.0
  ! B 0.5 0.5 0.5    0.0 0.0 0.0   0.0 0.0 0.0
  ! =====================================================
  allocate ( tmphess( 3 * ntype , 3 * ntype ) )
  allocate ( hessij ( 3 * ntype , 3 * ntype ) )

  ! ========================================
  !  generate kpath in structure kp
  ! ========================================
  allocate ( kp%kptx ( nkphon ) , kp%kpty ( nkphon ) , kp%kptz ( nkphon ) )
  CALL get_kpath ( ks , kf , nkphon , path , kp )  

  OPEN (UNIT = kunit_DOSKFF ,FILE = 'DOSKFF')

  lwork = 9
  jobz = 'V'
  uplo = 'U'

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  ! ==================================
  !   loop on kpath index
  ! ==================================
  do ik = 1 , kp%nk

    tmphess = 0.0_dp
    hessij  = 0.0_dp
    work    = 0.0_dp
    ww      = 0.0_dp

    ! direct to kartesian times 2 * pi 
    kx = tpi * ( kp%kptx(ik) * simu_cell%B(1,1) +  kp%kpty(ik) * simu_cell%B(1,2) + kp%kptz(ik) * simu_cell%B(1,3) ) * ncell
    ky = tpi * ( kp%kptx(ik) * simu_cell%B(2,1) +  kp%kpty(ik) * simu_cell%B(2,2) + kp%kptz(ik) * simu_cell%B(2,3) ) * ncell
    kz = tpi * ( kp%kptx(ik) * simu_cell%B(3,1) +  kp%kpty(ik) * simu_cell%B(3,2) + kp%kptz(ik) * simu_cell%B(3,3) ) * ncell

    do ia = 1 , natm

      rxi = rx(ia)
      ryi = ry(ia)
      rzi = rz(ia)

      do ja = 1 , natm

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
        p1 = itype(ia)
        p2 = itype(ja)

        if (rijsq .lt. rcutsq(p1,p2)) then
          kr =  kx * rxij + ky * ryij + kz * rzij
          sphase = SIN ( 0.5_dp * kr )
          sphase = sphase * sphase

          tmphess ( 1 , 1 ) = tmphess ( 1 , 1 ) + hess ( ja            , ia            ) * sphase 
          tmphess ( 2 , 2 ) = tmphess ( 2 , 2 ) + hess ( ja + natm     , ia + natm     ) * sphase
          tmphess ( 3 , 3 ) = tmphess ( 3 , 3 ) + hess ( ja + 2 * natm , ia + 2 * natm ) * sphase
          tmphess ( 1 , 2 ) = tmphess ( 1 , 2 ) + hess ( ja            , ia + natm     ) * sphase
          tmphess ( 1 , 3 ) = tmphess ( 1 , 3 ) + hess ( ja            , ia + 2 * natm ) * sphase
          tmphess ( 2 , 3 ) = tmphess ( 2 , 3 ) + hess ( ja + natm     , ia + 2 * natm ) * sphase

        endif
      enddo  ! j atom loop
    enddo ! i atom loop

    hessij ( 1 , 1 ) = - tmphess ( 1 , 1 ) * 2.0_dp / natm
    hessij ( 2 , 2 ) = - tmphess ( 2 , 2 ) * 2.0_dp / natm
    hessij ( 3 , 3 ) = - tmphess ( 3 , 3 ) * 2.0_dp / natm

    hessij ( 1 , 2 ) = - tmphess ( 1 , 2 ) * 2.0_dp / natm
    hessij ( 1 , 3 ) = - tmphess ( 1 , 3 ) * 2.0_dp / natm
    hessij ( 2 , 3 ) = - tmphess ( 2 , 3 ) * 2.0_dp / natm
 
    hessij ( 2 , 1 ) = hessij ( 1 , 2 )
    hessij ( 3 , 1 ) = hessij ( 1 , 3 )
    hessij ( 3 , 2 ) = hessij ( 2 , 3 )

    CALL DSYEV(jobz,uplo,3*ntype,hessij,3*ntype,ww,work,lwork,info)
    if ( info .ne. 0 ) then
      io_node WRITE ( stderr ,'(a,i5)') 'ERROR in band : improper termination. of DSYEV info = ',info
      STOP
    endif


    if ( MOD ( ik + 1 , 10 ) .eq. 0 ) &
    WRITE ( stdout ,'(i7,3f12.4,3f18.6)')        ik, kp%kptx ( ik ) , kp%kpty ( ik ) , kp%kptz ( ik ) , SQRT (ww(1)),SQRT (ww(2)), SQRT (ww(3))
    WRITE ( kunit_DOSKFF ,'(i7,3f12.4,3f18.6)')  ik, kp%kptx ( ik ) , kp%kpty ( ik ) , kp%kptz ( ik ) , SQRT (ww(1)),SQRT (ww(2)), SQRT (ww(3))
  enddo

  CLOSE ( kunit_DOSKFF )

  deallocate ( kp%kptx , kp%kpty , kp%kptz  )
  deallocate ( tmphess)
  deallocate ( hessij )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

#ifdef MPI
  ttt2 = MPI_WTIME ( ierr ) ! timing info
  bandtimetot = bandtimetot + ( ttt2 - ttt1) 
#endif

  return  

END SUBROUTINE band 

! *********************** SUBROUTINE doskpt ************************************
!
!> \brief
!! calculates the DOS for a given set of k-points
!
! ******************************************************************************
SUBROUTINE doskpt ( hess , eigenk , km )

  USE config,           ONLY :  natm , rx , ry , rz , itype , simu_cell , ntype
  USE control,          ONLY :  calc
  USE constants,        ONLY :  tpi , pi
  USE io,          ONLY :  ionode , stdout , stderr , kunit_IBZKPTFF , kunit_DKFF
  USE field,            ONLY :  rcutsq , sigsq , epsp , fc , uc 
  USE cell,             ONLY :  kardir , dirkar
  USE time,             ONLY :  doskpttimetot
  USE kspace,           ONLY :  kmesh
  USE cell,             ONLY :  dirkar_1

  implicit none

  ! global  
  TYPE(kmesh)                  :: km 
  real(kind=dp) , intent (in)  :: hess  (3 * natm, 3 * natm )
  real(kind=dp) , intent (out) :: eigenk( km%nk  , 3 * ntype)

  ! local
  integer          :: wi , ierr , si , sj , im
  integer          :: ik , ia , ja , p1 , p2
  real(kind=dp)    :: rxi , ryi , rzi
  real(kind=dp)    :: rxij , ryij , rzij , rijsq
  real(kind=dp)    :: sxij , syij , szij
  real(kind=dp)    :: sphase , ncell_tpi
  real(kind=dp)    :: kr 
  real(kind=dp)    :: kx , ky , kz
  real(kind=dp)    :: ttt1 , ttt2 
  real(kind=dp)    :: ttt1l , ttt2l ! inside loop
  real(kind=dp) , allocatable  :: tmphess( : , : ) , hessij ( : , : ) , ww(:)
  real(kind=dp)    :: work(9*ntype)
  character(len=1) :: jobz, uplo
  integer*4        :: info, lwork

  ! trash
  !integer          :: iiii

  ! =====================================================
  ! WARNING : 
  ! we assume ntype to be 
  ! the number of atoms in the primitiv cell 
  ! in case of several inequivalent position 
  ! of the same type, one could define a further 
  ! identical type .
  ! example :
  ! let say there is 3 atoms in the primitiv cell
  ! 2 differents types
  ! in POSFF :
  ! --------------------------------------------------
  ! 3
  ! cell  
  ! 3 
  ! A A B 
  ! 1 1 1
  ! A 0.0 0.0 0.0    0.0 0.0 0.0   0.0 0.0 0.0
  ! A 0.2 0.2 0.2    0.0 0.0 0.0   0.0 0.0 0.0
  ! B 0.5 0.5 0.5    0.0 0.0 0.0   0.0 0.0 0.0
  ! =====================================================
  allocate ( tmphess ( 3 * ntype , 3 * ntype ) )
  allocate ( hessij  ( 3 * ntype , 3 * ntype ) )
  allocate ( ww      ( 3 * ntype             ) )

#ifdef MPI
  ttt1 = MPI_WTIME(ierr) ! timing info
#endif

  lwork = 9*ntype
  jobz = 'N'
  uplo = 'U'

  ! some constants of the loop 
  ncell_tpi = REAL( ncell , kind = dp ) * tpi

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  wi = 1
        
#ifdef MPI
  ttt1l = MPI_WTIME(ierr) ! timing info
#endif

  do ik = km%kpt_dec%istart , km%kpt_dec%iend

    kx = km%kptx(ik) * REAL( ncell , kind = dp )
    ky = km%kpty(ik) * REAL( ncell , kind = dp )
    kz = km%kptz(ik) * REAL( ncell , kind = dp )

    tmphess = 0.0_dp
    hessij  = 0.0_dp
    work    = 0.0_dp
    ww      = 0.0_dp

    do ia = 1 , natm

      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )

      do ja = 1 , natm

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

        p1 = itype(ia)
        p2 = itype(ja)

        if ( rijsq .lt. rcutsq(p1,p2) ) then
          si = 3 * p1 - 3
          sj = 3 * p2 - 3
          kr =  kx * rxij + ky * ryij + kz * rzij
          sphase = SIN ( 0.5_dp * kr )
          sphase = sphase * sphase

          tmphess ( 1 + si , 1 + sj ) = tmphess ( 1 + si , 1 + sj ) + hess ( ja            , ia            ) * sphase
          tmphess ( 2 + si , 2 + sj ) = tmphess ( 2 + si , 2 + sj ) + hess ( ja + natm     , ia + natm     ) * sphase
          tmphess ( 3 + si , 3 + sj ) = tmphess ( 3 + si , 3 + sj ) + hess ( ja + 2 * natm , ia + 2 * natm ) * sphase
          tmphess ( 1 + si , 2 + sj ) = tmphess ( 1 + si , 2 + sj ) + hess ( ja            , ia + natm     ) * sphase
          tmphess ( 1 + si , 3 + sj ) = tmphess ( 1 + si , 3 + sj ) + hess ( ja            , ia + 2 * natm ) * sphase
          tmphess ( 2 + si , 3 + sj ) = tmphess ( 2 + si , 3 + sj ) + hess ( ja + natm     , ia + 2 * natm ) * sphase
          tmphess ( 2 + sj , 1 + si ) = tmphess ( 2 + sj , 1 + si ) + hess ( ia + natm     , ja            ) * sphase
          tmphess ( 3 + sj , 1 + si ) = tmphess ( 3 + sj , 1 + si ) + hess ( ia + 2 * natm , ja            ) * sphase
          tmphess ( 3 + sj , 2 + si ) = tmphess ( 3 + sj , 2 + si ) + hess ( ia + 2 * natm , ja + natm     ) * sphase

        endif
      enddo  ! j atom loop
    enddo ! i atom loop
   
#ifdef MPI
    if ( MOD ( ik + 1 , km%kpt_dec%dim_data / 20 ) .eq. 0 ) then
      ttt2l = MPI_WTIME(ierr) ! timing info
      io_node WRITE ( stdout , 110 ) ' kpoints : ',int(ik*numprocs),' DOSKPT  ', ttt2l - ttt1l
      ttt1l = MPI_WTIME(ierr) ! timing info
    endif
#endif

    hessij = - tmphess * 2.0_dp / natm
    
    ! diagonalisation of the 3x3 matrix
    CALL DSYEV(jobz,uplo,3*ntype,hessij,3*ntype,ww,work,lwork,info)
    if ( info .ne. 0 ) then
      io_node WRITE ( stderr ,'(a,i5)') 'ERROR in doskpt : improper termination. of DSYEV info = ', info
      STOP 
    endif


    ! remove for the parallel version it was probably wrong anyway
    !CALL dirkar_1 ( kx , ky , kz , simu_cell%B , ncell )
    !do i = 1 , wi
     ! if ( MOD ( ik , km%nk / 20 ) .eq. 0 ) then
     !   if ( ALL(ww.ge.0.0d0,3*ntype)) then
     !     WRITE ( stdout ,'(i7,3f12.4,<3*ntype>f18.6)')    ck, kx, ky, kz, ( SQRT ( ww ( im ) ) , im = 1 , 3 * ntype )
     !   else
     !     WRITE ( stdout ,'(a)')    'WARNING imaginary frequency'
     !     WRITE ( stdout ,'(i7,3f12.4,<3*ntype>f18.6)')    ck, kx, ky, kz, ( ww ( im )  , im = 1 , 3 * ntype )
     !     STOP
     !   endif
     !   
     ! endif
    !    if ( ALL(ww.ge.0.0d0,3*ntype)) then
    !      WRITE ( kunit_DKFF ,'(i7,3f12.4,<3*ntype>f18.6)')    ck, kx, ky, kz, ( SQRT ( ww ( im ) ) , im = 1 , 3 * ntype )
    !    else
    !      WRITE ( stdout ,'(a)')    'WARNING imaginary frequency'
    !      WRITE ( stdout ,'(i7,3f12.4,<3*ntype>f18.6)')    ck, kx, ky, kz, ( ww ( im )  , im = 1 , 3 * ntype )
    !      STOP
    !    endif
    !  ck = ck + 1
    !enddo

   eigenk ( ik , : ) = ww(:)

  enddo !ik loop

  do im = 1 , 3 * ntype
    CALL MPI_ALL_REDUCE_DOUBLE ( eigenk ( : , im ) , km%nk ) 
  enddo

!  CLOSE ( kunit_IBZKPTFF )
  !CLOSE ( kunit_DKFF )

  deallocate ( tmphess )
  deallocate ( hessij  )
  deallocate ( ww  )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

#ifdef MPI
  ttt2 = MPI_WTIME(ierr) ! timing info
  doskpttimetot = doskpttimetot + ( ttt2 - ttt1 )
#endif

110   FORMAT(2X,A8,I8,A20,' :  cpu time',F9.2)

  return  

END SUBROUTINE doskpt 

! *********************** SUBROUTINE gene_IBZKPTFF *****************************
!
!> \brief
!! calculates the DOS for a given set of k-points
!
! ******************************************************************************
SUBROUTINE write_IBZKPTFF ( km ) 

  USE constants,        ONLY :  pi
  USE io  ,             ONLY :  ionode, kunit_IBZKPTFF
  USE kspace,           ONLY :  kmesh

  implicit none
  ! global
  TYPE ( kmesh ), intent(in) :: km

  ! local
  integer :: ik 
  integer :: wi , nktot

  wi = 1
  nktot = ( ( nkphon + 1 ) *( nkphon + 1 ) * ( nkphon + 1 ) )

  OPEN (UNIT = kunit_IBZKPTFF ,FILE = 'IBZKPTFF')

  if ( ionode ) then
    WRITE (kunit_IBZKPTFF,*) 'Automatically generated mesh FF'
    WRITE (kunit_IBZKPTFF,*)  nktot
    WRITE (kunit_IBZKPTFF,*) 'Reciprocal lattice'
  endif

  do ik = 1 , nktot
    io_node WRITE (kunit_IBZKPTFF,'(3f16.12)') km%kptx(ik) , km%kpty(ik) , km%kptz(ik) 
  enddo

  CLOSE(kunit_IBZKPTFF)

  return

END SUBROUTINE write_IBZKPTFF


END MODULE vib

! ===== fmV =====
