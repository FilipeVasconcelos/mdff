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
MODULE pim

  USE constants,                        ONLY :  dp 
  USE config,                           ONLY :  ntypemax 
  USE io,                               ONLY :  ionode, stdin, stdout, ioprintnode
  USE mpimdff

  implicit none

  logical, PRIVATE :: symmetric_pot


  ! ============================================================  
  !                    force field type info 
  ! ============================================================  

  ! type dependent properties
  real(kind=dp)    :: poldip   ( ntypemax , 3 , 3 )    !< dipole     polarizability if ldip_polar( it ) = .true. 
  real(kind=dp)    :: poldip_iso ( ntypemax )          !< isotropic dipole polarizability if ldip_polar( it ) = .true.
  real(kind=dp)    :: polquad  ( ntypemax , 3 , 3 , 3 )!< quadrupole polarizability if lquad_polar( it ) = .true.
  real(kind=dp)    :: polquad_iso ( ntypemax )         !< isotropic quadrupole polarizability if ldip_polar( it ) = .true. 


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


  character(len=4) :: algo_ext_dipole                  !< set the algorithm used to get induced moments from polarization 
  character(len=4) :: algo_ext_dipole_allowed(2)
  data                algo_ext_dipole_allowed       / 'poly','aspc'/

  character(len=11) :: algo_moment_from_pola            !< set the algorithm used to get induced moments from polarization 
  character(len=11) :: algo_moment_from_pola_allowed(6) !< set the algorithm used to get induced moments from polarization 
  data                 algo_moment_from_pola_allowed / 'scf' , 'scf_kO_v1' , 'scf_kO_v2' , 'scf_kO_v3' , 'scf_kO_v4_1' , 'scf_kO_v4_2'/  !! scf ( self consistent ) 
  
  logical, SAVE     :: ldip_wfc          !< calculate electrostatic contribution from dipolar momemt coming from wfc
  logical, SAVE     :: lwrite_dip_wfc    !< write dipoles from wannier centers to file
  integer           :: lwfc     ( ntypemax )            !< moment from wannier centers 
  real(kind=dp)     :: rcut_wfc                         !< radius cut-off for WFs searching




CONTAINS

! *********************** SUBROUTINE field_default_tag *************************
!> \brief
!! set default values to field tags
! ******************************************************************************
SUBROUTINE pim_default_tag

  implicit none

  ! =================
  !  default values
  ! =================

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
  symmetric_pot = .true.

  omegakO               = 0.7_dp 

  ! wannier centers related
  lwfc           = 0
  ldip_wfc       = .true.            
  rcut_wfc       = 0.5_dp

  return

END SUBROUTINE pim_default_tag


! *********************** SUBROUTINE field_check_tag ***************************
!> \brief
!! check field tag values
! ******************************************************************************
SUBROUTINE pim_check_tag

  USE coulomb,                  ONLY :  task_coul, qch
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

  if ( any(qch .ne.0.0_dp) ) lqch =.true. !charge
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

  ! =================================================
  ! symetrization of input potentials !
  ! =================================================
  if ( symmetric_pot ) then
    do it = 1 , ntype
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

  ! =======================
  !  algo_moment_from_pola 
  ! =======================
  CALL check_allowed_tags( size( algo_moment_from_pola_allowed), algo_moment_from_pola_allowed, algo_moment_from_pola, 'in fieldtag','algo_moment_from_pola' ) 



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

END SUBROUTINE pim_check_tag

! *********************** SUBROUTINE field_init ********************************
!> \brief
!! force field initialisation
! ******************************************************************************
SUBROUTINE pim_init

  USE config,                   ONLY :  ntype
!  USE control,                  ONLY :  calc , lnmlj , lcoulomb , lmorse , longrange, lnon_bonded

  implicit none

  ! local
  character(len=132)    :: filename
  integer               :: it, ioerr
  logical               :: linduced, ldamp

  namelist /pimtag/      symmetric_pot , &
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
                         omegakO       , &
                         lwfc          , &            
                         lwrite_dip_wfc, &            
                         ldip_wfc      , &            
                         rcut_wfc      , &            
                         ldip_polar    , &
                         ldip_damping  , &
                         lquad_polar   
  
  ! ================================
  ! defaults values for field tags 
  ! ================================
  CALL pim_default_tag

  ! ================================
  ! read field tags
  ! ================================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
    READ ( stdin , pimtag, iostat=ioerr)
    if ( ioerr .lt. 0 )  then
      io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : coulombtag section is absent'
      STOP
    elseif ( ioerr .gt. 0 )  then
      io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : coulombtag wrong tag',ioerr
      STOP
    endif
  CLOSE ( stdin )

   linduced = .false.
   do it = 1 , ntype
     if ( ldip_polar(it) )  linduced = .true.
   enddo
   ldamp = .false.
   if ( any ( ldip_damping )  )  ldamp = .true.

   if ( .not. linduced ) return

  ! ================================
  ! check field tags values
  ! ================================
  CALL pim_check_tag

END SUBROUTINE pim_init

! *********************** SUBROUTINE field_print_info **************************
!> \brief
!! print force field information to standard output
! ******************************************************************************
SUBROUTINE pim_print_info ( kunit )

  USE config,           ONLY :  ntype, atypei

  implicit none

  ! global
  integer :: kunit
  integer :: it, it1, it2, i, j 
  logical :: linduced, ldamp


  linduced = .false.
  do it = 1 , ntype
    if ( ldip_polar(it) )  linduced = .true.
  enddo
  ldamp = .false.
  if ( any ( ldip_damping )  )  ldamp = .true.

  if ( ionode ) then
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
   
  endif

  return

END SUBROUTINE pim_print_info


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

END MODULE pim 
! ===== fmV =====
