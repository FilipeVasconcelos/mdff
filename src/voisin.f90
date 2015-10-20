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
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
! USA.
! ===== fmV =====
#include "symbol.h"
! ======= Hardware =======
!#define debug_vois1_fixed
!#define debug_vois1_sann
!#define debug_vois1_voronoi
!#define sort_voronoi
!#define further_info_voronoi    
! ======= Hardware =======

! ******************************************************************************
! voisin module : voisin ( vois1 )  means neighbour in french
! in construction 
!
! main reference :
![1] JChem_Phys_136_234107.pdf
!
! I follow [1] to implement three different method to indentify
! nearest neighbors.
!
! - fixed_distance cutoff
! - voronoi constrution 
! - SANN ( derived and presented in [1]

! ******************************************************************************
MODULE voisin

  USE constants,                        ONLY :  dp
  USE config,                           ONLY :  ntypemax
  USE mpimdff

  implicit none

  integer , PARAMETER ::  nmaxneigh = 500

  integer       :: nconf                            ! number of configurations readed for vois1 analysis (only when calc = 'vois1') 
  real(kind=dp) :: cutvois1 ( ntypemax , ntypemax ) ! cutoff of the first neighbour sphere for each pair of types

  ! ===============
  !  distributions
  ! ===============
  integer       , dimension(:,:)  , allocatable :: dib_nb 
  integer                                         :: nbmax
  integer       , dimension (:)   , allocatable   :: kk 

  ! =========================================
  !     neighbor algorithm
  ! =========================================
  character(len=60) :: vois1algo
  character(len=60) :: vois1algo_allowed(4)
  data                 vois1algo_allowed  / 'fixed' , 'sann' , 'sannsq' , 'voronoi' /

CONTAINS


! *********************** SUBROUTINE vois1_init ********************************
! ******************************************************************************
SUBROUTINE vois1_init

  USE io,                       ONLY :  ionode , stdin , stdout , stderr
  USE control,                  ONLY :  calc

  implicit none

  ! local
  integer            :: ioerr
  character(len=132) :: filename

  namelist /vois1tag/  vois1algo        , &
                       nconf            , &
                       nbmax            , & 
                       cutvois1
                      
  if ( calc .ne. 'vois1' ) return

  ! ====================
  !  set default values
  ! ====================
  CALL vois1_default_tag
  ! =====================
  !  reads vois1tag tags
  ! =====================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , vois1tag , iostat=ioerr)
  if ( ioerr .lt. 0 )  then
    io_node WRITE ( stderr, '(a)') 'ERROR reading input_file : vois1tag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node WRITE ( stderr, '(a,i8)') 'ERROR reading input_file : vois1tag wrong tag'
    STOP
  endif

  CLOSE ( stdin )
  ! ===============
  !  check vois1tag
  ! ===============
  CALL vois1_check_tag
  ! =============
  !  print info
  ! =============


END SUBROUTINE vois1_init

! *********************** SUBROUTINE vois1_default_tag *************************
! ******************************************************************************
SUBROUTINE vois1_default_tag

  implicit none

  ! =================
  !  default values
  ! =================
  nconf    = 0
  cutvois1 = 0.0_dp

  return

END SUBROUTINE vois1_default_tag

! *********************** SUBROUTINE vois1_check_tag ***************************
! ******************************************************************************
SUBROUTINE vois1_check_tag

  USE io,                  ONLY :  ionode , stdout 

  implicit none

  ! local
  logical :: allowed
  integer :: i

  allowed = .false.
  do i = 1 , size ( vois1algo_allowed )
    if ( TRIM( vois1algo ) == vois1algo_allowed(i) ) allowed = .TRUE.
  enddo
  if ( .not. allowed ) then
    io_node WRITE ( stdout ,'(a,a)') 'ERROR vois1tag: vois1algo should be ',vois1algo_allowed
    STOP
  endif

  return


END SUBROUTINE vois1_check_tag

! *********************** SUBROUTINE vois1_print_info **************************
! ******************************************************************************
SUBROUTINE vois1_print_info(kunit)

  USE config,                   ONLY :  ntype , atypei
  USE io,                       ONLY :  ionode

  implicit none

  ! global 
  integer :: kunit
  ! local
  integer :: it1,it2 

  if ( ionode ) then
      separator(kunit)
      blankline(kunit)
      WRITE ( kunit ,'(a)')       'neighbour analysis'
      WRITE ( kunit ,'(a)')       'read config from file : TRAJFF'
      WRITE ( kunit ,'(a,i10)')   'numbers of config read nconf = ',nconf
      if ( vois1algo .eq. 'fixed' ) then
        WRITE ( kunit ,'(a,f12.5)') 'fixed distance algorithm'
        WRITE ( kunit ,'(a,f12.5)') 'first neighbour sphere cutoff radius for each types ' 
        do it1 = 1 , ntype
          do it2 = it1 , ntype
            WRITE ( kunit ,'(2a,f12.5)') atypei(it1),atypei(it2),cutvois1(it1,it2) 
          enddo
        enddo
      else if ( vois1algo .eq. 'sann' .or. vois1algo .eq. 'sannsq') then
        WRITE ( kunit ,'(a,f12.5)') 'the solid-angle based nearest-neighbor algorithm (SANN)'
        WRITE ( kunit , '(a)'     ) 'J. Chem. Phys. 136, 234107 (2012); http://dx.doi.org/10.1063/1.4729313 (12 pages)'
        WRITE ( kunit , '(a)'     ) 'Jacobus A. van Meel, Laura Filion, Chantal Valeriani, and Daan Frenkel'
        if ( vois1algo .eq. 'sannsq' ) then
          WRITE ( kunit , '(a)'     ) '"square" (S^2) version : to get second nearest neighbour'
        endif
      endif
  endif

  return

END SUBROUTINE vois1_print_info

! *********************** SUBROUTINE vois1_driver ******************************
! ******************************************************************************
SUBROUTINE vois1_driver

  USE io,                       ONLY :  ionode , stdout , stderr , kunit_TRAJFF , kunit_DTNBFF , kunit_VOIS1FF
  USE config,                   ONLY :  system , natm , ntype , atype , rx , ry , rz , itype , & 
                                        atypei , natmi , rho , simu_cell , config_alloc , &
                                        config_print_info , coord_format_allowed , atom_dec , read_traj_header , read_traj 
  USE cell,                     ONLY :  lattice, dirkar
  USE control,                  ONLY :  itraj_format , trajff_data

  implicit none

  ! local 
  integer            :: i , inb , iconf ,it
  character(len=20) :: FMT
  dectime


  OPEN ( UNIT = kunit_VOIS1FF , FILE = 'VOIS1FF')

  ! get systems info from TRAJFF
  if ( itraj_format .ne. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' )
  if ( itraj_format .eq. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' , form = 'unformatted')
  CALL read_traj_header( kunit_TRAJFF , itraj_format )
  if ( itraj_format .ne. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' )
  if ( itraj_format .eq. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' , form = 'unformatted')

  CALL lattice ( simu_cell )
  rho = natm / simu_cell%omega
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc
  CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms')
  CALL vois1_alloc
  CALL typeinfo_init
  ! =============
  !  print info
  ! =============
  CALL vois1_print_info ( stdout )

  CALL config_print_info ( stdout )

  ! ============================================
  ! LOOP OVER CONFIGURATIONS 
  ! ============================================
  allocate ( kk ( 0:nbmax ) ) 
  kk = 0
  do iconf = 1, nconf
#ifdef MPI
    statime
#endif

    ! ====================
    !  read config iconf
    ! ====================
    CALL read_traj ( kunit_TRAJFF , itraj_format , trajff_data )

    CALL lattice ( simu_cell )

    !CALL distance_tab 

    ! ======================================
    ! fixed-distance method :
    !  cutvois1 ( it ) should be set
    ! ======================================
    if ( vois1algo .eq. 'fixed' ) then
      CALL fixed_distance
    ! ======================================
    ! 3D voronoi construction ( Delaunay triangulation ) 
    !  Allen-Tildesley implementation
    ! ======================================
    else if ( vois1algo .eq. 'voronoi' ) then
      CALL voronoi_construction 
    ! ======================================
    ! solid-angle based nearest-neighbor algorithm (SANN)
    ! J. Chem. Phys. 136, 234107 (2012)
    ! Jacobus A. van Meel et al.
    ! parameter-free method
    ! ======================================
    else if ( vois1algo .eq. 'sann' .or. vois1algo .eq. 'sannsq' ) then
      CALL sann 
    endif

#ifdef MPI
    stotime
    writime('config : ',' VOIS1  ',iconf)
#endif

  enddo
  do i=0,nbmax
    write(20000,*) i,kk(i)
  enddo
  deallocate ( kk ) 

  ! ===================
  !  write distrib Nb 
  ! ===================
  OPEN (UNIT = kunit_DTNBFF , FILE = 'DTNBFF')
  WRITE( kunit_DTNBFF , '(a)') '#nb    P(Nb)'
  do inb = 0 , nbmax
#ifdef GFORTRAN
     WRITE ( FMT , * ) ntype + 1
     WRITE( kunit_DTNBFF , '(i6,'// ADJUSTL(FMT) //'f12.6)' ) inb , ( REAL ( dib_nb ( it , inb ) ) / ( REAL ( nconf , kind = dp ) ) , it = 0 , ntype )
#else
     WRITE( kunit_DTNBFF , '(i6,<ntype+1>f12.6)' ) inb , ( REAL ( dib_nb ( it , inb ) ) / ( REAL ( nconf , kind = dp ) ) , it = 0 , ntype )
#endif
  enddo
  CLOSE ( kunit_DTNBFF )

  CLOSE(kunit_TRAJFF)
  CLOSE(kunit_VOIS1FF)

  CALL vois1_dealloc

  return 

END SUBROUTINE vois1_driver

! *********************** SUBROUTINE fixed_distance ****************************
! ******************************************************************************
SUBROUTINE fixed_distance

  USE config,                   ONLY :  natm , ntype , rx , ry, rz , simu_cell , atype , itype , atypei 
  USE cell,                     ONLY :  kardir , dirkar , dirkar_1
  USE io,                       ONLY :  ionode , stdout , stderr , kunit_VOIS1FF

  implicit none

  ! local
  integer :: ia , ja , it , jt , im , jnb , jtnb , kkkk
  real(kind=dp) :: dist 
  real(kind=dp) :: rxi ,ryi ,rzi
  real(kind=dp) :: rxij ,ryij ,rzij
  real(kind=dp) :: sxij ,syij ,szij
  integer       , dimension (:)   , allocatable  :: natmi_cluster 
  character(len=2), dimension (:) , allocatable  :: atypei_cluster 

  integer , dimension (:)   , allocatable  :: Nb 
  integer , dimension (:,:) , allocatable  :: selectedneighbors ! neighbour table
  integer , dimension (:,:) , allocatable  :: spec_vois1 ! neighbour table
  real(kind=dp)                                  :: xxx , yyy , zzz 
  character(len=20) :: FMT

  ! distributions
  integer :: knb
 
  ! ===============================================================================================
  ! definitions :
  !   - Nb (ia )                 : number of neighbours of ia (inside the cutoff defined for each type )
  !   - selectedneighbors ( ia , k )         : is the kth neighbour of ia
  !   - spec_vois1 ( ia , it )   : number of neighbours of type it of ia 
  !   - lab_count ( it , nb )    : number of local config with nb numbers of neighbours for each type it
  !   - label_loconfig ( it , nb , 
  ! ===============================================================================================
  allocate ( selectedneighbors ( natm , nmaxneigh ) )
  allocate ( Nb ( natm ) )
  selectedneighbors = 0

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

    do ia = 1 , natm
      Nb (ia ) = 1
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )

      do ja = 1 , natm 

        if ( ia .eq. ja ) cycle
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        dist = rxij * rxij + ryij * ryij + rzij * rzij
        if ( sqrt(dist) .lt. cutvois1(itype(ia),itype(ja)) ) then
           selectedneighbors ( ia , Nb ( ia ) ) = ja
           Nb ( ia ) = Nb ( ia ) + 1 
#ifdef debug_vois1_fixed
           WRITE ( stdout , '(a,i4,2x,a,i4,f16.8,2i8)' ) atype(ia),ia, atype(ja),ja , sqrt(dist) , selectedneighbors (ia , 0 ),Nb ( ia )
#endif
        endif
      enddo
    enddo

   ! Nb was used used as the index of selectedneighbors
   ! Nb of neighbours inside cutoff
   Nb = Nb - 1

#ifdef debug_vois1_fixed
  WRITE ( stdout , * )
  do ia = 1 , natm
    WRITE ( stdout , '(a,a,i6)' )       'atom      : ', atype(ia), ia 
    WRITE ( stdout , '(a,  i6)' )       '#nb vois1 = ', Nb ( ia ) 
    WRITE ( stdout , '(a,<Nb(ia)>i6)' ) 'neighb    : ', ( selectedneighbors (ia , im ), im = 1 , Nb ( ia ) )
  enddo
#endif 

  allocate ( natmi_cluster  ( 0: ntype ) )
  allocate ( atypei_cluster ( 0: ntype ) )
  do ia = 1 , natm
    it = itype ( ia )
    atypei_cluster(0) = 'C'//atype(ia)
    atypei_cluster(1) = 'A'
    atypei_cluster(2) = 'B'
    natmi_cluster     = 0
    do ja = 1 , Nb ( ia )
      jnb = selectedneighbors ( ia , ja )
      jtnb= itype ( jnb )
      natmi_cluster(jtnb) = natmi_cluster(jtnb) + 1
    enddo
    natmi_cluster(0)    = 1
    kkkk = 10000+Nb ( ia ) + 1
    kk(Nb ( ia ) + 1) = kk(Nb ( ia ) + 1) + 1
    WRITE ( kkkk , * ) Nb ( ia ) + 1
    WRITE ( kkkk , '(a)' ) 'cluster'
    WRITE ( kkkk , '(3f20.12)' )  1000.0_dp ,    0.0_dp ,    0.0_dp
    WRITE ( kkkk , '(3f20.12)' )     0.0_dp , 1000.0_dp ,    0.0_dp
    WRITE ( kkkk , '(3f20.12)' )     0.0_dp ,    0.0_dp , 1000.0_dp
    WRITE ( kkkk , '(i4)'      )  ntype+1
#ifdef GFORTRAN
    WRITE ( FMT , * ) ntype+1 
    WRITE ( kkkk , '('// ADJUSTL(FMT) //'a3 )') ( atypei_cluster(it) , it=0,ntype )
#else
    WRITE ( kkkk , '(<ntype+1>a3 )') ( atypei_cluster(it) , it=0,ntype )
#endif
    WRITE ( kkkk , *           ) ( natmi_cluster (it) , it=0,ntype )
    WRITE ( kkkk ,'(A)')         'Cartesian'
    WRITE ( kkkk ,'(a,3e20.12)') atypei_cluster(0) , 0.0_dp , 0.0_dp , 0.0_dp
    do it = 1 , ntype
      do ja = 1 , Nb ( ia )
        jnb = selectedneighbors ( ia , ja )
        jtnb= itype ( jnb )
        xxx = rx ( jnb ) - rx ( ia ) - NINT ( rx ( jnb ) - rx ( ia ) )
        yyy = ry ( jnb ) - ry ( ia ) - NINT ( ry ( jnb ) - ry ( ia ) )
        zzz = rz ( jnb ) - rz ( ia ) - NINT ( rz ( jnb ) - rz ( ia ) )
        CALL dirkar_1 ( xxx , yyy , zzz , simu_cell%A , 1 )
        if ( jtnb .eq. it ) WRITE ( kkkk ,'(a,3e20.12)') atype ( jnb ) , xxx , yyy , zzz
      enddo
    enddo
    !WRITE ( 10000 , '(a,i6,a,i6,a,<Nb(ia)>(i6,a))' ) 'neighbors of atom : ', ia
    !, ' #nb : ', Nb ( ia ), ' : ', ( selectedneighbors ( ia , ja ) , atype
    !(selectedneighbors ( ia , ja )) , ja = 1 , Nb ( ia ) )
  enddo
  deallocate ( natmi_cluster )
  deallocate ( atypei_cluster )



  ! =================================
  !         some statistic
  ! =================================
  ! on Nb
  do ia = 1 , natm
    it = itype  ( ia ) 
    knb = int( Nb(ia) )
    ! ====================== 
    !  test out of bound
    ! ====================== 
    if (knb.lt.0.or.knb.gt.nbmax) then
      io_node WRITE ( stderr , * ) 'ERROR: out of bound dist_nb'
      io_node WRITE ( stderr , '(4i6)') ia , knb , Nb(ia) , nbmax
      STOP
    endif
    dib_nb ( 0  , knb ) = dib_nb ( 0  , knb ) + 1
    dib_nb ( it , knb ) = dib_nb ( it , knb ) + 1
  enddo

  ! ==================================
  !       single speciation
  ! dimension of table speciation is :
  !         natm x ntype 
  ! ==================================
  allocate ( spec_vois1 ( natm , ntype ) )
  spec_vois1 = 0

  do ia = 1 , natm
    do im = 1 , Nb ( ia )
      ja = selectedneighbors ( ia , im )
      jt = itype ( ja )
      spec_vois1 ( ia , jt ) = spec_vois1 ( ia , jt ) + 1   
    enddo
  enddo

  WRITE ( kunit_VOIS1FF , '(a)' ) '#atom    nb_neighbors       ... of different types ... '
  do ia = 1 , natm
#ifdef GFORTRAN
   WRITE ( FMT ,  * ) ntype
   WRITE ( kunit_VOIS1FF , '(2i6,i6,'// ADJUSTL(FMT) //'(i6))' ) ia , Nb ( ia ), ( spec_vois1 ( ia , jt ) , jt = 1 ,ntype )
#else
   WRITE ( kunit_VOIS1FF , '(2i6,i6,<ntype>(i6))' ) ia , Nb ( ia ), ( spec_vois1 ( ia , jt ) , jt = 1 ,ntype )
#endif
  enddo
  deallocate ( spec_vois1 )


  deallocate ( selectedneighbors )
  deallocate ( Nb )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE fixed_distance

! *********************** SUBROUTINE SANN **************************************
!
!     based on sann.f :
!     Fortran implementation of the SANN algorithm                        
!     van Meel, Filion, Valeriani and Frenkel November (2011)             
!
! ******************************************************************************
SUBROUTINE sann

  USE config,           ONLY :  natm , simu_cell , rx , ry , rz , atype , atypei , ntype , itype , natmi
  USE control,          ONLY :  cutshortrange
  USE cell,             ONLY :  kardir , dirkar, dirkar_1
  USE io,               ONLY :  ionode , stderr , stdout , kunit_VOIS1FF
  
  implicit none
  
  ! local
  integer :: ia, ja , im , sm , it , jt  , jnb , jtnb , kkkk 
  integer :: m
  real(kind=dp) :: dist
  real(kind=dp) :: rxi ,ryi ,rzi
  real(kind=dp) :: rxij ,ryij ,rzij
  real(kind=dp) :: sxij ,syij ,szij
  real(kind=dp) :: rm , rm1  !     R(m) as in Eq.3 in the manuscript

  integer       , dimension (:)   , allocatable  :: Nb
  integer       , dimension (:)   , allocatable  :: natmi_cluster 
  character(len=2), dimension (:) , allocatable  :: atypei_cluster 
  integer       , dimension (:)   , allocatable  :: countneighbors 
  integer       , dimension (:,:) , allocatable  :: neighbor 
  integer       , dimension (:,:) , allocatable  :: sortneighbor 
  integer       , dimension (:,:) , allocatable  :: selectedneighbors 
  integer       , dimension (:)   , allocatable  :: labeltd 
  integer       , dimension (:)   , allocatable  :: labeld 
  integer       , dimension (:,:) , allocatable  :: spec_vois1
  real(kind=dp) , dimension (:,:) , allocatable  :: distance 
  real(kind=dp) , dimension (:)   , allocatable  :: tmpdist 
  real(kind=dp) , dimension (:)   , allocatable  :: td 
  real(kind=dp)                                  :: xxx , yyy , zzz 
  real(kind=dp) , dimension (:,:) , allocatable  :: distancesorted 
  character(len=20) :: FMT

  ! distributions
  integer :: knb
 
  logical :: notselected 

  ! ====================
  ! allocate main arrays
  ! ====================
  ! countneighbors = number of neighbours of particle ia
  allocate ( countneighbors ( natm ) ) 
  ! neighbor = list of neighbours of particles ia
  allocate ( neighbor ( natm , nmaxneigh ) )
  ! distance = list of distances between each neighbour of particle ia and particle ja 
  allocate (  distance ( natm , nmaxneigh ) )
  allocate (  distancesorted ( natm , nmaxneigh ) )
  ! working arrays
  allocate (  tmpdist ( nmaxneigh ) )
  allocate (  td ( ( nmaxneigh + 1 ) / 2 ) )
  allocate (  labeld ( nmaxneigh ) )
  allocate (  labeltd ( ( nmaxneigh + 1 ) / 2 ) )
  ! sortneighbor = sorted neighbours
  allocate ( sortneighbor ( natm , nmaxneigh ) )
  allocate ( selectedneighbors ( natm , nmaxneigh ) )
  ! Nb = final number of neighbours of particle ia
  allocate ( Nb ( natm ) )

  ! ===============
  !     init 
  ! ===============
  countneighbors    = 0 
  neighbor          = 0 
  distance          = 0.0_dp
  distancesorted    = 0.0_dp
  tmpdist           = 0.0_dp
  td                = 0.0_dp
  labeld            = 0
  labeltd           = 0
  sortneighbor      = 0
  selectedneighbors = 0
  Nb                = 0

  ! =======================================================================
  !     Step 1:
  !     first we identify the particles within a cutoff radius
  ! =======================================================================

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = 1 , natm
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    do ja = 1 , natm
      if ( ia .eq. ja ) cycle
      rxij = rxi - rx ( ja )
      ryij = ryi - ry ( ja )
      rzij = rzi - rz ( ja )
      sxij = rxij - nint ( rxij )
      syij = ryij - nint ( ryij )
      szij = rzij - nint ( rzij )
      rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
      ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
      rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
      dist = rxij * rxij + ryij * ryij + rzij * rzij
!      print*,ia,ja,sqrt(dist),cutshortrange*0.75d0
      if ( sqrt(dist) .lt. cutshortrange*0.75d0 ) then
        ! ja is a neighbour of ia
          countneighbors(ia) = countneighbors(ia) + 1
        ! build a list of neighbours
          if ( countneighbors(ia) .gt. nmaxneigh ) then
            io_node WRITE ( stderr , '(a,4i12)' ) 'ERROR: out of bound neighbor in sann',countneighbors(ia),nmaxneigh,ia,ja
          endif
          neighbor(ia,countneighbors(ia)) = ja
        ! create a list with the distance between ia and ja 
          distance(ia,countneighbors(ia)) = dist
       endif
    enddo
  enddo

  ! =============================================================== 
  ! Step 2:
  ! for every particle i sort all (countneighbors) 
  ! neighbours (neighbor) according to their 
  ! distances (distance) and create  a new list of 
  ! particle i's (sortneighbor)
  ! and a new sorted list of distances (distancesorted)
  ! =============================================================== 
  do ia = 1 , natm
    ! store distance to sort
    tmpdist = distance ( ia , : ) 
    ! store label of neighbor to sort  
    do im = 1 , countneighbors ( ia )
      labeld ( im ) = im
    enddo
    ! ===========================================
    !  arrays are sorted for increasing distance 
    !  (see tools.f90 for more details )
    !  the old labels are stored in labeld
    !  td , labeltd : working arrays used recursively in merge_sort
    ! ===========================================
    call merge_sort ( tmpdist , countneighbors ( ia ) , td , labeld , labeltd )
    distancesorted ( ia , : ) = tmpdist
    ! get sortneighbor form labeld  
    do im = 1 , countneighbors ( ia )
      sm = labeld ( im ) 
      sortneighbor ( ia , im )  = neighbor( ia , sm )
    enddo 
  enddo

#ifdef debug_vois1_sann
  WRITE ( stdout , '(a)' ) 'before merge_sort'
  !do ia = 1 , natm
    ia = 1
    do im = 1 , countneighbors ( ia )
      WRITE ( stdout ,*) distance ( ia , im ) , neighbor ( ia , im ) 
    enddo 
  !enddo

  WRITE ( stdout , '(a)' ) 'after merge_sort'
  !do ia = 1 , natm
    ia = 1
    do im = 1 , countneighbors ( ia )
      WRITE ( stdout , * ) distancesorted ( ia , im ) , sortneighbor ( ia , im ) 
    enddo 
  !enddo
#endif

  do ia = 1 , natm
    ! =================================== 
    ! Step 3: 
    ! start with 3 neighbours
    ! ===================================
    m = 3
    ! ===================================
    ! Step 4: 
    ! compute R(m) as in Eq.3 
    ! ===================================
    rm = 0
    do im=1,m
      rm = rm + distancesorted ( ia , im )
    enddo
    rm = rm/(m-2)
    ! ===================================
    ! compute r(m+1)
    ! ===================================
    do ja = 1 , countneighbors ( ia )      
      rm1 = 0
      do im = 1 , m
        rm1 = rm1 + distancesorted ( ia , im )
      enddo
      rm1 = rm1/(m-2)
      ! ===================================
      ! Step 5:  
      ! if rm > rm1     
      ! ===================================
      if ( rm .ge. rm1 ) then     
        rm = rm1
        ! ===================================
        ! increase m
        ! ===================================
        m = m + 1
      else
      ! ==================================================
      ! Step 6:
      ! if rm < rm1, m is the final number of neighbours
      ! ==================================================
        exit
      endif
    enddo
    ! =============================================
    ! the final number of neighbours is m = Nb(i) 
    ! and the neighbours are  selectedneighbors
    ! =============================================
    Nb(ia) = m
    do ja = 1 , Nb(ia)
      selectedneighbors ( ia , ja ) = sortneighbor ( ia , ja )
    enddo
  enddo

  ! output clusters 
  allocate ( natmi_cluster ( 0: ntype ) )
  allocate ( atypei_cluster ( 0: ntype ) )
  do ia = 1 , natm
    it = itype ( ia ) 
    atypei_cluster(0) = 'C'//atype(ia)
    atypei_cluster(1) = 'A'
    atypei_cluster(2) = 'B'
    natmi_cluster     = 0
    do ja = 1 , Nb ( ia )
      jnb = selectedneighbors ( ia , ja )
      jtnb= itype ( jnb ) 
      natmi_cluster(jtnb) = natmi_cluster(jtnb) + 1
    enddo
    natmi_cluster(0)    = 1
    kkkk = 10000+Nb ( ia ) + 1
    kk(Nb ( ia ) + 1) = kk(Nb ( ia ) + 1) + 1
    WRITE ( kkkk , * ) Nb ( ia ) + 1
    WRITE ( kkkk , '(a)' ) 'cluster'
    WRITE ( kkkk , '(3f20.12)' )  1000.0_dp ,    0.0_dp ,    0.0_dp
    WRITE ( kkkk , '(3f20.12)' )     0.0_dp , 1000.0_dp ,    0.0_dp
    WRITE ( kkkk , '(3f20.12)' )     0.0_dp ,    0.0_dp , 1000.0_dp
    WRITE ( kkkk , '(i4)'      )  ntype+1
#ifdef GFORTRAN
    WRITE ( FMT , * ) ntype+1
    WRITE ( kkkk , '('// ADJUSTL(FMT) //'a3 )') ( atypei_cluster(it) , it=0,ntype )
#else
    WRITE ( kkkk , '(<ntype+1>a3 )') ( atypei_cluster(it) , it=0,ntype )
#endif
    WRITE ( kkkk , *           ) ( natmi_cluster (it) , it=0,ntype )
    WRITE ( kkkk ,'(A)')         'Cartesian'
    WRITE ( kkkk ,'(a,3e20.12)') atypei_cluster(0) , 0.0_dp , 0.0_dp , 0.0_dp 
    do it = 1 , ntype
      do ja = 1 , Nb ( ia )
        jnb = selectedneighbors ( ia , ja )
        jtnb= itype ( jnb )
        xxx = rx ( jnb ) - rx ( ia ) - NINT ( rx ( jnb ) - rx ( ia ) )
        yyy = ry ( jnb ) - ry ( ia ) - NINT ( ry ( jnb ) - ry ( ia ) )
        zzz = rz ( jnb ) - rz ( ia ) - NINT ( rz ( jnb ) - rz ( ia ) )
        CALL dirkar_1 ( xxx , yyy , zzz , simu_cell%A , 1 )
        if ( jtnb .eq. it ) WRITE ( kkkk ,'(a,3e20.12)') atype ( jnb ) , xxx , yyy , zzz 
      enddo 
    enddo
    !WRITE ( 10000 , '(a,i6,a,i6,a,<Nb(ia)>(i6,a))' ) 'neighbors of atom : ', ia , ' #nb : ', Nb ( ia ), ' : ', ( selectedneighbors ( ia , ja ) , atype (selectedneighbors ( ia , ja )) , ja = 1 , Nb ( ia ) )
  enddo
  deallocate ( natmi_cluster )
  deallocate ( atypei_cluster )
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  ! ===================================================
  !  Square version : discard selectedneighbors atoms
  ! ===================================================
  if ( vois1algo .eq. 'sannsq' ) then
    ! ===============
    !     init 
    ! ===============
    countneighbors    = 0
    neighbor          = 0
    distance          = 0.0_dp
    distancesorted    = 0.0_dp
    tmpdist           = 0.0_dp
    td                = 0.0_dp
    labeld            = 0
    labeltd           = 0
    sortneighbor      = 0
    Nb                = 0

    do ia = 1 , natm
      rxi = rx ( ia )
      ryi = ry ( ia )
      rzi = rz ( ia )
      do ja = 1 , natm
        notselected = .FALSE.
        if ( ia .eq. ja ) cycle
        do im = 1 , Nb(ia)
          if ( ja .eq. selectedneighbors ( ia , im ) ) notselected = .TRUE.
        enddo
        if ( notselected ) cycle
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        dist = rxij * rxij + ryij * ryij + rzij * rzij
        if ( sqrt(dist) .lt. cutshortrange*0.75d0 ) then
          ! ja is a neighbour of ia
            countneighbors(ia) = countneighbors(ia) + 1
            if ( countneighbors(ia) .gt. nmaxneigh ) then
              io_node WRITE ( stderr , '(a,4i12)' ) 'ERROR: out of bound neighbor in sann',countneighbors(ia),nmaxneigh,ia,ja
            endif
            neighbor(ia,countneighbors(ia)) = ja
            distance(ia,countneighbors(ia)) = dist
         endif
      enddo
    enddo
    ! =============== 
    !     sort
    ! =============== 
    do ia = 1 , natm
      tmpdist = distance ( ia , : )
      do im = 1 , countneighbors ( ia )
        labeld ( im ) = im
      enddo
      call merge_sort ( tmpdist , countneighbors ( ia ) , td , labeld , labeltd )
      distancesorted ( ia , : ) = tmpdist
      do im = 1 , countneighbors ( ia )
        sm = labeld ( im )
        sortneighbor ( ia , im )  = neighbor( ia , sm )
      enddo
    enddo
    ! ================================
    !  Compute rm
    ! ================================
    do ia = 1 , natm
      m = 3
      rm = 0
      do im=1,m
        rm = rm + distancesorted ( ia , im )
      enddo
      rm = rm/(m-2)
      do ja = 1 , countneighbors ( ia )
        rm1 = 0
        do im = 1 , m
          rm1 = rm1 + distancesorted ( ia , im )
        enddo
        rm1 = rm1/(m-2)
        if ( rm .ge. rm1 ) then
          rm = rm1
          m = m + 1
        else
          exit
        endif
      enddo
      ! =============================================
      ! the final number of neighbours is m = Nb(i) 
      ! and the neighbours are  selectedneighbors
      ! =============================================
      Nb(ia) = m
      do ja = 1 , Nb(ia)
        selectedneighbors ( ia , ja ) = sortneighbor ( ia , ja )
      enddo
    enddo

  else if ( (vois1algo .eq. 'sannsq') .or. (vois1algo .eq. 'sann') ) then
    ! ==================================
    !        some statistics
    ! ==================================
    ! on Nb
    do ia = 1 , natm
      it = itype ( ia ) 
      knb = int( Nb(ia) ) 
      ! ====================== 
      !  test out of bound
      ! ====================== 
      if (knb.lt.0.or.knb.gt.nbmax) then
        io_node WRITE ( stderr , * ) 'ERROR: out of bound dist_nb'
        io_node WRITE ( stderr , '(4i6)') ia , knb , Nb(ia) , nbmax
        STOP
      endif
      dib_nb ( 0  , knb ) = dib_nb ( 0  , knb ) + 1
      dib_nb ( it , knb ) = dib_nb ( it , knb ) + 1
    enddo

    ! ==================================
    !       single speciation
    ! dimension of table speciation is :
    !         natm x ntype 
    ! ==================================
    allocate ( spec_vois1 ( natm , ntype ) )
    spec_vois1 = 0

    do ia = 1 , natm
      do im = 1 , Nb ( ia ) 
        ja = selectedneighbors ( ia , im )
        jt = itype ( ja )
        spec_vois1 ( ia , jt ) = spec_vois1 ( ia , jt ) + 1
      enddo
    enddo

    WRITE ( kunit_VOIS1FF , '(a)' ) '#atom    nb_neighbors       ... of different types ... '
    do ia = 1 , natm
#ifdef GFORTRAN
      WRITE ( FMT , *  ) ntype
      WRITE ( kunit_VOIS1FF , '(2i6,i6,'// ADJUSTL(FMT) // '(i6))' ) ia , Nb ( ia ), ( spec_vois1 ( ia , jt ) , jt = 1 ,ntype )
#else
      WRITE ( kunit_VOIS1FF , '(2i6,i6,<ntype>(i6))' ) ia , Nb ( ia ), ( spec_vois1 ( ia , jt ) , jt = 1 ,ntype )
#endif
    enddo
    deallocate ( spec_vois1 )

  endif

  ! =================
  !    deallocate
  ! =================
  deallocate ( countneighbors ) 
  deallocate ( neighbor )
  deallocate ( distance )
  deallocate ( distancesorted )
  deallocate ( tmpdist )
  deallocate ( td )
  deallocate ( labeld )
  deallocate ( labeltd )
  deallocate ( sortneighbor )
  deallocate ( selectedneighbors )
  deallocate ( Nb )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
      
  return

END SUBROUTINE sann 

! *********************** SUBROUTINE voronoi_construction **********************
!
! adapted from Allen and Tildesley
! we put the original description. But the actual implementation is slightly
! more general ( non cubic )
!    FICHE F.35  - PART B                                          
!    THE VORONOI CONSTRUCTION IN 3D.                               
!
!     CONSTRUCTION OF VORONOI POLYHEDRON IN 3D.                     
!                                                                   
!     THIS PROGRAM TAKES IN A CONFIGURATION IN A CUBIC BOX WITH     
!     CONVENTIONAL PERIODIC BOUNDARY CONDITIONS AND FOR EACH ATOM   
!     OBTAINS THE SURROUNDING VORONOI POLYHEDRON, DEFINED AS THAT   
!     REGION OF SPACE CLOSER TO THE CHOSEN ATOM THAN TO ANY OTHER.  
!     NEIGHBOURING POLYHEDRA DEFINE NEIGHBOURING ATOMS.             
!     THIS PROGRAM IS SLOW BUT ESSENTIALLY FOOLPROOF.               
!     WE USE THE MINIMUM IMAGE CONVENTION AND SET A CUTOFF BEYOND   
!     WHICH ATOMS ARE ASSUMED NOT TO BE NEIGHBOURS: BOTH OF THESE   
!     MEASURES ARE DANGEROUS FOR SMALL AND/OR RANDOM SYSTEMS.       
!     WE DELIBERATELY DO NOT USE PREVIOUSLY-FOUND NEIGHBOURS IN     
!     CONSTRUCTING NEIGHBOUR LISTS, SO THAT AN INDEPENDENT CHECK    
!     MAY BE MADE AT THE END.                                       
!     HERE WE SIMPLY PRINT OUT THE GEOMETRICAL INFORMATION AT THE   
!     END.  THE OUTPUT IS QUITE LENGTHY.  IN PRACTICE, IT WOULD     
!     PROBABLY BE ANALYZED DIRECTLY WITHOUT PRINTING IT OUT.        
!     NB: BEWARE DEGENERATE CONFIGURATIONS, I.E. ONES IN WHICH MORE 
!     THAN FOUR VORONOI DOMAINS SHARE A VERTEX. THE SIMPLE CUBIC    
!     AND FACE-CENTRED CUBIC LATTICES ARE EXAMPLES.               
!
! ******************************************************************************
SUBROUTINE voronoi_construction

  USE config,           ONLY :  natm , simu_cell , rx, ry , rz , itype , ntype 
  USE control,          ONLY :  cutshortrange
  USE cell,             ONLY :  kardir, dirkar
  USE io,          ONLY :  ionode, stdout , stderr , kunit_VOIS1FF

  implicit none

  integer       :: ia, ja , im , jm , it , jt 
  integer       :: can 
  integer       :: ncan , nver , nedge, nface , ncoord
  logical       :: ok
  real(kind=dp) :: dist , coord , rcut , rcutsq
  real(kind=dp) :: rxi ,ryi ,rzi
  real(kind=dp) :: rxij ,ryij ,rzij
  real(kind=dp) :: sxij ,syij ,szij

  integer       , dimension (:)     , allocatable :: tag , edges
  real(kind=dp) , dimension (:)     , allocatable :: px , py , pz , ps
  real(kind=dp) , dimension (:)     , allocatable :: rxver , ryver , rzver ! vertices position
  integer       , dimension (:)     , allocatable :: iver , jver , kver ! vertices indices
  integer       , dimension (:)     , allocatable :: nb 
  integer       , dimension (:,:)   , allocatable :: selectedneighbors 
  integer       , dimension (:,:)   , allocatable :: spec_vois1
  character(len=20) :: FMT

  integer :: knb
#ifdef further_info_voronoi
  integer :: ver
#endif

#ifdef sort_voronoi
  integer :: sm      
  integer       , dimension (:)   , allocatable  :: labeltd
  integer       , dimension (:)   , allocatable  :: labeld
  real(kind=dp) , dimension (:)   , allocatable  :: tmpdist
  real(kind=dp) , dimension (:)   , allocatable  :: td
  real(kind=dp) , dimension (:)   , allocatable  :: tmpx , tmpy , tmpz 
#endif

  
  ! ========================
  !        allocation
  ! ======================== 
  allocate ( px ( nmaxneigh ) , py ( nmaxneigh ) , pz ( nmaxneigh ) , ps ( nmaxneigh ) )
  allocate ( tag ( nmaxneigh ) , edges ( nmaxneigh ) )
  allocate ( rxver ( nmaxneigh ) , ryver ( nmaxneigh ) , rzver ( nmaxneigh ) ) 
  allocate ( iver ( nmaxneigh ) , jver ( nmaxneigh ) , kver ( nmaxneigh ) ) 
  allocate ( nb ( natm ) )
  allocate ( selectedneighbors ( natm , nmaxneigh ) )

  ! ========================
  !     some constants
  ! ========================
  rcut   = cutshortrange*0.75_dp
  rcutsq = rcut * rcut 

  ! ========================
  !   ZERO ACCUMULATORS 
  ! ========================
  nb = 0
  selectedneighbors = 0

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  ! ======================================
  !         MAIN LOOP STARTS  
  ! ======================================
  do ia = 1 , natm
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )

#ifdef debug_vois1_voronoi
    WRITE ( stdout ,'('' RESULTS FOR ATOM '',i6)') ia
#endif
    can = 0

    ! ======================================
    !         SELECT CANDIDATES 
    ! ======================================
    do ja = 1 , natm

      if ( ia .ne. ja ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        dist = rxij * rxij + ryij * ryij + rzij * rzij
        if ( dist .lt. rcutsq ) then
          can = can + 1
          if ( can .gt. nmaxneigh ) then
            io_node WRITE ( stderr , '(a)') 'ERROR out of bound nmaxneigh too small in voronoi_construction',can,nmaxneigh
            STOP
          endif
          px (can)  = rxij
          py (can)  = ryij
          pz (can)  = rzij
          ps (can)  = dist
          tag(can)  = ja 
        endif
      endif
    enddo
    ! ====================================
    !   CANDIDATES HAVE BEEN SELECTED 
    ! ====================================
    ncan = can

    ! =========================================
    !   SORT INTO ASCENDING ORDER OF DISTANCE 
    !   THIS SHOULD IMPROVE EFFICIENCY        
    ! =========================================

    !CALL SORT ( maxcan, PX, PY, PZ, PS, TAG, NCAN )
#ifdef sort_voronoi
    ! =============================================================== 
    ! Step 2:
    ! for every particle i sort all (countneighbors) 
    ! neighbours (neighbor) according to their 
    ! distances (distance) and create  a new list of 
    ! particle i's (sortneighbor)
    ! and a new sorted list of distances (distancesorted)
    ! =============================================================== 
    ! working arrays
    allocate ( tmpdist ( nmaxneigh ) )
    allocate ( tmpx ( nmaxneigh ) )
    allocate ( tmpy ( nmaxneigh ) )
    allocate ( tmpz ( nmaxneigh ) )
    allocate ( td ( ( nmaxneigh + 1 ) / 2 ) )
    allocate ( labeld ( nmaxneigh ) )
    allocate ( labeltd ( ( nmaxneigh + 1 ) / 2 ) )

    tmpdist=0.0_dp
    tmpx = 0.0_dp ; tmpy = 0.0_dp ; tmpz = 0.0_dp
    td = 0.0_dp
    labeld = 0
    labeltd = 0

    tmpdist =  ps 
    ! store label of neighbor to sort  
    do im = 1 , ncan 
      labeld ( im ) = im
    enddo
#ifdef debug_vois1_voronoi
  WRITE ( stdout , '(a)' ) 'before merge_sort'
  do im = 1 , ncan
    WRITE ( stdout , '(3f16.8)' ) ps ( im ) , px ( im ) , py ( im )
  enddo
#endif
    ! ===========================================
    !  arrays are sorted for increasing distance 
    !  (see tools.f90 for more details )
    !  the old labels are stored in labeld
    !  td , labeltd : working arrays used recursively in merge_sort
    ! ===========================================
    call merge_sort ( tmpdist , ncan , td , labeld , labeltd )
    ps = tmpdist
    ! get px py pz sorted form labeld  
    tmpx = px
    tmpy = py
    tmpz = pz
    do im = 1 , ncan 
      sm = labeld ( im )
      px ( im )  = tmpx ( sm )
      py ( im )  = tmpy ( sm )
      pz ( im )  = tmpz ( sm )
    enddo

#ifdef debug_vois1_voronoi
  WRITE ( stdout , '(a)' ) 'after merge_sort'
  do im = 1 , ncan
    WRITE ( stdout ,*) ps ( im ) , px ( im ) , py ( im )
  enddo
#endif
    deallocate ( td , labeld , labeltd , tmpdist , tmpx , tmpy , tmpz ) 
#endif 

    ! ====================================
    !    PERFORM VORONOI ANALYSIS
    ! ====================================
    CALL work_voronoi ( nmaxneigh , ncan , nver , nedge , nface , px , py , pz , ps , edges , rxver , ryver , rzver , iver , jver , kver )

    ! ====================================
    !       WRITE OUT RESULTS
    ! ====================================
#ifdef further_info_voronoi    
    WRITE(*,'(1x,''number of neighbours '',i6)') NFACE
    WRITE(*,'(1x,''neighbour list '')')
    WRITE(*,10001)
#endif
    do can = 1, ncan
      if ( edges ( can ) .ne. 0 ) then
        ps ( can ) = sqrt ( ps ( can ) )
#ifdef further_info_voronoi    
        WRITE(*,'(2i5,4f12.5)') tag ( can ), edges ( can ), px ( can ), py ( can ), pz ( can ), ps ( can )
#endif
        nb ( ia ) = nb ( ia ) + 1
        selectedneighbors ( ia, nb (ia) ) = tag ( can )
      endif
    enddo

#ifdef further_info_voronoi    
    WRITE(*,'(1x,''number of edges '',i6)')    nedge 
    WRITE(*,'(1x,''number of vertices '',i6)') nver
    WRITE(*,'(1x,''vertex list '')')
    WRITE(*,10002)

    do ver = 1, nver
      WRITE(*,'(1x,3i5,3x,3f12.5)') tag ( iver ( ver ) ), tag ( jver ( ver ) ), tag ( kver ( ver ) ) , &
                                    rxver ( ver ) , ryver ( ver ) , rzver ( ver ) 
    enddo
#endif

  enddo
  ! ================================
  !       MAIN LOOP ENDS   
  ! ================================
#ifdef further_info_voronoi    
  WRITE(*,'(''FINAL SUMMARY'')')
  WRITE(*,10003)
#endif

  ncoord = 0
  do ia = 1 , natm
    ncoord = ncoord + nb (ia)
    ! ========================================
    !  CHECK THAT IF I IS A NEIGHBOUR OF J 
    !  THEN J IS ALSO A NEIGHBOUR OF I     
    ! ========================================
    do im = 1, Nb ( ia )
      ja = selectedneighbors ( ia , im )
      ok = .FALSE.
      jm = 1
      do while ( jm .LE. Nb ( ja ) ) 
        ok = ( ia .EQ. selectedneighbors ( ja , im ) )
        jm = jm + 1
      enddo

      if ( .NOT. OK ) then
        WRITE(*,'(1x,i6,'' is not a neighbour of '',i6)') ia, ja
      endif 

    enddo
  enddo

  COORD = REAL ( NCOORD ) / REAL ( natm )

  WRITE(*,'(1X,'' average coordination number = '',f10.5)') COORD

  ! =================================
  !         some statistic
  ! =================================
  ! on Nb 
  do ia = 1 , natm
    it = itype  ( ia )
    knb = int( nb (ia) )
    ! ====================== 
    !  test out of bound
    ! ====================== 
    if (knb.lt.0.or.knb.gt.nbmax) then
      io_node WRITE ( stderr , * ) 'ERROR: out of bound dist_nb'
      io_node WRITE ( stderr , '(4i6)') ia , knb , nb(ia) , nbmax
      STOP
    endif
    dib_nb ( 0  , knb ) = dib_nb ( 0  , knb ) + 1
    dib_nb ( it , knb ) = dib_nb ( it , knb ) + 1
  enddo
  ! ==================================
  !       single speciation
  ! dimension of table speciation is :
  !         natm x ntype 
  ! ==================================
  allocate ( spec_vois1 ( natm , ntype ) )
  spec_vois1 = 0

  do ia = 1 , natm
    do im = 1 , nb ( ia )
      ja = selectedneighbors ( ia , im ) 
      jt = itype ( ja )
      spec_vois1 ( ia , jt ) = spec_vois1 ( ia , jt ) + 1
    enddo
  enddo

  OPEN ( UNIT = kunit_VOIS1FF , FILE = 'VOIS1FF')
  WRITE ( kunit_VOIS1FF , '(a)' ) '#atom    nb_neighbors       ... of different types ... '
  do ia = 1 , natm
#ifdef GFORTRAN
    WRITE ( FMT , * ) ntype
    WRITE ( kunit_VOIS1FF , '(2i6,i6,'// ADJUSTL(FMT) //'(i6))' ) ia , nb ( ia ), ( spec_vois1 ( ia , jt ) , jt = 1 ,ntype )
#else
    WRITE ( kunit_VOIS1FF , '(2i6,i6,<ntype>(i6))' ) ia , nb ( ia ), ( spec_vois1 ( ia , jt ) , jt = 1 ,ntype )
#endif
  enddo
  CLOSE( kunit_VOIS1FF )
  deallocate ( spec_vois1 )

  ! ========================
  !        deallocation
  ! ======================== 
  deallocate ( px , py , pz , ps )
  deallocate ( tag , edges )
  deallocate ( rxver , ryver , rzver )
  deallocate ( iver , jver , kver )
  deallocate ( nb )
  deallocate ( selectedneighbors )

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )


10001   FORMAT(1x,'atom ',3x,'face ',1x,'index ',3x,'edges ',3x,'            relative posisiton         ',3x,'  distance')
10002   FORMAT(1x,'      indices           relative position')
10003   FORMAT(1x,'index    nabs    ... neighbour indices ... ')

END SUBROUTINE voronoi_construction

! *********************** SUBROUTINE work_voronoi ******************************
!
!     ROUTINE TO PERFORM VORONOI ANALYSIS                           
!     FMV note : en toute rigueur c'est une triangulation de Delaunay (vraiment ?? )
!                                                                   
!     WE WORK INITIALLY ON DOUBLE THE CORRECT SCALE,                
!     I.E. THE FACES OF THE POLYHEDRON GO THROUGH THE POINTS.       
!
! ******************************************************************************
SUBROUTINE work_voronoi ( maxcan, nn , nv , ne , nf , px , py , pz , ps , edges , vx , vy , vz , iv , jv , kv )

  USE io,          ONLY :  stdout

  implicit none

  integer :: maxcan, nn, nv, ne, nf
  integer :: edges(maxcan)
  integer ::  iv(maxcan), jv(maxcan), kv(maxcan)
  real(kind=dp) :: px (maxcan), py (maxcan), pz (maxcan), ps (maxcan)
  real(kind=dp) :: vx (maxcan), vy (maxcan), vz (maxcan)

  logical :: ok
  integer ::     i, j, k, l, nn1, nn2, n, v
  real(kind=dp) :: AI, BI, CI, DI, AJ, BJ, CJ, DJ, AK, BK, CK, DK
  real(kind=dp) :: AB, BC, CA, DA, DB, DC, det , detinv 
  real(kind=dp) :: vxijk , vyijk , vzijk  
  real(kind=dp) , PARAMETER :: tol = 1.E-6_dp

  ! =============================================
  !  IF THERE ARE LESS THAN 4 POINTS GIVEN 
  !  WE CANNOT CONSTRUCT A POLYHEDRON      
  ! =============================================
  if ( nn .LT. 4 ) then
    WRITE(*,'('' LESS THAN 4 POINTS GIVEN TO WORK '',I5)') nn 
    STOP
  endif
  nn1 = nn - 1
  nn2 = nn - 2
  v = 0
  ! =============================================
  !  WE AIM TO EXAMINE EACH POSSIBLE VERTEX  
  !  DEFINED BY THE INTERSECTION OF 3 PLANES 
  !  EACH PLANE IS SPECIFIED BY px,py,pz,ps  
  ! =============================================
  do i = 1, nn2
    AI =  px (i)
    BI =  py (i)
    CI =  pz (i)
    DI = -ps (i)
    do j = i + 1, nn1
      AJ =  px (j)
      BJ =  py (j)
      CJ =  pz (j)
      DJ = -ps (j)
      AB = AI * BJ - AJ * BI
      BC = BI * CJ - BJ * CI
      CA = CI * AJ - CJ * AI
      DA = DI * AJ - DJ * AI
      DB = DI * BJ - DJ * BI
      DC = DI * CJ - DJ * CI

      do k = j + 1 , nn
        AK =  px (k)
        BK =  py (k)
        CK =  pz (k)
        DK = -ps (k)

        det = AK * BC + BK * CA + CK * AB

        if ( abs ( det ) .GT. tol ) then
        ! ===========================
        !   THE PLANES INTERSECT 
        ! ===========================
        detinv = 1.0_dp / DET
        vxijk  = ( - DK * BC + BK * DC - CK * DB ) * detinv
        vyijk  = ( - AK * DC - DK * CA + CK * DA ) * detinv
        vzijk  = (   AK * DB - BK * DA - DK * AB ) * detinv

        ! =================================
        !  NOW WE TAKE SHOTS AT THE VERTEX
        !  USING THE REMAINING PLANES ....
        ! =================================

        ok = .TRUE.
        l  = 1
        do while ( ok .AND. ( l .LE. NN ) ) 
          if ( ( l .NE. i ) .AND. ( l .NE. j ) .AND. ( l .NE. k )       ) then
            ok = ( ( px (l) * vxijk + py (l) * vyijk + pz (l) * vzijk  ) .le. ps (l) )
          endif
          l = l + 1
        enddo
        ! ============================
        !  IF THE VERTEX MADE IT      
        !  ADD IT TO THE HALL OF FAME 
        !  CONVERT TO CORRECT SCALE   
        ! ============================
        if ( ok ) then
          v = v + 1
          if ( v .gt. maxcan ) STOP 'TOO MANY VERTICES'
            iv (v)  = i
            jv (v)  = j
            kv (v)  = k
            vx (v) = 0.5_dp * vxijk
            vy (v) = 0.5_dp * vyijk
            vz (v) = 0.5_dp * vzijk
          endif
        endif 
      enddo
    enddo
  enddo

  nv = v 

  if ( nv .lt. 4 ) then
    WRITE ( stdout ,'('' ERROR : less than 4 vertices in work_voronoi '',i5)') nv
    STOP 
  endif

  ! ===================================
  !   IDENTIFY NEIGHBOURING POINTS 
  ! ===================================
  do n = 1, nn
    edges (n) = 0
  enddo

  do v = 1, nv
    edges ( iv (v) ) = edges ( iv (v) ) + 1
    edges ( jv (v) ) = edges ( jv (v) ) + 1
    edges ( kv (v) ) = edges ( kv (v) ) + 1
  enddo
  ! ============================================
  !  POINTS WITH NONZERO edges ARE NEIGHBOURS 
  ! ============================================

  ! =========================
  !  CHECK EULER RELATION 
  !  
  ! =========================
  nf = 0
  ne = 0

  do n = 1, nn
    if ( edges (n) .gt. 0 ) nf = nf + 1
    ne = ne + edges (n)
  enddo

  if ( MOD ( ne, 2 ) .ne. 0 ) then
    WRITE(*,'('' noninteger number of edges'',i5)') ne
    STOP
  endif

  ne = ne / 2

  if ( ( nv - ne + nf ) .ne. 2 ) then
    WRITE(*,'('' ERROR : euker relation degeneracy in work_voronoi '')')
    STOP
  endif

  return

END SUBROUTINE work_voronoi

SUBROUTINE vois1_alloc
 
  USE config,                   ONLY :  ntype 
  USE control,                  ONLY :  calc 

  implicit none

  if ( calc .ne. 'vois1' ) return

  allocate( dib_nb(0:ntype,0:nbmax) )
  dib_nb  = 0

  return

END SUBROUTINE vois1_alloc

SUBROUTINE vois1_dealloc

  USE control,                  ONLY :  calc

  implicit none

  if ( calc .ne. 'vois1' ) return

  deallocate( dib_nb )

  return

END SUBROUTINE vois1_dealloc

END MODULE voisin
