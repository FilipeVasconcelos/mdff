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
! ======= Hardware =======

! *********************** MODULE radial_distrib ********************************
!> \brief
!! Module related to radial function distribution calculation and/or static
!! factor structure
! ******************************************************************************
MODULE radial_distrib 

  USE io,               ONLY :  ionode, stdout, stdin, stderr
  USE constants,        ONLY :  dp
  USE mpimdff

  implicit none

  integer :: nbins            !< (internal) number of bins in g(r) distribution
  integer :: npairs            !< (internal) number of bins in g(r) distribution
  integer :: nconf            !< number of configurations used in g(r) calculation
  real(kind=dp) :: cutgr      !< radial cut-off 
  real(kind=dp) :: resg       !< resolution in g(r) distribution 
  !> g(r) function ( bin x ntype x ntype )
  integer, dimension(:,:,:), allocatable :: gr  

CONTAINS

! *********************** SUBROUTINE grcalc_init *******************************
!
!> \brief
!! initialize radial distribution calculation parameters
!
! ******************************************************************************
SUBROUTINE gr_init

  USE config,                   ONLY :  simu_cell
  USE control,                  ONLY :  calc

  implicit none

  ! local
  integer            :: ioerr 
!28/05/13  integer            :: npangr, i
  character(len=132) :: filename

  namelist /grtag/   nconf , &
                     cutgr , &
                     resg  

  if ( calc .ne. 'gr' ) return

  CALL gr_default_tag
  
  ! =================
  !  read grtag tags
  ! =================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , grtag , iostat=ioerr )
  if ( ioerr .lt. 0 )  then
   io_node WRITE ( stdout, '(a)') 'ERROR reading input_file : grtag section is absent'
   STOP
  elseif ( ioerr .gt. 0 )  then
   io_node WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : grtag wrong tag'
   STOP
  endif
  CLOSE ( stdin )

  ! ==========================================
  ! define a new resolution to be 2^N points
  ! ==========================================
 nbins=int(cutgr/resg)+1
 print*,cutgr,resg,nbins
!28/05/13 ! i = 1
!28/05/13 ! do while ( 2**i .lt. nbins )
!28/05/13 !    i = i + 1
!28/05/13 ! enddo
!28/05/13 ! npangr = i
!28/05/13 ! nbins = 2** npangr
!28/05/13 ! resg = cutgr / DBLE ( nbins  )

  CALL gr_print_info(stdout)

  return 
 
END SUBROUTINE gr_init

! *********************** SUBROUTINE gr_alloc **********************************
!
!> \brief
!! Allocate g(r) function 
!
! ******************************************************************************
SUBROUTINE gr_alloc

  USE control,                  ONLY :  calc
  USE config,                   ONLY :  ntype

  implicit none

  if ( calc .ne. 'gr' .and. calc .ne. 'rmc' ) return
  allocate(  gr ( 0 : nbins - 1 , 0 : ntype , 0 : ntype ) ) 
  gr = 0      

  return 
 
END SUBROUTINE gr_alloc


! *********************** SUBROUTINE gr_dealloc ********************************
!
!> \brief
!! Deallocate g(r) function 
!
! ******************************************************************************
SUBROUTINE gr_dealloc

  USE control,                  ONLY :  calc

  implicit none

  if ( calc .ne. 'gr' ) return
  
  deallocate( gr )

  return 
 
END SUBROUTINE gr_dealloc


! *********************** SUBROUTINE gr_default_tag ****************************
!
!> \brief
!! set default values to gr tag
!
! ******************************************************************************
SUBROUTINE gr_default_tag

  USE config,           ONLY : simu_cell

  implicit none

  ! ===============
  !  default value
  ! ===============
  resg = 0.1_dp
  nconf = 0
  cutgr=0.5_dp * MIN(simu_cell%WA,simu_cell%WB,simu_cell%WC)

  return

END SUBROUTINE gr_default_tag


! *********************** SUBROUTINE gr_print_info *****************************
!
!> \brief
!! print infog on g(r) calculation
!
! ******************************************************************************
SUBROUTINE gr_print_info(kunit)

  USe control,                  ONLY :  calc

  implicit none
 
  ! local
  integer :: kunit

   if ( ionode ) then
                  WRITE ( kunit ,'(a,f10.5,a)')         'resolution of g(r) function resg     = ',resg,' new value to have 2^N points in g(r)'
                  WRITE ( kunit ,'(a,i5)')              'number of points in g(r)             = ',nbins
                  WRITE ( kunit ,'(a)')                 'save radial_distribution in file     :   GRTFF' 
      if ( calc .eq. 'gr' )     then 
                  WRITE ( kunit ,'(a)')                 'read configuration from file         :   TRAJFF'
                  blankline(kunit)
                  WRITE ( kunit ,'(a,i5)')              'number of config. in TRAJFF          = ',nconf        
                  blankline(kunit)
      endif
   endif 
  return

END SUBROUTINE gr_print_info

! *********************** SUBROUTINE grcalc ************************************
!
!> \brief
!! main driver of radial distribution function calculation
!! this subroutine read the trajectory, Allocate, call the  
!
! ******************************************************************************
SUBROUTINE grcalc

  USE control,                  ONLY :  itraj_format , trajff_data
  USE config,                   ONLY :  system , natm , ntype , rx , ry , rz , atype , &
                                        rho , config_alloc , simu_cell , atypei , itype, natmi, &
                                        coord_format_allowed , atom_dec , read_traj , read_traj_header
  USE io,                       ONLY :  kunit_TRAJFF , kunit_GRTFF , kunit_NRTFF
  USE constants,                ONLY :  pi 
  USE cell,                     ONLY :  lattice , dirkar , periodicbc, kardir
  USE time,                     ONLY :  grtimetot_comm

  implicit none

  ! local 
  integer                                              :: ic , ngr , igr 
  integer                                              :: it1 , it2 , mp , ierr 
  real(kind=dp),     dimension ( : , : ) , allocatable :: grr 
  integer,           dimension ( : )     , allocatable :: nr 
  character(len=15), dimension ( : )     , allocatable :: cint
  real(kind=dp)                                        :: rr , vol
  real(kind=dp)                                        :: ttt1 , ttt2      
  real(kind=dp)                                        :: average_volume
  real(kind=dp)                                        :: rho_av 


  rho_av = 0.0_dp
  average_volume = 0.0_dp      

  OPEN ( kunit_GRTFF , file = 'GRTFF' )
  OPEN ( kunit_NRTFF , file = 'NRTFF' )

  if ( itraj_format .ne. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' )
  if ( itraj_format .eq. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' , form = 'unformatted')
  CALL read_traj_header ( kunit_TRAJFF , itraj_format )
  if ( itraj_format .ne. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' )
  if ( itraj_format .eq. 0 ) OPEN ( UNIT = kunit_TRAJFF  , FILE = 'TRAJFF' , form = 'unformatted')

  CALL lattice ( simu_cell ) 
  rho = DBLE ( natm )  / simu_cell%omega 

#ifdef debug
  write(*,*) simu_cell%A
  write(*,*) simu_cell%omega
  write(*,*) simu_cell%ANORM
#endif

  CALL gr_init

  CALL print_general_info( stdout )

#ifdef debug
  write(*,*)
  write(*,*) simu_cell%A
  write(*,*) simu_cell%omega
  write(*,*) simu_cell%ANORM
#endif
  ! ===================================
  !  here we know natm, then alloc 
  !  and decomposition can be applied 
  ! ================================== 
  CALL config_alloc 
  CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms' )
  CALL gr_alloc

  npairs =  ntype * ( ntype + 1 ) / 2
  allocate ( grr ( 0 : npairs , 0 : nbins-1 ) , nr ( 0 : npairs ) , cint ( 0 : npairs  ))
  grr  = 0.0_dp
  nr   = 0
  cint = ''

#ifdef debug
  if ( ionode ) then 
    WRITE ( stdout , '(a,2i6)' ) 'debug : atom decomposition istart, iend ', atom_dec%istart , atom_dec%iend
    WRITE ( stdout , '(a,2i6)' ) 'debug : number of type npairs ', npairs
  endif
#endif

  CALL typeinfo_init

  ngr = 0
  do ic = 1, nconf

    CALL read_traj ( kunit_TRAJFF , itraj_format , trajff_data ) 

    CALL lattice ( simu_cell )
    rho_av = rho_av + ( REAL ( natm ,kind=dp )  / simu_cell%omega )
    average_volume = average_volume + simu_cell%omega    
    io_node WRITE ( stdout , '(a,i6,a,i6,a,f12.3)' ) 'config : [ ',ic,' / ',nconf,' ]   vol : ',simu_cell%omega
#ifdef debug
    CALL distance_tab
!    print*,simu_cell%omega,average_volume/ REAL(ic,kind=dp)
#endif
    CALL kardir     ( natm , rx , ry , rz , simu_cell%B ) 
    CALL periodicbc ( natm , rx , ry , rz  )
    CALL dirkar     ( natm , rx , ry , rz , simu_cell%A ) 

#ifdef debug
    CALL distance_tab
#endif

    ngr=ngr+1 
    ! ==========================
    !  calc radial_distribution 
    ! ==========================  
    call gr_main 

  enddo !nconf 
  rho_av = rho_av / REAL(nconf,kind=dp)      
  average_volume = average_volume / REAL(nconf,kind=dp)
  if ( ionode .and. average_volume .ne. simu_cell%omega ) write(stdout,'(a,e16.8)') 'average volume : ',average_volume

  ! ===========================================
  !        merge results  
  ! ===========================================
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
  grtimetot_comm = grtimetot_comm + ( ttt2 - ttt1 )
#endif


#ifdef debug
  do igr=0, nbins
    io_node WRITE (stdout , '(a,5i6)') 'debug ( total ) : ',igr,gr(igr,1,1)
  enddo
#endif

  ! ======================================================= 
  !  write output files GRTFF , NRTFF
  ! ======================================================= 
  cint ( 0) = atypei ( 0 )//' - '//atypei ( 0 )
  mp = 1
  do it1 = 1 , ntype
    do it2 = it1 , ntype
      cint(mp) = atypei(it1)//' - '//atypei(it2)
      mp = mp + 1 
    enddo
  enddo

#ifdef GFORTRAN
  io_node WRITE ( kunit_GRTFF , '(8a)' )          '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
  io_node WRITE ( kunit_NRTFF , '(8a)' )          '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
#else
  io_node WRITE ( kunit_GRTFF , '(<npairs+2>a)' ) '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
  io_node WRITE ( kunit_NRTFF , '(<npairs+2>a)' ) '#       rr           ', ( cint ( mp ) , mp = 0 , npairs ) 
#endif

  nr ( 0 ) = 0 

  do igr = 0 , nbins-1
    rr  = ( REAL ( igr ,kind=dp )+0.5d0)*resg
    vol = 4.d0 * pi * ( resg * rr* rr + ( resg**3 ) / 12.d0 )
    vol = vol / average_volume  ! only  make sense for trajectory in NPT otherwise the average volume is the current volume ) 
    ! all - all pairs 
    grr ( 0 , igr ) = DBLE ( gr ( igr , 0 , 0 ) ) / ( ngr * vol * natm * natm )
!    grr ( 0 , igr ) = DBLE ( gr ( 0 , 0 , igr ) ) / ( ngr * vol * natm * natm )
    ! type pairs
    mp = 1
    do it1 = 1 , ntype
      do it2 = it1 , ntype
#ifdef debug
       io_node WRITE ( stdout , '(a,3i5)' ) 'debug ( pair ) : ', mp , it1 , it2
#endif        
        if ( mp .lt. 0 .and. mp .gt. npairs ) then
          WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
          STOP
        endif 
        nr ( mp ) = it1          
        grr ( mp , igr ) = DBLE ( gr ( igr , it1 , it2 ) ) / DBLE ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) )  
!        grr ( mp , igr ) = DBLE ( gr ( it1 , it2 , igr ) ) / DBLE ( ngr * vol * natmi ( it1 ) * natmi ( it2 ) )  
        mp = mp + 1
      enddo
    enddo
      if ( ionode ) then
#ifdef GFORTRAN
        WRITE ( kunit_GRTFF ,'(8e20.10)') rr , ( grr ( mp , igr ) , mp = 0 , npairs ) 
        WRITE ( kunit_NRTFF ,'(8e20.10)') rr , ( DBLE ( grr ( mp , igr ) ) * 4.0_dp * pi * rr * rr * DBLE ( natmi(nr(mp)) * vol ) , mp = 0 , npairs )
#else
        WRITE ( kunit_GRTFF ,'(<npairs+2>e20.10)') rr , ( grr ( mp , igr ) , mp = 0 , npairs ) 
        WRITE ( kunit_NRTFF ,'(<npairs+2>e20.10)') rr , ( DBLE ( grr ( mp , igr ) ) * 4.0_dp * pi * rr * rr * DBLE ( natmi(nr(mp)) * vol ) , mp = 0 , npairs )
#endif
      endif
  enddo


  CLOSE ( kunit_NRTFF )
  CLOSE ( kunit_GRTFF )

  CLOSE( kunit_TRAJFF )

  !CALL static_struc_fac ( grr , nbins , npairs ) 

  deallocate ( grr ,nr , cint)
  CALL gr_dealloc

  return

END SUBROUTINE grcalc

! *********************** SUBROUTINE gr_main ***********************************
!
!> \brief
!! based on Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE gr_main 

  USE control,                  ONLY :  myrank
  USE config,                   ONLY :  natm , natmi , rx , ry , rz , atype , simu_cell , ntype , itype, atom_dec
  USE time,                     ONLY :  grtimetot
  USE cell,                     ONLY :  kardir , dirkar

  implicit none

  ! local
  integer :: ia , ja , ierr , ita , jta 
  integer :: igr 
  real(kind=dp) :: cut2 , rijsq , rr 
  real(kind=dp) :: rxi , ryi , rzi
  real(kind=dp) :: rxij , ryij , rzij
  real(kind=dp) :: sxij , syij , szij
  real(kind=dp) :: ttt1 , ttt2      

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif
  ! =======================
  !  cut-off half box
  ! =======================
  cut2 = cutgr * cutgr 
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  do ia = atom_dec%istart , atom_dec%iend

    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    ita = itype( ia )
    do ja = 1, natm
      if ( ja .ne. ia ) then  
        jta = itype( ja )
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
        if ( rijsq.lt.cut2 ) then
          rr = SQRT ( rijsq )
          igr = INT ( rr / resg ) 
          if ( igr .lt. 0 .and. igr .gt. nbins-1 ) then
            WRITE ( stderr , '(a)' ) 'ERROR out of bound of gr in gr_main'
            STOP
          endif 
          ! all pairs
          gr ( igr , 0    , 0   ) = gr ( igr , 0   ,   0 ) + 1 
          gr ( igr , ita  , jta ) = gr ( igr , ita , jta  ) + 1
!          gr ( 0    , 0   , igr ) = gr ( 0   ,   0 , igr ) + 1 
!          gr ( ita  , jta , igr ) = gr ( ita , jta , igr ) + 1
        endif
      endif
    enddo
  enddo
  
#ifdef debug2
  do igr=0, nbins-1
    WRITE (stdout , '(a,5i12)') 'debug: ',myrank,igr,gr(igr,0,0)
  enddo
#endif 

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  grtimetot = grtimetot + ( ttt2 - ttt1 )
#endif

  return
 
END SUBROUTINE gr_main

! *********************** SUBROUTINE static_struc_fac **************************
!
!
! ******************************************************************************
!SUBROUTINE static_struc_fac ( gr , nbins , npairs )

!  USE io,                       ONLY :  ionode , kunit_STRFACFF , stdout 
!  USE config,                   ONLY :  rho
!  USE constants,                ONLY :  pi , tpi , imag

!  implicit none
!  ! global
!  integer :: nbins , npairs
!  real(kind=dp) :: gr ( 0 : nbins-1 , 0 : npairs ) , rr
!
!  ! local
!  integer :: i , j , is , NN
!  real(kind=dp) :: q , qj , ri , rip
!  complex(kind=dp) :: Sk
!  real(kind=dp) , dimension (:) , allocatable  :: stat_str
!  real(kind=dp) , dimension (:) , allocatable :: in 
!  real(kind=dp) , dimension (:,:) , allocatable :: Uji 
!  complex(kind=dp)   , dimension (:) , allocatable :: out
!  real(kind=dp) :: res , shift
!  real(kind=dp) :: x , k
!
!  io_node WRITE ( stdout , '(a)' ) 'in static_struc_fac'
!
!  allocate ( in ( nbins ) , out ( nbins /2 + 1 ) )
!
!  in  = ( 0.0,0.0)
!  out = ( 0.0,0.0)
!  ! ========
!  !   FFT
!  ! ========
!  do i=0,nbins-1
!    in ( i + 1 ) = gr ( i , 0 ) 
!  enddo
!
!!  CALL fft_1D_complex ( in , out , nbins )
!  CALL fft_1D_real(in,out,nbins)
!
!  do i= 1 , nbins/2+1
!    q = ( dble ( i )  + 0.5_dp ) / DBLE ( nbins ) / resg
!    Sk = 1.0_dp + rho * out( i + 1 )  
!    io_node WRITE ( 20000 , '(3e16.8)' )  q , Sk  
!  enddo
!
!  deallocate ( in , out )
!
!! other version
!  ! Uji (eq 12) J Phys Cndens Matter V17 2005 )
!!  allocate ( Uji ( nbins , nbins ) ) 
!
!  do i = 1 , nbins
!!    ri  = ( dble ( i )  + 0.5_dp )  * res
!    rip = ( dble ( i + 1 )  + 0.5_dp )  * res
!!    do j = 1 , nbins
!      qj = tpi * DBLE ( j ) + 0.5_dp / DBLE ( nbins / 2 + 1 ) / resg
!      Uji ( j , i ) = SR ( ri , qj ) - SR ( rip , qj ) 
! !     Uji ( j , i ) = Uji ( j , i ) / qj  
!    enddo
! ! enddo
!  Uji = 2.0_dp * tpi * Uji
!!
!  do i= 1 , nbins/2+1
!    q = tpi * ( dble ( i )  + 0.5_dp ) / DBLE ( nbins / 2 + 1) / resg
!!    do j = 1 , nbins
!      Sk =  Sk + Uji ( j , i ) * ( gr ( j , 0 ) -1.0_dp )  
!    enddo
!    io_node WRITE ( 30000 , '(3e16.8)' )  q , Sk
!  enddo
!
!
!  deallocate ( Uji )

! test purpose because I'm dumb
! I was enable to get k in "real life" 
! I did a simple discret case  ( N = 4 ) 
! to find the relation between k and q 

!  NN = 4
!  allocate ( in ( NN ) , out  ( NN ) )
!
!  in = ( 0.0,0.0)
!  is = 3 
!  in(is) = ( 1.0,0.0)
!  print*,'in in '
!  CALL fft_1D_complex ( in , out , NN )
!
!  res = 2.0_dp
!  shift= (is-1) * res
!  write( stdout , '(8a)' ) '            x       Re in(i)        Im in(i)            k          Re out(i)       Im out(i)        Re phase        Im phase'
!  do i=1,NN
!    x = DBLE(i-1)*res
!    k = DBLE(i-1) / DBLE ( NN ) 
!    q = DBLE(i-1) / DBLE ( NN ) / res
!    write( stdout , '(8e16.8)' ) x , in(i) , k , out(i) , exp ( -imag * tpi * k * (is-1) )
!    write( stdout , '(8e16.8)' ) x , in(i) , q , out(i) , exp ( -imag * tpi * q * shift )
!  enddo
!
!  deallocate ( in ,out )

!  return
!CONTAINS

!real(kind=dp) function SR(r,qj)
!  implicit none
!  real(kind=dp) :: r , qj  
!  SR = SIN ( qj * r ) / qj / qj 
!  SR = SR - r * COS ( qj * r ) / qj
!end function 

!END SUBROUTINE static_struc_fac


SUBROUTINE read_GRTFF( grr , filename)

  USE constants,         only : dp
  implicit none

  integer :: mp , bin
  real(kind=dp) , dimension ( 0:npairs , 0: nbins-1 ) :: grr
  real(kind=dp) :: rr
  character(*) :: filename

  OPEN(UNIT=1000,FILE=filename)
    read(1000,*)      
    do bin=0,nbins-1
      read(1000,*) rr, ( grr ( mp , bin ), mp = 0 , npairs )
    enddo 
  CLOSE(1000)

END SUBROUTINE

END MODULE radial_distrib 
! ===== fmV =====
