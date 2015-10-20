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
! ======= Hardware =======

MODULE polarizability

  USE constants,                ONLY :  dp
  USE mpimdff
  implicit none

  real(kind=dp) :: epsw
  real(kind=dp) :: applied_field(3)

  ! f = e + sum T * dmu 
  ! pola = f^-1 * dmu

  real(kind=dp), dimension ( : , : , : ), allocatable :: f
  real(kind=dp), dimension ( : , : , : ), allocatable :: invf 
  real(kind=dp)                                       :: e (3,3)
  real(kind=dp), dimension ( : , : , : ), allocatable :: dmu
  real(kind=dp), dimension ( : , : , : ), allocatable :: pola 

  ! =====================================================
  ! algorithm 
  ! =====================================================
  character(len=60), SAVE :: algo_pola      
  character(len=60), SAVE :: algo_pola_allowed(2)
  data algo_pola_allowed / 'from_dmu' , 'from_wfc' /



CONTAINS

! *********************** SUBROUTINE pola_init ***********************************
!> \brief
!! initialisation of main parameters
! ******************************************************************************
SUBROUTINE pola_init

  USE control,  ONLY :  lstatic , calc
  USE io,       ONLY :  ionode ,stdin, stdout

  implicit none

  integer            :: ioerr
  character(len=132) :: filename

  namelist /polatag/    applied_field,&
                        algo_pola,&
                        epsw  

  if ( calc .ne. 'polarizability' ) return
  ! ======================
  !  polatag default values
  ! ======================
  CALL pola_default_tag
  ! =====================
  !  read polatag namelist
  ! =====================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , polatag ,iostat =ioerr)
  if ( ioerr .lt. 0 )  then
    io_node &
    WRITE ( stdout, '(a)') 'ERROR reading input_file : polatag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node &
    WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : polatag wrong tag',ioerr
    STOP
  endif
  CLOSE  ( stdin )
  ! ======================
  !  check polatag namelist
  ! ======================
  CALL pola_check_tag

  ! ===================
  !  print polatag info
  ! ===================
  CALL pola_print_info(stdout)

  return

END SUBROUTINE pola_init

SUBROUTINE pola_print_info(kunit)

  implicit none

  ! global
  integer :: kunit

  return 

END SUBROUTINE pola_print_info

SUBROUTINE pola_check_tag

  implicit none

  e = 0.0_dp
  e(1,1) = applied_field(1)
  e(2,2) = applied_field(2)
  e(3,3) = applied_field(3)


  return 

END SUBROUTINE pola_check_tag

SUBROUTINE pola_default_tag

  implicit none

  return 

END SUBROUTINE pola_default_tag

SUBROUTINE pola_main

  USE control,		ONLY :	lvnlist, cutlongrange
  USE constants,       	ONLY : tpi, piroot
  USE config,		ONLY :	natm, verlet_coul, atom_dec
  USE cell
  USE field,		ONLY :	ewald_param, multipole_ES_dir , multipole_ES_rec , alphaES
  USE kspace 
  USE io,		ONLY :	stdout, kunit_DMU
  USE time,             ONLY :  fcoultimetot1 , fcoultimetot2

  implicit none

  ! local 
  integer :: ia 
  integer :: i, j, k  
  integer, parameter :: LWORK=1000
  real(kind=dp) :: WORK ( LWORK )
  integer :: ipiv ( 3 )
  integer :: ierr
  real(kind=dp), dimension ( : , : , : ), allocatable :: f_rec
  real(kind=dp), dimension ( : , : , : ), allocatable :: f_dir
  real(kind=dp), dimension ( : , : , : ), allocatable :: f_self
  real(kind=dp), dimension ( : , : )    , allocatable :: mu_dumb
  real(kind=dp), dimension ( : , : , : ), allocatable :: theta_dumb
  logical           :: do_efield , do_efg, do_forces, do_stress , damp_ind , task(6)
  real(kind=dp) :: ttt1, ttt2
  real(kind=dp) :: tau_dumb( 3 , 3 )
  real(kind=dp), dimension ( : )        , allocatable :: fx_dumb, fy_dumb, fz_dumb 
  real(kind=dp), dimension ( : , : )    , allocatable :: ef_dumb
  real(kind=dp), dimension ( : , : , :) , allocatable :: efg_dumb 
  real(kind=dp) :: u_dumb
  real(kind=dp) :: selfa,selfa2 , alpha2
  character(len=2) :: xxxx


  CALL read_pos

  CALL ewald_param
  CALL do_split ( natm , myrank , numprocs , atom_dec , 'atoms' )
  !CALL efg_alloc
  !CALL efg_mesh_alloc
  !if ( lcoulomb ) CALL do_split ( km_coul%nk , myrank , numprocs , km_coul%kpt_dec , 'k-pts' )


  verlet_coul%listname='coul'
  verlet_coul%cut=cutlongrange
  if ( lvnlist )    CALL vnlist_pbc !( verlet_coul )

  allocate ( dmu ( natm , 3 , 3 ) ) 
  allocate ( pola ( natm , 3 , 3 ) ) 
  allocate ( f_rec ( natm ,3 ,3 ) , f_dir ( natm ,3 ,3 ) , f_self ( natm ,3 ,3 ) , f ( natm ,3, 3 ) )
  allocate ( invf ( natm ,3 ,3 ) )
  dmu    = 0.0_dp
  pola   = 0.0_dp
  f_rec  = 0.0_dp
  f_dir  = 0.0_dp
  f_self = 0.0_dp
  f      = 0.0_dp
  invf   = 0.0_dp
 
  alpha2 = alphaES * alphaES
  selfa  = alphaES / piroot
  selfa2 = 2.0_dp * selfa * alpha2 / 3.0_dp



  if ( algo_pola .eq. 'from_dmu' ) then

  do_efield = .true.
  do_efg    = .false.
  do_forces = .false.
  do_stress = .false.
  damp_ind  = .false.
  task      = .false.
  task(3)   = .true. ! dipole-dipole
  
  allocate( mu_dumb ( natm , 3 ) )   
  allocate( theta_dumb ( natm , 3 , 3 ) )   
  allocate( ef_dumb  (natm,3) )
  allocate( efg_dumb(natm,3,3) )
  allocate( fx_dumb  (natm)    , fy_dumb  (natm)   , fz_dumb  (natm) )
  mu_dumb = 0.0_dp
  theta_dumb = 0.0_dp
  ef_dumb = 0.0_dp
  efg_dumb = 0.0_dp
  fx_dumb = 0.0_dp
  fy_dumb = 0.0_dp
  fz_dumb = 0.0_dp


  ! E field in x
  OPEN( kunit_DMU, FILE = 'DMUX' ) 
  do ia = 1 , natm 
    READ ( kunit_DMU, * ) xxxx , mu_dumb ( ia , 1 ) , mu_dumb ( ia , 2 ) , mu_dumb ( ia , 3 )
    print*,ia
  enddo 
  CLOSE( kunit_DMU )
  dmu(:,:,1) = mu_dumb


  ttt1 = MPI_WTIME(ierr)
  CALL multipole_ES_dir ( u_dumb , ef_dumb, efg_dumb, fx_dumb , fy_dumb , fz_dumb , tau_dumb , mu_dumb , theta_dumb , task , damp_ind , &
                          do_efield , do_efg , do_forces , do_stress )
  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 )
  f_dir (:,:,1) = ef_dumb

  ttt1 = MPI_WTIME(ierr)
  CALL multipole_ES_rec ( u_dumb , ef_dumb , efg_dumb , fx_dumb , fy_dumb , fz_dumb , tau_dumb , mu_dumb , theta_dumb , task , &
                          do_efield , do_efg , do_forces , do_stress )
  ttt2 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt2 - ttt1 )
  f_rec (:,:,1) = ef_dumb


  f_self( : , : , 1) = 2.0_dp * selfa2 * mu_dumb


  ! E field in y 
  OPEN( kunit_DMU, FILE = 'DMUY' ) 
  do ia = 1 , natm 
    READ ( kunit_DMU, * ) xxxx , mu_dumb ( ia , 1 ) , mu_dumb ( ia , 2 ) , mu_dumb ( ia , 3 )
    print*,ia
  enddo 
  CLOSE( kunit_DMU )
  dmu(:,:,2) = mu_dumb


  ttt1 = MPI_WTIME(ierr)
  CALL multipole_ES_dir ( u_dumb , ef_dumb, efg_dumb, fx_dumb , fy_dumb , fz_dumb , tau_dumb , mu_dumb , theta_dumb , task , damp_ind , &
                          do_efield , do_efg , do_forces , do_stress )
  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 )
  f_dir (:,:,2) = ef_dumb

  ttt1 = MPI_WTIME(ierr)
  CALL multipole_ES_rec ( u_dumb , ef_dumb , efg_dumb , fx_dumb , fy_dumb , fz_dumb , tau_dumb , mu_dumb , theta_dumb , task , &
                          do_efield , do_efg , do_forces , do_stress )
  ttt2 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt2 - ttt1 )
  f_rec (:,:,2) = ef_dumb

  f_self( : , : , 2) = 2.0_dp * selfa2 * mu_dumb

  ! E field in z
  OPEN( kunit_DMU, FILE = 'DMUZ' ) 
  do ia = 1 , natm 
    READ ( kunit_DMU, * ) xxxx , mu_dumb ( ia , 1 ) , mu_dumb ( ia , 2 ) , mu_dumb ( ia , 3 )
    print*,ia
  enddo 
  CLOSE( kunit_DMU )
  dmu(:,:,3) = mu_dumb


  ttt1 = MPI_WTIME(ierr)
  CALL multipole_ES_dir ( u_dumb , ef_dumb, efg_dumb, fx_dumb , fy_dumb , fz_dumb , tau_dumb , mu_dumb , theta_dumb , task , damp_ind , &
                          do_efield , do_efg , do_forces , do_stress )
  ttt2 = MPI_WTIME(ierr)
  fcoultimetot1 = fcoultimetot1 + ( ttt2 - ttt1 )
  f_dir (:,:,3) = ef_dumb

  ttt1 = MPI_WTIME(ierr)
  CALL multipole_ES_rec ( u_dumb , ef_dumb , efg_dumb , fx_dumb , fy_dumb , fz_dumb , tau_dumb , mu_dumb , theta_dumb , task , &
                          do_efield , do_efg , do_forces , do_stress )
  ttt2 = MPI_WTIME(ierr)
  fcoultimetot2 = fcoultimetot2 + ( ttt2 - ttt1 )
  f_rec (:,:,3) = ef_dumb

  f_self( : , : , 3) = 2.0_dp * selfa2 * mu_dumb

  endif

  deallocate( mu_dumb )   
  deallocate( theta_dumb )   
  deallocate( ef_dumb )
  deallocate( efg_dumb )
  deallocate( fx_dumb , fy_dumb , fz_dumb )


  do ia=1,natm
  f(ia,:,:) = f_dir(ia,:,:) + f_rec(ia,:,:) + f_self(ia,:,:) + e(:,:)
  enddo

  do ia = 1 , 1 
    WRITE( stdout,'(a)') '------------------------------'
    WRITE( stdout,'(a,i)') 'ion : ', ia
    WRITE( stdout,'(a)') '------------------------------'
    WRITE( stdout,'(a)') ''
    WRITE( stdout,'(a)') 'f =  '
    WRITE( stdout,'(a,3e16.8)') 'x',( f(ia,1,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'y',( f(ia,2,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'z',( f(ia,3,j) , j=1,3)
    WRITE( stdout,'(a)') ''
    WRITE( stdout,'(a)') 'f (rec)=  '
    WRITE( stdout,'(a,3e16.8)') 'x',( f_rec(ia,1,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'y',( f_rec(ia,2,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'z',( f_rec(ia,3,j) , j=1,3)
    WRITE( stdout,'(a)') ''
    WRITE( stdout,'(a)') 'f (dir)=  '
    WRITE( stdout,'(a,3e16.8)') 'x',( f_dir(ia,1,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'y',( f_dir(ia,2,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'z',( f_dir(ia,3,j) , j=1,3)
    WRITE( stdout,'(a)') ''
    WRITE( stdout,'(a)') 'f (self)=  '
    WRITE( stdout,'(a,3e16.8)') 'x',( f_self(ia,1,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'y',( f_self(ia,2,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'z',( f_self(ia,3,j) , j=1,3)


    invf ( ia , : , : )  = f ( ia , : , : )
    CALL DGETRF( 3, 3, invf(ia,:,:), 3, ipiv, ierr )
    if ( ierr.lt.0 ) then
      WRITE( 6 , '(a,i6)' ) 'ERROR call to DGETRF failed in induced_moment',ierr
      STOP
    endif
    CALL DGETRI( 3 , invf(ia,:,:) , 3 ,  ipiv , WORK, LWORK, ierr )
    if ( ierr.lt.0 ) then
      WRITE( 6, '(a,i6)' ) 'ERROR call to DGETRI failed in induced_moment',ierr
      STOP
    endif
    WRITE( stdout,'(a)') 'f^-1 =  '
    WRITE( stdout,'(a,3e16.8)') 'x',( invf(ia,1,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'y',( invf(ia,2,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'z',( invf(ia,3,j) , j=1,3)
    WRITE( stdout,'(a)') ''
    WRITE( stdout,'(a)') 'dmu =  '
    WRITE( stdout,'(a,3e16.8)') 'x',( dmu(ia,1,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'y',( dmu(ia,2,j) , j=1,3)
    WRITE( stdout,'(a,3e16.8)') 'z',( dmu(ia,3,j) , j=1,3)
    WRITE( stdout,'(a)') ''

  pola(ia,:,:)=0.0_dp
  call dgemm("N","N",3,3,3,1.0_dp,invf(ia,:,:),3,dmu(ia,:,:),3,0.0_dp,pola(ia,:,:),3)

  WRITE( stdout,'(a)') 'pola =  '
  WRITE( stdout,'(a,3e16.8)') 'x',( pola(ia,1,j) , j=1,3) 
  WRITE( stdout,'(a,3e16.8)') 'y',( pola(ia,2,j) , j=1,3) 
  WRITE( stdout,'(a,3e16.8)') 'z',( pola(ia,3,j) , j=1,3) 
  WRITE( stdout,'(a)') ''
  WRITE( stdout,'(a,e16.8)') 'alpha = ',pola(ia,1,1)+pola(ia,2,2)+pola(ia,3,3) / 3.0_dp

  enddo

  deallocate ( dmu ) 
  deallocate ( pola ) 
  deallocate ( f_rec , f_dir , f_self , f )
  deallocate ( invf )

  return

END SUBROUTINE pola_main

END MODULE polarizability
