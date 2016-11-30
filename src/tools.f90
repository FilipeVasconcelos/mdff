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
!#define debug_vnl   
!#define debug_vnl2  
! ======= Hardware =======

! *********************** SUBROUTINE estimate_alpha ****************************
!
!> \brief
!! This subroutine estimate the alpha parameter in the Ewald summation.
!! This estimation is based on the cell size
!
!> \param[in]  epsw relative error target
!> \param[in]  rcut radial cut-off in the real space summation 
!> \param[out] alpha alpha parameter in the Ewald summation  
!
!> \note
! The parameter should be anyway optimised for each problem
!
! ******************************************************************************
SUBROUTINE estimate_alpha(alpha,epsw,rcut)
  
  USE constants,                ONLY :  dp

  implicit none

  ! global
  real(kind=dp) , intent (in)  :: epsw , rcut
  real(kind=dp) , intent (out) :: alpha
  ! local 
  real(kind=dp) :: alpha2 , rcut2 , rcut3 , e

  alpha = 4.0_dp
  alpha2 = alpha*alpha
  rcut2=rcut*rcut
  rcut3=rcut2*rcut

  e = EXP ( -alpha2*rcut2 )
  e = e / alpha2 / rcut3
  e = e * 0.56_dp  
  do while ( e .le. epsw )
    alpha = alpha - 0.01_dp
    alpha2 = alpha*alpha
    e = EXP ( -alpha2*rcut2 )
    e = e / alpha2 / rcut3
    e = e * 0.56_dp
  enddo

  return

END SUBROUTINE estimate_alpha

! *********************** SUBROUTINE accur_ES_frenkel_smit *********************
!
!> \brief
!! This subroutine estimate the alpha parameter and the number of cell 
!! in reciprocal space to perform the Ewald summation.
!! This estimation is adapted from frenkel and smit consideration 
!
!> \param[in]  epsw relative error target
!> \param[in]  rc radial cut-off in the real space summation 
!> \param[out] alpha alpha parameter in the Ewald summation  
!> \param[out] nc number of cells in reciprocal space 3D  
!
!> \note
! The parameter should be anyway optimised for each problem
!
! ******************************************************************************
SUBROUTINE accur_ES_frenkel_smit ( epsw , alpha , rc , nc ) 

  USE constants,                ONLY :  dp , pi
  USE config,                   ONLY :  simu_cell 
  USE io,                       ONLY :  stderr

  implicit none

  ! global
  real(kind=dp) , intent (in)  :: epsw
  real(kind=dp) , intent (in)  :: rc
  integer       , intent (out) :: nc ( 3 )
  real(kind=dp) , intent (out) :: alpha
  
  ! local
  real(kind=dp) :: s , ss , ra
  integer :: i , k

  k = 0
  s = 10.0_dp
  do 
    s = s - 0.01_dp
    ss = s * s   
    ra = exp ( - ss )  / ss
    if ( abs ( ra - epsw ) .gt. epsw ) exit  
    if ( s .le. 0.0_dp .or. k .gt.1e6 ) then
      WRITE( stderr , * ) 'ERROR in accur_ES_frenkel_smit',s,k 
      stop
    endif
    k = k + 1
  enddo 

  alpha = s / rc
  do i = 1 , 3 
    nc ( i ) = int ( s * simu_cell%ANORM(i) * alpha / pi ) 
  enddo

  return

END SUBROUTINE accur_ES_frenkel_smit

! *********************** SUBROUTINE distance_tab ******************************
!
!> \brief
!! this subroutine calculates distance between atoms and 
!! check if the smallest distance is not too small ( i.e < sigmaAA*0.001_dp) 
!
!> \author
!!  FMV
!
!> \history
!! 01/03/13 : do not check distance to wannier centers
!
! ******************************************************************************
SUBROUTINE distance_tab 

  USE control,                  ONLY :  calc
  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , ry , rz , simu_cell , itype, atype
  USE field,                    ONLY :  sigmalj , lwfc
  USE io,                       ONLY :  ionode , stdout, stderr
  USE cell,                     ONLY :  kardir , dirkar 

  implicit none
  
  ! local
  integer                             :: ia , ja , PANdis, kdis , it , jt , bin
  real(kind=dp)                       :: rxi, ryi, rzi
  real(kind=dp)                       :: sxij, syij, szij
  real(kind=dp)                       :: rxij, ryij, rzij, rij, rijsq, norm
  real(kind=dp)                       :: resdis,mindis 
  integer, dimension (:) ,allocatable :: dist
  integer :: indxmin(2)

  resdis = 0.01_dp ! should be keep hardware no need to be controled

  PANdis = MAX(simu_cell%WA,simu_cell%WB,simu_cell%WC) / resdis

  allocate ( dist ( 0:PANdis ) )
  dist = 0
  mindis = 100000.0_dp

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  kdis =0
  do ia = 1 , natm 
    it = itype ( ia ) 
    !if ( lwfc ( it ) .eq. -1 ) cycle
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    do ja = 1 , natm
      if ( ia .ne. ja ) then
        jt = itype ( ja ) 
    !    if ( lwfc ( jt ) .eq. -1 ) cycle
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        sxij = rxij - nint ( rxij )
        syij = ryij - nint ( ryij )
        szij = rzij - nint ( rzij )
        rxij = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
        ryij = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
        rzij = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
        rijsq = rxij *rxij + ryij * ryij + rzij * rzij
        rij = SQRT ( rijsq ) 
        if ( rij .lt. mindis  ) then
          mindis = rij
          indxmin(1) = ia
          indxmin(2) = ja
        endif
        if ( rij .lt. sigmalj(1,1) * 0.001_dp ) then
          if ( ionode ) &
          WRITE ( stdout ,'(a,i5,a,i5,a,f12.6)') &
          'ERROR: DISTANCE between atoms', ia ,' and ', ja ,' is very small',rij
          STOP 
        endif
        kdis = INT ( rij / resdis )
        if ( kdis .lt. 0 .or. kdis .gt. PANdis ) then
          io_node WRITE ( stdout , '(a,2i12,f48.8,2i12,2a)' ) 'ERROR: out of bound dist in distance_tab',kdis,PANdis,rij,ia,ja,atype(ia),atype(ja)
        endif
        dist ( kdis ) = dist ( kdis ) + 1
      endif
    enddo 
  enddo

  norm = natm * ( natm - 1 ) / 2
 
  if ( ionode ) then
    WRITE ( stderr ,'(a)')       'distance check subroutine'
    WRITE ( stderr ,'(a,f13.5,2i6,2a6)') 'smallest distance = ',mindis,indxmin(1),indxmin(2),atype(indxmin(1)),atype(indxmin(2))
    if  ( any(lwfc .eq. -1) )  then
    WRITE ( stderr ,'(a)')       'WARNING : we do not check distance to wannier centers'
    endif
    blankline(stderr)
  endif

  if ( calc .eq. 'dist' ) then
    do bin=0,PANdis
      write(10000,'(f16.8,i8)') REAL(bin,kind=dp)*resdis,dist(bin)
    enddo
  endif

  
  deallocate(dist)

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE distance_tab


! *********************** SUBROUTINE vnlist_pbc ********************************
!
!> \brief
!! verlet list subroutine :
!! Periodic boundaries condition version ( minimum image convention ) 
!! Parallelized and tested
! 
!> \decription
!! 
!! output:
!!  
!!          point    : array of size natm+1 
!!                     gives the starting and finishing index of array list for a given atom i
!!                     jbegin = point(i)  jend = point(i+1) - 1
!!          list     : index list of neighboring atoms
!!
!! how to use it :
!!                   do ia = 1, natm
!!                     jbegin = point(i)
!!                     jend = point(i+1) - 1
!!                     do jvnl = jbegin , jend
!!                       ja = list ( jvnl ) 
!!                       then ia en ja are neighboors   
!!                     enddo
!!                   enddo
!
! ******************************************************************************
SUBROUTINE vnlist_pbc  

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , ry , rz , itype, ntype , simu_cell , vnlmax , atom_dec, verlet_vdw , verlet_coul
  USE control,                  ONLY :  skindiff 
  USE cell,                     ONLY :  kardir , dirkar
  USE io,                       ONLY :  ionode , stdout , stderr

  implicit none

  ! local
  integer :: ia , ja
  integer :: icount_vdw , icount_coul
  integer :: k_vdw , k_coul
  integer :: p1 , p2
  real(kind=dp)  :: rskinsq_vdw , rskinsq_coul 
  real(kind=dp) :: rxi , ryi , rzi , rxij , ryij , rzij , rijsq , sxij , syij , szij 
  !debug
  !integer :: iv

  verlet_vdw%list=0
  verlet_coul%list=0
  verlet_vdw%point=0
  verlet_coul%point=0

  rskinsq_vdw = verlet_vdw%cut + skindiff
  rskinsq_vdw = rskinsq_vdw * rskinsq_vdw 
  rskinsq_coul = verlet_coul%cut + skindiff
  rskinsq_coul = rskinsq_coul * rskinsq_coul 

  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  icount_vdw = 1
  icount_coul = 1

  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    k_vdw = 0
    k_coul = 0
    do ja = 1 , natm
      if ( ( ia .gt. ja .and. ( MOD ( ia + ja , 2 ) .eq. 0 ) ) .or. &
           ( ia .lt. ja .and. ( MOD ( ia + ja , 2 ) .ne. 0 ) ) ) then

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
        !vdw sphere
        if ( rijsq .le. rskinsq_vdw ) then
          icount_vdw = icount_vdw + 1
          k_vdw = k_vdw + 1
          if ( icount_vdw .lt. 1 .or. icount_vdw - 1 .gt. vnlmax*natm ) then
            io_node WRITE ( stderr , '(a,2i12,f48.8)' ) 'ERROR: out of bound list in vnlist_pbc',icount_vdw-1,vnlmax*natm
            STOP
          endif
          verlet_vdw%list(icount_vdw-1) = ja
        endif
        ! coul sphere
        if ( rijsq .le. rskinsq_coul ) then
          icount_coul = icount_coul + 1
          k_coul = k_coul + 1
          if ( icount_coul .lt. 1 .or. icount_coul - 1 .gt. vnlmax*natm ) then
            io_node WRITE ( stderr , '(a,2i12,f48.8)' ) 'ERROR: out of bound list in vnlist_pbc',icount_coul-1,vnlmax*natm
            STOP
          endif
          verlet_coul%list(icount_coul-1) = ja
        endif
      endif
    enddo
    verlet_vdw%point(ia)  = icount_vdw-k_vdw
    verlet_coul%point(ia) = icount_coul-k_coul
  enddo
  verlet_vdw%point(atom_dec%iend + 1 ) = icount_vdw
  verlet_coul%point(atom_dec%iend + 1 )= icount_coul

  ! ======================================
  !         direct to cartesian
  ! ======================================
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  return

END SUBROUTINE vnlist_pbc

! *********************** SUBROUTINE vnlist_nopbc ******************************
!
!> \brief
!! no periodic boundaries condition version of verlet list
!! same as vnlist_pbc but no periodic boundaries
!
! ******************************************************************************
SUBROUTINE vnlist_nopbc ( vlist )

  USE constants, ONLY : dp
  USE config,    ONLY :  natm ,  rx , ry , rz , itype , verlet_list , ntype , atom_dec
  USE control,   ONLY :  skindiff

  implicit none

  TYPE(verlet_list) :: vlist
  ! local
  integer :: icount , ia , ja , it , jt , k
  integer :: p1,p2
  real(kind=dp) :: rskinsq ( ntype , ntype ) , rcut ( ntype , ntype ) , rskin ( ntype , ntype )
  real(kind=dp) :: rxi,ryi,rzi,rxij,ryij,rzij,rijsq

  do jt = 1, ntype
    do it = 1, ntype
       rcut    ( it , jt ) = vlist%cut
       rskin   ( it , jt ) = rcut  ( it , jt ) + skindiff
       rskinsq ( it , jt ) = rskin ( it , jt ) * rskin ( it , jt )
    enddo
  enddo

  icount = 1
  do ia = atom_dec%istart , atom_dec%iend
    rxi = rx ( ia )
    ryi = ry ( ia )
    rzi = rz ( ia )
    k = 0
    do ja = 1 , natm
      if ( ( ia .gt. ja .and. ( MOD ( ia + ja , 2 ) .eq. 0 ) ) .or. &
           ( ia .lt. ja .and. ( MOD ( ia + ja , 2 ) .ne. 0 ) ) ) then
        rxij = rxi - rx ( ja )
        ryij = ryi - ry ( ja )
        rzij = rzi - rz ( ja )
        rijsq = rxij * rxij + ryij * ryij + rzij * rzij
        p1 = itype ( ia )
        p2 = itype ( ja )
        if (rijsq .le. rskinsq(p1,p2)) then
          icount = icount + 1
          k = k + 1
          vlist%list ( icount - 1 ) = ja
        endif
      endif
    enddo
    vlist%point (ia ) = icount-k
  enddo
  vlist%point( atom_dec%iend + 1 ) = icount


  !print*,'verlet list info'
  !print*,vlist%list
  !print*,vlist%point
  !print*,vlist%listname
  !if ( vlist%listname .eq. 'coul' ) stop
 
 
  return

END SUBROUTINE vnlist_nopbc

! *********************** SUBROUTINE vnlistcheck *******************************
!
!> \brief
!! check wether verlet list should be updated
!
! ******************************************************************************
SUBROUTINE vnlistcheck  

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , ry , rz , xs , ys , zs , verlet_vdw , verlet_coul , atom_dec , simu_cell
  USE control,                  ONLY :  skindiff
  USE md,                       ONLY :  updatevnl , itime
  USE time,                     ONLY :  vnlisttimetot
  USE io,                       ONLY :  stdout , ionode
  USE cell,                     ONLY :  kardir , dirkar
  USE mpimdff

  implicit none

  ! local
  integer :: ia , ierr
  real(kind=dp) :: drneimax , drneimax2 , drnei
  real(kind=dp) :: rxsi,rysi,rzsi
  real(kind=dp) :: sxij,syij,szij
  real(kind=dp) :: ttt1 , ttt2

!#ifdef debug_vnl
!  write(stdout,'(a,a)') 'debug : in vnlistcheck',vlist%listname
!#endif

#ifdef MPI
  ttt1 = MPI_WTIME(ierr)
#endif
  ! ======================================
  !         cartesian to direct 
  ! ======================================
  CALL kardir ( natm , rx , ry , rz , simu_cell%B )

  drneimax = 0.0_dp
  drneimax2 = 0.0_dp
  do ia = 1, natm
    rxsi = rx ( ia ) - xs ( ia ) 
    rysi = ry ( ia ) - ys ( ia ) 
    rzsi = rz ( ia ) - zs ( ia ) 
    sxij = rxsi - nint ( rxsi )
    syij = rysi - nint ( rysi )
    szij = rzsi - nint ( rzsi )
    rxsi = sxij * simu_cell%A(1,1) + syij * simu_cell%A(1,2) + szij * simu_cell%A(1,3)
    rysi = sxij * simu_cell%A(2,1) + syij * simu_cell%A(2,2) + szij * simu_cell%A(2,3)
    rzsi = sxij * simu_cell%A(3,1) + syij * simu_cell%A(3,2) + szij * simu_cell%A(3,3)
    drnei = SQRT ( rxsi * rxsi + rysi * rysi + rzsi * rzsi ) 
#ifdef debug_vnl
    write(stdout,*) rxsi,rx ( ia ),xs ( ia )
#endif
    if ( drnei .gt. drneimax ) then
      drneimax2 = drneimax
      drneimax  = drnei
    else
      if ( drnei .gt. drneimax2 ) then
        drneimax2 = drnei
      endif
    endif        
  enddo 

  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
        
#ifdef debug_vnl
  write(stdout,'(a,2e16.8)') 'debug : verlet list displacement =',drneimax + drneimax2,skindiff
#endif
  ! ========================================
  ! DISPL = 2.0 * SQRT  ( 3.0 * DISPL ** 2 ) 
  !  I don't know from where this come from, 
  !  It was in the very first version
  ! ========================================

  if ( drneimax + drneimax2 .gt. skindiff ) then
    updatevnl = updatevnl + 1
#ifdef debug_vnl
  if ( ionode ) write ( stdout , '(a,i6,5f12.8)' ) 'verlet list update frequency =',updatevnl,REAL(itime)/REAL(updatevnl),skindiff, drneimax,drneimax2
#endif
    CALL vnlist_pbc 
    ! save coordinates it in direct coordinates
    CALL kardir ( natm , rx , ry , rz , simu_cell%B )
    do ia = 1, natm 
      xs ( ia ) = rx ( ia )
      ys ( ia ) = ry ( ia )
      zs ( ia ) = rz ( ia )
    END do
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )
     
  endif

#ifdef MPI
  ttt2 = MPI_WTIME(ierr)
  vnlisttimetot = vnlisttimetot + ( ttt2 - ttt1 ) 
#endif

  return

END SUBROUTINE vnlistcheck

! *********************** SUBROUTINE print_tensor_n ******************************
!
!> \brief
!! subroutine which print an (3,3) array in a tensor format 
!! the trace is also given in output 
!
!> \param[in] tens tensor being printed
!> \param[in] key tensr label
!
! ******************************************************************************
SUBROUTINE print_tensor_n ( tens ) 

  USE constants, ONLY : dp
  USE config,   ONLY :  natm
  USE io,  ONLY :  ionode , stdout

  implicit none

  ! global
  real(kind=dp) :: tens(3,3)

  ! local
  integer :: i

  if ( ionode ) then
    do i = 1 , 3
      WRITE ( stdout ,'(3e16.8)') tens(i,1) , tens(i,2) , tens(i,3)
    enddo
    WRITE ( stdout ,'(a,e16.8,a,e16.8,a)') '  iso = ',(tens(1,1) + tens(2,2) + tens(3,3))/3.0_dp , '(',(tens(1,1) + tens(2,2) + tens(3,3)) / 3.0_dp / dble(natm),')'
    blankline(stdout)
  endif

  return

END SUBROUTINE print_tensor_n 

! *********************** SUBROUTINE print_tensor ******************************
!
!> \brief
!! subroutine which print an (3,3) array in a tensor format 
!! the trace is also given in output 
!
!> \param[in] tens tensor being printed
!> \param[in] key tensr label
!
! ******************************************************************************
SUBROUTINE print_tensor( tens , key )

  USE constants,        ONLY :  dp
  USE config,           ONLY :  natm
  USE io,               ONLY :  ionode , stdout

  implicit none

  ! global
  real(kind=dp) :: tens(3,3)
  character(len=8) :: key

  ! local
  integer :: i

  if ( ionode ) then
    blankline(stdout)
    WRITE ( stdout ,'(a)') key
    do i = 1 , 3
      WRITE ( stdout ,'(3e16.8)') tens(i,1) , tens(i,2) , tens(i,3)
    enddo
    WRITE ( stdout ,'(a,e16.8,a,e16.8,a)') '  iso = ',(tens(1,1) + tens(2,2) + tens(3,3))/3.0_dp , '(',(tens(1,1) + tens(2,2) + tens(3,3)) / 3.0_dp / dble(natm),')'
    blankline(stdout)
  endif

  return

END SUBROUTINE print_tensor

! *********************** SUBROUTINE print_tensor_6x6 **************************
!
! subroutine which print an (n,n) array in a tensor format 
! the trace is also given in output 
!
! ******************************************************************************
SUBROUTINE print_tensor_nxn ( tens , key , n )

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm
  USE io,                       ONLY :  ionode , stdout

  implicit none

  ! global
  integer          :: n 
  real(kind=dp)    :: tens(n,n) 
  character(len=8) :: key
  character(len=20) :: FMT

  ! local
  integer          :: i , j
  real(kind=dp)    :: trace
  

  if ( ionode ) then
    blankline(stdout)
    WRITE ( stdout ,'(a)') key
    do i = 1 , n
#ifdef GFORTRAN
      WRITE ( FMT    ,* ) n
      WRITE ( stdout ,'('// ADJUSTL(FMT) //'e16.8)') ( tens(i,j) , j=1,n)
#else
      WRITE ( stdout ,'(<n>e16.8)') ( tens(i,j) , j=1,n)
#endif
      trace = trace + tens ( i , i ) 
    enddo
    WRITE ( stdout ,'(a,e16.8,a,e16.8,a)') 'iso = ',( trace )/3.0_dp , '(', ( trace ) / 3.0_dp / dble(natm),')'
    blankline(stdout)
  endif

  return

END SUBROUTINE print_tensor_nxn

! *********************** SUBROUTINE merge_1 ***********************************
!
!  adapted from :
!  http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
! 
!  I changed the routine to keep the initial labels during the sort process
!
! ******************************************************************************
SUBROUTINE merge_1(A,NA,B,NB,C,NC,labela,labelb,labelc)
  
  USE constants, ONLY : dp
  implicit none
    
  ! global
  integer       , intent(in)       :: NA,NB,NC     ! Normal usage: NA+NB = NC
  real(kind=dp) , intent(inout)    :: A(NA)        ! B overlays C(NA+1:NC)
  real(kind=dp) , intent(in)       :: B(NB)
  real(kind=dp) , intent(inout)    :: C(NC)
  integer       , intent(inout)    :: labelA(NA)       
  integer       , intent(in)       :: labelB(NB)
  integer       , intent(inout)    :: labelC(NC)

  ! local
  integer :: i,j,k

  !dir$ vector aligned
  i = 1; j = 1; k = 1;
  do while(i .le. NA .and. j .le. NB)
    if (A ( i ) .le. B ( j )) then
      C ( k ) = A ( i )
      labelc ( k ) = labela ( i ) 
      i = i + 1
    else
      C ( k ) = B ( j )
      labelc ( k ) = labelb ( j )
      J = J+1
    endif
    k = k + 1
  enddo
  do while (i .le. NA)
    C ( k ) = A ( i )
    labelc ( k ) = labela ( i )
    i = i + 1
    k = k + 1
  enddo

  return
 
END SUBROUTINE merge_1
 
! *********************** SUBROUTINE merge_sort ********************************
!
!  adapted from :
!  http://rosettacode.org/wiki/Sorting_algorithms/Merge_sort#Fortran
! 
!  I changed the routine to keep the initial labels during the sort process
!
!  A : input array of size N
!  labela : 1 --- N
!
! !example merge_sort
!  a(1)=6.0
!  a(2)=2.0
!  a(3)=-2.0
!  a(4)=3.0
!  a(5)=1.0

!  do ia = 1 , 5
!  labela(ia) = ia
!  enddo
!  write(*,*) 'befo',a
!  call merge_sort ( a , 5 , t , labela , labelt )
!  write(*,*) 'after1 ',a

!
! ******************************************************************************
RECURSIVE SUBROUTINE merge_sort(A,N,T,labela,labelt)
  
  USE constants, ONLY : dp

  implicit none 

  ! global
  integer, intent(in)                                :: N
  real(kind=dp), dimension(N), intent(in out)        :: A
  real(kind=dp), dimension((N+1)/2), intent (out)    :: T
  integer, dimension(N), intent(in out)              :: labelA
  integer, dimension((N+1)/2), intent(out)           :: labelt

  ! local
  real(kind=dp) :: V
  integer :: NA,NB,labelv
 
  if (N .lt. 2) return
  if (N .eq. 2) then
    if (A(1) .gt. A(2)) then
      V = A(1)
      labelv=labelA(1)
      A(1) = A(2)
      labela(1)=labelA(2)
      A(2) = V
      labela(2)=labelv
    endif
    return
  endif      
  NA=(N+1)/2
  NB=N-NA

  call merge_sort(A,NA,T,labela,labelt)
  call merge_sort(A(NA+1),NB,T,labela(NA+1),labelt)
 
  if (A(NA) .gt. A(NA+1)) then
    T(1:NA)=A(1:NA)
    labelt(1:NA)=labelA(1:NA)
    call merge_1(T,NA,A(NA+1),NB,A,N,labelt,labela(NA+1),labela)
  endif

  return
 
END SUBROUTINE merge_sort

! *********************** SUBROUTINE expro *************************************
!
! EXPRO
! caclulates the x-product of two vectors
! adapted from VASP ;) 
!
! ******************************************************************************
SUBROUTINE expro (H,U1,U2)

  USE constants, ONLY : dp
  IMPLICIT none 
  real(kind=dp), dimension ( 3 ) :: H ,U1 ,U2

  H(1)=U1(2)*U2(3)-U1(3)*U2(2)
  H(2)=U1(3)*U2(1)-U1(1)*U2(3)
  H(3)=U1(1)*U2(2)-U1(2)*U2(1)

  RETURN

END SUBROUTINE
 
! *********************** SUBROUTINE print_config_sample ***********************
!
! print a sample of the configurations (pos , vel , force ) 
! essentially for debug purpose
!
! input : 
!          time ,rank :  just to have some information in output
!
! ******************************************************************************
SUBROUTINE print_config_sample ( time , rank )

  USE config,   ONLY :  natm , atype , itype , rx , vx , fx , qia , dipia , ipolar , massia
  USE field,    ONLY :  mu_t,theta_t
  USE mpimdff,  ONLY :  myrank
  USE io,       ONLY :  stdout

  implicit none

  ! global
  integer, intent(in) :: time

  ! local 
  integer :: ia , rank

  if ( myrank.eq.rank ) then
       bigseparator_noionode(stdout) 
       WRITE ( stdout ,'(a)') 'debug :  SAMPLE OF THE CONFIGIURATION '
       WRITE ( stdout ,'(a5,i10)') 'time = ',time
       WRITE ( stdout ,'(a5,i10)') 'rank = ',rank
       WRITE ( stdout ,'(a)') '     i    atype       itype      ipolar      q      mass    d_x    d_y    d_z             rx                vx                  fx                 mu_x         theta_xx'
    if ( natm .ge. 32)   &
       WRITE ( stdout ,'(i6,a10,i10,l10,4x,5f8.3,5f20.10)') &
       ( ia , atype ( ia ) , itype ( ia ) , ipolar ( ia ) , qia ( ia ) , massia(ia), dipia ( 1 , ia ), dipia ( 2 , ia ) ,dipia ( 3 , ia ), &
        rx ( ia ) , vx ( ia ) , fx ( ia ) , mu_t( 1 , ia ) , theta_t(1,1,ia) , ia = 1 , 16 )
    if ( natm .ge. 32)   &
       WRITE ( stdout ,'(i6,a10,i10,l10,4x,5f8.3,5f20.10)') &
       ( ia , atype ( ia ) , itype ( ia ) , ipolar ( ia ) , qia ( ia ) , massia(ia), dipia ( 1, ia  ), dipia ( 2, ia ) ,dipia ( 3, ia  ), &
        rx ( ia ) , vx ( ia ) , fx ( ia ) , mu_t( 1 , ia ), theta_t(1,1,ia) , ia = natm - 16  , natm )
    if ( natm .lt. 32)   &
       WRITE ( stdout ,'(i6,a10,i10,l10,4x,5f8.3,5f20.10)') &
       ( ia , atype ( ia ) , itype ( ia ) , ipolar ( ia ) , qia ( ia ) , massia(ia), dipia ( 1, ia ), dipia ( 2, ia ) ,dipia ( 3, ia  ), &
        rx ( ia ) , vx ( ia ) , fx ( ia ) , mu_t( 1 , ia ),  theta_t(1,1,ia) , ia = 1 , natm )
       blankline(stdout) 
       bigseparator_noionode(stdout) 
  endif

  return

END SUBROUTINE print_config_sample

 
! *********************** SUBROUTINE print_general_info ************************
!
! ******************************************************************************
SUBROUTINE print_general_info (kunit)

  USE io,  ONLY :  ionode 
  USE config,   ONLY : natm , ntype , rho , simu_cell

  implicit none

  ! global 
  integer :: kunit 
  ! local 
  integer :: i 

  if ( ionode ) then
    blankline(kunit)
    WRITE ( kunit ,'(a)')          'Remind some parameters of the system:'
    WRITE ( kunit ,'(a,i5)')       'natm            = ', natm
    WRITE ( kunit ,'(a,i5)')       'ntype           = ', ntype
    WRITE ( kunit ,'(a,f10.3)')    'density         = ', rho
    WRITE ( kunit ,'(a,3f10.3)')   'cell parameters = ', (simu_cell%ANORM(i),i=1,3)
    WRITE ( kunit ,'(a,f10.3)')    'volume          = ', simu_cell%omega
  endif

  return

END SUBROUTINE print_general_info 


! *********************** SUBROUTINE write_all_conf_proc ******************************
!
!>\brief
! write configuration (pos,vel) to CONTFF file
!
! ******************************************************************************
SUBROUTINE write_all_conf_proc

  USE constants,                ONLY :  dp
  USE mpimdff,                  ONLY :  myrank, numprocs
  USE io,                       ONLY :  kunit_bak_proc, ionode
  USE config,                   ONLY :  system , simu_cell, natm, atype, ntype, natmi , atypei , rx, ry ,rz , vx ,vy ,vz ,fx,fy,fz
  USE md,                       ONLY :  itime

  implicit none

  ! local
  integer :: ia , it
  character(len=2)  :: str1
  character(len=2)  :: str2
  character(len=60) :: str3
  character(len=60) :: filename_out

  write(str1,'(i2)' ) numprocs
  write(str2,'(i2)' ) myrank
  write(str3,'(i8)' ) itime
  filename_out=trim(adjustl( 'CONTFF.np'//trim(adjustl(str1))//'.rank'//trim(adjustl(str2))//'.step'//trim(adjustl(str3)) ) )
  !write(*,*) 'write conf to unit',kunit_bak_proc(myrank),filename_out
  OPEN ( kunit_bak_proc(myrank) ,file = filename_out , STATUS = 'UNKNOWN')
      WRITE (  kunit_bak_proc(myrank),'(i8)') natm
      WRITE (  kunit_bak_proc(myrank),'(a)') system
      WRITE (  kunit_bak_proc(myrank),'(3f20.12)') simu_cell%A ( 1 , 1 ) , simu_cell%A ( 2 , 1 ) , simu_cell%A ( 3 , 1 )
      WRITE (  kunit_bak_proc(myrank),'(3f20.12)') simu_cell%A ( 1 , 2 ) , simu_cell%A ( 2 , 2 ) , simu_cell%A ( 3 , 2 )
      WRITE (  kunit_bak_proc(myrank),'(3f20.12)') simu_cell%A ( 1 , 3 ) , simu_cell%A ( 2 , 3 ) , simu_cell%A ( 3 , 3 )
      WRITE (  kunit_bak_proc(myrank),'(i4)') ntype
      WRITE (  kunit_bak_proc(myrank),*) ( atypei(it) , it=1,ntype )
      WRITE (  kunit_bak_proc(myrank),*) ( natmi (it) , it=1,ntype )
      WRITE (  kunit_bak_proc(myrank),'(A)') 'Cartesian'
      WRITE (  kunit_bak_proc(myrank),'(a,9e20.12)') ( atype ( ia ) , rx ( ia ) , ry ( ia ) , rz ( ia ) , &
                                                                      vx ( ia ) , vy ( ia ) , vz ( ia ) , &
                                                                      fx ( ia ) , fy ( ia ) , fz ( ia ) , ia = 1 , natm )
  CLOSE (kunit_bak_proc(myrank))

  return

END SUBROUTINE write_all_conf_proc 




! *********************** SUBROUTINE dumb_guy **********************************
!
!> \brief
!!  This subroutine print out the dumb guy!
!
! ******************************************************************************

SUBROUTINE dumb_guy(kunit)

  USE io,  ONLY :  ionode 

  implicit none

  integer :: kunit

!  Here is the guy ...
  if ( ionode ) then
     WRITE ( kunit ,'(a)') '                            \\|//                            '
     WRITE ( kunit ,'(a)') '                           -(o o)-                           '
     WRITE ( kunit ,'(a)') '========================oOO==(_)==OOo========================'
  endif

  return

END SUBROUTINE dumb_guy


! *********************** SUBROUTINE check_allowed_tags **********************************
!
!> \brief
!!  This subroutine check if a specific character tag is allowed
!
! ******************************************************************************
SUBROUTINE check_allowed_tags( size_allowed , allowed_values , tag , tagsection , tagname )

  USE io,       ONLY :  ionode, stderr

  implicit none

  ! global
  integer      :: size_allowed
  character(*) :: allowed_values(size_allowed)
  character(*) :: tag
  character(*) :: tagsection
  character(*) :: tagname
  !
  ! local
  integer :: i
  logical :: allowed

  allowed = .false.
  do i = 1 , size_allowed 
    if ( trim ( tag ) .eq. allowed_values ( i ) )  allowed = .true.
  enddo
  if ( .not. allowed ) then
    if ( ionode )  WRITE ( stderr , '(a,a,a,a,a,<size_allowed>a)' ) 'ERROR ',tagsection,': ',tagname,' should be :'
    do i = 1 , size_allowed 
      if ( ionode )  WRITE ( stderr ,'(a)') allowed_values(i)
    enddo
    STOP
  endif

  return

END SUBROUTINE check_allowed_tags

! ===== fmV =====
