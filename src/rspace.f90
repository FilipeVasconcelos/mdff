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

! *********************** MODULE rspace ****************************************
!> @brief
!> module related to real space summation 
!> @author
!> FMV
! ******************************************************************************
MODULE rspace

  USE constants,                ONLY :  dp

  implicit none      

  ! ===================
  !  direct summation
  ! ===================
  TYPE rmesh
    integer                                     :: ncell     !< nb of cell in the direct calc. in each direction  
    integer                                     :: ncmax     !< (internal) number of cell in direct summation
    real(kind=dp), dimension(:,:) , allocatable :: boxxyz    !< vectors in real space
    real(kind=dp), dimension(:)   , allocatable :: rr        !< vectors module in real space
    integer      , dimension(:)   , allocatable :: lcell     !< lcell is 1 in the central box 
    character(len=15)                           :: meshlabel !< label name 
  END TYPE rmesh

CONTAINS

! *********************** SUBROUTINE direct_sum_init ***************************
!!
!> @brief 
!> generate vectors in real space of the neigboring cells 
!> in a cubic cutoff for the direct summation 
!
! cubic cutoff: -ncelldirect to ncelldirect in each direction
! vectors are stored in boxxyz(alpha,nc) alpha = 1,2,3 (x,y,z)
! and nc is the index of the given neighboring cell 
! lcell (nc) is set to 1 if nc is not the central box (ncell !=0)
!> @author 
!> FMV
!
! ******************************************************************************
SUBROUTINE direct_sum_init ( rm )

  USE io,       ONLY : ionode , stdout 

  implicit none
 
  ! global
  TYPE ( rmesh ) :: rm

  ! local
  integer :: nc , ncellx , ncelly , ncellz , ncelldirect

  io_node WRITE ( stdout      ,'(a,a,a)') 'generate real-space  (full) ',rm%meshlabel,' mesh'

  ncelldirect = rm%ncell

  nc = 0
  rm%lcell = 0
  do ncellx = -ncelldirect,ncelldirect
    do ncelly = -ncelldirect,ncelldirect
      do ncellz = -ncelldirect,ncelldirect
        nc = nc + 1
        rm%boxxyz(1,nc) = DBLE (ncellx)
        rm%boxxyz(2,nc) = DBLE (ncelly)
        rm%boxxyz(3,nc) = DBLE (ncellz)
        rm%rr(nc) = rm%boxxyz(1,nc) * rm%boxxyz(1,nc) + &
                rm%boxxyz(2,nc) * rm%boxxyz(2,nc) + &
                rm%boxxyz(3,nc) * rm%boxxyz(3,nc)  
        if (ncellx .ne. 0 .or. ncelly .ne. 0 .or. ncellz .ne. 0) then
          rm%lcell(nc) = 1
        endif
      enddo
    enddo
  enddo

  if ( nc .ne. rm%ncmax ) then
    io_node WRITE ( stdout ,'(a,3i7,a)') 'number of ncells do not match in direct_sum_init', nc , rm%ncmax , rm%meshlabel
    STOP
  endif

  ! ======================
  !  organized rpt arrays 
  ! ======================
  call reorder_rpt ( rm )
  io_node WRITE ( stdout      ,'(a)') '(full) real space arrays sorted'
  io_node WRITE ( stdout      ,'(i9)') nc 

  return

END SUBROUTINE direct_sum_init

! *********************** SUBROUTINE reorder_rpt *******************************
!
! this subroutine reorder kpt arrays (increasing k^2)
!
! ******************************************************************************

SUBROUTINE reorder_rpt ( rm )

  USE io,  ONLY :  stdout 

  implicit none

  !global
  TYPE(rmesh) :: rm

  !local
  integer :: ir, lr
  real(kind=dp), dimension (:), allocatable :: trpt
  real(kind=dp), dimension (:), allocatable :: tmprx , tmpry , tmprz

  integer, dimension (:), allocatable :: labelrpt, labelt , tmplc

  allocate ( trpt ((rm%ncmax+1)/2) , labelrpt(rm%ncmax), labelt((rm%ncmax+1)/2))
  allocate ( tmprx(rm%ncmax) , tmpry(rm%ncmax) , tmprz(rm%ncmax)  , tmplc(rm%ncmax) )

  ! ==============================
  !  set the initial array labels
  ! ==============================
  do ir=1,rm%ncmax
    labelrpt(ir)=ir
  enddo

  ! ===========================================
  !  arrays are sorted for increasing rr^2 
  !  (see tools.f90 for more details )
  !  the old labels are stored in labelkpt
  !  that will be used to reorganized kx,ky,kz
  ! ===========================================
  call merge_sort ( rm%rr , rm%ncmax , trpt , labelrpt , labelt )

  ! ==============================
  !  save previous order
  ! ==============================
  tmprx=rm%boxxyz(1,:)
  tmpry=rm%boxxyz(2,:)
  tmprz=rm%boxxyz(3,:)
  tmplc=rm%lcell(:)
  ! ===============================================
  !  change boxxyz, lcell following rr sort
  ! ===============================================
  do ir=1,rm%ncmax
    lr=labelrpt(ir)
    rm%boxxyz(1,ir) = tmprx(lr)
    rm%boxxyz(2,ir) = tmpry(lr)
    rm%boxxyz(3,ir) = tmprz(lr)
    rm%lcell( ir )  = tmplc(lr)
  enddo

#ifdef debug
  do ir=1,rm%ncmax
    WRITE (stdout , '(4f16.6,i6)' ) rm%boxxyz(1,ir),rm%boxxyz(2,ir),rm%boxxyz(3,ir),rm%rr(ir),rm%lcell(ir)
  enddo
#endif

  deallocate ( trpt , labelrpt , labelt )
  deallocate ( tmprx , tmpry , tmprz    )

  return

END SUBROUTINE reorder_rpt



END MODULE rspace
! ===== fmV =====
