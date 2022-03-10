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
!#define debug_TT
! ======= Hardware =======
MODULE tt_damp


  USE constants,        ONLY :  dp
  USE io,               ONLY :  stdout, ionode

  integer      , PARAMETER                             :: maximum_of_TT_expansion=8
  real(kind=dp), dimension (0:maximum_of_TT_expansion) :: E_TT 

CONTAINS

! *********************** SUBROUTINE get_TT_damp *******************************
!> \brief
!!  
!> \author
!! FMV
!> \date 
!! February 2014
! ******************************************************************************
SUBROUTINE get_TT_damp 

  implicit none

  integer :: k 

  E_TT(0) = 1.0_dp
  do k = 1 , maximum_of_TT_expansion
    E_TT(k) = E_TT(k-1) / REAL ( k ,kind=dp )
#ifdef debug_TT  
    io_node write( stdout , '(a,i0,e20.12)') 'TT_damp coeff =',k,E_TT(k)
#endif
  enddo

  return

END SUBROUTINE get_TT_damp

! *********************** SUBROUTINE TT_damping_functions **********************
!> \brief
!! Tang-Toennies damping function
!> \author
!! FMV
!> \date 
!! February 2014
! ******************************************************************************
SUBROUTINE TT_damping_functions(b,c,r,f,fd,order)

!  USE tt_damp,          ONLY : E_TT ! Tang-Toennies coefficients define once up to order 8 

  implicit none

  ! global 
  integer :: order
  real(kind=dp) :: b , c , r, f , fd  ! f damping function , fd first derivative

  ! local 
  integer :: k
  real(kind=dp) :: expbdr , br

  br = b * r
  expbdr = EXP(-br) * c

  f = E_TT(order)
  do k=order-1,1,-1
    f = f * br + E_TT(k)
  enddo
  f = f * br + E_TT(0)
  f = 1.0_dp - f * expbdr

  ! derivative of f checked on sage 18/02/14 worksheet BMHFTD
  fd = E_TT(order) * ( br ) ** ( order ) * expbdr * b

  return

END SUBROUTINE TT_damping_functions



END MODULE tt_damp
