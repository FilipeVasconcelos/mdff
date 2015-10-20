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
! ======= Hardware =======

!> \brief
!> complementary error function
!> \author
!> W. Press et al. : Numerical Recipes, p. 214
!> \note
!> built-in erfc() is buggy e.g. on some Linux distributions
  real(kind=dp) FUNCTION errfc(x)
  USE constants, ONLY : dp
  implicit none
  ! local
  real(kind=dp) :: x , z , t    
  z=ABS(x)
  t=1._dp/(1._dp+0.5_dp*z)
  errfc=t*EXP(-z*z-1.26551223_dp+t*(1.00002368_dp+t*(.37409196_dp+ &
  t*(.09678418_dp+t*(-.18628806_dp+t*(.27886807_dp+t*(-1.13520398_dp+ &
  t*(1.48851587_dp+t*(-.82215223_dp+t*.17087277_dp)))))))))
  if(x.lt.0._dp) errfc=2._dp-errfc
  return
  END FUNCTION

!> \brief
!> error function
  real(kind=dp) FUNCTION errf(x)
  USE constants, ONLY : dp
  implicit none
  ! local
  real(kind=dp) :: x , errfc   
  errf=1.0_dp-errfc(x)
  END FUNCTION
! ===== fmV =====
