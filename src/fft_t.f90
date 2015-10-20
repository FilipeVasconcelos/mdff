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

! *********************** SUBROUTINE fft_1D_complex ****************************
!
! driver for fftw routines
!
! ******************************************************************************

SUBROUTINE fft_1D_complex(in,out,N)

  USE constants , ONLY : dp 
  implicit none
  INCLUDE "fftw3.f"
!  INCLUDE "mpif.h"

  integer :: N
  complex(kind=dp), dimension(N) :: in, out
  integer plan

  CALL dfftw_plan_dft_1d ( plan , N  , in , out , FFTW_FORWARD , FFTW_ESTIMATE )
  CALL dfftw_execute_dft ( plan , in , out )
  CALL dfftw_destroy_plan( plan )

  return 

END SUBROUTINE fft_1D_complex

! *********************** SUBROUTINE fft_1D_complex ****************************
!
! driver for fftw routines
!
! ******************************************************************************

SUBROUTINE fft_1D_real(in,out,N)
  
  USE constants , ONLY : dp 
  implicit none
  INCLUDE "fftw3.f"

  integer :: N
  real(kind=dp), dimension(N)      :: in
  complex(kind=dp)  ,dimension(N/2 +1 ) :: out
  integer*8 plan

  call dfftw_plan_dft_r2c_1d(plan,N,in,out,FFTW_ESTIMATE)
  call dfftw_execute_dft_r2c(plan, in, out)
  call dfftw_destroy_plan(plan)

  return

END SUBROUTINE fft_1D_real
! ===== fmV =====
