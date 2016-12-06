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

! ======= Hardware =======
#include "symbol.h"
!#define debug 
!#define GFORTRAN
! ======= Hardware =======


MODULE restart

  USE io,  ONLY :  ionode, stdout

CONTAINS

SUBROUTINE restart_init ( MDFF ) 

  USE io,       ONLY :  kunit_RESTART
  USE cell
  USE control
  USE config
  USE non_bonded
  USE field
  USE tt_damp
  USE md

  implicit none

  ! global
  character(len=80)  :: MDFF

  ! local 
  integer :: i
  character(len=60) :: cpos
  character(len=20) :: FMT 

  CALL print_RESTART_info ( stdout )
 
  WRITE(stdout,'(a)') 'end of restart_init'

  return

END SUBROUTINE restart_init


SUBROUTINE write_RESTART

  USE io,  ONLY :  kunit_RESTART
  USE cell
  USE control
  USE config
  USE field
  USE non_bonded
  USE md

  implicit none

  return

END SUBROUTINE write_RESTART

SUBROUTINE print_RESTART_info ( kunit ) 

  implicit none

  ! local 
  integer :: kunit

  if ( ionode ) call dumb_guy(kunit)
  if ( ionode ) then
     WRITE ( kunit ,'(a)')       "       _____  ______  _____ _______       _____ _______ "
     WRITE ( kunit ,'(a)')       "      |  __ \|  ____|/ ____|__   __|/\   |  __ \__   __|"
     WRITE ( kunit ,'(a)')       "      | |__) | |__  | (___    | |  /  \  | |__) | | |   "
     WRITE ( kunit ,'(a)')       "      |  _  /|  __|  \___ \   | | / /\ \ |  _  /  | |   "
     WRITE ( kunit ,'(a)')       "      | | \ \| |____ ____) |  | |/ ____ \| | \ \  | |   "
     WRITE ( kunit ,'(a)')       "      |_|  \_\______|_____/   |_/_/    \_\_|  \_\ |_|   "
     blankline(kunit)
     separator(kunit)
     blankline(kunit)
     WRITE ( kunit ,'(a)')       "WARING : the full restart mode is only reading parameters in"
     WRITE ( kunit ,'(a)')       "the local RESTART file. The control file is not used at all" 
     WRITE ( kunit ,'(a)')       "(except full_restart tag;)"
     WRITE ( kunit ,'(a)')       "Be carefull !! "  
     WRITE ( kunit ,'(a)')       "Some of the most important parameters read from the RESTART "
     WRITE ( kunit ,'(a)')       "file are reminded below. But not all of them !!!"

     blankline(kunit)
  endif

  return

END SUBROUTINE print_RESTART_info

END MODULE restart
