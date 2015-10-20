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

! *********************** MODULE CONSTANTS *************************************
!> \author FMV
!> \brief Main constants and units definition of the code
!> \note
!!
!! Units of MDFF in input/output : 
!! length        : angstom 
!! energy        : eV 
!! mass          : atomic mass unit
!! time          : ps
!! temperature   : Kelvin
!! pressure      : GPa
!! charge        : atomic charge ( proton charge )
!! dipole moment : Debye
!! efg           : V  / angstom**2 
!! forces        : eV / angtrom
!!
!! Internal units:
!! length        : angstom 
!! energy        : eV 
!! mass          : atomic mass unit
!! time          : angstrom * ( atomicmassunit / eV ) ** 0.5 
!! temperature   : eV
!! pressure      : ( eV / angstrom**3)
!! charge        : atomic charge ( proton charge )
!! dipole moment : angstrom * atomiccharge 
!! 1/4piepsilon0 : eV  * angstrom / atomiccharge **2 == 1.0
!! efg           : atomiccharge / angstrom ** 3  

! ******************************************************************************
MODULE constants 

  implicit none

  save

  ! -------------------------------------------------------------------------------------------------------------------
  ! kind definitions 
 
  integer, PARAMETER, PUBLIC    :: dp          = selected_real_kind(15,300)            !< double precision definition 
  integer, PARAMETER, PUBLIC    :: sp          = selected_real_kind(6,30)              !< single precision definition 
  ! -------------------------------------------------------------------------------------------------------------------
  ! mathematical constants
  real(kind=dp),      PARAMETER :: pi          = 3.14159265358979323846264338327950_dp !< pi
  real(kind=dp),      PARAMETER :: dzero       = 0.0_dp                                !< zero 
  real(kind=dp),      PARAMETER :: done        = 1.0_dp                                !< one 
  real(kind=dp),      PARAMETER :: pisq        = pi * pi                               !< pi^2
  real(kind=dp),      PARAMETER :: tpi         = 2.0_dp*pi                             !< 2pi
  real(kind=dp),      PARAMETER :: fpi         = 2.0_dp*tpi                            !< 4pi
  real(kind=dp),      PARAMETER :: piroot      = SQRT ( pi )                           !< sqrt(pi)
  real(kind=dp),      PARAMETER :: radian      = 180.0_dp / pi                         !< radian 
  complex(kind=dp),   PARAMETER :: imag        = (0.0_dp, 1.0_dp)                      !< imaginary number i
  complex(kind=dp),   PARAMETER :: mimag       = (0.0_dp,-1.0_dp)                      !< negative imaginary number (-imag)
  complex(kind=dp),   PARAMETER :: citpi       = imag*tpi                              !< 2*i*pi
  ! -------------------------------------------------------------------------------------------------------------------
  ! physical constants
  real(kind=dp),      PARAMETER :: rytoev      = 13.6056919282781_dp                   !< Rydberg constant a.u. to eV
  real(kind=dp),      PARAMETER :: hart        = rytoev*2.0_dp                         !< Hartree energy 
  real(kind=dp),      PARAMETER :: bohr        = 0.52917721092_dp                      !< a.u in angstrom  
  real(kind=dp),      PARAMETER :: evtoj       = 1.602176487e-19_dp                    !< eV to J 
  real(kind=dp),      PARAMETER :: hplanck     = 6.62617636e-34_dp                     !< Planck constant 
  real(kind=dp),      PARAMETER :: eh          = evtoj / 1e14_dp / hplanck             !< e/h (electron charge / Planck constant )
  real(kind=dp),      PARAMETER :: cq_unit     = eh / 1000.0_dp                        !< Cq unit conversion 
  real(kind=dp),      PARAMETER :: g_to_am     = 1.660538782_dp                        !< gram to atomic mass unit 
  ! -------------------------------------------------------------------------------------------------------------------
  ! units related
  real(kind=dp)      :: press_unit  = 0.00624150964712042_dp                !< GPa = > internal unit of pressure ( eV / angstrom**3) 
  real(kind=dp)      :: boltz_unit  = 8.61734229648141e-05                  !< boltzmann constant ( energy in eV)
  real(kind=dp)      :: time_unit   = 98.2269514139276_dp                   !< unit of time picosecond => angstrom * ( atomicmassunit / eV ) ** 0.5
  real(kind=dp)      :: Debye_unit  = 2.54174618782479816355_dp / bohr      !< Debye unit for dipole moments in eA
  real(kind=dp)      :: coul_unit   = 14.3996441494161_dp                   !< 1 / 4pi epsilon0  in eV



CONTAINS 
!------------------------------------------------------------------------------!
!
!   Print information about the used data types.
!   from Quantum-espresso package
!
SUBROUTINE print_kind_info (stdout)
!
!------------------------------------------------------------------------------!
!
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: stdout
!
        WRITE( stdout,'(/,T2,A)') 'DATA TYPE INFORMATION:'
!
        WRITE( stdout,'(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8))') &
          'REAL: Data type name:', 'DP', '      Kind value:', kind(0.0_DP), &
          '      Precision:', precision(0.0_DP), &
          '      Smallest nonnegligible quantity relative to 1:', &
          epsilon(0.0_DP), '      Smallest positive number:', tiny(0.0_DP), &
          '      Largest representable number:', huge(0.0_DP)
        WRITE( stdout,'(/,T2,A,T78,A,2(/,T2,A,T75,I6),3(/,T2,A,T67,E15.8))') &
          '      Data type name:', 'sp', '      Kind value:', kind(0.0_SP), &
          '      Precision:', precision(0.0_SP), &
          '      Smallest nonnegligible quantity relative to 1:', &
          epsilon(0.0_SP), '      Smallest positive number:', tiny(0.0_SP), &
          '      Largest representable number:', huge(0.0_SP)
        WRITE( stdout,'(/,T2,A,T72,A,4(/,T2,A,T61,I20))') &
          'INTEGER: Data type name:', '(default)', '         Kind value:', &
          kind(0), '         Bit size:', bit_size(0), &
          '         Largest representable number:', huge(0)
        WRITE( stdout,'(/,T2,A,T72,A,/,T2,A,T75,I6,/)') 'LOGICAL: Data type name:', &
          '(default)', '         Kind value:', kind(.TRUE.)
        WRITE( stdout,'(/,T2,A,T72,A,/,T2,A,T75,I6,/)') &
          'CHARACTER: Data type name:', '(default)', '           Kind value:', &
          kind('C')
!
      END SUBROUTINE print_kind_info
!
!------------------------------------------------------------------------------!


   SUBROUTINE reduced_units

     implicit none

     time_unit   = 1.0_dp
     coul_unit   = 1.0_dp
     press_unit  = 1.0_dp
     boltz_unit  = 1.0_dp

     return

   END SUBROUTINE reduced_units


END MODULE constants
! ===== fmV =====
