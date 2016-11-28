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

! *********************** MODULE calc_thermo ***********************************
!> \brief
!! Module related to thermodynamic quantities
! ******************************************************************************
MODULE thermodynamic

  USE io,                       ONLY :  ionode , ioprint 
  USE constants,                ONLY :  dp
  USE block,                    ONLY :  accu
  USE md,                       ONLY :  nve_ensemble , nvt_ensemble , npt_ensemble , npe_ensemble

  implicit none

  real(kind=dp) :: h_tot          !< general conserved quantity (ensemble dependent) 
  real(kind=dp) :: e_tot          !< total energy  ( potential + kinetic )
  real(kind=dp) :: u_tot          ! total potential energy 
  
  real(kind=dp) :: e_kin          !< kinetic energy
  real(kind=dp) :: u_lj           !< potential energy from lennard_jones interaction
  real(kind=dp) :: u_morse        !< potential energy from morse interaction
  real(kind=dp) :: u_bmhft        !< potential energy from bmhft interaction
  real(kind=dp) :: u_coul
  real(kind=dp) :: u_pol 
  real(kind=dp) :: u_harm         !< harmonic oscilaltor lharm  
  real(kind=dp) :: e_nvt          !< Nose Hoover thermostat "energy" (conserved quantity) 
  real(kind=dp) :: e_npt          !< Nose-Hoover-Andersen (thermostat+barostat) "energy" (conserved quantity) 
  real(kind=dp) :: csvr_conint    !< CSVR Stochastic velocity rescaling (conserved quantity) 

  ! output (correct units and reduced or not )
  real(kind=dp) :: e_kin_r        !< kinetic energy
  real(kind=dp) :: temp_r         !< temperature ( from kinetic energy )  
  real(kind=dp) :: u_lj_r         !< potential energy from lennard_jones interaction
  real(kind=dp) :: u_morse_r      !< potential energy from morse interaction
  real(kind=dp) :: u_bmhft_r      !< potential energy from bmhft interaction
  real(kind=dp) :: u_coul_r       !< potential energy from coulombic interaction 
  real(kind=dp) :: e_nvt_r        !< further Nose Hoover energy for the conserved quantity 
  real(kind=dp) :: e_npt_r        !< further Nose Hoover energy for the conserved quantity 
  real(kind=dp) :: csvr_conint_r  !< CSVR Stochastic velocity rescaling (conserved quantity) 


  real(kind=dp) :: vir_tot        !< total virial
  real(kind=dp) :: vir_lj         !< virial of lj interaction
  real(kind=dp) :: vir_morse      !< virial of morse interaction
  real(kind=dp) :: vir_bmhft      !< virial of bmhft interaction
  real(kind=dp) :: vir_coul_tot   !< virial of coulombic interaction
!  real(kind=dp) :: vir_coul_qq    !< virial of coulombic interaction
!  real(kind=dp) :: vir_coul_dd    !< virial of coulombic interaction
!  real(kind=dp) :: vir_coul_qd    !< virial of coulombic interaction

  real(kind=dp) :: pvirial_tot    !< virial correction to the pressure
  real(kind=dp) :: pvirial_lj     !< virial correction to the pressure
  real(kind=dp) :: pvirial_morse  !< virial correction to the pressure
  real(kind=dp) :: pvirial_bmhft  !< virial correction to the pressure
  real(kind=dp) :: pvirial_coul   !< virial correction to the pressure
  real(kind=dp) :: pvirial_tot_r  !< virial correction to the pressure
  real(kind=dp) :: pvirial_lj_r   !< virial correction to the pressure
  real(kind=dp) :: pvirial_morse_r!< virial correction to the pressure
  real(kind=dp) :: pvirial_bmhft_r!< virial correction to the pressure
  real(kind=dp) :: pvirial_coul_r !< virial correction to the pressure

  real(kind=dp) :: pressure_tot   !< total pressure
  real(kind=dp) :: pressure_lj    !< pressure from lj interactions
  real(kind=dp) :: pressure_morse !< pressure from lj interactions
  real(kind=dp) :: pressure_bmhft !< pressure from lj interactions
  real(kind=dp) :: pressure_coul  !< pressure from coulombic interactions
  real(kind=dp) :: pressure_tot_r !< total pressure
  real(kind=dp) :: pressure_lj_r  !< pressure from lj interactions
  real(kind=dp) :: pressure_morse_r  !< pressure from lj interactions
  real(kind=dp) :: pressure_bmhft_r  !< pressure from lj interactions
  real(kind=dp) :: pressure_coul_r!< pressure from coulombic interactions

  TYPE (accu) acc_e_tot          !< total energy  ( potential + kinetic ) 
  TYPE (accu) acc_u_tot          !< total potential energy 
  
  TYPE (accu) acc_e_kin_r        !< kinetic energy
  TYPE (accu) acc_u_lj_r         !< potential energy from lennard_jones interaction
  TYPE (accu) acc_u_coul_r       !< potential energy from coulombic interaction 
  TYPE (accu) acc_temp_r         !< temperature ( from kinetic energy )  

  TYPE (accu) acc_vir_tot        !< total virial
  TYPE (accu) acc_vir_lj         !< virial of lj interaction
  TYPE (accu) acc_vir_coul       !< virial of coulombic interaction

  TYPE (accu) acc_pressure_tot   !< total pressure
  TYPE (accu) acc_pressure_lj    !< pressure from lj interactions
  TYPE (accu) acc_pressure_coul  !< pressure from coulombic interactions

CONTAINS

! *********************** SUBROUTINE calc_thermo *******************************
!> \brief
!! calculation of the main thermodynamic quantities 
! ******************************************************************************
SUBROUTINE calc_thermo

  USE constants,                ONLY :  boltz_unit , press_unit
  USE control,                  ONLY :  lreducedN , lcsvr
  USE config,                   ONLY :  natm ,  simu_cell 
  USE md,                       ONLY :  integrator, press!, itime

  implicit none
  real(kind=dp) :: omega, pv

  omega = simu_cell%omega
 !!!! WARNING virials are wrong

  if (lreducedN) then
    u_lj_r   = u_lj         / REAL ( natm , kind = dp )
    u_morse_r= u_morse      / REAL ( natm , kind = dp )
    u_bmhft_r= u_bmhft      / REAL ( natm , kind = dp )
    u_coul_r = u_coul       / REAL ( natm , kind = dp )
    e_kin_r  = e_kin        / REAL ( natm , kind = dp )
    e_nvt_r  = e_nvt        / REAL ( natm , kind = dp )
    e_npt_r  = e_npt        / REAL ( natm , kind = dp )
    csvr_conint_r = csvr_conint / REAL ( natm , kind = dp )
  else
    u_lj_r   = u_lj  
    u_morse_r= u_morse  
    u_bmhft_r= u_bmhft
    u_coul_r = u_coul
    e_kin_r  = e_kin     
    e_nvt_r  = e_nvt    
    e_npt_r  = e_npt   
    csvr_conint_r = csvr_conint 
  endif
 
  u_tot    = u_lj_r + u_coul_r + u_morse_r + u_harm + u_bmhft_r
  e_tot    = u_tot  + e_kin_r

  vir_tot  = vir_lj + vir_coul_tot + vir_morse + vir_bmhft

  pvirial_lj    = vir_lj       / omega 
  pvirial_morse = vir_morse    / omega 
  pvirial_bmhft = vir_bmhft    / omega 
  pvirial_tot   = pvirial_lj + pvirial_coul + pvirial_morse + pvirial_bmhft
  pressure_tot  = pvirial_tot +  temp_r * boltz_unit / omega
  pressure_lj   = pvirial_lj   
  pressure_morse= pvirial_morse
  pressure_bmhft= pvirial_bmhft
  pressure_coul = pvirial_coul 
  if (lreducedN) then
    pvirial_lj_r    = pvirial_lj    / REAL ( natm , kind = dp )
    pvirial_morse_r = pvirial_morse / REAL ( natm , kind = dp ) 
    pvirial_coul_r  = pvirial_coul  / REAL ( natm , kind = dp ) 
    pvirial_bmhft_r = pvirial_bmhft / REAL ( natm , kind = dp ) 
  else
    pvirial_lj_r    = pvirial_lj   
    pvirial_morse_r = pvirial_morse 
    pvirial_coul_r  = pvirial_coul  
    pvirial_bmhft_r = pvirial_bmhft 
  endif
  pvirial_tot_r   = pvirial_lj_r + pvirial_coul_r + pvirial_morse_r + pvirial_bmhft_r 
  pressure_tot_r  = ( pvirial_tot_r + temp_r * boltz_unit / omega ) / press_unit
  pressure_lj_r   = pvirial_lj_r    / press_unit
  pressure_morse_r= pvirial_morse_r / press_unit   
  pressure_bmhft_r= pvirial_bmhft_r / press_unit  
  pressure_coul_r = pvirial_coul_r  / press_unit

  ! conserved quantity extended Hamiltonian
  if ( any ( integrator .eq. nve_ensemble )) then 
    h_tot = e_tot
    if ( lcsvr ) h_tot = e_tot + csvr_conint_r  ! special case of stochastic velocity rescaling with a specific conserved quantity 
  endif
  if ( any ( integrator .eq. npe_ensemble )) then
    pv = press * omega 
    h_tot = e_tot + pv
    io_print print*,'PV',pv
  endif
  if ( any ( integrator .eq. nvt_ensemble )) then
    !print*,'NVT ensemble'
    if ( integrator .eq. 'nvt-and' ) continue ! no conserved quantity for andersen thermostat
    h_tot = e_tot + e_nvt_r
  endif
  if ( any ( integrator .eq. npt_ensemble )) then
  !  print*,itime,'NPT ensemble',e_npt_r
    h_tot = e_tot + e_npt_r
  endif

  return

END SUBROUTINE calc_thermo


! *********************** SUBROUTINE init_general_accumulator ******************
!> \brief
!! initializationj of accumulators
! ******************************************************************************
SUBROUTINE init_general_accumulator

  implicit none

  acc_e_tot%accval           = 0.0_dp 
  acc_e_tot%accvalsq         = 0.0_dp
  acc_e_tot%counter          = 0

  acc_u_tot%accval           = 0.0_dp
  acc_u_tot%accvalsq         = 0.0_dp
  acc_u_tot%counter          = 0

  acc_e_kin_r%accval         = 0.0_dp
  acc_e_kin_r%accvalsq       = 0.0_dp
  acc_e_kin_r%counter        = 0

  acc_u_lj_r%accval          = 0.0_dp
  acc_u_lj_r%accvalsq        = 0.0_dp
  acc_u_lj_r%counter         = 0
 
  acc_u_coul_r%accval        = 0.0_dp
  acc_u_coul_r%accvalsq      = 0.0_dp
  acc_u_coul_r%counter       = 0

  acc_temp_r%accval          = 0.0_dp
  acc_temp_r%accvalsq        = 0.0_dp
  acc_temp_r%counter         = 0

  acc_vir_tot%accval         = 0.0_dp
  acc_vir_tot%accvalsq       = 0.0_dp
  acc_vir_tot%counter        = 0

  acc_vir_lj%accval          = 0.0_dp
  acc_vir_lj%accvalsq        = 0.0_dp
  acc_vir_lj%counter         = 0

  acc_vir_coul%accval        = 0.0_dp
  acc_vir_coul%accvalsq      = 0.0_dp
  acc_vir_coul%counter       = 0

  acc_pressure_tot%accval    = 0.0_dp
  acc_pressure_tot%accvalsq  = 0.0_dp
  acc_pressure_tot%counter   = 0

  acc_pressure_lj%accval     = 0.0_dp
  acc_pressure_lj%accvalsq   = 0.0_dp
  acc_pressure_lj%counter    = 0

  acc_pressure_coul%accval   = 0.0_dp
  acc_pressure_coul%accvalsq = 0.0_dp
  acc_pressure_coul%counter  = 0

  return

END SUBROUTINE init_general_accumulator

! *********************** SUBROUTINE general_accumulator ***********************
!> \brief
!! add values to accumulators
!> \note
!! not really used 
! ******************************************************************************
SUBROUTINE general_accumulator

  implicit none

  acc_e_tot%accval           = acc_e_tot%accval           + e_tot
  acc_e_tot%accvalsq         = acc_e_tot%accvalsq         + e_tot * e_tot
  acc_e_tot%counter          = acc_e_tot%counter          + 1

  acc_u_tot%accval           = acc_u_tot%accval           + u_tot
  acc_u_tot%accvalsq         = acc_u_tot%accvalsq         + u_tot * u_tot
  acc_u_tot%counter          = acc_u_tot%counter          + 1

  acc_e_kin_r%accval         = acc_e_kin_r%accval         + e_kin_r
  acc_e_kin_r%accvalsq       = acc_e_kin_r%accvalsq       + e_kin_r * e_kin_r
  acc_e_kin_r%counter        = acc_e_kin_r%counter        + 1

  acc_u_lj_r%accval          = acc_u_lj_r%accval          + u_lj_r
  acc_u_lj_r%accvalsq        = acc_u_lj_r%accvalsq        + u_lj_r * u_lj_r
  acc_u_lj_r%counter         = acc_u_lj_r%counter         + 1
 
  acc_u_coul_r%accval        = acc_u_coul_r%accval        + u_coul_r
  acc_u_coul_r%accvalsq      = acc_u_coul_r%accvalsq      + u_coul_r * u_coul_r
  acc_u_coul_r%counter       = acc_u_coul_r%counter       + 1

  acc_temp_r%accval          = acc_temp_r%accval          + temp_r
  acc_temp_r%accvalsq        = acc_temp_r%accvalsq        + temp_r * temp_r
  acc_temp_r%counter         = acc_temp_r%counter         + 1

  acc_vir_tot%accval         = acc_vir_tot%accval         + vir_tot
  acc_vir_tot%accvalsq       = acc_vir_tot%accvalsq       + vir_tot * vir_tot 
  acc_vir_tot%counter        = acc_vir_tot%counter        + 1

  acc_vir_lj%accval          = acc_vir_lj%accval          + vir_lj
  acc_vir_lj%accvalsq        = acc_vir_lj%accvalsq        + vir_lj * vir_lj
  acc_vir_lj%counter         = acc_vir_lj%counter         + 1

  acc_vir_coul%accval        = acc_vir_coul%accval        + vir_coul_tot
  acc_vir_coul%accvalsq      = acc_vir_coul%accvalsq      + vir_coul_tot * vir_coul_tot
  acc_vir_coul%counter       = acc_vir_coul%counter       + 1

  acc_pressure_tot%accval    = acc_pressure_tot%accval    + pressure_tot
  acc_pressure_tot%accvalsq  = acc_pressure_tot%accvalsq  + pressure_tot * pressure_tot
  acc_pressure_tot%counter   = acc_pressure_tot%counter   + 1

  acc_pressure_lj%accval     = acc_pressure_lj%accval     + pressure_lj
  acc_pressure_lj%accvalsq   = acc_pressure_lj%accvalsq   + pressure_lj * pressure_lj
  acc_pressure_lj%counter    = acc_pressure_lj%counter    + 1

  acc_pressure_coul%accval   = acc_pressure_coul%accval   + pressure_coul
  acc_pressure_coul%accvalsq = acc_pressure_coul%accvalsq + pressure_coul * pressure_coul
  acc_pressure_coul%counter  = acc_pressure_coul%counter  + 1

  return

END SUBROUTINE general_accumulator

! *********************** SUBROUTINE write_thermo ******************************
!> \brief
!! write thermodynamic quantities to file OSZIFF or print to standard output 
! ******************************************************************************
SUBROUTINE write_thermo ( step , kunit , key )

  USE constants,        ONLY :  time_unit
  USE config,           ONLY :  simu_cell 
  USE md,               ONLY :  dt

  implicit none

  ! global
  integer          , intent(in)           :: kunit , step
  character(len=3) , intent (in)          :: key 

  ! local 
  real(kind=dp) :: omega , acell ,bcell , ccell
  real(kind=dp) :: u_vdw_r , pvirial_vdw_r
  
  omega = simu_cell%omega
  acell = simu_cell%ANORM(1)
  bcell = simu_cell%ANORM(2)
  ccell = simu_cell%ANORM(3)

  u_vdw_r = u_lj_r + u_bmhft_r + u_morse_r
  pvirial_vdw_r   = pvirial_lj_r + pvirial_morse_r + pvirial_bmhft_r 

  if ( key .eq. 'osz' ) then
    if ( ionode ) then
        WRITE ( kunit , 200 ) &
        step , REAL ( step * dt / time_unit , kind = dp ) , e_tot   , e_kin_r      , u_tot        , u_vdw_r      , u_coul_r   
        WRITE ( kunit , 201 ) &
        step , REAL ( step * dt / time_unit , kind = dp ) , temp_r  , pressure_tot_r , pvirial_vdw_r , pvirial_coul_r , omega , h_tot  
    endif
  endif
  if ( key .eq. 'std' ) then
        io_node WRITE ( kunit, 300)  step , REAL ( step * dt / time_unit , kind = dp ) , &
                                     e_kin_r , temp_r , u_tot  , u_vdw_r , u_coul_r  , &
                                     pressure_tot_r , pvirial_vdw_r , pvirial_coul_r , pvirial_tot_r, omega ,    &
                                     acell , bcell, ccell , e_tot, h_tot
  endif

 200 FORMAT(' step = ',I9,2X,' Time = 'E15.8,'  Etot = ',E15.8,'  Ekin  = ',E15.8,'  Utot      = ',&
                E15.8,'  U_vdw     = ',E15.8,'  U_coul   = ',E15.8)                     
 201 FORMAT(' step = ',I9,2X,' Time = 'E15.8,'  Temp = ',E15.8,'  Press = ',E15.8,'  Pvir_vdw  = ',&
                E15.8,'  Pvir_coul = ',E15.8,'  Volume   = ',E15.8,'  Htot = ',E15.8)

 300 FORMAT(/ &
              '  Thermodynamic information '/ &
     &        '  ---------------------------------------------'/ &
     &        '  step                  = ',I9/ &
     &        '  time                  = ',E19.12/ &
     &        '  Ekin                  = ',E19.12/ &
     &        '  Temp                  = ',E19.12/ &
     &        '  Utot                  = ',E19.12/ &
     &        '  U_vdw                 = ',E19.12/ &
     &        '  U_coul                = ',E19.12/ &
     &        '  Pressure              = ',E19.12/ &
     &        '  Pvir_vdw              = ',E19.12/ &
     &        '  Pvir_coul             = ',E19.12/ &
     &        '  Pvir_tot              = ',E19.12/ &
     &        '  volume                = ',E19.12/ &
     &        '  a cell                = ',E19.12/ &
     &        '  b cell                = ',E19.12/ &
     &        '  c cell                = ',E19.12/ &
     &        '  ---------------------------------------------'/ &
     &        '  Etot                  = ',E19.12/ &
     &        '  Htot                  = ',E19.12)

  return

END SUBROUTINE write_thermo

! *********************** SUBROUTINE write_average_thermo **********************
!> \brief
!! write time average thermodynamic quantities to file OSZIFF or print to standard output 
! ******************************************************************************

SUBROUTINE write_average_thermo ( kunit )

  implicit none

  ! global
  integer, intent(in) :: kunit

  !local
  real(kind=dp) :: e_tot_av , e_kin_r_av , u_tot_av , u_lj_r_av
  real(kind=dp) :: u_coul_r_av , temp_r_av  , pressure_tot_av , pressure_lj_av  , pressure_coul_av               

  real(kind=dp) :: e_tot_avsq , e_kin_r_avsq , u_tot_avsq , u_lj_r_avsq
  real(kind=dp) :: u_coul_r_avsq , temp_r_avsq  , pressure_tot_avsq , pressure_lj_avsq  , pressure_coul_avsq

  real(kind=dp) :: e_tot_sig , e_kin_r_sig , u_tot_sig , u_lj_r_sig
  real(kind=dp) :: u_coul_r_sig , temp_r_sig  , pressure_tot_sig , pressure_lj_sig  , pressure_coul_sig
 
  if ( acc_e_tot%counter         .eq. 0 ) return
  if ( acc_e_kin_r%counter       .eq. 0 ) return
  if ( acc_u_tot%counter         .eq. 0 ) return 
  if ( acc_u_lj_r%counter        .eq. 0 ) return
  if ( acc_u_coul_r%counter      .eq. 0 ) return
  if ( acc_temp_r%counter        .eq. 0 ) return
  if ( acc_pressure_tot%counter  .eq. 0 ) return
  if ( acc_pressure_lj%counter   .eq. 0 ) return
  if ( acc_pressure_coul%counter .eq. 0 ) return
 
  ! ========
  !  < A >
  ! ========
  e_tot_av           = acc_e_tot%accval           / REAL ( acc_e_tot%counter , kind = dp )
  e_kin_r_av         = acc_e_kin_r%accval         / REAL ( acc_e_kin_r%counter , kind = dp ) 
  u_tot_av           = acc_u_tot%accval           / REAL ( acc_u_tot%counter , kind = dp )
  u_lj_r_av          = acc_u_lj_r%accval          / REAL ( acc_u_lj_r%counter , kind = dp )
  u_coul_r_av        = acc_u_coul_r%accval        / REAL ( acc_u_coul_r%counter , kind = dp )
  temp_r_av          = acc_temp_r%accval          / REAL ( acc_temp_r%counter , kind = dp )
  pressure_tot_av    = acc_pressure_tot%accval    / REAL ( acc_pressure_tot%counter , kind = dp )
  pressure_lj_av     = acc_pressure_lj%accval     / REAL ( acc_pressure_lj%counter , kind = dp )
  pressure_coul_av   = acc_pressure_coul%accval   / REAL ( acc_pressure_coul%counter , kind = dp )

  ! ========
  !  < A² >
  ! ========
  e_tot_avsq         = acc_e_tot%accvalsq         / REAL ( acc_e_tot%counter , kind = dp )
  e_kin_r_avsq       = acc_e_kin_r%accvalsq       / REAL ( acc_e_kin_r%counter , kind = dp ) 
  u_tot_avsq         = acc_u_tot%accvalsq         / REAL ( acc_u_tot%counter , kind = dp )
  u_lj_r_avsq        = acc_u_lj_r%accvalsq        / REAL ( acc_u_lj_r%counter , kind = dp )
  u_coul_r_avsq      = acc_u_coul_r%accvalsq      / REAL ( acc_u_coul_r%counter , kind = dp )
  temp_r_avsq        = acc_temp_r%accvalsq        / REAL ( acc_temp_r%counter , kind = dp )
  pressure_tot_avsq  = acc_pressure_tot%accvalsq  / REAL ( acc_pressure_tot%counter , kind = dp )
  pressure_lj_avsq   = acc_pressure_lj%accvalsq   / REAL ( acc_pressure_lj%counter , kind = dp )
  pressure_coul_avsq = acc_pressure_coul%accvalsq / REAL ( acc_pressure_coul%counter , kind = dp )

  ! =============================
  !  sqrt ( < A² > - < A > ²)
  ! =============================
  e_tot_sig         = SQRT ( e_tot_avsq         - e_tot_av         * e_tot_av         )
  e_kin_r_sig       = SQRT ( e_kin_r_avsq       - e_kin_r_av       * e_kin_r_av       )
  u_tot_sig         = SQRT ( u_tot_avsq         - u_tot_av         * u_tot_av         )
  u_lj_r_sig        = SQRT ( u_lj_r_avsq        - u_lj_r_av        * u_lj_r_av        )
  u_coul_r_sig      = SQRT ( u_coul_r_avsq      - u_coul_r_av      * u_coul_r_av      )
  temp_r_sig        = SQRT ( temp_r_avsq        - temp_r_av        * temp_r_av        )
  pressure_tot_sig  = SQRT ( pressure_tot_avsq  - pressure_tot_av  * pressure_tot_av  )
  pressure_lj_sig   = SQRT ( pressure_lj_avsq   - pressure_lj_av   * pressure_lj_av   )
  pressure_coul_sig = SQRT ( pressure_coul_avsq - pressure_coul_av * pressure_coul_av )


  if ( ionode ) then
  !   WRITE ( kunit , '(30i6)' )  acc_e_tot%counter , acc_e_kin_r%counter, acc_u_tot%counter , acc_u_lj_r%counter , acc_u_coul_r%counter , acc_temp_r%counter , acc_pressure_tot%counter , acc_pressure_lj%counter , acc_pressure_coul%counter 

      WRITE ( kunit , 100 ) e_tot_av    , e_kin_r_av       , u_tot_av         , u_lj_r_av         , u_coul_r_av        
      WRITE ( kunit , 102 ) e_tot_sig   , e_kin_r_sig      , u_tot_sig        , u_lj_r_sig        , u_coul_r_sig
      WRITE ( kunit , 101 ) temp_r_av   , pressure_tot_av  , pressure_lj_av   , pressure_coul_av  
      WRITE ( kunit , 103 ) temp_r_sig  , pressure_tot_sig , pressure_lj_sig  , pressure_coul_sig  
  endif

 100 FORMAT(2X,'Aver. values:  <Etot>     = ',E15.8,'  <Ekin>      = ',&
                     E15.8,'  <Utot>      = ',E15.8,'  <U_lj>      = ',E15.8,'  <U_coul>       = ',E15.8)
 101 FORMAT(2X,'Aver. values:  <Temp>     = ',E15.8,'  <Press>     = ',&
                     E15.8,'  <P_lj>      = ',E15.8,'  <P_coul>    = ',E15.8)
 102 FORMAT(2X,'Var.  estim :  sig2(Etot) = ',E15.8,'  sig2(Ekin)  = ',&
                     E15.8,'  sig2(Utot)  = ',E15.8,'  sig2(U_lj)  = ',E15.8,'  sig2(U_coul)   = ',E15.8)
 103 FORMAT(2X,'Var.  estim :  sig2(Temp) = ',E15.8,'  sig2(Press) = ',&
                     E15.8,'  sig2(P_lj)  = ',E15.8,'  sig2(P_coul)= ',E15.8)

  return

END SUBROUTINE write_average_thermo

END MODULE thermodynamic
! ===== fmV =====
