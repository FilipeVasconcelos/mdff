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

! *********************** SUBROUTINE init_velocities ***************************
!> \brief 
!!  This routine initialize the velocities. 
! ******************************************************************************
SUBROUTINE init_velocities

  USE constants,        ONLY :  dp 
  USE config,           ONLY :  vx , vy , vz , natm , ntype , ntypemax , atypei , center_of_mass
  USE md,               ONLY :  nequil , setvel , temp
  USE io,               ONLY :  ionode , stdout
  USE control,          ONLY :  full_restart

  implicit none

  ! local
  real(kind=dp) :: T, ekin 
  real(kind=dp) :: com ( 0:ntypemax , 3 )
  integer :: key , ia , it 

  ! =======================================
  !  set key: 
  !   key = 0 if all velocities  are null
  !   key = 1 if at least one is not null
  ! =======================================
  key = 0
  do ia = 1 , natm
    if ( vx (ia) .ne. 0.0_dp ) key = 1
    if ( vy (ia) .ne. 0.0_dp ) key = 1
    if ( vz (ia) .ne. 0.0_dp ) key = 1
  enddo

  ! ===============================================
  !  generate velocities from a given distribution
  ! ===============================================
  if ( (key .eq. 0 .or. .not. full_restart) .and. temp .eq. 0.0_dp ) then

    separator(stdout)    
    io_node blankline(stdout)    
    io_node WRITE ( stdout ,'(a)') 'generate velocities'

    ! ================================
    !  Maxwell-Boltzmann distribution
    ! ================================
    if (setvel .eq. 'MaxwBoltz') CALL maxwellboltzmann_velocities

    ! ============================================
    !  Uniform distribution for test purpose only
    ! ============================================
    if (setvel .eq. 'Uniform')   CALL uniform_random_velocities

    ! ====================
    !  rescale velocities 
    ! ====================
    if ( setvel.ne.'Uniform' )   CALL rescale_velocities(1)

  ! ===========
  ! non null velocities in input from POSFF 
  ! ===========
  elseif ( key .eq. 1 ) then

    ! =======================
    !  input temperature  
    ! =======================
    CALL calc_temp(T, ekin)

    if ( ionode )  WRITE ( stdout ,'(a,f10.4)') 'input temperature                    = ',T

  else 

    ! ============================
    !  no initial kinetic energy
    ! ============================
    vx = 0.0_dp
    vy = 0.0_dp
    vz = 0.0_dp

    io_node  WRITE ( stdout , '(a)' ) 'no initial kinetic energy' 

  endif

  CALL center_of_mass ( vx , vy , vz , com )
  io_node WRITE ( stdout ,'(a,4e16.6)') 'c.o.m. vel. ALL ',com( 0 , :)
  do it = 1 , ntype
    io_node WRITE ( stdout ,'(a,a,a,4e16.6)') 'c.o.m. vel. ', atypei( it ),' ',com( it , :)
  enddo
  io_node blankline(stdout)

  return

END SUBROUTINE init_velocities


! *********************** SUBROUTINE rescale_velocities ************************
!> \brief
!! this subroutine rescale velocities with beredsen thermostat.
!! If tauTberendsen = dt , this becomes a simple rescale procedure
!> \param[in] quite : make the subroutine quite
! ******************************************************************************
SUBROUTINE rescale_velocities (quite)

  USE constants,                ONLY :  dp , boltz_unit
  USE control,                  ONLY :  lcsvr
  USE config,                   ONLY :  natm , vx , ntype, vy , vz, ntypemax, atypei ,center_of_mass
  USE md,                       ONLY :  dt , temp , tauTberendsen, taucsvr, annealing
  USE io,                       ONLY :  ionode , stdout
  USE thermodynamic,            ONLY :  csvr_conint

  implicit none

  ! global
  integer, intent(in) :: quite

  ! local
  integer :: ia , ndeg
  real(kind=dp) :: T, lambda, ekin , sigma , SUMX,SUMY,SUMZ
  real(kind=dp) , external :: resamplekin

  CALL calc_temp(T,ekin)

  lambda = ( 1.0_dp + (dt / tauTberendsen) * (  ( temp / T / boltz_unit ) - 1.0_dp) ) ** 0.5_dp
  
  if ( lcsvr ) then 
    ndeg        = 3 * natm - 3 
    sigma       = ndeg * temp * 0.5d0 
    lambda      = resamplekin(ekin,sigma,ndeg,taucsvr)
    csvr_conint = csvr_conint - (lambda-ekin)
    lambda      = sqrt(lambda/ekin)
#ifdef debug
    write(*,'(a,3e20.8)') 'resamplekin',ekin,sigma,lambda
#endif
  endif

  if ( annealing .ne. 1.0_dp ) lambda = annealing

  do ia = 1 , natm
     vx ( ia ) = vx ( ia ) * lambda
     vy ( ia ) = vy ( ia ) * lambda
     vz ( ia ) = vz ( ia ) * lambda
  enddo

!  if ( lcsvr ) then 
    SUMX = 0.0_dp
    SUMY = 0.0_dp
    SUMZ = 0.0_dp
    do ia = 1 , natm
      SUMX = SUMX + vx ( ia )
      SUMY = SUMY + vy ( ia )
      SUMZ = SUMZ + vz ( ia )
    enddo
    SUMX = SUMX / REAL ( natm , kind = dp )
    SUMY = SUMY / REAL ( natm , kind = dp )
    SUMZ = SUMZ / REAL ( natm , kind = dp )
    do ia = 1 , natm
       vx ( ia ) = vx ( ia ) - SUMX
       vy ( ia ) = vy ( ia ) - SUMY
       vz ( ia ) = vz ( ia ) - SUMZ
    enddo
!  endif

  if ( ionode .and. quite .eq. 1) then
    WRITE ( stdout ,'(a,f10.4)') 'effective temperature      T        = ',T
    WRITE ( stdout ,'(a,f10.4)') 'wanted temperature         T0       = ',temp / boltz_unit
    WRITE ( stdout ,'(a,2e16.8)') 'velocities rescaled by              = ',lambda,(dt / tauTberendsen) * (  ( temp / T / boltz_unit ) - 1.0_dp)
#ifdef debug    
    CALL calc_temp(T,ekin)
    WRITE ( stdout ,'(a,f10.4)') 'debug : rescaled temperature       = ',T
#endif   
    blankline(stdout) 
    
  endif


  return

END SUBROUTINE rescale_velocities

! *********************** SUBROUTINE rescale_volume ************************
!> \brief
!! this subroutine rescale the volume with a beredsen barostat.
!! If tauPberendsen = dt , this becomes a simple rescale procedure
!> \param[in] quite make the subroutine quite
! ******************************************************************************
SUBROUTINE rescale_volume (quite)

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx, ry, rz , simu_cell
  USE md,                       ONLY :  dt , press , tauPberendsen
  USE io,                       ONLY :  ionode , stdout
  USE cell,                     ONLY :  lattice, dirkar, kardir
  USE thermodynamic,            ONLY :  pressure_tot!, calc_thermo

  implicit none

  ! global
  integer, intent(in) :: quite

  ! local
  integer :: ia
  real(kind=dp) :: P, lambda, lambda3

  !CALL calc_thermo()
  P = pressure_tot

  lambda  = ( 1.0_dp - (dt / tauPberendsen) * ( press - P ) ) 
  lambda3 = lambda ** ( 1.0d0 / 3.0d0 )
  !print*,lambda3,P,press

  simu_cell%A(:,1) = simu_cell%A(:,1) * lambda3 ! A
  simu_cell%A(:,2) = simu_cell%A(:,2) * lambda3 ! B
  simu_cell%A(:,3) = simu_cell%A(:,3) * lambda3 ! C

  CALL kardir ( natm , rx , ry , rz , simu_cell%B ) 
  CALL lattice ( simu_cell )
  do ia = 1, natm
    rx(ia) = rx(ia) * lambda3
    ry(ia) = ry(ia) * lambda3
    rz(ia) = rz(ia) * lambda3
  enddo
  CALL dirkar ( natm , rx , ry , rz , simu_cell%A )

  if ( ionode .and. quite .eq. 1) then
    WRITE ( stdout ,'(a,f10.4)') 'Berendsen barostat'
    WRITE ( stdout ,'(a,f10.4)') 'current pressure           P        = ',P
    WRITE ( stdout ,'(a,f10.4)') 'wanted  pressure           P0       = ',press
    WRITE ( stdout ,'(a,f10.4)') 'volume rescaled by                  = ',lambda
    blankline(stdout)
  endif


  return

END SUBROUTINE rescale_volume


! *********************** SUBROUTINE andersen_velocities ***********************
!!> \brief
!! andersen thermostat 
!> \note
!! based on algorithm 15 by Frenkel and Smit
! ******************************************************************************
SUBROUTINE andersen_velocities

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz
  USE md,                       ONLY :  dt , temp , nuandersen 
  USE control,                  ONLY :  dgauss

  implicit none

  ! local
  integer :: ia , iseed
  real(kind=dp) :: T, ekin, sigma, U, G
  
  
  CALL calc_temp(T,ekin)
  sigma = sqrt( temp ) 
 
  do ia = 1, natm
  CALL RANDOM_SEED(SIZE = ISEED)
  CALL RANDOM_NUMBER(HARVEST = U)
     if ( U .lt. nuandersen * dt) then  
       if ( dgauss .eq. 'boxmuller_basic' ) then
         CALL boxmuller_basic(G,0.0_dp,sigma)
         vx ( ia ) = G
         CALL boxmuller_basic(G,0.0_dp,sigma)
         vy ( ia ) = G
         CALL boxmuller_basic(G,0.0_dp,sigma)
         vz ( ia ) = G
       elseif ( dgauss .eq. 'boxmuller_polar' ) then
         CALL boxmuller_polar(G,0.0_dp,sigma)
         vx ( ia ) = G
         CALL boxmuller_polar(G,0.0_dp,sigma)
         vy ( ia ) = G
         CALL boxmuller_polar(G,0.0_dp,sigma)
         vz ( ia ) = G
       elseif ( dgauss .eq. 'knuth' ) then
         CALL knuth(G,0.0_dp,sigma)
         vx ( ia ) = G
         CALL knuth(G,0.0_dp,sigma)
         vy ( ia ) = G
         CALL knuth(G,0.0_dp,sigma)
         vz ( ia ) = G
       endif
     endif
  enddo

  return

END SUBROUTINE andersen_velocities


! *********************** SUBROUTINE uniform_random_velocities *****************
!> \brief
!! uniform random velocities distribution 
!> \author
!! Frenkel and Smit
!> \note
!! only used for test
! ******************************************************************************
SUBROUTINE uniform_random_velocities

  USE constants,        ONLY :  dp 
  USE config,   ONLY :  natm , vx , vy , vz
  USE md,       ONLY :  dt , temp  
  USE control,  ONLY :  dgauss
  USE io,  ONLY :  ionode , stdout
  USE mpimdff

  implicit none

  ! local
  integer :: i
  real(kind=dp) :: G
  integer :: iseed
  DOUBLE PRECISION v2, vx0, vy0, vz0, Vx0t, Vy0t, Vz0t, f

  if ( ionode ) then
    WRITE ( stdout ,'(a)') 'Velocities from uniform distribution' 
    WRITE ( stdout ,'(a)') 'routine from Frenkel Smit Case Study 4'
    WRITE ( stdout ,'(a)') 'only used to test the code'
  endif

  ! ===========================
  !  give particle a velocity
  ! ===========================
  vx0 = 0.0_dp
  vy0 = 0.0_dp
  vz0 = 0.0_dp
  v2 = 0.0_dp
  DO i = 1, natm
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VX(i) = G - 0.5_dp
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VY(i) = G - 0.5_dp
    CALL RANDOM_SEED(SIZE = iseed)
    CALL RANDOM_NUMBER(HARVEST = G)
    VZ(i) = G - 0.5_dp
    vx0 = vx0 + VX(i)
    vy0 = vy0 + VY(i)
    vz0 = vz0 + VZ(i)
    v2 = v2 + VX(i) ** 2 + VY(i) ** 2 + VZ(i) ** 2
  ENDDO
  ! ======================================
  !   set centre of mass movement to zero
  ! ======================================
  vx0 = vx0/natm
  vy0 = vy0/natm
  vz0 = vz0/natm
  Vx0t = 0.0_dp
  Vy0t = 0.0_dp
  Vz0t = 0.0_dp
  f = SQRT ( 3 * REAL ( natm , kind = dp ) * temp / v2 )
  v2 = 0.0_dp
  DO i = 1, natm
    VX(i) = (VX(i)-vx0) * f
    VY(i) = (VY(i)-vy0) * f
    VZ(i) = (VZ(i)-vz0) * f
    Vx0t = Vx0t + VX(i)
    Vy0t = Vy0t + VY(i)
    Vz0t = Vz0t + VZ(i)
    v2 = v2 + VX(i) ** 2 + VY(i) ** 2 + VZ(i) ** 2
  ENDDO
  v2 = v2 / REAL (3 * natm , kind = dp )
  Vx0t = Vx0t/natm
  Vy0t = Vy0t/natm
  Vz0t = Vz0t/natm
  Temp = v2
  io_node WRITE ( stdout , 99001) v2
  io_node WRITE ( stdout , 99002) Vx0t, Vy0t, Vz0t

  return

99001 FORMAT('Initial temperature     : ',f6.3)
99002 FORMAT('Velocity centre of mass : ',/,'          x = ',e9.2,/, '          y = ',e9.2,/, '          z = ',e9.2)

END SUBROUTINE uniform_random_velocities

! *********************** SUBROUTINE maxwellboltzmann_velocities ***************
!
!> \brief
!! maxwell-boltzmann velocities
!! \f$ f(v) = \sqrt{ \frac{m}{2\pi kT}} \exp{-\frac{mv^2}{2kT}}\f$
!
!> \author
!! Allen-Tildsley
!
!> \note
!! adapted from F.24 
!
! ******************************************************************************
SUBROUTINE maxwellboltzmann_velocities

  USE constants,        ONLY :  dp , boltz_unit
  USE config,           ONLY :  natm , vx , vy , vz, massia
  USE md,               ONLY :  dt , temp  
  USE control,          ONLY :  dgauss
  USE io,               ONLY :  ionode , stdout 
  USE mpimdff

  implicit none

  ! local
  integer :: ia
  real(kind=dp) :: RTEMP, SUMX, SUMY, SUMZ
  real(kind=dp) :: G

#ifdef fun
  if ( ionode ) then
    WRITE ( stdout ,'(a)') 'Heat:Hot, as _____:Cold'
    blankline(stdout) 
    WRITE ( stdout ,'(a)') 'a poem by Roald Hoffman'
    WRITE ( stdout ,'(a)') 'from: Chemistry Imagined, Reflections on Science'
    blankline(stdout) 
    WRITE ( stdout ,'(a)') 'Deep in,'
    WRITE ( stdout ,'(a)') "they're there, they're"
    WRITE ( stdout ,'(a)') "at it all the time, it's jai"
    WRITE ( stdout ,'(a)') 'alai on the hot molecular fronton-'
    WRITE ( stdout ,'(a)') 'a bounce off walls onto the packed aleatory'
    WRITE ( stdout ,'(a)') 'dance floor where sideswipes are medium of exchange,'
    WRITE ( stdout ,'(a)') 'momentum trades sealed in swift carom sequences,'
    WRITE ( stdout ,'(a)') 'or just that quick kick in the rear, the haphaz-'
    WRITE ( stdout ,'(a)') 'ard locomotion of the warm, warm world.'
    WRITE ( stdout ,'(a)') 'But spring nights grow cold in Ithaca;'
    WRITE ( stdout ,'(a)') 'the containing walls, glass or metal,'
    WRITE ( stdout ,'(a)') 'are a jagged rough rut of tethered'
    WRITE ( stdout ,'(a)') 'masses, still vibrant, but now'
    WRITE ( stdout ,'(a)') 'retarding, in each collision,'
    WRITE ( stdout ,'(a)') 'the cooling molecules.'
    WRITE ( stdout ,'(a)') "There, they're there,"
    WRITE ( stdout ,'(a)') 'still there,'
    WRITE ( stdout ,'(a)') 'in deep,'
    WRITE ( stdout ,'(a)') 'slow'
    WRITE ( stdout ,'(a)') '.'
    blankline(stdout) 
  endif
#endif

  if ( ionode ) then
    WRITE ( stdout      ,'(a)')     'Velocities from Maxwell-Boltzmann distribution'
    WRITE ( stdout      ,'(a,a20)') 'normal distribution method = ', dgauss  
  endif      


  do ia = 1 , natm
! added 28/01/13
     RTEMP = SQRT ( temp / massia(ia) )
! added 28/01/13
     if ( dgauss .eq. 'boxmuller_basic' ) then
       CALL boxmuller_basic(G,0.0_dp,1.0_dp)
       vx ( ia ) = RTEMP * G
       CALL boxmuller_basic(G,0.0_dp,1.0_dp)
       vy ( ia ) = RTEMP * G
       CALL boxmuller_basic(G,0.0_dp,1.0_dp)
       vz ( ia ) = RTEMP * G
     elseif ( dgauss .eq. 'boxmuller_polar' ) then
       CALL boxmuller_polar(G,0.0_dp,1.0_dp)
       vx ( ia ) = RTEMP * G
       CALL boxmuller_polar(G,0.0_dp,1.0_dp)
       vy ( ia ) = RTEMP * G
       CALL boxmuller_polar(G,0.0_dp,1.0_dp)
       vz ( ia ) = RTEMP * G
     elseif ( dgauss .eq. 'knuth' ) then
       CALL knuth(G,0.0_dp,1.0_dp)
       vx ( ia ) = RTEMP * G 
       CALL knuth(G,0.0_dp,1.0_dp)
       vy ( ia ) = RTEMP * G
       CALL knuth(G,0.0_dp,1.0_dp)
       vz ( ia ) = RTEMP * G
     endif
  enddo

  SUMX = 0.0_dp
  SUMY = 0.0_dp
  SUMZ = 0.0_dp

  do ia = 1 , natm
    SUMX = SUMX + vx ( ia )
    SUMY = SUMY + vy ( ia )
    SUMZ = SUMZ + vz ( ia )
  enddo

  SUMX = SUMX / REAL ( natm ,kind=dp)
  SUMY = SUMY / REAL ( natm ,kind=dp)
  SUMZ = SUMZ / REAL ( natm ,kind=dp)

  do ia  = 1 , natm
     vx ( ia ) = vx ( ia ) - SUMX
     vy ( ia ) = vy ( ia ) - SUMY
     vz ( ia ) = vz ( ia ) - SUMZ
  enddo
  
  return

END SUBROUTINE maxwellboltzmann_velocities

! *********************** SUBROUTINE calc_temp *********************************
!> \brief
!! calculate temperature and kinetic enegy from velocities
!> \param[out] T temperature
!> \param[out] ekin kinetic energy
!> \note
!! units : 
! ******************************************************************************
SUBROUTINE calc_temp (T, ekin)

  USE constants,        ONLY :  dp, boltz_unit
  USE config,           ONLY :  natm , massia, vx , vy , vz 

  implicit none

  ! global
  real(kind=dp), intent(out) :: ekin , T

  ! local
  integer :: ia, L

  L = 3.0_dp * REAL(natm, kind=dp) 

  ekin = 0.0_dp
  do ia = 1 , natm
    ekin =  ekin + ( vx ( ia ) ** 2 + vy ( ia ) ** 2 + vz ( ia ) ** 2 ) * massia(ia) 
  enddo
  ekin = ekin * 0.5_dp
  T = 2.0_dp * ekin / boltz_unit
  T = T / REAL ( L ,kind = dp ) 

  return

END SUBROUTINE calc_temp

! *********************** FUNCTION resamplekin ******************************
! Stochastic velocity rescale, as described in
! Bussi, Donadio and Parrinello, J. Chem. Phys. 126, 014101 (2007)
!
! This subroutine implements Eq.(A7) and returns the new value for the 
! kinetic energy, which can be used to rescale the velocities.
! The procedure can be applied to all atoms or to smaller groups.
! If it is applied to intersecting groups in sequence, the kinetic energy
! that is given as an input (kk) has to be up-to-date with respect to the 
! previous rescalings.
!
! When applied to the entire system, and when performing standard molecular 
! dynamics (fixed c.o.m. (center of mass))
! the degrees of freedom of the c.o.m. have to be discarded in the calculation 
! of ndeg, and the c.o.m. momentum HAS TO BE SET TO ZERO.
! When applied to subgroups, one can chose to:
!       (a) calculate the subgroup kinetic energy in the usual reference frame, 
!           and count the c.o.m. in ndeg
!       (b) calculate the subgroup kinetic energy with respect to its c.o.m. 
!           motion, discard the c.o.m. in ndeg
!           and apply the rescale factor with respect to the subgroup c.o.m. velocity.
! They should be almost equivalent.
! If the subgroups are expected to move one respect to the other, the choice 
! (b) should be better.
!
! If a null relaxation time is required (taut=0.0), the procedure reduces 
! to an istantaneous randomization of the kinetic energy, as described in paragraph IIA.
!
! HOW TO CALCULATE THE EFFECTIVE-ENERGY DRIFT
! The effective-energy (htilde) drift can be used to check the integrator against discretization errors.
! The easiest recipe is:
! htilde = h + conint
! where h is the total energy (kinetic + potential)
! and conint is a quantity accumulated along the trajectory as minus the sum of all the increments of kinetic
! energy due to the thermostat.
! ******************************************************************************
function resamplekin(kk,sigma,ndeg,taut)
  USE constants, ONLY : dp
  implicit none
  real(kind=dp)               :: resamplekin
  real(kind=dp) ,  intent(in)  :: kk    ! present value of the kinetic energy of the atoms to be thermalized (in arbitrary units)
  real(kind=dp) ,  intent(in)  :: sigma ! target average value of the kinetic energy (ndeg k_b T/2)  (in the same units as kk)
  integer, intent(in)  :: ndeg          ! number of degrees of freedom of the atoms to be thermalized
  real(kind=dp) ,  intent(in)  :: taut  ! relaxation time of the thermostat, in units of 'how often this routine is called'
  real(kind=dp)  :: factor,rr,G2
  if(taut>0.1) then
    factor=exp(-1.0/taut)
  else
    factor=0.0
  end if
  CALL boxmuller_polar (G2, 0.0d0, 1.0d0)
  rr = G2 
  resamplekin = kk + ( 1.0 - factor ) * (sigma*(sumnoises(ndeg-1)+rr**2)/ndeg-kk) &
                   +   2.0 * rr * sqrt( kk*sigma/ndeg*(1.0-factor)*factor )

contains

! returns the sum of n independent gaussian noises squared
! (i.e. equivalent to summing the square of the return values of nn calls to gasdev)
double precision function sumnoises(nn)

  USE constants, ONLY : dp
  implicit none
  integer, intent(in) :: nn
  real(kind=dp) :: G1, G2
  
  if(nn==0) then
    sumnoises=0.0
  else if(nn==1) then
    CALL boxmuller_polar (G2, 0.0d0, 1.0d0)
    sumnoises=G2**2
  else if(modulo(nn,2)==0) then
    CALL gammadev(G1,nn/2)
    sumnoises=2.0*G1
  else
    CALL gammadev(G1,(nn-1)/2)
    CALL boxmuller_polar (G2, 0.0d0, 1.0d0)
    sumnoises=2.0*G1 + G2**2
  end if
end function sumnoises

end function resamplekin
! ===== fmV =====
