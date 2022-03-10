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
!#define debug2
!#define debug_nvt_nhc2
!#define debug_nvt_nhcn
!#define debug_npt_nhcnp
! ======= Hardware =======

!================================================================================
!
!> \file
!! all routines related to dynamical integration of phase-space 
!> \brief
!! this subroutines are not clear at all. need more documentation and references.
!! example: leap-frog and verlet are not clearly distinguised
!! to be rigorous I should test all of them and compare to existing codes
! 
!================================================================================


! *********************** SUBROUTINE prop_leap_frog ****************************
!
!> \brief
! leap-frog algorithm  ( or verlet algorithm ) 
!
!> \todo
!! leap-frog and verlet are not clearly distinguished !!
!
! ******************************************************************************
SUBROUTINE prop_leap_frog 

  USE constants,                ONLY :  dp 
  USE md,                       ONLY :  dt
  USE config,                   ONLY :  natm , massia , rx , ry , rz , fx , fy , fz , vx , vy , vz , rxs , rys , rzs
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE engforce_driver,          ONLY :  engforce

  implicit none

  ! local
  integer                                   :: ia 
  real(kind=dp), dimension (:), allocatable :: urx , ury , urz
  real(kind=dp)                             :: idt , dtsq
  real(kind=dp)                             :: tempi , kin
 
  allocate ( urx ( natm ) , ury ( natm ) , urz ( natm ) )
  urx = 0.0_dp
  ury = 0.0_dp
  urz = 0.0_dp

  idt = 0.5_dp / dt 
  dtsq = dt * dt

  ! ==========================
  ! force + potential f(t)
  ! ==========================
  CALL engforce

  ! ================================================= 
  !  r(t+dt) = 2 r(t) - r (t-dt) + f(t)/m dt*dt
  ! ================================================= 
  do ia = 1 , natm
    urx ( ia ) = 2.0_dp * rx ( ia ) - rxs ( ia ) + fx ( ia ) * dtsq / massia(ia)
    ury ( ia ) = 2.0_dp * ry ( ia ) - rys ( ia ) + fy ( ia ) * dtsq / massia(ia)
    urz ( ia ) = 2.0_dp * rz ( ia ) - rzs ( ia ) + fz ( ia ) * dtsq / massia(ia)
  enddo

  ! =====================================
  !  v(t) = ( r(t+dt) - r(t) ) / ( 2 dt) 
  ! =====================================
  do ia = 1, natm
    vx ( ia ) = idt * ( urx ( ia ) - rxs ( ia ) )
    vy ( ia ) = idt * ( ury ( ia ) - rys ( ia ) )
    vz ( ia ) = idt * ( urz ( ia ) - rzs ( ia ) )       
  enddo

  ! ==========================================================
  ! calculate kinetic energy of the full time step velocities
  ! ==========================================================
  CALL calc_temp( tempi , kin ) 
  temp_r = tempi
  e_kin  = kin      

  ! =========================================================
  !  updated positions r(t-dt) <= r(t)  and r(t) <= r (t+dt)
  ! =========================================================
  do ia = 1, natm
    rxs ( ia ) = rx  ( ia )
    rys ( ia ) = ry  ( ia )
    rzs ( ia ) = rz  ( ia )
    rx  ( ia ) = urx ( ia )
    ry  ( ia ) = ury ( ia )
    rz  ( ia ) = urz ( ia )
  enddo

  deallocate ( urx , ury , urz )

  return

END SUBROUTINE prop_leap_frog


! *********************** SUBROUTINE prop_velocity_verlet **********************
!
!> \brief
!! propagation for velocity-verlet algotithm (found a paper)
!
!> \todo
!! test it relatively to prop_leap_frog
!
! ******************************************************************************
SUBROUTINE prop_velocity_verlet 

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm, massia, rx, ry, rz, vx, vy, vz, fx, fy, fz
  USE md,                       ONLY :  dt, integrator
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE io,                       ONLY :  ionode
  USE engforce_driver,          ONLY :  engforce

  implicit none

  ! local
  integer :: ia
  real(kind=dp) :: dtsq2 , dt2
  real(kind=dp), dimension (:), allocatable :: fsx , fsy , fsz
  real(kind=dp) :: tempi , kin

  ! save previous forces
  allocate (fsx(natm),fsy(natm),fsz(natm))
  fsx = 0.0_dp
  fsy = 0.0_dp
  fsz = 0.0_dp
  fsx = fx
  fsy = fy
  fsz = fz

  dtsq2 = dt * dt * 0.5_dp
  dt2 = dt * 0.5_dp

  ! =================================================
  !  r(t+dt) = r(t) + v(t)*dt + f(t)/m * dt*dt/2 
  !  store forces of the previous step  
  ! =================================================
  do ia = 1 , natm
    rx  ( ia ) = rx ( ia ) + vx ( ia ) * dt + (fx ( ia ) * dtsq2 ) / massia(ia)
    ry  ( ia ) = ry ( ia ) + vy ( ia ) * dt + (fy ( ia ) * dtsq2 ) / massia(ia)
    rz  ( ia ) = rz ( ia ) + vz ( ia ) * dt + (fz ( ia ) * dtsq2 ) / massia(ia)
  enddo

  ! ==========================
  ! force + potential f(t+dt)
  ! ==========================
  CALL engforce

  ! ==============================================
  !  v(t+dt) = v(t) + ( f(t-dt) + f(t) ) * dt / 2
  ! ==============================================
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + ( fsx ( ia ) + fx ( ia ) ) * dt2 / massia(ia)
    vy ( ia ) = vy ( ia ) + ( fsy ( ia ) + fy ( ia ) ) * dt2 / massia(ia)
    vz ( ia ) = vz ( ia ) + ( fsz ( ia ) + fz ( ia ) ) * dt2 / massia(ia)
  enddo

  ! full t+dt kinetic energy
  CALL calc_temp(tempi, kin)
  temp_r = tempi      
  e_kin  = kin

  deallocate ( fsx , fsy , fsz )

  return

END SUBROUTINE prop_velocity_verlet

! *********************** SUBROUTINE nose_hoover *******************************
! ******************************************************************************
SUBROUTINE nose_hoover

  USE io,                       ONLY :  stdout,ioprintnode 
  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm,rx,ry,rz,vx,vy,vz,fx,fy,fz,massia
  USE md,                       ONLY :  timesca_thermo, temp,dt,vxi,xi
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt 
  USE engforce_driver,          ONLY :  engforce

  implicit none
  ! local
  integer                :: ia
  real(kind=dp)          :: kin , tempi, L 
  real(kind=dp)          :: invQ, Q
  real(kind=dp)          :: dt2,dt4,s
  
  dt2=dt*0.5_dp
  dt4=dt2*0.5_dp
  ! degrees of freedom
  L = 3.0_dp * ( REAL ( natm , kind=dp) )
  ! thermostat mass
  Q = timesca_thermo**2.0_dp * temp *  L 

  CALL calc_temp ( tempi , kin )
  vxi(1) = vxi(1) + ( 2.0_dp*kin - L * temp ) * dt4 / Q
  s = EXP ( - vxi(1) * dt2 ) 
  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
  enddo
  vxi(1) = vxi(1) + ( 2.0_dp*kin*s*s - L * temp ) * dt4 / Q
  xi(1) = xi(1) + vxi(1) * dt2 
  
  CALL prop_pos_vel_verlet ( kin )
  
  vxi(1) = vxi(1) + ( 2.0_dp*kin - L * temp ) * dt4 / Q
  s = EXP ( - vxi(1) * dt2 ) 
  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
  enddo
  vxi(1) = vxi(1) + ( 2.0_dp*kin*s*s - L * temp ) * dt4 / Q
  xi(1) = xi(1) + vxi(1) * dt2 
  
  CALL calc_temp ( tempi , kin )
  e_kin = kin
  temp_r = tempi

  e_nvt = 0.0_dp
  e_nvt = e_nvt + L * temp * xi(1)
  e_nvt = e_nvt + vxi(1) * vxi(1) * 0.5_dp * Q
  io_printnode blankline(stdout) 
  io_printnode write(stdout,'(a,5f16.8)') '  (extended system contribution) e_nvt',e_nvt,xi(1),vxi(1)
 
  return 

END SUBROUTINE

! *********************** SUBROUTINE nose_hoover_chain2 ************************
!
!> \brief
!! Nose-Hoover two chains
!
!> \note
!! adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE nose_hoover_chain2 

  USE io,                       ONLY :  ioprintnode, ioprint, stdout
  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz , center_of_mass, ntypemax
  USE md,                       ONLY :  dt, vxi, xi , timesca_thermo, temp
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt , e_tot,u_lj_r,h_tot

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi, L 
  real(kind=dp)          :: Q(2)

  ! degrees of freedom
  L = 3.0_dp * ( REAL ( natm , kind=dp) )

  ! thermostat mass
  Q(1) = timesca_thermo**2.0_dp * temp * L 
  Q(2) = timesca_thermo**2.0_dp * temp
 
  CALL calc_temp ( tempi , kin )
  CALL chain_nhc2 ( kin , vxi , xi , Q , L )

  CALL prop_pos_vel_verlet ( kin )

  CALL chain_nhc2( kin, vxi, xi , Q , L ) 

  CALL calc_temp ( tempi , kin )
  e_kin = kin
  temp_r = tempi

  e_nvt = 0.0_dp
  e_nvt = e_nvt + L * temp * xi(1)
  do inhc=1,2
    e_nvt = e_nvt + vxi(inhc) * vxi(inhc) * 0.5_dp / Q(inhc)
  enddo
  e_nvt = e_nvt + temp * xi(2)
  io_printnode blankline(stdout) 
  io_printnode write(stdout,'(a,5f16.8)') '  (extended system contribution) e_nvt',e_nvt
 
  return

END SUBROUTINE nose_hoover_chain2


! *********************** SUBROUTINE chain_nhc2 ********************************
!
!> \brief
!! intermediate routine used by nose_hoover_chain2
!
!> \note
!! adapted from Frenkel and Smit
!
! ******************************************************************************
SUBROUTINE chain_nhc2 ( kin, vxi, xi , Q , L )

  USE io,                       ONLY :  ioprint
  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz 
  USE md,                       ONLY :  temp , dt , nhc_n

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin
  real(kind=dp), intent (inout) :: vxi(2), xi(2) 
  real(kind=dp), intent (in)    :: L, Q(2)

  ! local
  integer       :: ia
  real(kind=dp) :: G1, G2  
  real(kind=dp) :: dt2, dt4, dt8
  real(kind=dp) :: s
 
  ! some constants related integrator step dt
  s   = 1.0_dp
  dt2 = dt  * 0.5_dp
  dt4 = dt2 * 0.5_dp
  dt8 = dt4 * 0.5_dp

  G2     = ( Q(1) * vxi(1) * vxi(1) - temp ) / Q(2)
#ifdef debug_nvt_nhc2
  write(*,'(a,e16.8)') "G2",G2
#endif
  vxi(2) = vxi(2) + G2 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G1     = ( 2.0_dp * kin - L * temp ) / Q(1)
#ifdef debug_nvt_nhc2
  write(*,'(a,3e16.8)') "G1",G1,Q(1),L
#endif
  vxi(1) = vxi(1) + G1 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
#ifdef debug_nvt_nhc2
    io_print write(stdout,'(a,2e16.8)') "vxi(1),vxi(2)",vxi
    io_print write(stdout,'(a,2e16.8)') "G1,G2",G1,G2
#endif
  s   = s * EXP ( - vxi(1) * dt2 )
  !kin = s * s * kin
  xi     = xi + vxi * dt2
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G1     = ( 2.0_dp * kin * s * s - L * temp) / Q(1)
  vxi(1) = vxi(1) + G1 * dt4
  vxi(1) = vxi(1) * EXP ( - vxi(2) * dt8 )
  G2     = ( Q(1) * vxi(1) * vxi(1) - temp ) / Q(2)
  vxi(2) = vxi(2) + G2 * dt4
#ifdef debug_nvt_nhc2
    io_print write(stdout,'(a,2e16.8)') "xi(1),xi(2)",xi
    io_print write(stdout,'(a,2e16.8)') "G1,G2",G1,G2
#endif

  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
  enddo

  return
 
END SUBROUTINE chain_nhc2


! *********************** SUBROUTINE prop_pos_vel_verlet ***********************
!
!> \brief
!! propagates position and position in the velet algorithm
!
! ******************************************************************************
SUBROUTINE prop_pos_vel_verlet ( kin )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , massia , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz 
  USE md,                       ONLY :  dt
  USE engforce_driver,          ONLY :  engforce

  implicit none

  ! global
  real(kind=dp), intent (out) :: kin

  ! local
  integer :: ia
  real(kind=dp) :: dt2

  dt2 = dt * 0.5_dp

  ! =========================================
  !  v(t+dt2) = v(t) + f(t) * dt2 / m
  !  r(t+dt)  = r(t) + v(t+dt2) * dt  
  ! note : dt2 = dt / 2
  ! =========================================
  do ia = 1 , natm  
    vx ( ia ) = vx ( ia ) + fx ( ia ) * dt2 / massia(ia) 
    vy ( ia ) = vy ( ia ) + fy ( ia ) * dt2 / massia(ia)
    vz ( ia ) = vz ( ia ) + fz ( ia ) * dt2 / massia(ia)
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt 
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt  
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt  
  enddo

  ! ==========================
  ! f(t+dt)
  ! ==========================
  CALL engforce
  ! =========================================
  !   v(t) = v(t+dt2) + f(t+dt) * dt2 / m
  !   v(t) = v(t) + dt2 * ( f(t) + f(t+ft) ) 
  ! =========================================
  kin  = 0.0_dp
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + fx ( ia ) * dt2 / massia(ia)
    vy ( ia ) = vy ( ia ) + fy ( ia ) * dt2 / massia(ia)
    vz ( ia ) = vz ( ia ) + fz ( ia ) * dt2 / massia(ia)
    kin = kin + ( vx ( ia ) * vx ( ia ) +  vy ( ia ) * vy ( ia ) + vz ( ia ) * vz ( ia ) ) * massia (ia)
  enddo      
  kin = kin * 0.5_dp  

  return

END SUBROUTINE prop_pos_vel_verlet

! *********************** SUBROUTINE beeman ************************************
!
!> \brief
!! D. Beeman "Some multistep methods for use in molecular dynamics calculations",
!
!> \note
!! Journal of Computational Physics 20 pp. 130-139 (1976)
!
! ******************************************************************************
SUBROUTINE beeman 

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , massia , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz , fxs , fys , fzs 
  USE thermodynamic,            ONLY :  temp_r , e_kin
  USE md,                       ONLY :  dt
  USE engforce_driver,          ONLY :  engforce

  implicit none

  ! local
  integer :: ia
  real(kind=dp) :: onesix , twothree , onethree , fivesix , dtsq
  real(kind=dp) , dimension (:) , allocatable :: fxtmp , fytmp, fztmp
  real(kind=dp) :: kin , tempi


  allocate ( fxtmp (natm) , fytmp (natm) , fztmp (natm) ) 

  ! ================
  !  some constants
  ! ================
  onesix = 1.0_dp / 6.0_dp
  onethree = 2.0_dp * onesix
  twothree = 4.0_dp * onesix
  fivesix  = 5.0_dp * onesix
  dtsq = dt * dt
  
  ! ==================================================================
  ! r (t+dt) = r (t) + v (t) dt  + ( 2/3 f (t) - 1/6 f (t-dt) ) / m dt^2
  ! ==================================================================
  do ia = 1 , natm
    rx ( ia ) = rx ( ia ) + vx ( ia ) * dt + ( twothree * fx ( ia ) - onesix * fxs ( ia ) ) * dtsq / massia(ia)
    ry ( ia ) = ry ( ia ) + vy ( ia ) * dt + ( twothree * fy ( ia ) - onesix * fys ( ia ) ) * dtsq / massia(ia)
    rz ( ia ) = rz ( ia ) + vz ( ia ) * dt + ( twothree * fz ( ia ) - onesix * fzs ( ia ) ) * dtsq / massia(ia)
  enddo

  ! ========================
  ! save current force f(t)
  ! ========================
  fxtmp = fx
  fytmp = fy
  fztmp = fz

  ! ===========================
  ! force + potential  f(t+dt)
  ! ===========================
  CALL engforce

  ! ==================================================================
  ! v (t+dt) = v (t) + ( 1/3 f(t+dt) + 5/6 f(t) - 1/6  f(t-dt) )  dt
  ! ==================================================================
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) + ( onethree * fx ( ia ) + fivesix * fxtmp ( ia ) - onesix * fxs ( ia ) ) * dt / massia(ia)    
    vy ( ia ) = vy ( ia ) + ( onethree * fy ( ia ) + fivesix * fytmp ( ia ) - onesix * fys ( ia ) ) * dt / massia(ia)
    vz ( ia ) = vz ( ia ) + ( onethree * fz ( ia ) + fivesix * fztmp ( ia ) - onesix * fzs ( ia ) ) * dt / massia(ia)
  enddo

  CALL calc_temp(tempi, kin)
  temp_r = tempi      
  e_kin  = kin

  ! ==============
  ! store f(t-dt)
  ! ==============
  fxs = fxtmp
  fys = fytmp
  fzs = fztmp

  deallocate ( fxtmp , fytmp , fztmp ) 

  return

END SUBROUTINE beeman

! *********************** SUBROUTINE nose_hoover_chain_n ***********************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE nose_hoover_chain_n

  USE io,                       ONLY :  ioprintnode, stdout
  USE constants,                ONLY :  dp  
  USE config,                   ONLY :  natm , rx , ry , rz , vx , vy , vz , fx , fy , fz , center_of_mass, ntypemax
  USE md,                       ONLY :  dt, vxi, xi , timesca_thermo, nhc_n,temp
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_nvt , e_tot,u_lj_r,h_tot
  USE time,                     ONLY :  integratimetot
  USE mpimdff

  implicit none

  ! local
  integer                :: inhc , ierr
  real(kind=dp)          :: kin , tempi, L 
  real(kind=dp), dimension(:), allocatable :: Q
  real(kind=dp)          :: ttt1
  real(kind=dp) :: nvt1, nvt2, nvt3, nvt4

#ifdef MPI
  ttt1 = MPI_WTIME(ierr) 
#endif
  ! degrees of freedom
  L = 3.0_dp * ( REAL ( natm , kind=dp) )

  ! thermostat mass
  allocate( Q (nhc_n) )
  Q = timesca_thermo**2.0_dp * temp 
  Q(1) = Q(1) * L 

  CALL calc_temp ( tempi , kin )
  CALL chain_nhcn ( kin , vxi , xi , Q , L )

  CALL prop_pos_vel_verlet ( kin )

  CALL chain_nhcn( kin, vxi, xi , Q , L )

  ! full t + dt kinetic energy 
  CALL calc_temp ( tempi , kin )
  e_kin = kin
  temp_r = tempi

  ! ===========================================
  !  conserved quantity energy in NVT ensemble
  ! ===========================================
  ! note on units :
  ! [ vxi ] = [ eV ] [ T ]
  ! [ Q ]   = [ eV ] [ T ]**2 
  ! vxi * vxi / Q =  [ eV ]
  ! [ temp ] * [xi] = [ eV ]  ( [xi] = sans unité )
  e_nvt = 0.0_dp
  nvt1 = L * temp * xi(1)
  nvt2 = vxi(1) * vxi(1) * 0.5_dp / Q(1)
  nvt3 = 0.0_dp
  nvt4 = 0.0_dp
  do inhc = 2 , nhc_n
    nvt3 = nvt3 + vxi(inhc) * vxi(inhc) * 0.5_dp / Q(inhc)
    nvt4 = nvt4 + xi(inhc)    
  enddo
  e_nvt = nvt1 + nvt2 + nvt3 + temp * nvt4
  io_printnode blankline(stdout) 
  io_printnode write(stdout,'(a,5f16.8)') '  e_nvt',e_nvt,nvt1,nvt2,nvt3,nvt4

  deallocate( Q ) 
#ifdef MPI
  integratimetot = integratimetot + ( MPI_WTIME(ierr) - ttt1 )
#endif

  return

END SUBROUTINE nose_hoover_chain_n 

! *********************** SUBROUTINE chain_nhcn ***********************
!
!> \brief
! ref : 
! [1] Molecular Physics, (1996), v87, n5, p1117 Martyna and al.
! [2] Phys Lett. A, (1190), v150 n5,6,7, p262, Yoshida
! [3] https://files.nyu.edu/mt33/public/abstracts/a6_19_s18.pdf
!
! ******************************************************************************
SUBROUTINE chain_nhcn ( kin , vxi , xi , Q , L )

  USE constants,                ONLY :  dp 
  USE config,                   ONLY :  natm , vx , vy , vz, massia 
  USE md,                       ONLY :  temp , dt , nhc_n , nhc_yosh_order,nhc_mults
  USE io,                       ONLY :  ionode, stdout , ioprint , ioprintnode

  implicit none

  ! global
  real(kind=dp), intent (inout) :: kin
  real(kind=dp), intent (inout) :: vxi(nhc_n), xi(nhc_n) , Q(nhc_n)

  ! local
  integer :: ia , j , k, inh
  real(kind=dp), dimension ( : ) , allocatable :: G
  real(kind=dp) :: s , dts , dts2 , dts4 , dts8 , L

  real(kind=dp) , dimension ( : ), allocatable :: yosh_w  ! integrator order as in [2] YOSHIDA 
  real(kind=dp) , dimension ( : ), allocatable :: dt_yosh ! yoshida time 

  allocate ( yosh_w  ( nhc_yosh_order) )       
  allocate ( dt_yosh ( nhc_yosh_order) )       
       
  CALL get_yoshida_weigth(nhc_yosh_order,yosh_w)

  allocate ( G ( nhc_n) )
  G = 0.0_dp
  dt_yosh =  yosh_w * dt

  s = 1.0_dp ! scale
  G(1) =  2.0_dp*kin - L * temp 

  msloop : do k=1,nhc_mults

    yoshloop:  do j=1, nhc_yosh_order 

    dts = dt_yosh( j ) / REAL(nhc_mults,kind=dp)
    dts2 = dts  * 0.5d0
    dts4 = dts2 * 0.5d0
    dts8 = dts4 * 0.5d0 

    G(nhc_n) = ( vxi(nhc_n-1) * vxi(nhc_n-1) / Q(nhc_n-1) - temp) 
    ! exp1 
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4 
    do inh=nhc_n-1,1,-1
      !exp2 : scale thermo momentum
      vxi ( inh )= vxi ( inh ) * EXP ( - vxi ( inh + 1 ) * dts8 / Q ( inh+1 ) )
      !exp3 : propagate thermo momentum
      vxi ( inh )= vxi ( inh ) + G ( inh ) * dts4 
      !exp4 : scale thermo momentum
      vxi ( inh )= vxi ( inh ) * EXP ( - vxi ( inh + 1 ) * dts8 / Q ( inh+1 ) )
    enddo
    ! exp5: scale velocities
    s = s * EXP ( - vxi(1) * dts2 / Q ( 1 ) ) 

    ! exp6 : propagating xi 
    xi = xi + vxi * dts2 / Q !!! minus in [3] seems to be wrong ???? 
    
    G(1) = 2.0_dp*kin*s*s - L * temp 

    ! same as below
    do inh = 1 , nhc_n - 1
      vxi ( inh )     =   vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 / Q ( inh + 1 ) )  
      vxi ( inh )     =   vxi ( inh ) + G(inh) * dts4
      vxi ( inh )     =   vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 / Q ( inh + 1 ) )  
      G   ( inh + 1 ) = ( vxi ( inh ) * vxi(inh) / Q(inh) - temp)
    enddo
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4

  enddo yoshloop

enddo msloop 

  kin = 0.0_dp
  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
    kin =  kin + ( vx ( ia ) ** 2 + vy ( ia ) ** 2 + vz ( ia ) ** 2 ) * massia(ia)
  enddo
  kin = kin * 0.5_dp

  deallocate ( yosh_w  )       
  deallocate ( dt_yosh )       
  deallocate ( G )

  return

END SUBROUTINE chain_nhcn

! *********************** SUBROUTINE chain_nhcnp ***********************
!
!> \brief
! ref : 
! [1] Molecular Physics, (1996), v87, n5, p1117 Martyna and al.
! [2] Phys Lett. A, (1190), v150 n5,6,7, p262, Yoshida
!
! ******************************************************************************
SUBROUTINE chain_nhcnp ( kin , vxi , xi , vxib , xib , ve , Q , Qb , W , L , trotter_order )

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  simu_cell, massia , natm , vx , vy , vz
  USE md,                       ONLY :  temp , press, dt , nhc_n , nhc_yosh_order, nhc_mults
  USE thermodynamic,            ONLY :  pvirial_tot
  USE io,                       ONLY :  ionode , stdout, stderr, ioprint, ioprintnode
  USE mpimdff,                  ONLY :  myrank

  implicit none

  ! global
  integer :: trotter_order
  real(kind=dp), intent (inout) :: kin , W , L 
  real(kind=dp), intent (inout) :: vxi(nhc_n), xi(nhc_n) , Q(nhc_n) , vxib(nhc_n), xib(nhc_n) , Qb(nhc_n) , ve 

  ! local
  integer :: ia , j , inh, k 
  real(kind=dp), dimension ( : ) , allocatable :: G , Gb

  real(kind=dp) , dimension ( : ), allocatable :: yosh_w  ! integrator order as in [2] YOSHIDA 
  real(kind=dp) , dimension ( : ), allocatable :: dt_yosh ! yoshida time 
  real(kind=dp) :: s , sb , dt2 , dts , dts2 , dts4 , dts8 , Ge , odnf , P_kin

  !print*,'in'
  allocate ( yosh_w  ( nhc_yosh_order) )
  allocate ( dt_yosh ( nhc_yosh_order) )

  CALL get_yoshida_weigth(nhc_yosh_order,yosh_w)

  odnf = 1.0_dp + 3.0_dp / L

  allocate ( G ( nhc_n) , Gb ( nhc_n) )
  G  = 0.0_dp
  Gb = 0.0_dp
  dt_yosh =  yosh_w * dt
  dt2 = dt * 0.5_dp

  s     = 1.0_dp ! scale particule velocities
  sb    = 1.0_dp ! scale "piston" velocity
  Ge    = odnf * kin + 3.0_dp * simu_cell%omega * ( pvirial_tot - press ) 
  ve = ve + Ge * dt2 ! pe 
  !G(1)  = 2.0_dp*kin - L * temp
  Gb(1) = 0.5_dp * ve * ve / W - temp
  ! barostat
#ifdef debug
  write(stderr,'(a,i,<nhc_n>e60.48,5f60.48)') 'debug : ',myrank,G,kin,L,temp,2.0_dp*kin,L * temp
#endif  


  msloop : do k=1,nhc_mults

    yoshloop:  do j=1, nhc_yosh_order

    dts = dt_yosh( j ) / REAL(nhc_mults,kind=dp)
    dts2 = dts  * 0.5_dp
    dts4 = dts2 * 0.5_dp
    dts8 = dts4 * 0.5_dp

! test order of trotter expansion
if ( trotter_order == 1 ) then
    ! ==============
    !  thermo-baro
    ! ==============
    Gb(nhc_n) = ( 0.5_dp * vxib(nhc_n-1) * vxib(nhc_n-1) / Qb(nhc_n-1) - temp)
    vxib ( nhc_n ) = vxib ( nhc_n ) + Gb ( nhc_n ) * dts4
    do inh=nhc_n-1,1,-1
      vxib ( inh ) = vxib ( inh ) * EXP ( - vxib ( inh + 1 ) * dts8 / Qb ( inh+1 ) )
      vxib ( inh ) = vxib ( inh ) + Gb ( inh ) * dts4
      vxib ( inh ) = vxib ( inh ) * EXP ( - vxib ( inh + 1 ) * dts8 / Qb ( inh+1 ) )
    enddo
    sb = sb * EXP ( - vxib(1) * dts2 / Qb ( 1 ) )
    do inh=1,nhc_n
      xib(inh) = xib(inh) + vxib(inh) * dts2 / Qb(inh) 
    enddo
    Gb(1) = 0.5_dp * ve * ve * sb * sb / W - temp
    do inh = 1 , nhc_n - 1
      vxib ( inh )     =   vxib ( inh ) * EXP ( - vxib ( inh + 1) * dts8 / Qb ( inh + 1 ) )
      vxib ( inh )     =   vxib ( inh ) + Gb(inh) * dts4
      vxib ( inh )     =   vxib ( inh ) * EXP ( - vxib ( inh + 1) * dts8 / Qb ( inh + 1 ) )
      Gb   ( inh + 1 ) = ( 0.5_dp * vxib ( inh ) * vxib(inh) / Qb(inh) - temp)
    enddo
    vxib ( nhc_n ) = vxib ( nhc_n ) + Gb ( nhc_n ) * dts4
#ifdef debug2
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp i G    ',(G(inh)   , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp i xi   ',(xi(inh)  , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp i vxi  ',(vxi(inh) , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp i Gb   ',(Gb(inh)  , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp i xib  ',(xib(inh) , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp i vxib ',(vxib(inh), inh=1,nhc_n)
  write(stdout,'(2i,a,e60.48)')        trotter_order,myrank,'mp i ve   ',ve
#endif
    ! ===================
    !  thermo-particules
    ! ===================
    G(nhc_n) = ( 0.5_dp * vxi(nhc_n-1) * vxi(nhc_n-1) / Q(nhc_n-1) - temp)
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4
    do inh=nhc_n-1,1,-1
      vxi ( inh ) = vxi ( inh ) * EXP ( - vxi ( inh + 1 ) * dts8 / Q ( inh+1 ) )
      vxi ( inh ) = vxi ( inh ) + G ( inh ) * dts4
      vxi ( inh ) = vxi ( inh ) * EXP ( - vxi ( inh + 1 ) * dts8 / Q ( inh+1 ) )
    enddo
    s = s * EXP ( - vxi(1) * dts2 / Q ( 1 ) )
    do inh=1,nhc_n
      xi(inh) = xi(inh) + vxi(inh) * dts2 / Q(inh) 
    enddo
    G(1) = 2.0_dp*kin *s*s- L * temp
    do inh = 1 , nhc_n - 1
      vxi ( inh )     =   vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 / Q ( inh + 1 ) )
      vxi ( inh )     =   vxi ( inh ) + G(inh) * dts4
      vxi ( inh )     =   vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 / Q ( inh + 1 ) )
      G   ( inh + 1 ) = ( 0.5_dp * vxi ( inh ) * vxi(inh) / Q(inh) - temp)
    enddo
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4

! trotter_order == 2
else

    ! ===================
    !  thermo-particules
    ! ===================
    G(nhc_n) = ( 0.5_dp * vxi(nhc_n-1) * vxi(nhc_n-1) / Q(nhc_n-1) - temp)
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4
    do inh=nhc_n-1,1,-1
      vxi ( inh ) = vxi ( inh ) * EXP ( - vxi ( inh + 1 ) * dts8 / Q ( inh+1 ) )
      vxi ( inh ) = vxi ( inh ) + G ( inh ) * dts4
      vxi ( inh ) = vxi ( inh ) * EXP ( - vxi ( inh + 1 ) * dts8 / Q ( inh+1 ) )
    enddo
    s = s * EXP ( - vxi(1) * dts2 / Q ( 1 ) )
    do inh=1,nhc_n
      xi(inh) = xi(inh) + vxi(inh) * dts2 / Q(inh) 
    enddo
    G(1) = 2.0_dp*kin*s*s - L * temp
    do inh = 1 , nhc_n - 1
      vxi ( inh )     =   vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 / Q ( inh + 1 ) )
      vxi ( inh )     =   vxi ( inh ) + G(inh) * dts4
      vxi ( inh )     =   vxi ( inh ) * EXP ( - vxi ( inh + 1) * dts8 / Q ( inh + 1 ) )
      G   ( inh + 1 ) = ( 0.5_dp * vxi ( inh ) * vxi(inh) / Q(inh) - temp)
    enddo
    vxi ( nhc_n ) = vxi ( nhc_n ) + G ( nhc_n ) * dts4
    ! ==============
    !  thermo-baro
    ! ==============
    Gb(nhc_n) = ( 0.5_dp * vxib(nhc_n-1) * vxib(nhc_n-1) / Qb(nhc_n-1) - temp)
    vxib ( nhc_n ) = vxib ( nhc_n ) + Gb ( nhc_n ) * dts4
    do inh=nhc_n-1,1,-1
      vxib ( inh ) = vxib ( inh ) * EXP ( - vxib ( inh + 1 ) * dts8 / Qb ( inh+1 ) )
      vxib ( inh ) = vxib ( inh ) + Gb ( inh ) * dts4
      vxib ( inh ) = vxib ( inh ) * EXP ( - vxib ( inh + 1 ) * dts8 / Qb ( inh+1 ) )
    enddo
    sb = sb * EXP ( - vxib(1) * dts2 / Qb ( 1 ) )
    do inh=1,nhc_n
      xib(inh) = xib(inh) + vxib(inh) * dts2 / Qb(inh) 
    enddo
     Gb(1) = 0.5_dp * ve * ve * sb * sb / W - temp
    do inh = 1 , nhc_n - 1
      vxib ( inh )     =   vxib ( inh ) * EXP ( - vxib ( inh + 1) * dts8 / Qb ( inh + 1 ) )
      vxib ( inh )     =   vxib ( inh ) + Gb(inh) * dts4
      vxib ( inh )     =   vxib ( inh ) * EXP ( - vxib ( inh + 1) * dts8 / Qb ( inh + 1 ) )
      Gb   ( inh + 1 ) = ( 0.5_dp * vxib ( inh ) * vxib(inh) / Qb(inh) - temp)
    enddo
    vxib ( nhc_n ) = vxib ( nhc_n ) + Gb ( nhc_n ) * dts4

endif
    ! 
  enddo yoshloop

enddo msloop

#ifdef debug
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp e xi   ',(xi(inh)  , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp e vxi  ',(vxi(inh) , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp e xib  ',(xib(inh) , inh=1,nhc_n)
  write(stdout,'(2i,a,<nhc_n>e60.48)') trotter_order,myrank,'mp e vxib ',(vxib(inh), inh=1,nhc_n)
  write(stdout,'(2i,a,e60.48)')        trotter_order,myrank,'mp e ve   ',ve
  write(stdout,'(2i,a,<nhc_n>e60.48)')        trotter_order,myrank,'mp e Q   ',(Q(inh)    ,inh=1,nhc_n)
#endif


  ! thermostat
  kin = 0.0_dp
  do ia = 1, natm
    vx ( ia ) = s * vx ( ia )
    vy ( ia ) = s * vy ( ia )
    vz ( ia ) = s * vz ( ia )
    kin =  kin + ( vx ( ia ) ** 2 + vy ( ia ) ** 2 + vz ( ia ) ** 2 ) * massia(ia)
  enddo
  kin = kin * 0.5_dp

!  ve = ve + Ge * dt2 ! pe 
! barostat
  ve = ve * sb 
  if ( trotter_order == 1 ) then
    !P_kin = 2.0_dp * odnf * kin  / ( 3.0_dp * simu_cell%omega )  
    P_kin = odnf * kin  / ( 3.0_dp * simu_cell%omega )  
    Ge    = 3.0_dp * simu_cell%omega * ( P_kin + pvirial_tot - press ) 
  !  Ge    = odnf * kin + 3.0_dp * simu_cell%omega * ( pvirial_tot - press ) 
    ve = ve + Ge * dt2 ! pe 
#ifdef debug_npt_nhcnp
  io_printnode write(stdout,'(i,a,7e16.8)') trotter_order,'ve out ',s,sb,ve/W,Ge,pvirial_tot,press
#endif
  endif


#ifdef debug_nvt_nhcnp
  io_printnode write(stdout,'(a,5e16.8)') 've sb',ve/W,sb,Ge,pvirial_tot,press 
#endif

  deallocate ( G , Gb )
  deallocate ( yosh_w  )
  deallocate ( dt_yosh )

  return

END SUBROUTINE chain_nhcnp

! *********************** SUBROUTINE nose_hoover_chain_n_p ***********************
!
!> \brief
!
! ******************************************************************************
SUBROUTINE nose_hoover_chain_n_p

  USE io,                       ONLY :  stdout, ioprintnode
  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , simu_cell, rx , ry , rz , vx , vy , vz , fx , fy , fz
  USE md,                       ONLY :  dt, vxi, xi , vxib, xib , xe , ve , xe0, timesca_thermo , timesca_baro, nhc_n,temp , press
  USE thermodynamic,            ONLY :  temp_r , e_kin , e_npt, pvirial_tot , h_tot , pressure_tot , calc_thermo
  USE mpimdff,                  ONLY :  myrank

  implicit none

  ! local
  integer                :: inhc
  real(kind=dp)          :: kin , tempi , W, L 
  real(kind=dp), dimension ( : ) , allocatable :: Q, Qb
  character(len=20) :: FMT

  L = 3.0_dp * REAL(natm, kind=dp) 
  allocate ( Q(nhc_n) , Qb (nhc_n) )
  ! thermostat/particules mass coupled to velocities (vx...)
  Q    = timesca_thermo**2.0_dp * temp
  Q(1) = L * Q(1) 
  ! thermostat/barostat mass coupled to ve
  Qb   = timesca_baro**2.0_dp * temp
  Qb(1)= Qb(1) * 9.0_dp
  ! barostat "mass"
  W    = (L + 3.0_dp) * timesca_baro**2.0_dp * temp

  CALL calc_temp ( tempi , kin )
#ifdef debug2
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 1 xi   ',(xi(inhc)  , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 1 vxi  ',(vxi(inhc) , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 1 xib  ',(xib(inhc) , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 1 vxib ',(vxib(inhc), inhc=1,nhc_n)
  write(stdout,'(i,a,e60.48)')        myrank,'mp 1 xe   ',xe
  write(stdout,'(i,a,e60.48)')        myrank,'mp 1 xe0  ',xe0
  write(stdout,'(i,a,e60.48)')        myrank,'mp 1 ve   ',ve
#endif

  CALL chain_nhcnp( kin , vxi , xi , vxib , xib , ve , Q , Qb , W , L , 1 )
#ifdef debug2
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 2 xi   ',(xi(inhc)  , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 2 vxi  ',(vxi(inhc) , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 2 xib  ',(xib(inhc) , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp 2 vxib ',(vxib(inhc), inhc=1,nhc_n)
  write(stdout,'(i,a,e60.48)')        myrank,'mp 2 xe   ',xe
  write(stdout,'(i,a,e60.48)')        myrank,'mp 2 xe0  ',xe0
  write(stdout,'(i,a,e60.48)')        myrank,'mp 2 ve   ',ve
#endif
  CALL prop_pos_vel_verlet_npt ( kin , xe , ve , xe0 , L , W ) 
  CALL chain_nhcnp( kin , vxi , xi , vxib , xib , ve , Q , Qb , W , L , 1 )
  CALL calc_temp ( tempi , kin )
  e_kin = kin
  temp_r = tempi

#ifdef debug2
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp xi   ',(xi(inhc)  , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp vxi  ',(vxi(inhc) , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp xib  ',(xib(inhc) , inhc=1,nhc_n)
  write(stdout,'(i,a,<nhc_n>e60.48)') myrank,'mp vxib ',(vxib(inhc), inhc=1,nhc_n)
  write(stdout,'(i,a,e60.48)')        myrank,'mp xe   ',xe
  write(stdout,'(i,a,e60.48)')        myrank,'mp xe0  ',xe0
  write(stdout,'(i,a,e60.48)')        myrank,'mp ve   ',ve
#endif

  ! ==============================================
  !  conserved quantity of NPT ensemble
  ! ==============================================
  ! note on units :
  ! [PV] = [A]**3 * [eV/A**3] = [ eV ]
  ! [ ve ] = [ eV ] [ T ] 
  ! [ W ]   = [ eV ] [ T ]**2 
  ! ve * ve / W = [ eV ]
  ! [ xe  ] = sans unité 
  ! [ xi  ] = sans unité
  ! [ xib ] = sans unité
  ! [ vxi ] = [ eV ] [ T ]
  ! [ vxib ] = [ eV ] [ T ]
  e_npt = 0.0_dp
  e_npt = e_npt + press    * simu_cell%omega                     ! PV              # barostat potential energy of barostat
  e_npt = e_npt + ve * ve  * 0.5_dp / W                          ! pe^2 / 2W       # barostat kinetic energy
  e_npt = e_npt + L * temp * xi(1)                               ! Nf kB T xi(1)  
  e_npt = e_npt +     temp * xib(1)                        !    kB T xib
  e_npt = e_npt + vxi(1)  * vxi(1)  * 0.5_dp / Q (1)   ! pxi^2  / 2 Q
  e_npt = e_npt + vxib(1) * vxib(1) * 0.5_dp / Qb(1)   ! pxib^2 / 2 Qb
  do inhc = 2 , nhc_n
    e_npt = e_npt + vxi(inhc)  * vxi(inhc)  * 0.5_dp / Q (inhc)   ! pxi^2  / 2 Q
    e_npt = e_npt + vxib(inhc) * vxib(inhc) * 0.5_dp / Qb(inhc)   ! pxib^2 / 2 Qb
    e_npt = e_npt +    temp * xi(inhc)                            !    kB T xi 
    e_npt = e_npt +    temp * xib(inhc)                           !    kB T xib 
  enddo

  io_printnode blankline(stdout)
  io_printnode write(stdout,'(a,6e16.8)')         'e_npt1 ',  press     * simu_cell%omega                  , &
                                                              ve * ve   * 0.5_dp / W                       , &
                                                              L  * temp * xi(1)                            , &
                                                              temp * xib(1)                                , &
                                                              vxi(1)   * vxi(1)  * 0.5_dp / Q (1)          , &
                                                              vxib(1)  * vxib(1) * 0.5_dp / Qb(1)          
#ifdef GFORTRAN
  write ( FMT , * ) 4 * nhc_n
  io_printnode write(stdout,'(a,'// ADJUSTL(FMT) //'e16.8)') 'e_npt2 ',( vxi(inhc)  * vxi(inhc)  * 0.5_dp / Q (inhc)  , &
#else
  io_printnode write(stdout,'(a,<4*nhc_n>e16.8)') 'e_npt2 ',( vxi(inhc)  * vxi(inhc)  * 0.5_dp / Q (inhc)  , &
#endif
                                                              vxib(inhc) * vxib(inhc) * 0.5_dp / Qb(inhc)  , &
                                                              temp * xi(inhc)                              , & 
                                                              temp * xib(inhc)             , inhc=2,nhc_n)
  io_printnode write(stdout,'(a,e16.8)')          'e_npt  ', e_npt
  io_printnode blankline(stdout)
             
  deallocate(Q , Qb )
  CALL calc_thermo

  return

END SUBROUTINE nose_hoover_chain_n_p 

! *********************** SUBROUTINE prop_pos_vel_verlet ***********************
!
!> \brief
!! propagates position and position in the velet algorithm
!
! ******************************************************************************
SUBROUTINE prop_pos_vel_verlet_npt ( kin , xe , ve , xe0 , L , W )

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , massia , rx , ry , rz , ry , vx , vy , vz , fx , fy , fz, simu_cell, rhoN
  USE md,                       ONLY :  dt , nhc_n , temp , press , first_time_xe0
  USE engforce_driver,          ONLY :  engforce
  USE cell,                     ONLY :  lattice , kardir, dirkar, periodicbc
  USE thermodynamic,            ONLY :  calc_thermo, pvirial_tot
  USE io,                       ONLY :  stdout, ioprintnode

  implicit none
  ! global
  real(kind=dp) :: kin, xe, ve , xe0, L, W
  ! local
  integer :: ia
  real(kind=dp) :: dt2 , dt4 , e2, e4, e6, e8
  real(kind=dp) :: AA , AA2 , BB , poly , ARG , ARG2, odnf , SINHA

  ! useful constants 
  odnf = 1.0_dp + 3.0_dp / L 
  dt2 = dt  * 0.5_dp
  dt4 = dt2 * 0.5_dp
  e2 = 1.0_dp / 6.0_dp
  e4 = e2     / 20.0_dp
  e6 = e4     / 42.0_dp
  e8 = e6     / 72.0_dp

  ! =========================================
  !  v(t+dt2) = v(t) + f(t) * dt2
  !  r(t+dt)  = r(t) + v(t+dt2) * dt / m 
  ! note : dt2 = dt / 2
  ! =========================================
  ARG = odnf * ve * dt4 / W
  ARG2 = ARG*ARG
  AA  = EXP ( - ARG ) 
  AA2 = AA * AA
  poly = ( ( (e8*ARG2+E6)*ARG2 + E4) *ARG2 + E2 ) * ARG2 + 1.0_dp ! poly = sinh(x)/x = 1.0 + e2 * x^2 + e4 * x^4 + e6 * x^6 + ... 
  BB = AA * dt2 * poly !SINHA
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) * AA2  + fx ( ia ) * BB / massia(ia)
    vy ( ia ) = vy ( ia ) * AA2  + fy ( ia ) * BB / massia(ia)
    vz ( ia ) = vz ( ia ) * AA2  + fz ( ia ) * BB / massia(ia)
  enddo

  ! xe = log( V / V0 ) / 3
  ! 3 xe = log (V) - log (V0)
  ! V = exp ( 3xe + log (V0) ) 
  ! lambda = V**(1/3) = exp ( xe + log (V0)/3  ) 
  ! xe0 =  log (V0)/3
  ! lambda =  exp ( xe + xe0 ) 

  ! EXP ( 3 xe ) = V / V0
  ! EXP ( xe )   = lambda / lambda0

  xe = xe + ve * dt / W
  simu_cell%A(1,:) = simu_cell%A(1,:) * EXP(xe-xe0) 
  simu_cell%A(2,:) = simu_cell%A(2,:) * EXP(xe-xe0) 
  simu_cell%A(3,:) = simu_cell%A(3,:) * EXP(xe-xe0) 
  xe0 = xe 

  CALL lattice ( simu_cell )
  rhoN = REAL ( natm , kind =dp) / simu_cell%omega

  ARG = ve*dt2 / W
  ARG2 = ARG*ARG 
  AA  = EXP ( ARG ) 
  AA2 = AA * AA
  poly = ( ( (e8*ARG2+E6)*ARG2 + E4) *ARG2 + E2 ) * ARG2 + 1.0_dp
  BB = AA * poly * dt
  do ia = 1,natm
    rx ( ia ) = rx ( ia ) * AA2 + vx (ia ) * BB 
    ry ( ia ) = ry ( ia ) * AA2 + vy (ia ) * BB 
    rz ( ia ) = rz ( ia ) * AA2 + vz (ia ) * BB 
  enddo

  ! ==========================
  ! f(t+dt)
  ! ==========================
  CALL engforce

  ! =========================================
  !   v(t) = v(t+dt2) + f(t+dt) * dt2
  !   v(t) = v(t) + dt2 * ( f(t) + f(t+ft) ) 
  ! =========================================
  ARG = odnf * ve * dt4 / W
  ARG2 = ARG*ARG
  AA  = EXP ( - ARG ) 
  AA2 = AA * AA
  poly = ( ( (e8*ARG2+E6)*ARG2 + E4) *ARG2 + E2 ) * ARG2 + 1.0_dp
  BB = AA * poly * dt2
  kin  = 0.0_dp
  do ia = 1 , natm
    vx ( ia ) = vx ( ia ) * AA2  + fx ( ia ) * BB / massia(ia)
    vy ( ia ) = vy ( ia ) * AA2  + fy ( ia ) * BB / massia(ia)
    vz ( ia ) = vz ( ia ) * AA2  + fz ( ia ) * BB / massia(ia)
    kin = kin + ( vx ( ia ) * vx ( ia ) +  vy ( ia ) * vy ( ia ) + vz ( ia ) * vz ( ia ) ) * massia (ia)
  enddo
  kin = kin * 0.5_dp

  return

END SUBROUTINE prop_pos_vel_verlet_npt

SUBROUTINE get_yoshida_weigth(yosh_order,yosh_w)

  USE constants,        ONLY :  dp
  USE io ,              ONLY :  ionode, stdout
  USE md,               ONLY :  yosh_allowed
  
  implicit none
  !global
  integer,       intent (in)    :: yosh_order
  real(kind=dp), intent (inout) :: yosh_w(yosh_order) 
  
  SELECT CASE ( yosh_order )
  CASE DEFAULT
    io_node WRITE(stdout,'(a,5i3)') 'value of yoshida order not available try :',yosh_allowed
    STOP
  CASE (1)
     yosh_w(1) = 1.0_dp
  CASE (3)
     yosh_w(1) = 1.0_dp/(2.0_dp-(2.0_dp)**(1.0_dp/3.0_dp))
     yosh_w(2) = 1.0_dp - 2.0_dp*yosh_w(1)
     yosh_w(3) = yosh_w(1)
  CASE (5)
     yosh_w(1) = 1.0_dp/(4.0_dp-(4.0_dp)**(1.0_dp/3.0_dp))
     yosh_w(2) = yosh_w(1)
     yosh_w(3) = yosh_w(1)
     yosh_w(4) = yosh_w(1)
     yosh_w(5) = 1.0_dp - 4.0_dp*yosh_w(1)
  CASE (7)
     yosh_w(1) = .78451361047756_dp
     yosh_w(2) = .235573213359357_dp
     yosh_w(3) = -1.17767998417887_dp
     yosh_w(4) = 1.0_dp - 2.0_dp*(yosh_w(1)+yosh_w(2)+yosh_w(3))
     yosh_w(5) = yosh_w(3)
     yosh_w(6) = yosh_w(2)
     yosh_w(7) = yosh_w(1)
   CASE (9)
     yosh_w(1) = 0.192_dp
     yosh_w(2) = 0.554910818409783619692725006662999_dp
     yosh_w(3) = 0.124659619941888644216504240951585_dp
     yosh_w(4) = -0.843182063596933505315033808282941_dp
     yosh_w(5) = 1.0_dp - 2.0_dp*(yosh_w(1)+yosh_w(2)+yosh_w(3)+yosh_w(4))
     yosh_w(6) = yosh_w(4)
     yosh_w(7) = yosh_w(3)
     yosh_w(8) = yosh_w(2)
     yosh_w(9) = yosh_w(1)
  END SELECT
END SUBROUTINE
! ===== fmV =====
