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

!> \brief
!! cell related module
MODULE cell
  
  USE constants,                ONLY :  dp 
  USE mpimdff

  implicit none      

  !< 
  TYPE celltype
    real(kind=dp) :: A(3,3)               !< direct basis vector  
    real(kind=dp) :: B(3,3)               !< reciprocal basis vectors  
    real(kind=dp) :: G(3,3)               !< metric tensor (A^T A)  
    real(kind=dp) :: ANORM(3)             !< norm of direct basis vectors
    real(kind=dp) :: BNORM(3)             !< norm of reciprocal basis vectors
    real(kind=dp) :: OMEGA                !< volume ( direct )
    real(kind=dp) :: ROMEGA               !< volume ( reciprocal )
    real(kind=dp) :: WA , WB , WC         !< perpendicular width (direct) 
    real(kind=dp) :: ALPH , BET , GAMM    !< angles ( direct )
    real(kind=dp) :: RWA , RWB , RWC      !< perpendicular width (reciprocal)
    real(kind=dp) :: RALPH , RBET , RGAMM !< angles ( reciprocal )
  END TYPE

CONTAINS

! *********************** SUBROUTINE lattice ***********************************
!> \brief
!! subroutine for calculating the reciprocal lattice from the direct 
!! lattice in addition the norm of the lattice-vectors and the volume of 
!! the basis-cell are calculated
!! \note
!! adapted from vasp
!! \author 
!! gK (VASP)
!! \param[in,out] Mylatt
! ****************************************************************************** 
SUBROUTINE lattice ( Mylatt )

  USE constants,                ONLY : radian

  implicit none 
 
  ! global
  TYPE(celltype) Mylatt
  ! local 
  real(kind=dp) :: omega , romega
  integer :: i, j 
  intrinsic SUM 
  real(kind=dp) :: WA , WB , WC 
  real(kind=dp) :: RWA , RWB , RWC 
  real(kind=dp) :: alph, bet , gamm 
  real(kind=dp) :: ralph, rbet , rgamm 

  ! metric tensor
  ! A^TA
  Mylatt%G(1,1) = Mylatt%A(1,1)**2.0d0 + Mylatt%A(2,1)**2.0d0 + Mylatt%A(3,1)**2.0d0
  Mylatt%G(2,1) = Mylatt%A(2,1)*Mylatt%A(2,2) + Mylatt%A(3,1)*Mylatt%A(3,2)
  Mylatt%G(3,1) = Mylatt%A(3,1)*Mylatt%A(3,3) 
  Mylatt%G(1,2) = Mylatt%G(2,1) 
  Mylatt%G(2,2) = Mylatt%A(2,2)**2.0d0 + Mylatt%A(3,2)**2.0d0
  Mylatt%G(3,2) = Mylatt%A(3,2)*Mylatt%A(3,3)
  Mylatt%G(1,3) = Mylatt%G(3,1)
  Mylatt%G(2,3) = Mylatt%G(3,2)
  Mylatt%G(3,3) = Mylatt%A(3,3)**2.0d0


  CALL EXPRO(Mylatt%B(1:3,1),Mylatt%A(1:3,2),Mylatt%A(1:3,3))  ! B x C
  CALL EXPRO(Mylatt%B(1:3,2),Mylatt%A(1:3,3),Mylatt%A(1:3,1))  ! C x A
  CALL EXPRO(Mylatt%B(1:3,3),Mylatt%A(1:3,1),Mylatt%A(1:3,2))  ! A x B


  ! volume ( direct ) 
  omega = Mylatt%B(1,1)*Mylatt%A(1,1)+Mylatt%B(2,1)*Mylatt%A(2,1) + Mylatt%B(3,1)*Mylatt%A(3,1)
  Mylatt%omega=omega

  ! shortest distance between opposite faces
  WA = omega / SQRT ( Mylatt%B(1,1) * Mylatt%B(1,1) + Mylatt%B(2,1) * Mylatt%B(2,1) + Mylatt%B(3,1) * Mylatt%B(3,1) ) 
  WB = omega / SQRT ( Mylatt%B(1,2) * Mylatt%B(1,2) + Mylatt%B(2,2) * Mylatt%B(2,2) + Mylatt%B(3,2) * Mylatt%B(3,2) ) 
  WC = omega / SQRT ( Mylatt%B(1,3) * Mylatt%B(1,3) + Mylatt%B(2,3) * Mylatt%B(2,3) + Mylatt%B(3,3) * Mylatt%B(3,3) ) 

  Mylatt%WA=WA 
  Mylatt%WB=WB 
  Mylatt%WC=WC

  do i=1,3
    do j=1,3
      Mylatt%B(i,j)=Mylatt%B(i,j)/omega
    enddo
  enddo

  ! volume ( reciprocal ) 
  romega =  Mylatt%A(1,1)*Mylatt%B(1,1)+Mylatt%A(2,1)*Mylatt%B(2,1) + Mylatt%A(3,1)*Mylatt%B(3,1)
  Mylatt%Romega=romega

  ! shortest distance between opposite faces ( reicprocal )
  RWA = Romega / SQRT ( Mylatt%A(1,1) * Mylatt%A(1,1) + Mylatt%A(2,1) * Mylatt%A(2,1) + Mylatt%A(3,1) * Mylatt%A(3,1) )
  RWB = Romega / SQRT ( Mylatt%A(1,2) * Mylatt%A(1,2) + Mylatt%A(2,2) * Mylatt%A(2,2) + Mylatt%A(3,2) * Mylatt%A(3,2) )
  RWC = Romega / SQRT ( Mylatt%A(1,3) * Mylatt%A(1,3) + Mylatt%A(2,3) * Mylatt%A(2,3) + Mylatt%A(3,3) * Mylatt%A(3,3) )

  Mylatt%RWA=RWA
  Mylatt%RWB=RWB
  Mylatt%RWC=RWC

  do i=1,3
    Mylatt%ANORM(i)=SQRT(SUM(Mylatt%A(:,i)*Mylatt%A(:,i)))
    Mylatt%BNORM(i)=SQRT(SUM(Mylatt%B(:,i)*Mylatt%B(:,i)))
  enddo

  ! angles ( direct )
  alph = Mylatt%A(1,3) * Mylatt%A(1,2) + Mylatt%A(2,3) * Mylatt%A(2,2) + Mylatt%A(3,3) * Mylatt%A(3,2)    ! C . B 
  bet  = Mylatt%A(1,1) * Mylatt%A(1,3) + Mylatt%A(2,1) * Mylatt%A(2,3) + Mylatt%A(3,1) * Mylatt%A(3,3)    ! A . C
  gamm = Mylatt%A(1,1) * Mylatt%A(1,2) + Mylatt%A(2,1) * Mylatt%A(2,2) + Mylatt%A(3,1) * Mylatt%A(3,2)    ! A . B 

  alph = alph / ( Mylatt%ANORM(3) * Mylatt%ANORM(2) ) 
  bet  = bet  / ( Mylatt%ANORM(1) * Mylatt%ANORM(3) ) 
  gamm = gamm / ( Mylatt%ANORM(1) * Mylatt%ANORM(2) ) 
 
  alph = acos ( alph ) * radian
  bet  = acos ( bet  ) * radian
  gamm = acos ( gamm ) * radian

  Mylatt%ALPH  = alph
  Mylatt%BET   = bet
  Mylatt%GAMM  = gamm

  ! angles ( reciprocal )
  ralph = Mylatt%B(1,3) * Mylatt%B(1,2) + Mylatt%B(2,3) * Mylatt%B(2,2) + Mylatt%B(3,3) * Mylatt%B(3,2)    ! C* . B* 
  rbet  = Mylatt%B(1,1) * Mylatt%B(1,3) + Mylatt%B(2,1) * Mylatt%B(2,3) + Mylatt%B(3,1) * Mylatt%B(3,3)    ! A* . C*
  rgamm = Mylatt%B(1,1) * Mylatt%B(1,2) + Mylatt%B(2,1) * Mylatt%B(2,2) + Mylatt%B(3,1) * Mylatt%B(3,2)    ! A* . B* 

  ralph = ralph / ( Mylatt%BNORM(3) * Mylatt%BNORM(2) ) 
  rbet  = rbet  / ( Mylatt%BNORM(1) * Mylatt%BNORM(3) ) 
  rgamm = rgamm / ( Mylatt%BNORM(1) * Mylatt%BNORM(2) ) 
 
  ralph = acos ( ralph ) * radian
  rbet  = acos ( rbet  ) * radian
  rgamm = acos ( rgamm ) * radian

  Mylatt%RALPH  = ralph
  Mylatt%RBET   = rbet
  Mylatt%RGAMM  = rgamm

  RETURN

END SUBROUTINE lattice

! *********************** SUBROUTINE kardir ************************************
!> \brief
!! transform a set of vectors from cartesian coordinates to
!! ) direct lattice      (BASIS must be equal to B reciprocal lattice)
!! ) reciprocal lattice  (BASIS must be equal to A direct lattice)
!! \author 
!! gK (VASP)
!! \note 
!! adapted from VASP
!! \param[in] NMAX dimension of vectors VX , VY , VZ 
!! \param[in,out] VX , VY , VZ vectors being transformed
!! \param[in] BASIS basis vector ( direct or reciprocal lattice )
! ******************************************************************************
SUBROUTINE kardir ( NMAX , VX , VY , VZ , BASIS )

  USE time,     ONLY : kardirtottime

  implicit none 

  ! global
  integer , intent(in) :: NMAX
  real(kind=dp) :: VX(NMAX), VY(NMAX),VZ(NMAX), BASIS(3,3)

  ! local 
  integer :: N
  real(kind=dp) :: V1 , V2 , V3
#ifdef MPI
  dectime

  statime  
#endif
  do N=1,NMAX
    V1=VX(N)*BASIS(1,1)+VY(N)*BASIS(2,1)+VZ(N)*BASIS(3,1)
    V2=VX(N)*BASIS(1,2)+VY(N)*BASIS(2,2)+VZ(N)*BASIS(3,2)
    V3=VX(N)*BASIS(1,3)+VY(N)*BASIS(2,3)+VZ(N)*BASIS(3,3)
    VX(N)=V1
    VY(N)=V2
    VZ(N)=V3
  enddo

#ifdef MPI
  stotime
  addtime(kardirtottime)
#endif
  return

END SUBROUTINE

! *********************** SUBROUTINE dirkar ************************************
!> \brief
!! transform a set of vectors from
!! ) direct lattice      (BASIS must be equal to A direct lattice)
!! ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
!! to cartesian coordinates
!! \author 
!! gK (VASP)
!! \note 
!! adapted from VASP
!! \param[in] NMAX dimension of vectors VX , VY , VZ 
!! \param[in,out] VX , VY , VZ vectors being transformed
!! \param[in] BASIS basis vector ( direct or reciprocal lattice )
! ******************************************************************************
SUBROUTINE dirkar ( NMAX , VX , VY , VZ , BASIS )

  USE time,     ONLY : dirkartottime 
  implicit none

  ! global
  integer :: NMAX
  real(kind=dp) :: VX ( NMAX ) , VY ( NMAX ) , VZ ( NMAX ) , BASIS(3,3)

  ! local 
  integer :: N
  real(kind=dp) :: V1 , V2 , V3
#ifdef MPI
  dectime

  statime  
#endif
  do N=1,NMAX
    V1=VX(N)*BASIS(1,1)+VY(N)*BASIS(1,2)+VZ(N)*BASIS(1,3)
    V2=VX(N)*BASIS(2,1)+VY(N)*BASIS(2,2)+VZ(N)*BASIS(2,3)
    V3=VX(N)*BASIS(3,1)+VY(N)*BASIS(3,2)+VZ(N)*BASIS(3,3)
    VX(N)=V1
    VY(N)=V2
    VZ(N)=V3
  enddo

#ifdef MPI
  stotime
  addtime(dirkartottime)
#endif
  return

END SUBROUTINE dirkar 

! *********************** SUBROUTINE periodicpbc *******************************
!> \brief
!! replace atoms inside the MD cell 
!! \note
!! pbc are used in calculation of most properties but the positions are not
!! effectively reajust to prevent big jump when using verlet list. 
!! \param[in] natm number of atoms
!! \param[in] latt lattice type
!! \param[in,out] xxx , yyy , zzz position vectors 
! ******************************************************************************
SUBROUTINE periodicbc ( natm , xxx , yyy , zzz )

  implicit none

  ! global
  integer :: natm
  real(kind=dp) :: xxx ( natm ) , yyy ( natm ) , zzz ( natm )

  ! local
  integer :: ia

  do ia = 1 , natm
     xxx ( ia ) = xxx ( ia ) - NINT ( xxx ( ia )  ) 
     yyy ( ia ) = yyy ( ia ) - NINT ( yyy ( ia )  ) 
     zzz ( ia ) = zzz ( ia ) - NINT ( zzz ( ia )  )  
  enddo

  return

END SUBROUTINE periodicbc


! *********************** SUBROUTINE dirkar ************************************
!> \brief
!! transform a set of scalars from
!! ) direct lattice      (BASIS must be equal to A direct lattice)
!! ) reciprocal lattice  (BASIS must be equal to B reciprocal lattice)
!! to cartesian coordinates
!! \author 
!! gK (VASP)
!! \note 
!! adapted from VASP
!! \param[in] NMAX dimension of vectors VX , VY , VZ 
!! \param[in,out] VX , VY , VZ vectors being transformed
!! \param[in] BASIS basis vector ( direct or reciprocal lattice )
! ******************************************************************************
SUBROUTINE dirkar_1 ( VX , VY , VZ , BASIS , ncell )

  USE time,     ONLY : dirkartottime
  implicit none

  ! global
  integer :: ncell
  real(kind=dp) :: VX , VY , VZ , BASIS(3,3)

  ! local 
  real(kind=dp) :: V1 , V2 , V3
#ifdef MPI
  dectime

  statime
#endif
  V1=VX*BASIS(1,1)+VY*BASIS(1,2)+VZ*BASIS(1,3)
  V2=VX*BASIS(2,1)+VY*BASIS(2,2)+VZ*BASIS(2,3)
  V3=VX*BASIS(3,1)+VY*BASIS(3,2)+VZ*BASIS(3,3)
  VX=V1*REAL(ncell,kind=dp)
  VY=V2*REAL(ncell,kind=dp)
  VZ=V3*REAL(ncell,kind=dp)

#ifdef MPI
  stotime
  addtime(dirkartottime)
#endif
  return

END SUBROUTINE dirkar_1


END MODULE cell

! ===== fmV =====
