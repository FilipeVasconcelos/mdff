!#define debug_TT
MODULE tt_damp

#include "symbol.h"
  USE constants,        ONLY :  dp
  USE io,               ONLY :  stdout, ionode

  integer      , PARAMETER                             :: maximum_of_TT_expansion=8
  real(kind=dp), dimension (0:maximum_of_TT_expansion) :: E_TT 

CONTAINS

SUBROUTINE get_TT_damp 

  implicit none

  integer :: k 

  E_TT(0) = 1.0_dp
  do k = 1 , maximum_of_TT_expansion
    E_TT(k) = E_TT(k-1) / REAL ( k ,kind=dp )
#ifdef debug_TT  
    io_node write( stdout , '(a,i,e20.12)') 'TT_damp coeff =',k,E_TT(k)
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
