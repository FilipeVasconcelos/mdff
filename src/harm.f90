SUBROUTINE HARM

END SUBROUTINE

! ==========================
!  Deterime the force:
!  Potential U = 0.5 x^2
! ==========================
SUBROUTINE engforce_harm

  USE constants,                ONLY :  dp
  USE config,                   ONLY :  natm , rx , fx
  USE thermodynamic,            ONLY :  u_harm

  implicit none
  u_harm = 0.5_dp * rx(1) * rx(1) 
  fx(1)  = - rx(1)    
  return

END SUBROUTINE
