MODULE tensors_rk

  USE constants,        ONLY :  dp

  implicit none

  TYPE :: tensor_rank0
    real( kind=dp ) :: sca 
    !real( kind=dp ) :: sca_damp      
    !real( kind=dp ) :: sca_damp2      
  END TYPE
  TYPE :: tensor_rank1
    real( kind=dp ) :: a(3)      
    real( kind=dp ) :: a_damp(3)      
    real( kind=dp ) :: a_damp2(3)      
  END TYPE
  TYPE :: tensor_rank2
    real( kind=dp ) :: ab(3,3)      
    real( kind=dp ) :: ab_damp(3,3)      
    real( kind=dp ) :: ab_damp2(3,3)      
    real( kind=dp ) :: ab_thole(3,3)      
  END TYPE
  TYPE :: tensor_rank3
    real( kind=dp ) :: abc(3,3,3)      
    !real( kind=dp ) :: abc_damp(3,3,3)      
    !real( kind=dp ) :: abc_damp2(3,3,3)      
  END TYPE
 
  TYPE :: interaction_dd
    TYPE ( tensor_rank2 ) :: T2
    TYPE ( tensor_rank3 ) :: T3
  END TYPE


END MODULE tensors_rk
