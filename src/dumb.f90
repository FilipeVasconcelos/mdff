module dumb

  USE constants,        ONLY : dp
  real(kind=dp) :: t12,t23,t34,t45,t56,t67,t78,t89,t90

CONTAINS

SUBROUTINE init_dumb
  
  implicit none 

  t12 = 0.0_dp 
  t23 = 0.0_dp 
  t34 = 0.0_dp 
  t45 = 0.0_dp 
  t56 = 0.0_dp 
  t67 = 0.0_dp 
  t78 = 0.0_dp 
  t89 = 0.0_dp 
  t90 = 0.0_dp 

END SUBROUTINE init_dumb

SUBROUTINE print_dumb

  USE io,       ONLY :  stdout, ionode 

  implicit none

  real( kind=dp ) :: somme 

  if ( .not. ionode ) return
 
  somme = t12+t23+t34+t45+t56+t67+t78+t89+t90

  if ( somme .ne. 0.0_dp ) then
    write ( stdout , '(a)' ) '' 
    write ( stdout , '(a)' ) 'SUPPLEMENTARTY TIMING INFO (seed dumb.f90) ' 
  endif

  if ( t12 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t12 = ",t12
  if ( t23 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t23 = ",t23
  if ( t34 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t34 = ",t34
  if ( t45 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t45 = ",t45
  if ( t56 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t56 = ",t56
  if ( t67 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t67 = ",t67
  if ( t78 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t78 = ",t78
  if ( t89 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t89 = ",t89
  if ( t90 .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "t90 = ",t90
  if ( somme .ne. 0.0_dp ) write ( stdout,'(a,f8.2)') "sum = ",somme

  write ( stdout , '(a)' ) '' 

  return

end subroutine

end module dumb
