!#Example 1: LJ fcc structure at low temperature

! ======================================================================
! main parameters
! ======================================================================
&controltag
      calc           = 'md'      ! md calculation
      lstatic        = .true.    ! static
      cutshortrange  = 100.0d0   
      lnmlj          = .true.    ! lennard-jones
&end

! ======================================================================
! molecular dynamics parameters 
! ======================================================================
&mdtag
        dt=0.001
&end
! ======================================================================
! force-field parameters
! ======================================================================
&fieldtag
        ctrunc='linear'
        epslj(1,1)        = 1.0d0     ! espilon in the lennard jones 
        sigmalj(1,1)      = 1.0d0     ! sigma in the lennard jones
&end
! ======================================================================
! properties 
! ======================================================================
&proptag 
&end
