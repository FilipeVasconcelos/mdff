!#Example 0: Minimum settings
! ======================================================================
! main parameters
! ======================================================================
&controltag
        lnmlj   = .true.
        lstatic = .true.
        cutshortrange  = 10.d0
        lvnlist =.false.
&end
! ======================================================================
! parameters for the molecular dynamics 
! ======================================================================
&mdtag
&end
! ======================================================================
! parameters for the force-field 
! ======================================================================
&fieldtag
        mass = 1.0d0
&end
&non_bondedtag
        ctrunc='notrunc'
        epslj(1,1)=0.010323576 !119.8d0 K in eV 
        sigmalj(1,1)=3.405d0
&end
