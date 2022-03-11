&controltag
        calc      = 'vib+band'
        cutshortrange    = 8.5d0
        lnmlj     = .true.
        lreduced  = .false.
&end

&vibtag
  nconf        = 1
  resdos       = 0.05d0    
  omegamax     = 0.80d0
  nkphon       = 100
  path='G-X'
  ks = 0.0d0 0.0d0 0.0d0,
  ks = 0.0d0 1.0d0 1.0d0,
  ncell=4
&end

&mdtag
&end

&proptag 
&end
&fieldtag
&end
&non_bondedtag
        epslj(1,1)    = 0.010323576d0 ! lennard-jones potential in [eV] = 119.8 K
        sigmalj(1,1)  = 3.405d0       ! lennard-jones potential in [A]
&end
