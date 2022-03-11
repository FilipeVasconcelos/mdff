&controltag
        calc    = 'vib'
        lnmlj   = .true.
        lreduced= .false.
        cutshortrange=8.5d0
&end

&vibtag
        lwrite_vectff = .true.
        nconf        = 1
        resdos       = 0.1d0    
        omegamax     = 1.0d0
        ncell        = 4
&end

&mdtag
&end

&fieldtag
&end
&non_bondedtag
        epslj(1,1)    = 0.010323576d0 ! lennard-jones potential in [eV] = 119.8 K
        sigmalj(1,1)  = 3.405d0       ! lennard-jones potential in [A]
&end
