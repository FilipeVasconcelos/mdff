&controltag
        calc     = 'vib+dos'
        lnmlj    = .true.
        lreduced = .false.
        cutshortrange=8.5d0
&end

&vibtag
        nconf        = 1
        resdos       = 0.01d0    
        omegamax     = 0.80d0
        nkphon       = 50 ! (nkphon+1) * (nkphon+1) * (nkphon+1) generated in IBZKPTFF
        ncell        = 4
&end

&mdtag
&end

&fieldtag
        epslj(1,1)    = 0.010323576d0 ! lennard-jones potential in [eV] = 119.8 K
        sigmalj(1,1)  = 3.405d0       ! lennard-jones potential in [A]

&end
