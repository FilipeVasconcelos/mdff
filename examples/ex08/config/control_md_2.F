&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        posff_data    = 'rnn'
        lreduced      = .false.
        lwrite_restart= .true.
        lnmlj         = .true.
        lvnlist       = .false.
        ltraj         = .true.
        itraj_start   = 0 
        itraj_period  = 10 
&end
&mdtag
        npas          = 40000
        nprint        = 2000
        fprint        = 10

        integrator    = 'npt-nhcnp'
        dt            = 0.001
        temp          = 85.0d0
        press         = 0.1d0
        nhc_n         = 2
        nhc_yosh_order= 1 
        nhc_mults     = 1
        timesca_thermo= 1.0 
        timesca_baro  = 1.0 
&end
&fieldtag
        mass=39.95d0
&end
&non_bondedtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 0.010323576 
        sigmalj(1,1)= 3.405d0
&end
