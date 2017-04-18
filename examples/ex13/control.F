&controltag
        calc          = 'md'
        cutshortrange = 8.0d0
        skindiff      = 1.0d0
        posff_data  = 'rvf'
        lreduced      = .false.
        lrestart      = .true.
        lnmlj         = .true.
        lvnlist       = .false.
        ltraj         = .true.
        itraj_start   = 0 
        itraj_period  = 10 
&end
&mdtag
        npas          = 200000
        nprint        = 1000
        fprint        = 10

        !nequil =100000        
        !nequil_period =10

        integrator    = 'npt-nhcnp'
        dt            = 0.0005
        temp          = 119.8d0
        press         = 10.0d0 
        nhc_n         = 3
        nhc_yosh_order= 3 
        nhc_mults     = 3
        timesca_thermo= 2.0 
        timesca_baro  = 2.0

&end
&fieldtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 0.010323576 
        sigmalj(1,1)= 3.405d0
        mass=39.95d0
&end
