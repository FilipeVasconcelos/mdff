&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        posff_data    = 'rnn'
        lnmlj         = .true.
        lvnlist       = .false.
        ltraj         = .true.
        itraj_start   = 0 
        itraj_period  = 5
&end
&mdtag
        npas          = 40000 
        nprint        = 100
        fprint        = 100

        integrator    = 'nvt-nh'
        dt            = 0.006
        temp          = 119.8d0
        timesca_thermo= 2.0 
&end
&fieldtag
        mass=39.95d0
&end
&non_bondedtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 0.010323576 
        sigmalj(1,1)= 3.405d0
&end
