&controltag
        calc          = 'md'
        cutshortrange = 8.5d0
        skindiff      = 1.0d0
        posff_data    = 'rnn'
        lreduced      = .false.
        lwrite_restart= .true.
        full_restart  = .false.
        lnmlj         = .true.
        lvnlist       = .true.
        ltraj         = .true.
        itraj_start   = 0
        itraj_period  = 10
&end
&mdtag
        npas          = 4000
        nequil        = 2000
        nequil_period = 10
        nprint        = 100
        fprint        = 10
        spas          = 10000
        
        dt            = 0.001 
        temp          = 86.56d0 
&end
&fieldtag
        mass=39.95d0
&end
&non_bondedtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 0.010323576 
        sigmalj(1,1)= 3.405d0
&end
