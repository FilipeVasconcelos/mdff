&controltag
        calc          = 'md'
        cutshortrange = 2.5d0 
        skindiff      = 0.3d0
        posff_data    = 'rvf'
        !posff_data    = 'rnn'
        lnmlj         = .true.
        lvnlist       = .true.
        ltraj         = .true.
        itraj_start   = 100000 
        itraj_period  = 1000
        lreduced      = .true.
        lreducedN     = .true.
&end
&mdtag
        npas          = 1100000
        nprint        =    1000
        fprint        =     100

        integrator    = 'nve-vv'
        nequil        =  100000
        nequil_period =      10
        spas          =   10000


        !integrator    = 'nvt-nhcn'
        dt            = 0.003
        temp          = 1.0
        !nhc_n         = 2 
        !nhc_yosh_order= 1 
        !nhc_mults     = 1
        !timesca_thermo= 2.0 
&end
&fieldtag
        mass=1.0d0
&end
&non_bondedtag
        ctrunc = 'quadratic'
        epslj(1,1)  =  1.0d0
        epslj(1,2)  =  1.5d0
        epslj(2,2)  =  0.5d0

        sigmalj(1,1)= 1.0d0 
        sigmalj(1,2)= 0.8d0
        sigmalj(2,2)= 0.88d0
&end
