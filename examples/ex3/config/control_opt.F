&controltag
        calc          = 'opt'
        cutshortrange = 50.0d0
        skindiff      = 1.0d0
        posff_data  = 'rnn'
        lnmlj         = .true.
        lvnlist       = .false.
        ltraj         = .true.
        itraj_start   = 0 
        itraj_period  = 5
        lreduced = .true.
&end
&opttag
        nconf =1
        nmaxopt=10
        optalgo='m1qn3'
        nopt_print=1
        epsrel_m1qn3 = 1e-8 
        nopt_print=2

&end
&fieldtag
        ctrunc = 'notrunc'
        epslj(1,1)  = 1.0d0
        sigmalj(1,1)= 3.405d0

&end
