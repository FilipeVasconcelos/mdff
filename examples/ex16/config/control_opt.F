&controltag
        calc          = 'opt'
        cutshortrange = 2.5d0 
        skindiff      = 0.3d0
        !posff_data    = 'rvf'
        !posff_data    = 'rnn'
        lnmlj         = .true.
        !lvnlist       = .true.
        !ltraj         = .true.
        !itraj_start   = 100000 
        !itraj_period  = 1000
        lreduced      = .true.
        lreducedN     = .true.
&end
&opttag
        nconf =1000
        nmaxopt=1000
        optalgo='m1qn3'
        !optalgo='sastry'
        !optalgo='lbfgs'
        epsrel_m1qn3 = 1e-6
        !epsrel_lbfgs=1e-3
        nopt_print=16
        optvar='pos'
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
