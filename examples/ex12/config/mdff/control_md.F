&controltag
        calc          = 'md'
        cutshortrange = 5.0
        cutlongrange  = 5.0
        skindiff      = 1.0d0
        lcsvr         = .true.
        lvnlist       = .true.
        lwrite_restart= .false.
        posff_data    = 'rnn' 
        lbmhftd       = .true.
        lcoulomb      = .true.
        ltraj         = .true. 
        itraj_start   = 0 
        itraj_period  = 2
        !lstatic=.true.
&end

&mdtag
        integrator    = 'nvt-nhcn' 
        !integrator    = 'nve-vv' 
        npas          = 4000
        dt            = 0.0005
        nprint        = 10 
        fprint        = 10

        nhc_n         = 3 
        nhc_yosh_order= 3 
        nhc_mults     = 2
        timesca_thermo= 1.0 
        !taucsvr       = 1.0
        !INCLUDE <barostat>
        temp          = 2500.0
        !nequil        = 500
        !nequil_period = 1
&end
&coulombtag
        ! ewald summation
        lautoES       = .true.
        epsw          = 1.e-5
&end

&fieldtag
        !         O      Si 
        mass=15.99994d0 28.084d0
        !         O       Si
        qch= -2.0d0      4.0d0
&end
&pimtag
        algo_moment_from_pola = 'scf_kO_v1'
        omegakO       = 0.7 
        conv_tol_ind  = 1e-6
        algo_ext_dipole = 'aspc'
        extrapolate_order = 4
        ! ===========================================
        !               polar on O
        ldip_polar(1)=.true.
        poldip(1,1,1)=1.59150
        poldip(1,2,2)=1.59150
        poldip(1,3,3)=1.59150
        ! O - O
        ldip_damping(1,1,1) = .true.
        pol_damp_b(1,1,1) = 4.74888
        pol_damp_c(1,1,1) = 2.227
        pol_damp_k(1,1,1) = 4
        ! O - Si
        ldip_damping(1,1,2) = .true.
        pol_damp_b(1,1,2) = 3.66480
        pol_damp_c(1,1,2) = 1.44589
        pol_damp_k(1,1,2) = 4
&end
&non_bondedtag
        ! BMHFTD potential
        !O-O               eV               hartree
        Abmhftd(1,1)  = 7902.1862      !    290.4
        Bbmhftd(1,1)  = 4.54668
        Cbmhftd(1,1)  = 13.145548      !    0.48309
        Dbmhftd(1,1)  = 71.279951      !    2.61949
        BDbmhftd(1,1) = 2.64562 
        !O-Si               eV           hartree
        Abmhftd(1,2)  = 1311.5888      ! 48.2
        Bbmhftd(1,2)  = 3.18568
        Cbmhftd(1,2)  = 0.0
        Dbmhftd(1,2)  = 0.0
        BDbmhftd(1,2) = 0.0

&end
