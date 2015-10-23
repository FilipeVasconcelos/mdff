&controltag
        calc          = 'md'
        cutshortrange = 6.8793 
        cutlongrange  = 6.8793
        skindiff      = 1.0d0
        lvnlist       = .false.
        lrestart      = .false.
        restart_data  = 'rnn' 
        lbmhftd       = .true.
        lcoulomb      = .true.
        ltraj         = .false. 
        itraj_start   = 0 
        itraj_period  = 2000
        lstatic       = .false.
&end

&mdtag
        integrator    = 'nve-vv' 
        npas          = 0 
        dt            = 0.0005
        nprint        = 20
        fprint        = 20

        !INCLUDE <thermostat>
        !INCLUDE <barostat>
        temp          = 2500.0
        nequil        = 0 
        nequil_period = 0

&end

&fieldtag
        ! ewald summation
        lautoES       = .true.
        epsw          = 1.e-6
        algo_moment_from_pola='scf_kO_v4_1'
!        algo_moment_from_pola='scf_kO_v2'
        omegakO       = 0.7 
        conv_tol_ind  = 1e-6

        !         O      Si        Ca
        mass=15.9994d0 28.0855d0  40.078d0
        !         O       Si
        qch= -2.0d0      4.0d0     2.0d0

        ! ===========================================
        !                   DIPOLES
        ! ===========================================

        ! ===========================================
        !               dip polar on O
        ldip_polar(1)=.true.
        poldip_iso(1)=1.299150166912089
        ! O - O 
        ldip_damping(1,1,1) = .false.
        ! O - Si
        ldip_damping(1,1,2) = .true.
        pol_damp_b(1,1,2) = 4.015668032381992
        pol_damp_c(1,1,2) = 1.5933
        pol_damp_k(1,1,2) = 4
        ! O - Ca 
        ldip_damping(1,1,3) = .true.
        pol_damp_b(1,1,3) = 3.457775435124017
        pol_damp_c(1,1,3) = 2.35919526229
        pol_damp_k(1,1,3) = 4
        ! ===========================================
        !               dip polar on Si
        ldip_polar(2)=.false.
        ! ===========================================
        !               dip polar on Ca 
        ldip_polar(3)=.true.
        poldip_iso(3)=0.5186409672943927
        ! Ca - O
        ldip_damping(3,1,3) = .true.
        pol_damp_b(3,1,3) = 6.613550738353049
        pol_damp_c(3,1,3) = 1.003717874670000
        pol_damp_k(3,1,3) = 4
        ! Ca - Si
        ldip_damping(3,2,3) = .true.
        pol_damp_b(3,2,3) = 6.613550738353049 
        pol_damp_c(3,2,3) = 1.003717874670000
        pol_damp_k(3,2,3) = 4
        ! O - O 
        ldip_damping(3,3,3) = .false.

        ! ===========================================
        !                QUADRUPOLES 
        ! ===========================================

        ! ===========================================
        !               quad polar on O
        lquad_polar(1)=.true.
        polquad_iso(1)=0.4777179127647640
        ! ===========================================
        !               quad polar on Si 
        lquad_polar(2)=.false.
        ! ===========================================
        !               quad polar on Ca 
        lquad_polar(3)=.true.
        polquad_iso(3)=0.2066049444100642




        ! ===========================================
        !                   BMHFTD 
        ! ===========================================
        ! BMHFTD potential
        !O-O               eV               hartree
        Abmhftd(1,1)  = 29061.759      !    1068.0
        Bbmhftd(1,1)  = 5.0376 
        Cbmhftd(1,1)  = 26.512052      !    0.9743
        Dbmhftd(1,1)  = 142.77542      !    5.2469
        BDbmhftd(1,1) = 2.7184 
        !O-Si               eV           hartree
        Abmhftd(1,2)  = 1177.6271      !  43.2770
        Bbmhftd(1,2)  = 2.9136 
        Cbmhftd(1,2)  = 1.3034253      !  0.0479  
        Dbmhftd(1,2)  = 4.2340915      !  0.1556
        BDbmhftd(1,2) = 0.0            
        !O-Ca               eV           hartree
        Abmhftd(1,3)  = 1090.0391      ! 40.0582
        Bbmhftd(1,3)  = 2.8413
        Cbmhftd(1,3)  = 1.3034253      ! 0.0479
        Dbmhftd(1,3)  = 4.2340915      ! 0.1556
        BDbmhftd(1,3) = 0.0
        !Si-Si               eV           hartree
        Abmhftd(2,2)  = 0.0    
        Bbmhftd(2,2)  = 2.8170
        Cbmhftd(2,2)  = 0.12789351      ! 0.0047
        Dbmhftd(2,2)  = 0.14422034      ! 0.0053
        BDbmhftd(2,2) = 1889728.5785
        !Si-Ca               eV           hartree
        Abmhftd(2,3)  = 0.0   
        Bbmhftd(2,3)  = 2.8170
        Cbmhftd(2,3)  = 0.12789351     ! 0.0047
        Dbmhftd(2,3)  = 0.14422034     ! 0.0053
        BDbmhftd(2,3) = 1889728.5785
        !Ca-Ca               eV           hartree
        Abmhftd(3,3)  = 0.0  
        Bbmhftd(3,3)  = 2.8170
        Cbmhftd(3,3)  = 0.12789351      ! 0.0047
        Dbmhftd(3,3)  = 0.14422034      ! 0.0053
        BDbmhftd(3,3) = 1889728.5785

&end
