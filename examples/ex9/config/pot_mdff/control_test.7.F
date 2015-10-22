&controltag
        calc          = 'md'
        cutshortrange = 7.2 
        cutlongrange  = 7.2
        skindiff      = 1.0d0
        lvnlist       = .true.
        restart_data  = 'rnn'
        lreduced      = .false.
        lbmhft        = .true.
        lbmhftd       = .false.
        lcoulomb      = .true.
        longrange     = 'ewald'
        lstatic       = .true.
&end

&mdtag
&end

&fieldtag

        ctrunc = 'notrunc'
        lautoES = .true.
        epsw = 1e-12 

        min_scf_pol_iter = 6
        max_scf_pol_iter = 100
        conv_tol_ind = 1e-6 
        !algo_moment_from_pola = 'scf_kO_v4_1'
        algo_moment_from_pola = 'scf_kO_v1'

        !         O       Si         B
        mass=15.99994d0 28.0855d0 10.811d0
        !         O       Si         B
        qch= -2.0d0      4.0d0     3.0d0
        !INCLUDE <PIM.POT>
       ! polar on O
        lpolar(1)=.true.
        pol(1,1,1)=0.93971813676926158987
        pol(1,2,2)=0.93971813676926158987
        pol(1,3,3)=0.93971813676926158987

        ! O - O
        ldip_damping(1,1,1) = .true.
        pol_damp_b(1,1,1) = 5.58406251008641721012
        pol_damp_c(1,1,1) = 2.8904937083E+00
        pol_damp_k(1,1,1) = 4
        ! O - Si
        ldip_damping(1,1,2) = .true.
        pol_damp_b(1,1,2) = 3.66418041600447487325
        pol_damp_c(1,1,2) = 1.4460000000E+00
        pol_damp_k(1,1,2) = 4
        ! O - B 
        ldip_damping(1,1,3) = .true.
        pol_damp_b(1,1,3) = 4.48243215408077070621
        pol_damp_c(1,1,3) = 1.3650000000E+00
        pol_damp_k(1,1,3) = 4

        !O-O               eV                  hartree
        Abmhftd(1,1)  = 1.26277351E+04 ! 4.6406074534E+02
        Bbmhftd(1,1)  = 5.02221846716694036211
        Cbmhftd(1,1)  = 15.177181     ! 0.55774986759237149837 
        Dbmhftd(1,1)  = 82.256749     ! 3.02287057536699255464
        BDbmhftd(1,1) = 1.889
        !O-Si               eV                  hartree
        Abmhftd(1,2)  = 1.31294927E+03 !4.8250000000E+01
        Bbmhftd(1,2)  = 3.1860795159
        Cbmhftd(1,2)  = 1.1950537 !0.04391731240884814948
        Dbmhftd(1,2)  = 4.1831138 !0.153726127713925
        BDbmhftd(1,2) = 4.15739913110358159935
        !O-B               eV                  hartree
        Abmhftd(1,3)  = 6.51712643E+02 !2.3950000000E+01
        Bbmhftd(1,3)  = 3.251153569410612
        Cbmhftd(1,3)  = 0.0
        Dbmhftd(1,3)  = 0.0
        BDbmhftd(1,3) = 0.0


&end
