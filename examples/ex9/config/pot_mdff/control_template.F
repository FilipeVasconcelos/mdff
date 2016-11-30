&controltag
        calc          = 'md'
        cutshortrange = __RCUT__ 
        cutlongrange  = __RCUT__
        skindiff      = 1.0d0
        lvnlist       = .true.
        posff_data    = 'rnn'
        lreduced      = .false.
        lbmhft        = __BMH__
        lbmhftd       = __FT_ON__
        lcoulomb      = __COUL__
        longrange     = 'ewald'
        lstatic       = .true.
&end

&mdtag
&end

&fieldtag

        ctrunc = 'notrunc'
        lautoES = .true.
        epsw = __EPSW__ 

        min_scf_pol_iter = 6
        max_scf_pol_iter = 100
        conv_tol_ind = __CONV_SCF__ 
        !algo_moment_from_pola = 'scf_kO_v4_1'
        algo_moment_from_pola = 'scf_kO_v1'

        !INCLUDE <IONS.POT>
        !INCLUDE <PIM.POT>
        !INCLUDE <PIMD.POT>
        !INCLUDE <BMHFTD.POT>  

&end
