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

        !INCLUDE <IONS.POT>

&end

&non_bondedtag
        ctrunc = 'notrunc'
        !INCLUDE <BMHFTD.POT>  
        !INCLUDE <BMHFT.POT>  
&end


&coulombtag
        lautoES = .true.
        epsw = __EPSW__ 
&end

&pimtag
       ! ewald summation
       algo_moment_from_pola='scf'
       omegakO       = 0.7
       conv_tol_ind  = 1e-11
       algo_ext_dipole = 'aspc'
       extrapolate_order = 2

       !INCLUDE <PIM.POT>
       !INCLUDE <PIMD.POT>
&end
