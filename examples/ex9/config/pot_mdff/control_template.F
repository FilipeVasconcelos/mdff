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
&end


&coulombtag
        lautoES = .true.
        epsw = __EPSW__ 
&end

&pimtag
       !INCLUDE <PIM.POT>
&end
