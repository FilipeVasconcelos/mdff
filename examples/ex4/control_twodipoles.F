! field calculation
&controltag
        calc='md'
        lcoulomb  = .true.
        longrange ='ewald'
	cutlongrange = 5.0d0
        lstatic=.true.
        restart_data='rvf'
	lreduced=.true.            ! reduced units 
	lreducedN=.false.          ! NOT REDUCED BY N
&end
&fieldtag
        dip(1,1) =  1.0 
        dip(2,1) =  1.0 
	kES = 11 11 11,            ! Parameters from [1]
	alphaES = 0.8              ! from [1]

        doefg=.true.               ! do EFG
        doefield=.true.            ! do Electric Field
        lwrite_efg=.true.          ! and write ...
        lwrite_ef=.true.
&end
&mdtag
&end
