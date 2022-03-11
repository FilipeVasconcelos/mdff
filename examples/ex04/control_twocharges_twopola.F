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
	lsurf=.true.
&end
&fieldtag
	qch =  1.0 -1.0 , 
	kES = 11 11 11,            ! Parameters from [1]
	alphaES = 0.8              ! from [1]
        min_scf_pol_iter = 3
        max_scf_pol_iter = 3
        conv_tol_ind = 1e-6
	!algo_moment_from_pola = 'scf'
	algo_moment_from_pola = 'scf_kO_v2'


	lpolar(3) =.true.
	pol(3,1,1) = 0.1d0

	lpolar(4) =.true.
	pol(4,1,1) = 0.1d0
	pol(4,2,2) = 0.1d0
	pol(4,3,3) = 0.1d0
	ldip_damping = .false.

        doefg=.true.               ! do EFG
        doefield=.true.            ! do Electric Field
        lwrite_efg=.true.          ! and write ...
        lwrite_ef=.true.
&end
&mdtag
&end
