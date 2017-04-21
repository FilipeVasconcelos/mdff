&controltag
        calc          = 'md'
        cutshortrange = __CUTOFF_SR__
        cutlongrange  = __CUTOFF_LR__
        skindiff      = 1.0d0
        lvnlist       = .true.
        lwrite_restart= __RESTART__
        posff_data  = __DATA__ 
        lbmhftd       = .true.
        lcoulomb      = .true.
        ltraj         = __LTRAJ__ 
        itraj_start   = 0 
        itraj_period  = __ITRAJ_PERIOD__
&end

&mdtag
        integrator    = __INTEGRATOR__ 
        npas          = __NPAS__
        dt            = __DT__
        nprint        = __NPRINT__
        fprint        = __FPRINT__

        !INCLUDE <thermostat>
        !INCLUDE <barostat>
        !INCLUDE <rescaling>

&end

&fieldtag
        !INCLUDE <IONS.POT>
&end
&coulombtag
        lautoES       = .true.
        epsw          = 1.e-4
&end

&pimtag        
        ! ewald summation
        algo_moment_from_pola='scf_kO_v4_1'
        omegakO       = 0.7 
        conv_tol_ind  = 1e-6
        algo_ext_dipole = 'aspc'
        extrapolate_order = 4
        !INCLUDE <PIM.POT>
&end
&non_bondedtag        
        !INCLUDE <BMHFTD.POT>  
&end

