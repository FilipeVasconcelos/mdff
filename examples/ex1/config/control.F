!#Example 1: LJ fcc structure at low temperature
! ======================================================================
! main parameters
! ======================================================================
&controltag
        calc          = 'md'      ! md calculation
        lnmlj         = .true.    ! switch on lj potential
        cutshortrange = 8.5d0     ! lj short range cut off
        lvnlist       = .false.
&end
! ======================================================================
! parameters for the molecular dynamics 
! ======================================================================
&mdtag
        integrator    = 'nve-vv'  ! NVE ensemble velocity verlet
        npas          = 5000      ! number of md steps
        dt            = 0.001     ! time step    [ps]
        temp          = 90.0d0    ! temperature  [K]
        nprint        = 100       ! printing info each nprint
        fprint        = 1         ! printing info to OSZIFF
        nequil        = 1000      ! equilibration for nequil steps
        nequil_period = 1         ! period of rescaling 
&end
! ======================================================================
! parameters for the force-field 
! ======================================================================
&fieldtag
        epslj(1,1)    = 0.010323576d0 ! lennard-jones potential in [eV] = 119.8 K
        sigmalj(1,1)  = 3.405d0       ! lennard-jones potential in [A]
        mass          = 39.95         ! in atomic mass unit
&end
