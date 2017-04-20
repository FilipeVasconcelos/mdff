! MDFF parallel Molecular Dynamics ... For Fun
! Copyright (C) 2011  F. Vasconcelos
!
! This program is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
! ===== fmV =====

! ======= Hardware =======
#include "symbol.h"
! ======= Hardware =======

! *********************** MODULE md ********************************************
!> \brief 
!! module related to molecular dynamics calculation ( calc = 'md' )
! ******************************************************************************
MODULE md

  USE io,                       ONLY :  ionode ,stdin, stdout
  USE constants,                ONLY :  dp

  implicit none

  logical, SAVE :: lleapequi               !< leap-frog used in the equilibration part together with verlet -> vv + lf 

  integer :: npas                          !< number of time steps
  integer :: itime                         !< current time step
  integer :: itime0                        !< starting time step
  integer :: itime1                        !< finish time step
  integer :: nequil                        !< number of equilibration steps
  integer :: nequil_period                 !< equilibration period
  integer :: spas                          !< save configuration each spas step 
  integer :: npropr                        !< period to calculate on-the fly property 
  integer :: npropr_start                  !<        starting from npropr_start
  integer :: nprint                        !< print thermo info to standard output
  integer :: fprint                        !< print thermo info to file OSZIFF
  integer :: updatevnl                     !< number of verlet list update  

  real(kind=dp) :: dt                   !< time step
  real(kind=dp) :: temp                 !< temperature in input in Kelvin but it's kBT everywhere else
  real(kind=dp) :: press                !< pressure
  real(kind=dp) :: timesca_thermo       !< Nose-Hoover Chain : time scale of thermostat 
  real(kind=dp) :: timesca_baro         !< Nose-Hoover Chain : time scale of barostat 
  integer       :: nhc_yosh_order       !< Nose-Hoover Chain : order of the yoshida integrator 
  integer       :: nhc_mults            !< Nose-Hoover Chain : number of multiple timesteps 
  integer       :: nhc_n                !< Nose-Hoover Chain : length of the Nose-Hoover chain
  real(kind=dp) :: tauTberendsen        !< characteristic time in berendsen thermostat (simple rescale if tauTberendsen = dt )
  real(kind=dp) :: tauPberendsen        !< characteristic time in berendsen barostat   (simple rescale if tauPberendsen = dt )
  real(kind=dp) :: taucsvr              !< characteristic time in Stochastic velocity rescaling (simple rescale if taucsvr = 0.0 )
  real(kind=dp) :: nuandersen           !< characteristic frequency in andersen thermostat ( to be merged with timesca_thermo )
  real(kind=dp) :: annealing            !< velocity rescaling
  logical       :: first_time_xe0
  real(kind=dp), dimension(:)    , allocatable :: vxi , xi           !< thermostat coordinates coupled to the particules (Nose-Hoover Chain : nhcn )
  real(kind=dp), dimension(:)    , allocatable :: vxib, xib          !< thermostat coordinates coupled to the volume (Nose-Hoover Chain : nhcnp )
  real(kind=dp)                                :: ve, xe, xe0        !< coordinates of the barostat (Andersen) 

  ! ================================================
  !     algorithm for dynamic integration
  ! ================================================
  integer           :: yosh_allowed(5)
  data                 yosh_allowed / 1, 3 , 5 , 7 , 9 / 

  character(len=60) :: integrator               !< integration method   
  character(len=60) :: integrator_allowed(9)    
  data                 integrator_allowed / 'nve-vv'  , 'nve-lf'   , 'nve-be' ,  'nve-lfq',  &
                                            'nvt-and' , 'nvt-nhc2' , 'nvt-nhcn', & 
                                            'npe-vv'  , 'npt-nhcnp' /

  character(len=60) :: nve_ensemble(4)
  data                 nve_ensemble / 'nve-vv' , 'nve-lf', 'nve-be' , 'nve-lfq' /
  character(len=60) :: nvt_ensemble(4)
  data                 nvt_ensemble / 'nvt-and' , 'nvt-nh' , 'nvt-nhc2' , 'nvt-nhcn' /
  character(len=60) :: npe_ensemble(1) 
  data                 npe_ensemble / 'npe-vv' /
  character(len=60) :: npt_ensemble(1) 
  data                 npt_ensemble / 'npt-nhcnp' /

  ! ================================================
  !  velocity distribution (Maxwell-Boltzmann or uniform)
  ! ================================================
  character(len=60) :: setvel                   !< velocity distribution 
  character(len=60) :: setvel_allowed(2) 
  data setvel_allowed / 'MaxwBoltz', 'Uniform' /

CONTAINS


! *********************** SUBROUTINE md_init ***********************************
!> \brief
!! molecular dynamics initialisation of main parameters
! ******************************************************************************
SUBROUTINE md_init

  USE control,  ONLY :  calc, full_restart

  implicit none

  integer            :: ioerr
  character(len=132) :: filename

  namelist /mdtag/    integrator    , & 
                      setvel        , & 
                      npas          , & 
                      nequil        , &
                      nequil_period , & 
                      annealing     , &
                      npropr        , & 
                      npropr_start  , & 
                      nprint        , & 
                      nprint        , & 
                      fprint        , & 
                      spas          , & 
                      dt            , &
                      temp          , & 
                      press         , & 
                      nuandersen    , & 
                      taucsvr       , &
                      tauTberendsen , &
                      tauPberendsen , &
                      nhc_yosh_order, & 
                      nhc_mults     , & 
                      nhc_n         , & 
                      timesca_thermo, &      
                      timesca_baro        

  if ( calc .ne. 'md' ) return
  ! ======================
  !  mdtag default values
  ! ======================
  CALL md_default_tag 
  ! =====================
  !  read mdtag namelist
  ! =====================
  CALL getarg (1,filename)
  OPEN ( stdin , file = filename)
  READ ( stdin , mdtag ,iostat =ioerr)
  if ( ioerr .lt. 0 )  then
    io_node &
    WRITE ( stdout, '(a)') 'ERROR reading input_file : mdtag section is absent'
    STOP
  elseif ( ioerr .gt. 0 )  then
    io_node &
    WRITE ( stdout, '(a,i8)') 'ERROR reading input_file : mdtag wrong tag',ioerr
    STOP
  endif
  CLOSE  ( stdin )
  ! ======================
  !  check mdtag namelist
  ! ======================
  CALL md_check_tag

  if ( full_restart ) return
  CALL extended_coordinates_alloc
  ! ===================
  !  print mdtag info
  ! ===================
  CALL md_print_info(stdout)

  return

END SUBROUTINE md_init 


! *********************** SUBROUTINE md_default_tag ****************************
!> \brief
!! set default values to md tag
! ******************************************************************************
SUBROUTINE md_default_tag

  implicit none
  
  ! =================
  !  default values
  ! =================
  itime         = 1 
  lleapequi     = .false.
  integrator    = 'nve-vv'
  setvel        = 'MaxwBoltz'
  npas          = 10
  nequil        = 0
  nequil_period = 1
  nprint        = 1             
  fprint        = 1
  spas          = 1000           
  dt            = 0.0_dp
  temp          = 1.0_dp
  press         = 0.0_dp
  taucsvr       = 0.0_dp
  tauTberendsen = 0.0_dp
  tauPberendsen = 0.0_dp
  nhc_yosh_order= 3 
  nhc_mults     = 2 
  nhc_n         = 4
  annealing     = 1.0_dp
  npropr        = 1
  npropr_start  = 0
  timesca_thermo= 1.0_dp
  timesca_baro  = 1.0_dp

  first_time_xe0 = .true.

  return

END SUBROUTINE md_default_tag


! *********************** SUBROUTINE md_check_tag ******************************
!> \brief
!! check md tag values
! ******************************************************************************
SUBROUTINE md_check_tag

  USE constants,        ONLY :  boltz_unit , time_unit, press_unit
  USE control,          ONLY :  lstatic, lmsd, lvacf

  implicit none

  ! local
  integer :: i

  ! ===========
  ! properties on-the-fly
  ! ===========
  if ( ( lmsd .or. lvacf) .and. npropr .eq. 0) then
    if ( ionode ) WRITE ( stdout , '(a)' ) 'ERROR controltag: npropr need to be defined for lmsd, lvacf ... properties on-the-fly '
  endif

  if (dt.eq.0.0_dp .and. .not. lstatic ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR mdtag: timestep dt is zero'
    STOP
  endif

  ! =========================================
  !  scaling velocities berendsen, 
  !  if not defined simple velocity/volume rescale 
  ! =========================================
  if (tauTberendsen.eq.0.0_dp ) tauTberendsen = dt
  if (tauPberendsen.eq.0.0_dp ) tauPberendsen = dt

  ! ====================
  !  check integrator
  ! ====================
  CALL check_allowed_tags( size( integrator_allowed ), integrator_allowed, integrator, 'in mdtag', 'integrator' ) 
  ! ===============
  !  check setvel
  ! ===============
  CALL check_allowed_tags( size( setvel_allowed ), setvel_allowed, setvel, 'in mdtag', 'setvel' ) 

  ! ===================
  !  check timesca_thermo/baro 
  ! ===================
  if ( any ( integrator .eq. nvt_ensemble ) .and. timesca_thermo .eq. 0.0_dp ) then
     if ( ionode )  WRITE ( stdout ,'(a,f10.5)') 'ERROR mdtag: with integrator in nvt_ensemble timesca_thermo should be set : ',timesca_thermo
    STOP
  endif
  if ( any ( integrator .eq. npt_ensemble ) .and. timesca_baro .eq. 0.0_dp ) then
     if ( ionode )  WRITE ( stdout ,'(a,f10.5)') 'ERROR mdtag: with integrator in npt_ensemble timesca_baro should be set : ',timesca_baro
    STOP
  endif

  if (integrator.eq.'nvt-nh' ) then
    if ( ionode )  WRITE ( stdout ,'(a)') 'ERROR mdtag: integrator = "nvt-nh" not yet implemented try nvt-nhc2 '
    if ( ionode )  WRITE ( stdout ,'(a)') integrator
    STOP 
  endif

  ! ============================
  !  check if leap frog is 
  !  wanted with equilibration 
  !  then we use nve-vv for the 
  !  equilibration period
  ! ============================
  if ( integrator.eq.'nve-lfq' ) then
    lleapequi = .true.
    integrator = 'nve-vv' 
  endif


  !if ( nhc_n == 1 ) then
  !  write(stdout,'(a)') 'ERROR : nhc_n = 1 is not allowed with Nose Hoover chains thermostats'
  !  stop
  !endif

  ! annealing 
  if ( annealing .ne. 1.0_dp ) then
    nequil=npas
    nequil_period=1
  endif


  ! units          

  !  eV             K     eV/K
  temp           = temp * boltz_unit ! temp = kB * T
  ! eV / A**3    =  GPa     eV/A**3/GPa
  press          = press * press_unit
  ! angstrom*(atomicmassunit/eV)** 0.5  <= ps
  dt             = dt             * time_unit 
  timesca_thermo = timesca_thermo * time_unit
  timesca_baro   = timesca_baro   * time_unit
  tauTberendsen  = tauTberendsen  * time_unit
  tauPberendsen  = tauPberendsen  * time_unit

  itime0=itime
  itime1=(itime0-1)+npas

  return

END SUBROUTINE md_check_tag

SUBROUTINE extended_coordinates_alloc

  USE control,  ONLY :  calc

  implicit none

  if ( calc .ne. 'md' ) return

  ! allocation of thermostat coordinates
  if ( integrator .eq. 'nvt-nhc2' ) then
    allocate ( vxi(2) , xi(2) )
    vxi  = 0.0_dp
    xi   = 0.0_dp
    nhc_n= 2
  else if ( integrator .eq. 'nvt-nhcn' .or. integrator .eq. 'npt-nhcnp' ) then
    allocate ( vxi(nhc_n) , xi(nhc_n) )
    vxi = 0.0_dp
    xi  = 0.0_dp
    ve   = 0.0_dp
    xe   = 0.0_dp
    allocate ( vxib(nhc_n) , xib(nhc_n) )
    vxib = 0.0_dp
    xib  = 0.0_dp
  else
    !restart purpose if ensemble is NVT
    allocate ( vxi(1) , xi(1) )
    allocate ( vxib(1) , xib(1) )
    vxi  = 0.0_dp
    xi   = 0.0_dp
    vxib = 0.0_dp
    xib  = 0.0_dp
  endif


  return

END SUBROUTINE extended_coordinates_alloc

SUBROUTINE extended_coordinates_dealloc

  USE control,  ONLY :  calc

  implicit none

  if ( calc .ne. 'md' ) return

  ! Deallocation of thermostat coordinates
  if ( integrator .eq. 'nvt-nhc2' ) then
    deallocate ( vxi , xi )
  else if (  integrator .eq. 'nvt-nhcn' .or. integrator .eq. 'npt-nhcnp' ) then
    deallocate ( vxi , xi )
    if ( integrator .eq. 'npt-nhcnp' ) then
      deallocate ( vxib , xib )
    endif
  else
    deallocate (vxi,xi)
    deallocate (vxib,xib)
  endif


  return

END SUBROUTINE extended_coordinates_dealloc

! *********************** SUBROUTINE md_print_info *****************************
!> \brief
!! print general information to standard output for md control tag
! ******************************************************************************
SUBROUTINE md_print_info(kunit)

  USE constants,        ONLY :  boltz_unit , time_unit, press_unit
  USE control,          ONLY :  ltraj , lstatic , lvnlist , lreducedN , lreduced , lcsvr , itraj_start , itraj_period , itraj_format  

  implicit none
  
  !local
  integer :: kunit

  if ( ionode ) then 
                                          separator(kunit)    
                                          blankline(kunit)    
                                          WRITE ( kunit ,'(a)')       'MD MODULE ... WELCOME'
                                          blankline(kunit)    
      if (lstatic) then
                                          WRITE ( kunit ,'(a)')       'static  calculation ....boring                 '
      else
                                          WRITE ( kunit ,'(a)')       'periodic boundary conditions  '  
                                          WRITE ( kunit ,'(a)')       'using minimum image convention                 '  
        if ( lvnlist )                    WRITE ( kunit ,'(a)')       'verlet list used '  
        if ( lreduced )                   WRITE ( kunit ,'(a)')       'unit constant set to one == reduced units'
        if ( lreducedN )                  WRITE ( kunit ,'(a)')       'units reduced by the number of atom'
        if ( .not.lstatic)                WRITE ( kunit ,'(a)')       'dynamic calculation'
        if ( integrator .eq. 'nve-lf')    WRITE ( kunit ,'(a)')       'NVE ensemble --- leap-frog integrator          '
        if ( integrator .eq. 'nve-be')    WRITE ( kunit ,'(a)')       'NVE ensemble --- beeman integrator             '

        if ( integrator .eq. 'nve-vv' .and. lleapequi )   &
                                          WRITE ( kunit ,'(a)')       'NVE ensemble --- velocity verlet (equil) + leap-frog integrator        '
        if ( integrator .eq. 'nve-vv'.and.  .not. lleapequi )  &
                                          WRITE ( kunit ,'(a)')       'NVE ensemble --- velocity verlet integrator    '
        if ( integrator .eq. 'nvt-and' )  WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator    ' 
        if ( integrator .eq. 'nvt-and' )  WRITE ( kunit ,'(a)')       ' + Andersen thermostat'
        if ( integrator .eq. 'nvt-and' )  WRITE ( kunit ,'(a,f12.5)') 'nuandersen                         = ',nuandersen  
        if ( integrator .eq. 'nvt-nh' )   WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator'
        if ( integrator .eq. 'nvt-nh' )   WRITE ( kunit ,'(a)')       ' + Nose Hoover thermostat'
        if ( integrator .eq. 'nvt-nhc2' .or. &
             integrator .eq. 'nvt-nhcn' ) WRITE ( kunit ,'(a)')       'NVT ensemble --- velocity verlet integrator'
        if ( integrator .eq. 'nvt-nhc2' ) WRITE ( kunit ,'(a)')       ' + Nose Hoover chain 2 thermostat  (see Frenkel and Smit)'
        if ( integrator .eq. 'nvt-nhcn' ) WRITE ( kunit ,'(a)')       ' + Nose Hoover chain N thermostat  (see Martyna et al.)'
        if ( integrator .eq. 'nvt-nhc2' .or. & 
             integrator .eq. 'nvt-nhcn' ) WRITE ( kunit ,'(a,f12.5,a)') 'time scale thermostat: timesca_thermo = ',timesca_thermo/time_unit,' ps'
        if ( integrator .eq. 'npt-nhcnp') WRITE ( kunit ,'(a,f12.5,a)') 'time scale barostat  : timesca_baro   = ',timesca_baro/time_unit  ,' ps'
        if ( ( integrator .ne. 'nvt-and' )  .and. &
             ( integrator .ne. 'nvt-nhc2' ) .and. &
             ( integrator .ne. 'nvt-nhcn' ) .and. &
             ( integrator .ne. 'nvt-nh' ) ) then
              if ( nequil .eq. 0 ) then
                                          WRITE ( kunit ,'(a)')       'with no equilibration          '
              else
                                          WRITE ( kunit ,'(a)')       'with equilibration:             '
                                          WRITE ( kunit ,'(a)')       'berendsen scaling ( is not canonical ...and so NVE)'
               if ( integrator .eq. 'nve-lf' ) &
                                          WRITE ( kunit ,'(a)')       'WARNING WITH nve-lf no equilibration possible'
              endif !nequil
        endif !integrator
      endif !static
                                          WRITE ( kunit ,'(a,i12)')     'number of steps                       = ',npas 
                                          WRITE ( kunit ,'(a,e12.5,a)') 'timestep                              = ',dt / time_unit         ,'  ps'
                                          WRITE ( kunit ,'(a,e12.5,a)') 'time range                            = ',dt * npas / time_unit  ,'  ps'
                                          WRITE ( kunit ,'(a,f12.5,a)') 'temperature                           = ',temp / boltz_unit      ,'   K'
                                          WRITE ( kunit ,'(a,f12.5,a)') 'pressure                              = ',press / press_unit     ,' GPa'
      if ( integrator .eq. 'nve-vv' .and. nequil .ne. 0 ) then           
                                          WRITE ( kunit ,'(a,i12)')     'number of equilibration steps         = ',nequil
                                          WRITE ( kunit ,'(a,i12)')     'equilibration period                  = ',nequil_period
      endif 
      if ( nequil .ne. 0 .and.      lcsvr)WRITE ( kunit ,'(a,e12.5,a)') 'Stochastic velocity resc. time scale  = ',taucsvr/time_unit      ,'  ps'
      if ( nequil .ne. 0 .and.      lcsvr)WRITE ( kunit ,'(a,e12.5)')   'taucsvr = 0.0 -> simple rescale'
      if ( nequil .ne. 0 .and. .not.lcsvr)WRITE ( kunit ,'(a,e12.5,a)') 'Berendsen thermostat time scale       = ',tauTberendsen/time_unit,'  ps'
      if ( nequil .ne. 0 )                WRITE ( kunit ,'(a,e12.5,a)') 'Berendsen barostat time scale         = ',tauPberendsen/time_unit,'  ps'
      if ( nequil .ne. 0 .and. tauTberendsen .eq. dt )   &
                                          WRITE ( kunit ,'(a)')         'tau[T-P]berendsen = dt -> simple rescale'
                                          WRITE ( kunit ,'(a,i12)')     'print thermo  periodicity             = ',nprint
      if ( ltraj )                   then    
                                          WRITE ( kunit ,'(a,i12)')     'save trajectory from step             = ',itraj_start
                                          WRITE ( kunit ,'(a,i12)')     'saved trajectory periodicity          = ',itraj_period
      if ( itraj_format .eq. 0 )          WRITE ( kunit ,'(a,I7)')      'trajectory format                     : BINARY'
      if ( itraj_format .ne. 0 )          WRITE ( kunit ,'(a,I7)')      'trajectory format                     : FORMATTED'
                                          WRITE ( kunit ,'(a,I7)')      'starting from step                    :',itime0
                                          WRITE ( kunit ,'(a,I7)')      '           to step                    :',itime1
      endif       
                                          blankline(kunit)    
  endif !ionode

  return

END SUBROUTINE md_print_info


END MODULE md 
! ===== fmV =====
