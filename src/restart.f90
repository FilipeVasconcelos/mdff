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
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
! USA.

! ===== fmV =====

! ======= Hardware =======
#include "symbol.h"
!#define debug 
!#define GFORTRAN
! ======= Hardware =======


MODULE restart

  USE io,  ONLY :  ionode, stdout

CONTAINS

SUBROUTINE restart_init ( MDFF ) 

  USE io,       ONLY :  kunit_RESTART
  USE cell
  USE control
  USE config
  USE field
  USE tt_damp
  USE md

  implicit none

  ! global
  character(len=80)  :: MDFF

  ! local 
  integer :: i
  logical           :: allowed
  character(len=60) :: cpos
  character(len=20) :: FMT 

  CALL print_RESTART_info ( stdout )
 
  ! default values
  itime = itime0

  OPEN(kunit_RESTART, FILE='RESTART', form ='unformatted')
  ! reading RESTART (control parameters)
  READ( kunit_RESTART   ) lnmlj          
  READ( kunit_RESTART   ) lcoulomb       
  READ( kunit_RESTART   ) lmorse         
  READ( kunit_RESTART   ) lbmhft         
  READ( kunit_RESTART   ) lbmhftd        
  READ( kunit_RESTART   ) lsurf          
  READ( kunit_RESTART   ) lcsvr          
  READ( kunit_RESTART   ) lharm          
  READ( kunit_RESTART   ) ltraj          
  READ( kunit_RESTART   ) lvnlist        
  READ( kunit_RESTART   ) lstatic        
  READ( kunit_RESTART   ) lreduced       
  READ( kunit_RESTART   ) lreducedN       
  READ( kunit_RESTART   ) ltest          
  READ( kunit_RESTART   ) lmsd           
  READ( kunit_RESTART   ) lvacf          
  READ( kunit_RESTART   ) lrestart       
  READ( kunit_RESTART   ) cutlongrange   
  READ( kunit_RESTART   ) cutshortrange  
  READ( kunit_RESTART   ) calc           
  READ( kunit_RESTART   ) dgauss         
  READ( kunit_RESTART   ) longrange      
  READ( kunit_RESTART   ) itraj_start    
  READ( kunit_RESTART   ) itraj_period   
  READ( kunit_RESTART   ) itraj_format   
  READ( kunit_RESTART   ) trajff_data    
  READ( kunit_RESTART   ) iscff_format   
  READ( kunit_RESTART   ) iscff_data     
  READ( kunit_RESTART   ) iefall_format  
  READ( kunit_RESTART   ) iefgall_format 
  READ( kunit_RESTART   ) idipall_format 
  READ( kunit_RESTART   ) restart_data   
  READ( kunit_RESTART   ) skindiff

   CALL control_print_info( stdout , MDFF )

  ! ... continue reading RESTART (config info)  
  read( kunit_RESTART ) natm
  read( kunit_RESTART ) system
  read( kunit_RESTART ) simu_cell%A 
  read( kunit_RESTART ) ntype 
  READ( kunit_RESTART   ) ntype
  READ( kunit_RESTART   ) atypei 
  IF ( ionode ) WRITE ( stdout      ,'(A,20A3)' ) 'found type information on RESTART : ', atypei ( 1:ntype )
  READ( kunit_RESTART   ) natmi
  IF ( ionode ) WRITE ( stdout      ,'(A,20i4)' ) 'found type information on RESTART : ', natmi ( 1:ntype )

  CALL lattice ( simu_cell )
  rho = DBLE ( natm ) / simu_cell%omega

  CALL config_alloc

  ! ===============================
  ! config info
  ! ===============================
  CALL config_print_info(stdout)

  ! ... continue reading restart (config parameters)
  READ( kunit_RESTART   ) rx,ry,rz
  READ( kunit_RESTART   ) vx,vy,vz 
  READ( kunit_RESTART   ) fx,fy,fz


  READ( kunit_RESTART   ) fxs,fys,fzs
  READ( kunit_RESTART   ) rxs,rys,rzs
  READ( kunit_RESTART   ) xs,ys,zs
  READ( kunit_RESTART   ) atype 
  READ( kunit_RESTART   ) itype 
  READ( kunit_RESTART   ) verlet_vdw%list
  READ( kunit_RESTART   ) verlet_vdw%point
  READ( kunit_RESTART   ) verlet_coul%list
  READ( kunit_RESTART   ) verlet_coul%point
  READ( kunit_RESTART   ) verlet_vdw%cut
  READ( kunit_RESTART   ) verlet_coul%cut
  READ( kunit_RESTART   ) qia
  READ( kunit_RESTART   ) massia
  READ( kunit_RESTART   ) quadia_nuc
  READ( kunit_RESTART   ) dipia
  READ( kunit_RESTART   ) quadia 
  READ( kunit_RESTART   ) poldipia
  READ( kunit_RESTART   ) polquadia
  READ( kunit_RESTART   ) invpoldipia
  READ( kunit_RESTART   ) ipolar

  ! ... continue reading restart (md parameters)
  READ( kunit_RESTART   ) dt 
  READ( kunit_RESTART   ) npas
  READ( kunit_RESTART   ) itime
  READ( kunit_RESTART   ) nequil
  READ( kunit_RESTART   ) nequil_period
  READ( kunit_RESTART   ) spas
  READ( kunit_RESTART   ) npropr
  READ( kunit_RESTART   ) npropr_start
  READ( kunit_RESTART   ) nprint
  READ( kunit_RESTART   ) fprint
  READ( kunit_RESTART   ) updatevnl
  READ( kunit_RESTART   ) setvel 
  READ( kunit_RESTART   ) temp
  READ( kunit_RESTART   ) press 
  READ( kunit_RESTART   ) timesca_thermo 
  READ( kunit_RESTART   ) timesca_baro
  READ( kunit_RESTART   ) tauTberendsen
  READ( kunit_RESTART   ) tauPberendsen 
  READ( kunit_RESTART   ) taucsvr 
  READ( kunit_RESTART   ) nuandersen 
  READ( kunit_RESTART   ) annealing 
  READ( kunit_RESTART   ) first_time_xe0 
  ! thermostat
  READ( kunit_RESTART   ) integrator
  READ( kunit_RESTART   ) nhc_n
  READ( kunit_RESTART   ) nhc_yosh_order 
  READ( kunit_RESTART   ) nhc_mults
  print*,'here b'
  CALL extended_coordinates_alloc
  READ( kunit_RESTART   ) vxi
  print*,'here a'
  READ( kunit_RESTART   ) xi
  READ( kunit_RESTART   ) vxib
  READ( kunit_RESTART   ) xib
  READ( kunit_RESTART   ) ve
  READ( kunit_RESTART   ) xe
  READ( kunit_RESTART   ) xe0
  if ( ionode ) then
#ifdef GFORTRAN
    write(FMT,*) nhc_n
    write(*,'(a,' // ADJUSTL(FMT) //'f)') 'restart read ',vxi
    write(*,'(a,' // ADJUSTL(FMT) //'f)') 'restart read ',xi
#else
    write(*,'(a,<nhc_n>f)') 'restart read ',vxi
    write(*,'(a,<nhc_n>f)') 'restart read ',xi
#endif
  endif
  print*,'here'
  print*,'here',vxi
  itime0=itime
  itime1=itime0+npas
  
  ! ===================
  !  print mdtag info
  ! ===================
  CALL md_print_info(stdout) 

  ! ... continue reading RESTART (field parameters)

  READ( kunit_RESTART   ) lKA 
  READ( kunit_RESTART   ) ctrunc 
  READ( kunit_RESTART   ) symmetric_pot 
  READ( kunit_RESTART   ) ncelldirect 
  READ( kunit_RESTART   ) kES 
  READ( kunit_RESTART   ) alphaES 
  READ( kunit_RESTART   ) qlj 
  READ( kunit_RESTART   ) plj 
  READ( kunit_RESTART   ) sigmalj 
  READ( kunit_RESTART   ) epslj 
  READ( kunit_RESTART   ) sigmamor 
  READ( kunit_RESTART   ) epsmor 
  READ( kunit_RESTART   ) rhomor 
  READ( kunit_RESTART   ) mass 
  READ( kunit_RESTART   ) doefield 
  READ( kunit_RESTART   ) doefg
  READ( kunit_RESTART   ) qch           
  READ( kunit_RESTART   ) quad_nuc      
  READ( kunit_RESTART   ) dip           
  READ( kunit_RESTART   ) quad          
  READ( kunit_RESTART   ) poldip        
  READ( kunit_RESTART   ) poldip_iso    
  READ( kunit_RESTART   ) polquad       
  READ( kunit_RESTART   ) polquad_iso   
  READ( kunit_RESTART   ) pol_damp_b 
  READ( kunit_RESTART   ) pol_damp_c 
  READ( kunit_RESTART   ) pol_damp_k
  READ( kunit_RESTART   ) extrapolate_order
  READ( kunit_RESTART   ) conv_tol_ind
  READ( kunit_RESTART   ) min_scf_pol_iter
  READ( kunit_RESTART   ) max_scf_pol_iter
  READ( kunit_RESTART   ) algo_ext_dipole 
  READ( kunit_RESTART   ) thole_functions
  READ( kunit_RESTART   ) thole_function_type
  READ( kunit_RESTART   ) thole_param
  READ( kunit_RESTART   ) omegakO 
  READ( kunit_RESTART   ) epsw
  READ( kunit_RESTART   ) lautoES 
  READ( kunit_RESTART   ) lwrite_dip_wfc
  READ( kunit_RESTART   ) lwrite_dip
  READ( kunit_RESTART   ) lwrite_quad
  READ( kunit_RESTART   ) lwrite_ef
  READ( kunit_RESTART   ) lwrite_efg
  READ( kunit_RESTART   ) ldip_wfc
  READ( kunit_RESTART   ) rcut_wfc
  READ( kunit_RESTART   ) ldip_polar 
  READ( kunit_RESTART   ) ldip_damping
  READ( kunit_RESTART   ) lquad_polar 
  READ( kunit_RESTART   ) Abmhftd 
  READ( kunit_RESTART   ) Bbmhftd 
  READ( kunit_RESTART   ) Cbmhftd 
  READ( kunit_RESTART   ) Dbmhftd 
  READ( kunit_RESTART   ) BDbmhftd 
  READ( kunit_RESTART   ) task_coul
  READ( kunit_RESTART   ) algo_moment_from_pola 
   
  ! ================================
  ! initialize constant parameters
  ! ================================
  if ( non_bonded )    then
    CALL initialize_param_non_bonded
  endif
  if ( lcoulomb ) then
    CALL initialize_coulomb
  endif
  CALL get_TT_damp
  READ( kunit_RESTART   ) mu_t 

  ! ================================
  !  pint field info
  ! ================================
  CALL field_print_info(stdout,quiet=.false.)

  !allocate ( dipia_ind_t ( extrapolate_order+1, natm , 3 ) )

  CLOSE(kunit_RESTART)
#ifdef debug
   WRITE(stdout,'(a)') 'print config inside read RESTART'
   CALL print_config_sample(0,0)
#endif

  WRITE(stdout,'(a)') 'end of restart_init'

  return

END SUBROUTINE restart_init


SUBROUTINE write_RESTART

  USE io,  ONLY :  kunit_RESTART
  USE cell
  USE control
  USE config
  USE field
  USE md

  implicit none

  OPEN(kunit_RESTART, FILE='RESTART', form ='unformatted')
  ! reading RESTART (control parameters)
  WRITE( kunit_RESTART   ) lnmlj          
  WRITE( kunit_RESTART   ) lcoulomb       
  WRITE( kunit_RESTART   ) lmorse         
  WRITE( kunit_RESTART   ) lbmhft         
  WRITE( kunit_RESTART   ) lbmhftd        
  WRITE( kunit_RESTART   ) lsurf          
  WRITE( kunit_RESTART   ) lcsvr          
  WRITE( kunit_RESTART   ) lharm          
  WRITE( kunit_RESTART   ) ltraj          
  WRITE( kunit_RESTART   ) lvnlist        
  WRITE( kunit_RESTART   ) lstatic        
  WRITE( kunit_RESTART   ) lreduced       
  WRITE( kunit_RESTART   ) lreducedN       
  WRITE( kunit_RESTART   ) ltest          
  WRITE( kunit_RESTART   ) lmsd           
  WRITE( kunit_RESTART   ) lvacf          
  WRITE( kunit_RESTART   ) lrestart       
  WRITE( kunit_RESTART   ) cutlongrange   
  WRITE( kunit_RESTART   ) cutshortrange  
  WRITE( kunit_RESTART   ) calc           
  WRITE( kunit_RESTART   ) dgauss         
  WRITE( kunit_RESTART   ) longrange      
  WRITE( kunit_RESTART   ) itraj_start    
  WRITE( kunit_RESTART   ) itraj_period   
  WRITE( kunit_RESTART   ) itraj_format   
  WRITE( kunit_RESTART   ) trajff_data    
  WRITE( kunit_RESTART   ) iscff_format   
  WRITE( kunit_RESTART   ) iscff_data     
  WRITE( kunit_RESTART   ) iefall_format  
  WRITE( kunit_RESTART   ) iefgall_format 
  WRITE( kunit_RESTART   ) idipall_format 
  WRITE( kunit_RESTART   ) restart_data   
  WRITE( kunit_RESTART   ) skindiff
  ! ... continue reading RESTART (config info)  
  WRITE( kunit_RESTART ) natm
  WRITE( kunit_RESTART ) system
  WRITE( kunit_RESTART ) simu_cell%A 
  WRITE( kunit_RESTART ) ntype 
  WRITE( kunit_RESTART   ) ntype
  WRITE( kunit_RESTART   ) atypei 
  WRITE( kunit_RESTART   ) natmi
  ! ... continue reading restart (config parameters)
  WRITE( kunit_RESTART   ) rx,ry,rz
  WRITE( kunit_RESTART   ) vx,vy,vz 
  WRITE( kunit_RESTART   ) fx,fy,fz
  WRITE( kunit_RESTART   ) fxs,fys,fzs
  WRITE( kunit_RESTART   ) rxs,rys,rzs
  WRITE( kunit_RESTART   ) xs,ys,zs
  WRITE( kunit_RESTART   ) atype 
  WRITE( kunit_RESTART   ) itype 
  WRITE( kunit_RESTART   ) verlet_vdw%list
  WRITE( kunit_RESTART   ) verlet_vdw%point
  WRITE( kunit_RESTART   ) verlet_coul%list
  WRITE( kunit_RESTART   ) verlet_coul%point
  WRITE( kunit_RESTART   ) verlet_vdw%cut
  WRITE( kunit_RESTART   ) verlet_coul%cut
  WRITE( kunit_RESTART   ) qia
  WRITE( kunit_RESTART   ) massia
  WRITE( kunit_RESTART   ) quadia_nuc
  WRITE( kunit_RESTART   ) dipia
  WRITE( kunit_RESTART   ) quadia 
  WRITE( kunit_RESTART   ) poldipia
  WRITE( kunit_RESTART   ) polquadia
  WRITE( kunit_RESTART   ) invpoldipia
  WRITE( kunit_RESTART   ) ipolar
  ! ... continue reading restart (md parameters)
  WRITE( kunit_RESTART   ) dt 
  WRITE( kunit_RESTART   ) npas
  WRITE( kunit_RESTART   ) itime
  WRITE( kunit_RESTART   ) nequil
  WRITE( kunit_RESTART   ) nequil_period
  WRITE( kunit_RESTART   ) spas
  WRITE( kunit_RESTART   ) npropr
  WRITE( kunit_RESTART   ) npropr_start
  WRITE( kunit_RESTART   ) nprint
  WRITE( kunit_RESTART   ) fprint
  WRITE( kunit_RESTART   ) updatevnl
  WRITE( kunit_RESTART   ) setvel 
  WRITE( kunit_RESTART   ) temp
  WRITE( kunit_RESTART   ) press 
  WRITE( kunit_RESTART   ) timesca_thermo 
  WRITE( kunit_RESTART   ) timesca_baro
  WRITE( kunit_RESTART   ) tauTberendsen
  WRITE( kunit_RESTART   ) tauPberendsen 
  WRITE( kunit_RESTART   ) taucsvr 
  WRITE( kunit_RESTART   ) nuandersen 
  WRITE( kunit_RESTART   ) annealing 
  WRITE( kunit_RESTART   ) first_time_xe0 
  ! thermostat
  WRITE( kunit_RESTART   ) integrator
  WRITE( kunit_RESTART   ) nhc_n
  WRITE( kunit_RESTART   ) nhc_yosh_order 
  WRITE( kunit_RESTART   ) nhc_mults
  WRITE( kunit_RESTART   ) vxi
  WRITE( kunit_RESTART   ) xi
  WRITE( kunit_RESTART   ) vxib
  WRITE( kunit_RESTART   ) xib
  WRITE( kunit_RESTART   ) ve
  WRITE( kunit_RESTART   ) xe
  WRITE( kunit_RESTART   ) xe0
  ! ... continue reading RESTART (field parameters)
  WRITE( kunit_RESTART   ) lKA 
  WRITE( kunit_RESTART   ) ctrunc 
  WRITE( kunit_RESTART   ) symmetric_pot 
  WRITE( kunit_RESTART   ) ncelldirect 
  WRITE( kunit_RESTART   ) kES 
  WRITE( kunit_RESTART   ) alphaES 
  WRITE( kunit_RESTART   ) qlj 
  WRITE( kunit_RESTART   ) plj 
  WRITE( kunit_RESTART   ) sigmalj 
  WRITE( kunit_RESTART   ) epslj 
  WRITE( kunit_RESTART   ) sigmamor 
  WRITE( kunit_RESTART   ) epsmor 
  WRITE( kunit_RESTART   ) rhomor 
  WRITE( kunit_RESTART   ) mass 
  WRITE( kunit_RESTART   ) doefield 
  WRITE( kunit_RESTART   ) doefg
  WRITE( kunit_RESTART   ) qch           
  WRITE( kunit_RESTART   ) quad_nuc      
  WRITE( kunit_RESTART   ) dip           
  WRITE( kunit_RESTART   ) quad          
  WRITE( kunit_RESTART   ) poldip        
  WRITE( kunit_RESTART   ) poldip_iso    
  WRITE( kunit_RESTART   ) polquad       
  WRITE( kunit_RESTART   ) polquad_iso   
  WRITE( kunit_RESTART   ) pol_damp_b 
  WRITE( kunit_RESTART   ) pol_damp_c 
  WRITE( kunit_RESTART   ) pol_damp_k
  WRITE( kunit_RESTART   ) extrapolate_order
  WRITE( kunit_RESTART   ) conv_tol_ind
  WRITE( kunit_RESTART   ) min_scf_pol_iter
  WRITE( kunit_RESTART   ) max_scf_pol_iter
  WRITE( kunit_RESTART   ) algo_ext_dipole 
  WRITE( kunit_RESTART   ) thole_functions
  WRITE( kunit_RESTART   ) thole_function_type
  WRITE( kunit_RESTART   ) thole_param
  WRITE( kunit_RESTART   ) omegakO 
  WRITE( kunit_RESTART   ) epsw
  WRITE( kunit_RESTART   ) lautoES 
  WRITE( kunit_RESTART   ) lwrite_dip_wfc
  WRITE( kunit_RESTART   ) lwrite_dip
  WRITE( kunit_RESTART   ) lwrite_quad
  WRITE( kunit_RESTART   ) lwrite_ef
  WRITE( kunit_RESTART   ) lwrite_efg
  WRITE( kunit_RESTART   ) ldip_wfc
  WRITE( kunit_RESTART   ) rcut_wfc
  WRITE( kunit_RESTART   ) ldip_polar 
  WRITE( kunit_RESTART   ) ldip_damping
  WRITE( kunit_RESTART   ) lquad_polar 
  WRITE( kunit_RESTART   ) Abmhftd 
  WRITE( kunit_RESTART   ) Bbmhftd 
  WRITE( kunit_RESTART   ) Cbmhftd 
  WRITE( kunit_RESTART   ) Dbmhftd 
  WRITE( kunit_RESTART   ) BDbmhftd 
  WRITE( kunit_RESTART   ) task_coul
  WRITE( kunit_RESTART   ) algo_moment_from_pola 
  WRITE( kunit_RESTART   ) mu_t 
  CLOSE(kunit_RESTART)



!  if ( ionode ) then
!    write(*,'(a,<nhc_n>f)') 'restart write ',vxi
!    write(*,'(a,<nhc_n>f)') 'restart write ',xi
!  endif

  return

END SUBROUTINE write_RESTART

SUBROUTINE print_RESTART_info ( kunit ) 

  implicit none

  ! local 
  integer :: kunit

  if ( ionode ) call dumb_guy(kunit)
  if ( ionode ) then
     WRITE ( kunit ,'(a)')       "       _____  ______  _____ _______       _____ _______ "
     WRITE ( kunit ,'(a)')       "      |  __ \|  ____|/ ____|__   __|/\   |  __ \__   __|"
     WRITE ( kunit ,'(a)')       "      | |__) | |__  | (___    | |  /  \  | |__) | | |   "
     WRITE ( kunit ,'(a)')       "      |  _  /|  __|  \___ \   | | / /\ \ |  _  /  | |   "
     WRITE ( kunit ,'(a)')       "      | | \ \| |____ ____) |  | |/ ____ \| | \ \  | |   "
     WRITE ( kunit ,'(a)')       "      |_|  \_\______|_____/   |_/_/    \_\_|  \_\ |_|   "
     blankline(kunit)
     separator(kunit)
     blankline(kunit)
     WRITE ( kunit ,'(a)')       "WARING : the full restart mode is only reading parameters in"
     WRITE ( kunit ,'(a)')       "the local RESTART file. The control file is not used at all" 
     WRITE ( kunit ,'(a)')       "(except full_restart tag;)"
     WRITE ( kunit ,'(a)')       "Be carefull !! "  
     WRITE ( kunit ,'(a)')       "Some of the most important parameters read from the RESTART "
     WRITE ( kunit ,'(a)')       "file are reminded below. But not all of them !!!"

     blankline(kunit)
  endif

  return

END SUBROUTINE print_RESTART_info

END MODULE restart
