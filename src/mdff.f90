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

!> \mainpage
!! ===============================================================\n
!!                                                                \n
!> \brief PARALLEL MOLECULAR DYNAMICS ...for fun                  \n
!! This code perform molecular dynamics simulation with simple    \n
!! non-bonding potentials ( Lennard-jones, Buckingham and Morse). \n
!! It permits to calculate electric-field gradient (EFG) in a     \n
!! point charge and dipole approximations using direct and Ewald  \n
!! summation methods. More details in the manual svn/doc/mdff.tex \n
!! (not complete yet) 
!> \authors
!! email : filipe.manuel.vasconcelos@gmail.com                    \n
!! Based on different codes or books by:                          \n
!! F. Affouard University of Lille                                \n
!! S. Sastry   JNCASR                                             \n 
!! D. Frenkel and B. Smit                                         \n  
!! Allen-Tildesley                                                \n
!! vasp: some features ;)                                         \n
!> \date
!! - sep2005-jun2006                                        
!! - dec2007                                                        
!! - jan2011                                                
!! - jan2013                                                      \n
!! - jan/mars2014                                                 \n
!!
!  ===============================================================


! ======= Hardware =======
#include "symbol.h"
!#define MULTRUN
!#define block
! ======= Hardware =======


PROGRAM main_MDFF

  USE constants
  USE io
  USE control 
  USE config
  USE time
  USE md 
  USE field 
  USE vib
  USE opt
  USE efg
  USE radial_distrib 
  USE rmc
  USE voisin
  USE msd
  USE vacf 
  USE stochio
  USE restart 
  USE block
  USE mpimdff
  USE dumb

  implicit none

  ! local
  integer       :: argc               ! input argument of the executable
  dectime                             ! timing declaration ( ierr also declared ) 
  !character(len=80) :: SVNREV         ! svn revision 

#ifdef MULTRUN
  integer :: irun          
  real(kind=dp) :: tempi
#endif
  ! date time version infog as from vasp 5.3 ;)
  CHARACTER (LEN=80),PARAMETER :: MDFF = &
        'mdff.'//BRANCH// ' ' // &
        ' (build ' // __DATE__// ' ' //__TIME__// ') ' // &
        'parallel'

#ifdef MPI
  ! ====================================================
  ! MPI initialisation
  ! ====================================================
  CALL MPI_INIT ( ierr )                           
#endif

  ! ====================================================
  ! time information init 
  ! ====================================================
  CALL time_init
  CALL init_dumb
#ifdef MPI
  statime 
#endif

  ! ====================================================      
  ! random number generator init
  ! MD is too sensitive to the initial velocities
  ! for test purpose keep it commented
  ! uncomment it if you want uncorrelated runs  
  ! ====================================================
  !CALL init_random_seed 

#ifdef MPI
  ! ====================================================
  ! gives rank index (myrank) 
  ! ====================================================
  CALL MPI_COMM_RANK ( MPI_COMM_WORLD , myrank   , ierr )  

  ! ====================================================
  ! gives total proc number  
  ! ====================================================
  CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , numprocs , ierr )  
#endif

  ! ====================================================
  ! input/output init 
  ! ====================================================
  CALL io_init                                            

  !CALL print_kind_info ( stdout )
  
  ! ========================================
  ! reads command line
  ! ========================================
  argc = iargc ( )
  if ( argc .lt. 1 ) then
    if (  ionode  ) WRITE ( stdout , * ) 'usage : mdff.x <input_file>'
    STOP 
  endif

  ! ========================================
  ! set some default values to control tags 
  ! parameters before to read them 
  ! and check consistency
  ! ========================================
  CALL control_init ( MDFF ) 

  ! =====================================
  ! md_init is the main subroutine  
  ! only binary mixture lennard-jones
  ! particules are avaiable different 
  ! properties EFG, GR ... 
  ! are calculated inside this routine  
  ! =====================================
  if ( .not. ltest )  then

    ! =============================================
    ! configuration initialization + force_field
    ! =============================================
    if ( full_restart ) then
      CALL restart_init ( MDFF ) 
    else
      CALL config_init 
    endif

    ! =============================
    ! md initialization
    ! =============================
    if ( .not. full_restart ) CALL md_init
 
    ! =============================
    ! properties initialization	 
    ! ( removed feb 2013)
    ! =============================
    !CALL prop_init

    ! ===========================================================
    ! parallelization: atom decomposition split atom into
    ! different procs. dec%istart , dec%iend are the index of atoms 
    ! for rank = myrank.  only for calc='md'. 
    ! for the other type of calc this is done when natm is known
    ! 
    ! k-point parallelisation for lcoulomb calculation ( electric 
    ! field, electric field gradient ...
    ! TODO make bounds global in a      !done october 2013 
    ! specific module or where it used.
    ! ===========================================================
    if ( calc .eq. 'md' )                 CALL do_split ( natm       , myrank , numprocs , atom_dec , 'atoms' )
    if ( calc .eq. 'md' .and. lcoulomb )  CALL do_split ( km_coul%nk , myrank , numprocs , km_coul%kpt_dec  , 'k-pts' )
    ! =====================================
    ! verlet list generation
    ! mmmmh only md ? 
    ! vnlist should be test for calc='opt' 
    ! =====================================
    if ( calc .eq. 'md' ) then
      verlet_vdw%listname ='vdw ' 
      verlet_coul%listname='coul' 
      if ( lvnlist .and. non_bonded )    CALL vnlist_pbc!( verlet_vdw  )    
      if ( lvnlist .and. lcoulomb   )    CALL vnlist_pbc!( verlet_coul )    
    endif  

    !IF MD
    MDLOOP: if ( calc .eq. 'md') then

    ! ===============================================
    ! some initialization of on-the-fly properties
    ! ===============================================

    ! =====================================
    ! mean square displacement
    ! =====================================
    CALL msd_init
    CALL msd_alloc 

    ! =====================================
    ! velocity auto-correlation function
    ! =====================================
    CALL vacf_init
    CALL vacf_alloc

    ! ==================================================
    ! OUTSIDE LOOP (not used yet)
    ! can be useful for some applications
    ! temperature variations or any other parameter
    ! as only one parameters will variate with time
    ! try to find a condensate way to write it
    ! ==================================================

#ifdef MULTRUN
    RUN: DO irun = 1 , 1
    tempi = temp
    RUN: DO irun = 0 , 2
        temp = ( 0.1_dp * irun  ) + tempi
        io_node blankline(stdout)   
        io_node WRITE ( stdout , '(a,f8.5)' ) 'T = ' , temp
#endif
  
         ! =======================
         !        static 
         ! =======================
           if ( lstatic ) then  
             npas = 0
             integrator = 'nve-vv' 
             CALL md_run 
         ! =======================
         ! .... or dynamic   
         ! =======================
           else               
             ! =========================
             ! velocities intialization
             ! =========================
             CALL init_velocities 
             CALL md_run 
           endif
   
#ifdef MULTRUN
        ENDDO RUN
#endif
  
    ENDIF MDLOOP
    ! ==============================================
    !          end of calc = 'md' loop
    ! ==============================================

    ! ==============================================
    !                calc != 'md' 
    ! ==============================================

    if ( calc .eq. 'dist' ) then
      CALL distance_tab
      CALL write_CONTFF
    endif 

    ! ==============================================
    ! IF OPT :
    ! - optimisation of structures read in TRAJFF 
    ! - write the optimize strutures to ISCFF
    ! ==============================================
    if ( calc .eq. 'opt' ) then  
      CALL opt_init
      CALL opt_main 
    endif
    ! ==============================================
    ! IF VIB : 
    ! - phonon related calculations
    ! - structures read in ISCFF ( optimized structures )
    ! ==============================================
    if ( any( calc .eq. vib_allowed ) ) then       
      CALL vib_init
      CALL vib_main 
    endif 
    ! ==============================================
    ! IF EFG : 
    ! - electric field gradient calculation 
    ! - structures read in TRAJFF   
    ! - main output EFGALL file
    ! - calculate from point charges and/or dipoles
    ! - polarisation or wannier function centres (wfc)
    ! ==============================================
    if ( calc .eq. 'efg' ) then
      CALL efg_init
      CALL efgcalc 
    endif
    ! ==============================================
    ! IF EFG+ACF : 
    ! - calculate EFG auto-correlation functions 
    ! - data in EFGALL    
    ! ==============================================
    if ( calc .eq. 'efg+acf' ) then
      CALL efg_init
      CALL efg_acf 
    endif
    ! ==============================================
    ! IF GR : 
    ! - radial distribution function  of structures 
    ! - structures read in TRAJFF   
    ! ==============================================
    if ( calc .eq. 'gr' ) then
      CALL grcalc 
    endif
    if ( calc .eq. 'rmc' ) then 
      CALL rmc_init
      CALL rmc_main
    endif    
    ! ==============================================
    ! IF VOIS1 : 
    ! - first neighbour sphere analysis  
    ! - ( in construction )
    ! ==============================================
    if ( calc.eq. 'vois1' ) then
      CALL vois1_init
      CALL vois1_driver
    endif
    ! ========================================================
    ! write final config pos and vel (always) only for md run
    ! ========================================================
    if ( calc .eq. 'md' ) then 
      CALL write_CONTFF
      if ( lwrite_dip ) CALL write_DIPFF ( mu_t ) 
      !CALL write_RESTART
    endif

    if ( calc .eq. 'stochio' ) then
      CALL stochio_init    
      CALL stochio_main
    endif
 
    ! ==============================================
    ! Deallocate prop quantities
    ! ==============================================
    CALL msd_dealloc
    CALL vacf_dealloc

    ! ==============================================
    ! Deallocate coulombic related quantities
    ! ==============================================
    CALL finalize_coulomb
    ! ==============================================
    ! Deallocate principal quantites
    ! ==============================================
    CALL config_dealloc
    CALL extended_coordinates_dealloc

! test
#ifdef block
    CALL block_
#endif
! test 

  else
  ! =======================
  !  only for test purpose
  ! =======================

    separator(stdout)    
    io_node WRITE ( stdout , '(a)' ) '                           TEST                              '
    separator(stdout)    
    ! 
    !    WRITE HERE YOUR CODE TO BE TESTED 
    ! 
  endif

#ifdef MPI
  stotime
  addtime(timetot) 
#endif
  CALL print_time_info ( stdout ) 
  CALL print_dumb
  CALL io_end
#ifdef MPI
  CALL MPI_FINALIZE(ierr)
#endif

END PROGRAM main_MDFF
! ===== fmV =====
