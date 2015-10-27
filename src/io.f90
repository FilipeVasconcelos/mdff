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

! *********************** MODULE IO_FILE ***************************************
!
!> \brief definition of the units for output files
!
! ******************************************************************************
MODULE io

  implicit none

  !> rank for output (if true myrank.eq.0)
  logical :: ionode        
  logical :: ioprint
  logical :: ioprintnode        

  !> standard output
  integer, PARAMETER :: stdout          = 6 
  !> standard error output
  integer, PARAMETER :: stderr          = 0  

  !> standard input file ( code argument control.F )
  integer, PARAMETER :: stdin           = 1001

  !> tmp working file
  integer, PARAMETER :: kunit_tmp       =  8

  !> thermodynamic info
  integer, PARAMETER :: kunit_OSZIFF    = 10

  !> MD, RMC trajectory (positions, velocities , forces )
  integer, PARAMETER :: kunit_TRAJFF    = 20

  !> input configuration
  integer, PARAMETER :: kunit_POSFF     = 30

  !> end configuration 
  integer, PARAMETER :: kunit_CONTFF    = 40

  !> all efg for each atoms (if lefgprintall)
  integer, PARAMETER :: kunit_EFGALL    = 60 
  integer, PARAMETER :: kunit_EFGALLIT1 = 61
  integer, PARAMETER :: kunit_EFGALLIT2 = 62 

  !> all efg for each atoms (if lefgprintall)
  integer, PARAMETER :: kunit_NMRFF     = 65
  integer, PARAMETER :: kunit_NMRFFIT1  = 66
  integer, PARAMETER :: kunit_NMRFFIT2  = 67 
  
  !> all Efield for each atoms (if doefield)
  integer, PARAMETER :: kunit_EFALL    = 68 

  !> EFG eta distribution (average) 
  integer, PARAMETER :: kunit_DTETAFF   = 70
  integer, PARAMETER :: kunit_DTETAFFIT = 71

  !> EFG vzz distribution (average)
  integer, PARAMETER :: kunit_DTVZZFF   = 80
  integer, PARAMETER :: kunit_DTVZZFFIT = 81

  !> EFG tensor component distribution Ui  (average)
  integer, PARAMETER :: kunit_DTIBUFF   = 90
  integer, PARAMETER :: kunit_DTIBUFFIT = 91
  !> EFG tensor distribution  S
  integer, PARAMETER :: kunit_DTIBSFF   = 95

  !> GR radial distribution (average)
  integer, PARAMETER :: kunit_GRTFF     = 100

  !> OPT output configuration after optimisation
  integer, PARAMETER :: kunit_ISCFF     = 110

  !> OPT thermodynamics properties of optimized structure
  integer, PARAMETER :: kunit_ISTHFF    = 120 

  !> VIB eigenvalues (frequencies) of the hessian matrix
  integer, PARAMETER :: kunit_EIGFF     = 130

  !> VIB eigenvector (normal modes) of the hessian matrix 
  integer, PARAMETER :: kunit_VECTFF    = 140

  !> VIB density of states 
  integer, PARAMETER :: kunit_DOSFF     = 150

  !> VIB generated configuration of a given mode
  integer, PARAMETER :: kunit_MODFF     = 160

  !> VIB kpoint mesh for the complete dos
  integer, PARAMETER :: kunit_IBZKPTFF  = 170

  !> VIB "complete" density of states
  integer, PARAMETER :: kunit_DOSKFF    = 180
  integer, PARAMETER :: kunit_DKFF      = 181

  !> VIB fvibcalc output
  integer, PARAMETER :: kunit_VIBFF     = 190

  !> MSD output file
  integer, PARAMETER :: kunit_MSDFF     = 200

  !> stress tensor output file
  integer, PARAMETER :: kunit_STRESSFF  = 210

  !> velocity auto-correlation output file
  integer, PARAMETER :: kunit_VACFFF    = 220
 
  !> EFG auto-correlation output file
  integer, PARAMETER :: kunit_EFGACFFF  = 230
  integer, PARAMETER :: kunit_NMRACFFF  = 231
  integer, PARAMETER :: kunit_UIACFFF   = 232

  !> static structure factor
  integer, PARAMETER :: kunit_STRFACFF  = 240

  !> ??
  integer, PARAMETER :: kunit_EQUILFF   = 250

  !> mean number of atoms in a shell of width at distance r
  integer, PARAMETER :: kunit_NRTFF     = 260

  !> dipole moments on atoms from Wannier centers
  integer, PARAMETER :: kunit_DIPWFC    = 270
  !> total dipole moments on atoms
  integer, PARAMETER :: kunit_DIPFF     = 271
  integer, PARAMETER :: kunit_QUADFF    = 272

  !> neighbor info
  integer, PARAMETER :: kunit_VOIS1FF   = 280

  !> number of neighbors distribution
  integer, PARAMETER :: kunit_DTNBFF    = 290

  !> input config info for rmc calculation 
  integer, PARAMETER :: kunit_RMCFF     = 300
  !> output chi rmc 
  integer, PARAMETER :: kunit_RMCLOG    = 310

  integer, PARAMETER :: kunit_conf_proc = 400
  integer, dimension(:) , allocatable :: kunit_bak_proc

  integer, PARAMETER :: kunit_RESTART   = 600

CONTAINS

! *********************** SUBROUTINE io_init ***********************************
!
!> \brief
!! initialize the ionode logical variable
!
!> \note
!! usage : io_node WRITE ( unit , * )   (symbol.h)
!
! ******************************************************************************
SUBROUTINE io_init

  implicit none
#ifdef MPI
  include "mpif.h"
#endif

  ! local 
  integer :: myrank , numprocs , ierr

  myrank = 0 
#ifdef MPI 
  CALL MPI_COMM_RANK ( MPI_COMM_WORLD , myrank , ierr )    ! numero de processus (output myrank .eq. 0 )
  CALL MPI_COMM_SIZE ( MPI_COMM_WORLD , numprocs , ierr )  ! nombre de processus
#endif


  if ( myrank .eq. 0 ) then
    ionode = .true.
  else
    ionode = .false.
  endif

  allocate( kunit_bak_proc(0:numprocs) )
  kunit_bak_proc(myrank) = kunit_conf_proc + myrank

  return

END SUBROUTINE io_init

SUBROUTINE io_end
  implicit none

  deallocate( kunit_bak_proc )

  return

END SUBROUTINE io_end

! *********************** SUBROUTINE io_open ***********************************
!
!> \brief
!! open file and check status
!
! ******************************************************************************

SUBROUTINE io_open ( kunit , cunit , iostatus )

  implicit none

  ! global
  integer,      intent (in) :: kunit 
  character(*), intent (in) :: cunit 
  character(*), intent (in) :: iostatus

  ! local
  integer                   :: ioerr

  OPEN(UNIT=kunit,FILE=trim(sweep_blanks( cunit )),IOSTAT=ioerr,STATUS=iostatus)

  if ( ioerr .lt. 0 )  then
    io_node WRITE ( stdout, '(a,a,i4)') 'ERROR opening file : ', trim(sweep_blanks( cunit )) , kunit
    STOP
  endif

  return

END SUBROUTINE io_open

! *********************** SUBROUTINE io_close **********************************
!
!> \brief
!! close file and check status
!
! ******************************************************************************

SUBROUTINE io_close ( kunit ) 
  
  implicit none

  ! global
  integer, intent (in)    :: kunit 

  ! local
  integer                 :: ioerr

  CLOSE(UNIT=kunit,IOSTAT=ioerr)

  if ( ioerr .lt. 0 )  then
    io_node WRITE ( stdout, '(a,i4)') 'ERROR closing file : ', kunit
    STOP
  endif

  return

END SUBROUTINE io_close

! *********************** FUNCTION sweep_blanks ********************************
!
!> \brief
!! to remove leading and trailing spaces
!
! ******************************************************************************

character(len=30) FUNCTION sweep_blanks ( in_str )

  implicit none

  character(*), intent(in)  :: in_str
  character(len=30)         :: out_str
  character                 :: ch
  integer                   :: j

  out_str = " "

  do j=1, len_trim(in_str)

  ! get j-th char
    ch = in_str(j:j)
    if (ch .ne. " ") then
      out_str = trim(out_str) // ch
    endif
    sweep_blanks = out_str
  enddo
        
END FUNCTION sweep_blanks

END MODULE io
! ===== fmV =====
