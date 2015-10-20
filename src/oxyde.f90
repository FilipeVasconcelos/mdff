MODULE oxyde

  USE constants,        ONLY :  dp 
! 
! pour le moment il n'y a que 6 oxydes 
! mais l'ajout de nouveaux oxydes ne devraient pas poser de probleme.
! à modifier (dans le cas d'ajout d'un ou plusieurs oxydes) :
!        - le nombre d'oxyde = "noxyde"
!        - creer une variable reel pour etre lu en entré : real(kind=dp) :: sio2 , na2o , b2o3 ...
!        - ajouter le nouvel element à la liste des elements "ele"
!        - ajouter le degree d'oxydation "numoxyd"
!        - ajouter les infos pour les oxydes aux listes lox, ele_ox et nel_ox
!        - attention lox, ele_ox et nel_ox sont definies dans le meme ordre

! note:
! il y aurai certainement un moyen de reduire le nombre d'info 
! necessaire... juste à partir de la liste des elements et du numero
! d'oxydation par exemple.

  integer , PARAMETER :: noxyde = 10   
  integer , PARAMETER :: nelem  = 57

  real(kind=dp)       :: sio2 , na2o , b2o3 , cao , p2o5 , al2o3, geo2 , la2o3, sro , moo3 ! valeur en entré

  ! type structur pour l'element chimique 
  TYPE element
    character(len=2)    :: elename         !< nom de l'element
    real(kind=dp)       :: massele         !< mass atomic (amu)
    integer             :: valence         !< nombre d'electrons de valence (pour le caclcul du nombre de bandes )
    integer             :: numoxyd         !< correspond à la charge pour un champ de force "rigid ion"
  END TYPE

  ! type structur pour l'oxyde 
  TYPE oxy
    character(len=5)    :: nameox          !< label de l'oxyde ex: SiO2, Na2O
    character(len=2)    :: ele_ox  ( 2 )   !< element de l'oxyde ex:Si , Na . le second est toujours l'oxygene
    integer             :: nel_ox  ( 2 )   !< nombre d'ions par oxyde  
    real(kind=dp)       :: relcon          !< relative concentration
  END TYPE

  ! tableau de type element = tableau periodique  ( nelem_max = 57 )
  TYPE(element), dimension (:) , allocatable :: tabper
  ! tableau de type oxydes                        ( noxyde_max = 9 ) 
  TYPE(oxy)    , dimension (:) , allocatable :: oxydes 

CONTAINS

SUBROUTINE gen_tab_period

  USE io,               ONLY :  stdout 

  implicit none

  allocate ( tabper(nelem) )

  ! # 1
  ! # 2
  ! # 3
  ! # 4
  ! # 5
  tabper(5)%elename='B'
  tabper(5)%numoxyd=3
  tabper(5)%valence=3
  tabper(5)%massele=10.811_dp
  ! # 6
  ! # 7
  ! # 8
  tabper(8)%elename='O'
  tabper(8)%numoxyd=-2
  tabper(8)%valence=6
  tabper(8)%massele=15.9994_dp
  ! # 9
  ! # 10
  ! # 11
  tabper(11)%elename='Na'
  tabper(11)%numoxyd=1
  tabper(11)%valence=7
  tabper(11)%massele=22.9898_dp
  ! # 12
  ! # 13
  tabper(13)%elename='Al'
  tabper(13)%numoxyd=3
  tabper(13)%valence=3
  tabper(13)%massele=26.9815_dp
  ! # 14
  tabper(14)%elename='Si'
  tabper(14)%numoxyd=4
  tabper(14)%valence=4
  tabper(14)%massele=28.086_dp
  ! # 15
  tabper(15)%elename='P'
  tabper(15)%numoxyd=5
  tabper(15)%valence=5
  tabper(15)%massele=30.9738_dp
  ! # 16
  ! # 17
  ! # 18
  tabper(18)%elename='Ar'
  tabper(18)%numoxyd=0
  tabper(18)%valence=8
  tabper(18)%massele=39.948_dp
  ! # 19
  ! # 20
  tabper(20)%elename='Ca'
  tabper(20)%numoxyd=2
  tabper(20)%valence=8
  tabper(20)%massele=40.08_dp
  ! # 21
  ! # 22
  ! # 23
  ! # 24
  ! # 25
  ! # 26
  ! # 27
  ! # 28
  ! # 29
  ! # 30
  ! # 31
  ! # 32
  tabper(32)%elename='Ge'
  tabper(32)%numoxyd=4
  tabper(32)%valence=4
  tabper(32)%massele=72.92_dp
  ! # 33
  ! # 34
  ! # 35
  ! # 36
  ! # 37
  ! # 38
  tabper(38)%elename='Sr'
  tabper(38)%numoxyd=2
  tabper(38)%valence=2
  tabper(38)%massele=87.62_dp
  ! # 39
  ! # 40
  ! # 41
  ! # 42
  tabper(42)%elename='Mo'
  tabper(42)%numoxyd=6
  tabper(42)%valence=14
  tabper(42)%massele=95.95_dp
  ! # 43
  ! # 44
  ! # 45
  ! # 46
  ! # 47
  ! # 48
  ! # 49
  ! # 50
  ! # 51
  ! # 52
  ! # 53
  ! # 54
  ! # 55
  ! # 56
  ! # 57
  tabper(57)%elename='La'
  tabper(57)%numoxyd=3
  tabper(57)%valence=9
  tabper(57)%massele=138.90547_dp

  WRITE ( stdout , '(a)' ) 'periodic table generated'

  return

END SUBROUTINE 

SUBROUTINE gen_oxydes

  implicit none

  allocate ( oxydes(noxyde) )

  ! # 1 
  oxydes(1)%nameox='B2O3'
  oxydes(1)%ele_ox(1)='B'
  oxydes(1)%ele_ox(2)='O'
  oxydes(1)%nel_ox(1)=2
  oxydes(1)%nel_ox(2)=3
  oxydes(1)%relcon = b2o3
  ! # 2
  oxydes(2)%nameox='Na2O'
  oxydes(2)%ele_ox(1)='Na'
  oxydes(2)%ele_ox(2)='O'
  oxydes(2)%nel_ox(1)=2
  oxydes(2)%nel_ox(2)=1
  oxydes(2)%relcon = na2o
  ! # 3
  oxydes(3)%nameox='Al2O3'
  oxydes(3)%ele_ox(1)='Al'
  oxydes(3)%ele_ox(2)='O'
  oxydes(3)%nel_ox(1)=2
  oxydes(3)%nel_ox(2)=3
  oxydes(3)%relcon = al2o3
  ! # 4
  oxydes(4)%nameox='SiO2'
  oxydes(4)%ele_ox(1)='Si'
  oxydes(4)%ele_ox(2)='O'
  oxydes(4)%nel_ox(1)=1
  oxydes(4)%nel_ox(2)=2
  oxydes(4)%relcon = sio2
  ! # 5
  oxydes(5)%nameox='P2O5'
  oxydes(5)%ele_ox(1)='P'
  oxydes(5)%ele_ox(2)='O'
  oxydes(5)%nel_ox(1)=2
  oxydes(5)%nel_ox(2)=5
  oxydes(5)%relcon = p2o5
  ! # 6
  oxydes(6)%nameox='CaO'
  oxydes(6)%ele_ox(1)='Ca'
  oxydes(6)%ele_ox(2)='O'
  oxydes(6)%nel_ox(1)=1
  oxydes(6)%nel_ox(2)=1
  oxydes(6)%relcon = cao
  ! # 7
  oxydes(7)%nameox='GeO2'
  oxydes(7)%ele_ox(1)='Ge'
  oxydes(7)%ele_ox(2)='O'
  oxydes(7)%nel_ox(1)=1
  oxydes(7)%nel_ox(2)=2
  oxydes(7)%relcon = geo2 
  ! # 8 
  oxydes(8)%nameox='La2O3'
  oxydes(8)%ele_ox(1)='La'
  oxydes(8)%ele_ox(2)='O'
  oxydes(8)%nel_ox(1)=2
  oxydes(8)%nel_ox(2)=3
  oxydes(8)%relcon = la2o3 
  ! # 9 
  oxydes(9)%nameox='SrO'
  oxydes(9)%ele_ox(1)='Sr'
  oxydes(9)%ele_ox(2)='O'
  oxydes(9)%nel_ox(1)=1
  oxydes(9)%nel_ox(2)=1
  oxydes(9)%relcon = sro 
  ! # 10 
  oxydes(10)%nameox='MoO3'
  oxydes(10)%ele_ox(1)='Mo'
  oxydes(10)%ele_ox(2)='O'
  oxydes(10)%nel_ox(1)=1
  oxydes(10)%nel_ox(2)=3
  oxydes(10)%relcon = moo3 

  return

END SUBROUTINE

SUBROUTINE tabper_oxydes_dalloc

  implicit none

  IF (ALLOCATED(oxydes)) THEN
    deallocate(oxydes)
  ENDIF
  IF (ALLOCATED(tabper)) THEN
    deallocate(tabper)
  ENDIF

  return
END SUBROUTINE tabper_oxydes_dalloc


END MODULE oxyde
