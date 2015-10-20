#ifdef timing
! ====================
!  timing symbols
! ====================
! dectime : declare
#define dectime real(kind=dp) :: ttt0 , ttt , ierr 
! statime : start
#define statime ttt0 = MPI_WTIME(ierr)
! stotime : stop
#define stotime ttt  = MPI_WTIME(ierr)
! addtime(X) : add to variable X
#define addtime(X) X = X + ( ttt - ttt0 ) 
! writime(X1,X2,IC) : write
#define writime(X1,X2,IC) if ( ionode ) WRITE ( stdout , '(1X,A8,I12,A10,A,F9.2)' ) X1, IC ,X2,' : cpu time ',ttt-ttt0
#else
#define dectime
#define statime
#define stotime
#define addtime(X)
#define writime(X1,X2,IC)
#endif

! separator '='
#define separator(X)     if ( ionode ) WRITE ( X , '(a)' ) repeat('=',61) 
#define separator_noionode(X)          WRITE ( X , '(a)' ) repeat('=',61) 
#define bigseparator(X)     if ( ionode ) WRITE ( X , '(a)' ) repeat('=',100) 
#define bigseparator_noionode(X)          WRITE ( X , '(a)' ) repeat('=',100) 

! separator '-'
#define lseparator(X)     if ( ionode ) WRITE ( X , '(a)' ) repeat('-',61) 
#define lseparator_noionode(X)          WRITE ( X , '(a)' ) repeat('-',61) 
#define lseparator_ioprintnode(X)  if ( ioprintnode )        WRITE ( X , '(a)' ) repeat('-',61) 

! blank line
#define blankline(X)            WRITE ( X , '(a)' ) ''

! node condition io_node=if ( ionode )
#define io_node if ( ionode )

! md print condition io_print=if ( ioprint )
#define io_print if ( ioprint )

! md print condition io_printnode=if ( ioprintnode )
#define io_printnode if ( ioprintnode )
