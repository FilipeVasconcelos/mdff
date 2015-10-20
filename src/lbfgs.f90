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

#include "symbol.h"

! ************************************************************************************
!> \file 
!! This file contains the LBFGS algorithm and supporting routines
! ************************************************************************************
!C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!C                          JORGE NOCEDAL
!C                        *** July 1990 ***
!C
!C 
!> \brief
!!     This subroutine solves the unconstrained minimization problem
!! 
!!                      min F(x),    x= (x1,x2,...,xN),
!!
!!      using the limited memory BFGS method. The routine is especially
!!      effective on problems involving a large number of variables. In
!!      a typical iteration of this method an approximation Hk to the
!!      inverse of the Hessian is obtained by applying M BFGS updates to
!!      a diagonal matrix Hk0, using information from the previous M steps.
!!      The user specifies the number M, which determines the amount of
!!      storage required by the routine. The user may also provide the
!!      diagonal matrices Hk0 if not satisfied with the default choice.
!!      The algorithm is described in "On the limited memory BFGS method
!!      for large scale optimization", by D. Liu and J. Nocedal,
!!      Mathematical Programming B 45 (1989) 503-528.
!! 
!!      The user is required to calculate the function value F and its
!!      gradient G. In order to allow the user complete control over
!!      these computations, reverse  communication is used. The routine
!!      must be called repeatedly under the control of the parameter
!!      IFLAG. 
!!
!!      The steplength is determined at each iteration by means of the
!!      line search routine MCVSRCH, which is a slight modification of
!!      the routine CSRCH written by More' and Thuente.
!! 
!C      The calling statement is 
!C 
!C          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)
!C 
!C      where
!C 
!> \param[in]  
!!     N       is an INTEGER variable that must be set by the user to the
!!             number of variables. It is not altered by the routine.
!!             Restriction: N>0.
!C 
!> \param[in]     
!!     M       is an INTEGER variable that must be set by the user to
!!             the number of corrections used in the BFGS update. It
!!             is not altered by the routine. Values of M less than 3 are
!!             not recommended; large values of M will result in excessive
!!             computing time. 3<= M <=7 is recommended. Restriction: M>0.
!C 
!> \param[in]     
!!     X       is a real(kind=dp) array of length N. On initial entry
!!             it must be set by the user to the values of the initial
!!             estimate of the solution vector. On exit with IFLAG=0, it
!!             contains the values of the variables at the best point
!!             found (usually a solution).
!C 
!> \param[in]     
!!     F       is a real(kind=dp) variable. Before initial entry and on
!!             a re-entry with IFLAG=1, it must be set by the user to
!!             contain the value of the function F at the point X.
!C 
!> \param[in]     
!!     G       is a real(kind=dp) array of length N. Before initial
!!             entry and on a re-entry with IFLAG=1, it must be set by
!!             the user to contain the components of the gradient G at
!!             the point X.
!C 
!> \param[in]     
!!     DIAGCO  is a LOGICAL variable that must be set to .TRUE. if the
!!             user  wishes to provide the diagonal matrix Hk0 at each
!!             iteration. Otherwise it should be set to .FALSE., in which
!!             case  LBFGS will use a default value described below. If
!!             DIAGCO is set to .TRUE. the routine will return at each
!!             iteration of the algorithm with IFLAG=2, and the diagonal
!!              matrix Hk0  must be provided in the array DIAG.
!C 
!C 
!> \param[in]     
!!     DIAG    is a real(kind=dp) array of length N. If DIAGCO=.TRUE.,
!!             then on initial entry or on re-entry with IFLAG=2, DIAG
!!             it must be set by the user to contain the values of the 
!!             diagonal matrix Hk0.  Restriction: all elements of DIAG
!!             must be positive.
!C 
!> \param[in]     
!!     IPRINT  is an INTEGER array of length two which must be set by the
!!             user.
!C 
!!             IPRINT(1) specifies the frequency of the output:
!!                IPRINT(1) < 0 : no output is generated,
!!                IPRINT(1) = 0 : output only at first and last iteration,
!!                IPRINT(1) > 0 : output every IPRINT(1) iterations.
!! 
!!             IPRINT(2) specifies the type of output generated:
!!                IPRINT(2) = 0 : iteration count, number of function 
!!                                evaluations, function value, norm of the
!!                                gradient, and steplength,
!!                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of
!!                                variables and  gradient vector at the
!!                                initial point,
!!                IPRINT(2) = 2 : same as IPRINT(2)=1, plus vector of
!!                                variables,
!!                IPRINT(2) = 3 : same as IPRINT(2)=2, plus gradient vector.
!! 
!! 
!> \param[in]     
!!     EPS     is a positive real(kind=dp) variable that must be set by
!!             the user, and determines the accuracy with which the solution
!!             is to be found. The subroutine terminates when
!C
!!                         ||G|| < EPS max(1,||X||),
!!
!!             where ||.|| denotes the Euclidean norm.
!C 
!> \param[in]     
!!     XTOL    is a  positive real(kind=dp) variable that must be set by
!!             the user to an estimate of the machine precision (e.g.
!!             10**(-16) on a SUN station 3/60). The line search routine will
!!C            terminate if the relative width of the interval of uncertainty
!!             is less than XTOL.
!C 
!> \param[in]     
!!     W       is a real(kind=dp) array of length N(2M+1)+2M used as
!!             workspace for LBFGS. This array must not be altered by the
!!             user.
!C 
!> \param[in]     
!!     IFLAG   is an INTEGER variable that must be set to 0 on initial entry
!!             to the subroutine. A return with IFLAG<0 indicates an error,
!!             and IFLAG=0 indicates that the routine has terminated without
!!             detecting errors. On a return with IFLAG=1, the user must
!!             evaluate the function F and gradient G. On a return with
!!             IFLAG=2, the user must provide the diagonal matrix Hk0.
!C 
!!             The following negative values of IFLAG, detecting an error,
!!             are possible:
!C 
!!              IFLAG=-1  The line search routine MCSRCH failed. The
!!                        parameter INFO provides more detailed information
!!                        (see also the documentation of MCSRCH):
!C
!!                       INFO = 0  IMPROPER INPUT PARAMETERS.
!C
!!                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF
!C                                 UNCERTAINTY IS AT MOST XTOL.
!C
!!                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE
!C                                 REQUIRED AT THE PRESENT ITERATION.
!C
!!                       INFO = 4  THE STEP IS TOO SMALL.
!C
!!                       INFO = 5  THE STEP IS TOO LARGE.
!C
!!                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. 
!C                                 THERE MAY NOT BE A STEP WHICH SATISFIES
!C                                 THE SUFFICIENT DECREASE AND CURVATURE
!C                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL.
!C
!C 
!!              IFLAG=-2  The i-th diagonal element of the diagonal inverse
!C                        Hessian approximation, given in DIAG, is not
!C                        positive.
!C           
!!              IFLAG=-3  Improper input parameters for LBFGS (N or M are
!C                        not positive).
!C    COMMON:
!C 
!C     The subroutine contains one common area, which the user may wish to
!C    reference:
!C
!> \note      
!!    GTOL is a real(kind=dp) variable with default value 0.9, which
!!        controls the accuracy of the line search routine MCSRCH. If the
!!        function and gradient evaluations are inexpensive with respect
!!        to the cost of the iteration (which is sometimes the case when
!!        solving very large problems) it may be advantageous to set GTOL
!!        to a small value. A typical small value is 0.1.  Restriction:
!!        GTOL should be greater than 1.D-04.
!C 
!> \note
!!    STPMIN , STPMAX are non-negative real(kind=dp) variables which
!!        specify lower and uper bounds for the step in the line search.
!!        Their default values are 1.D-20 and 1.D+20, respectively. These
!!        values need not be modified unless the exponents are too large
!!        for the machine being used, or unless the problem is extremely
!!        badly scaled (in which case the exponents should be increased).
!C 
!C
!> \note
!! ON THE DRIVER:
!!
!!    The program that calls LBFGS must contain the declaration:
!!
!!                       EXTERNAL LB2
!!
!!    LB2 is a BLOCK DATA that defines the default values of several
!!    parameters described in the COMMON section. 
!!
!> \note
!!  MACHINE DEPENDENCIES
!!
!!        The only variables that are machine-dependent are XTOL,
!!        STPMIN and STPMAX.
!C 
!C
!> \note
!!  GENERAL INFORMATION
!! 
!!    Other routines called directly:  DAXPY, DDOT, LB1, MCSRCH
!C 
!> \note
!!    Input/Output  :  FMV : stdout and ionode 
!C 
!> \author
!! JORGE NOCEDAL ( July 1990 )
!> \note
!! adapted to be included in mdff by fmv october 2011
! ************************************************************************************
      SUBROUTINE LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG)

      USE constants , ONLY : dp 
      USE io,      ONLY :  ionode , stdout

      IMPLICIT NONE
!C
      INTEGER N,M,IPRINT(2),IFLAG
      real(kind=dp) X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M)
      real(kind=dp) F,EPS,XTOL
      LOGICAL DIAGCO
      COMMON /LB3/GTOL,STPMIN,STPMAX
!C
      real(kind=dp) GTOL,ONE,ZERO,GNORM,DDOT,STP1,FTOL,STPMIN,STPMAX,STP,YS,YY,SQ,YR,BETA,XNORM
      INTEGER ITER,NFUN,POINT,ISPT,IYPT,MAXFEV,INFO,BOUND,NPT,CP,I,NFEV,INMC,IYCN,ISCN
      LOGICAL FINISH
!C
      SAVE
      DATA ONE,ZERO/1.0D+0,0.0D+0/
!C
!C     INITIALIZE
!C     ----------
!C
      IF(IFLAG.EQ.0) GO TO 10
      GO TO (172,100) IFLAG
  10  ITER= 0
      IF(N.LE.0.OR.M.LE.0) GO TO 196
      IF(GTOL.LE.1.D-04) THEN
        io_node WRITE ( stdout , 245 )
        GTOL=9.D-01
      ENDIF
      NFUN= 1
      POINT= 0
      FINISH= .FALSE.
      IF(DIAGCO) THEN
         DO 30 I=1,N
 30      IF (DIAG(I).LE.ZERO) GO TO 195
      ELSE
         DO 40 I=1,N
 40      DIAG(I)= 1.0_dp
      ENDIF
!C
!C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
!C     ---------------------------------------
!C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
!C         OTHER TEMPORARY INFORMATION.
!C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
!C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
!C         IN THE FORMULA THAT COMPUTES H*G.
!C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
!C         STEPS.
!C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
!C         GRADIENT DIFFERENCES.
!C
!C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
!C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
!C
      ISPT= N+2*M
      IYPT= ISPT+N*M     
      DO 50 I=1,N
 50   W(ISPT+I)= -G(I)*DIAG(I)
      GNORM= SQRT (DDOT(N,G,1,G,1))
      STP1= ONE/GNORM
!C
!C     PARAMETERS FOR LINE SEARCH ROUTINE
!C     
! original     
      FTOL= 1.0D-04
      MAXFEV= 100
!      FTOL= 1.0D-2
!      MAXFEV= 100
!C
      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
!C
!C    --------------------
!C     MAIN ITERATION LOOP
!C    --------------------
!C
 80   ITER= ITER+1
      INFO=0
      BOUND=ITER-1
      IF(ITER.EQ.1) GO TO 165
      IF (ITER .GT. M)BOUND=M
!C
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
      IF(.NOT.DIAGCO) THEN
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         DO 90 I=1,N
   90    DIAG(I)= YS/YY
      ELSE
         IFLAG=2
         RETURN
      ENDIF
 100  CONTINUE
      IF(DIAGCO) THEN
        DO 110 I=1,N
 110    IF (DIAG(I).LE.ZERO) GO TO 195
      ENDIF
!C
!C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!C     "Updating quasi-Newton matrices with limited storage",
!C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!C     ---------------------------------------------------------
!C
      CP= POINT
      IF (POINT.EQ.0) CP=M
      W(N+CP)= ONE/YS
      DO 112 I=1,N
 112  W(I)= -G(I)
      CP= POINT
      DO 125 I= 1,BOUND
         CP=CP-1
         IF (CP.EQ. -1)CP=M-1
         SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
         INMC=N+M+CP+1
         IYCN=IYPT+CP*N
         W(INMC)= W(N+CP+1)*SQ
         CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
 125  CONTINUE
!C
      DO 130 I=1,N
 130  W(I)=DIAG(I)*W(I)
!C
      DO 145 I=1,BOUND
         YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
         BETA= W(N+CP+1)*YR
         INMC=N+M+CP+1
         BETA= W(INMC)-BETA
         ISCN=ISPT+CP*N
         CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
         CP=CP+1
         IF (CP.EQ.M)CP=0
 145  CONTINUE
!C
!C     STORE THE NEW SEARCH DIRECTION
!C     ------------------------------
!C
       DO 160 I=1,N
 160   W(ISPT+POINT*N+I)= W(I)
!C
!C     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION 
!C     BY USING THE LINE SEARCH ROUTINE MCSRCH
!C     ----------------------------------------------------
 165  NFEV=0
      STP=ONE
      IF (ITER.EQ.1) STP=STP1
      DO 170 I=1,N
 170  W(I)=G(I)
 172  CONTINUE
      CALL MCSRCH(N,X,F,G,W(ISPT+POINT*N+1),STP,FTOL,XTOL,MAXFEV,INFO,NFEV,DIAG)
      IF (INFO .EQ. -1) THEN
        IFLAG=1
        RETURN
      ENDIF
      IF (INFO .NE. 1) GO TO 190
      NFUN= NFUN + NFEV
!C
!C     COMPUTE THE NEW STEP AND GRADIENT CHANGE 
!C     -----------------------------------------
!C
      NPT=POINT*N
      DO 175 I=1,N
      W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
 175  W(IYPT+NPT+I)= G(I)-W(I)
      POINT=POINT+1
      IF (POINT.EQ.M)POINT=0
!C
!C     TERMINATION TEST
!C     ----------------
!C
      GNORM= SQRT (DDOT(N,G,1,G,1))
      XNORM= SQRT (DDOT(N,X,1,X,1))
      XNORM= DMAX1(1.0_dp,XNORM)
      IF (GNORM/XNORM .LE. EPS) FINISH=.TRUE.
!C
! fmv moved to the driver  
      IF(IPRINT(1).GE.0) CALL LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
      IF (FINISH) THEN
         IFLAG=0
         RETURN
      ENDIF
      GO TO 80
!C
!C     ------------------------------------------------------------
!C     END OF MAIN ITERATION LOOP. ERROR EXITS.
!C     ------------------------------------------------------------
!C
 190  IFLAG=-1
      IF ( ionode ) WRITE ( stdout , 200 ) INFO
      RETURN
 195  IFLAG=-2
      IF ( ionode ) WRITE ( stdout , 235 ) I
      RETURN
 196  IFLAG= -3
      IF ( ionode ) WRITE ( stdout , 240 )
!C
!C     FORMATS
!C     -------
!C
 200  FORMAT(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE', &
               ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN', &
               ' OF LINE SEARCH: INFO= ',I2,/ &
               ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/,&
               ' OR INCORRECT TOLERANCES')
 235  FORMAT(/' IFLAG= -2',/' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/, &
            ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
 240  FORMAT(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M',&
            ' ARE NOT POSITIVE)')
 245  FORMAT(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04', &
            / ' IT HAS BEEN RESET TO 9.D-01')
      RETURN
      END
!C
!C     LAST LINE OF SUBROUTINE LBFGS
!C
!C
      SUBROUTINE LB1(IPRINT,ITER,NFUN,GNORM,N,M,X,F,G,STP,FINISH)
 
      USE constants , ONLY : dp 
      USE io,    ONLY :    stdout !, ionode 

      IMPLICIT NONE
!C
!C     -------------------------------------------------------------
!C     THIS ROUTINE PRINTS MONITORING INFORMATION. THE FREQUENCY AND
!C     AMOUNT OF OUTPUT ARE CONTROLLED BY IPRINT.
!C     -------------------------------------------------------------
!C
      INTEGER IPRINT(2),ITER,NFUN,N,M,I
      real(kind=dp) X(N),G(N),F,GNORM,STP,GTOL,STPMIN,STPMAX
      LOGICAL FINISH
      COMMON /LB3/GTOL,STPMIN,STPMAX
!C
      IF (ITER.EQ.0)THEN
           WRITE ( stdout , 10 )
           WRITE ( stdout , 20 ) N , M
           WRITE ( stdout , 30 ) F  ,GNORM
                 IF (IPRINT(2).GE.1)THEN
                     WRITE ( stdout ,40 )
                     WRITE ( stdout ,50 ) (X(I),I=1,N)
                     WRITE ( stdout ,60 )
                     WRITE ( stdout ,50 ) (G(I),I=1,N)
                  ENDIF
           WRITE ( stdout , 10 )
           WRITE ( stdout , 70 )
      ELSE
          IF ((IPRINT(1).EQ.0).AND.(ITER.NE.1.AND..NOT.FINISH))RETURN
              IF (IPRINT(1).NE.0)THEN
                   IF(MOD(ITER-1,IPRINT(1)).EQ.0.OR.FINISH)THEN
                         IF(IPRINT(2).GT.1.AND.ITER.GT.1) WRITE ( stdout , 70 )
                         WRITE ( stdout , 80 ) ITER , NFUN , F ,GNORM , STP
                   ELSE
                         RETURN
                   ENDIF
              ELSE
                   IF( IPRINT(2).GT.1.AND.FINISH) WRITE ( stdout , 70 )
                   WRITE ( stdout , 80 ) ITER,NFUN,F,GNORM,STP
              ENDIF
              IF (IPRINT(2).EQ.2.OR.IPRINT(2).EQ.3)THEN
                    IF (FINISH)THEN
                        WRITE ( stdout , 90 )
                    ELSE
                        WRITE ( stdout , 40 )
                    ENDIF
                      WRITE ( stdout , 50 ) ( X ( I ) , I = 1 , N )
                  IF (IPRINT(2).EQ.3)THEN
                      WRITE ( stdout , 60 )
                      WRITE ( stdout , 50 ) ( G ( I ) , I = 1 , N )
                  ENDIF
              ENDIF
            IF ( FINISH ) WRITE ( stdout , 100 )
      ENDIF
!C
 10   FORMAT('*************************************************')
 20   FORMAT('  N=',I5,'   NUMBER OF CORRECTIONS=',I2,/,  '       INITIAL VALUES')
 30   FORMAT(' F= ',1PD10.3,'   GNORM= ',1PD10.3)
 40   FORMAT(' VECTOR X= ')
 50   FORMAT(6(2X,1PE16.8))
 60   FORMAT(' GRADIENT VECTOR G= ')
 70   FORMAT(/'   I   NFN',4X,'FUNC',8X,'GNORM',7X,'STEPLENGTH'/)
 80   FORMAT(2(I4,1X),3X,3(1PE16.8,2X))
 90   FORMAT(' FINAL POINT X= ')
 100  FORMAT(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.',/' IFLAG = 0')
!C
      RETURN
      END
!C     ******
!C
!C
!C   ----------------------------------------------------------
!C     DATA 
!C   ----------------------------------------------------------
!C
      BLOCK DATA LB2
      USE constants , ONLY : dp 
      real(kind=dp) GTOL,STPMIN,STPMAX
      COMMON /LB3/GTOL,STPMIN,STPMAX
      DATA GTOL,STPMIN,STPMAX/9.0D-01,1.0D-20,1.0D+20/
      !DATA GTOL,STPMIN,STPMAX/1.0D-03,1.0D-20,1.0D+20/
      END
!C
!C
!C   ----------------------------------------------------------
!C
      subroutine daxpy(n,da,dx,incx,dy,incy)
      USE constants , ONLY : dp 
!c
!c     constant times a vector plus a vector.
!c     uses unrolled loops for increments equal to one.
!c     jack dongarra, linpack, 3/11/78.
!c
      real(kind=dp) dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
!c
      if (n.le.0)return
      if (da .eq. 0.0_dp) return
      if (incx.eq.1.and.incy.eq.1)go to 20
!c
!c        code for unequal increments or equal increments
!c          not equal to 1
!c
      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!c
!c        code for both increments equal to 1
!c
!c
!c        clean-up loop
!c
   20 m = mod(n,4)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if ( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
!C
!C
!C   ----------------------------------------------------------
!C
      real(kind=dp) FUNCTION DDOT(n,dx,incx,dy,incy)
      USE constants , ONLY : dp 
      IMPLICIT NONE  
!c
!c     forms the dot product of two vectors.
!c     uses unrolled loops for increments equal to one.
!c     jack dongarra, linpack, 3/11/78.
!c
      real(kind=dp) dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!c
      ddot = 0.0_dp
      dtemp = 0.0_dp
      if (n.le.0)return
      if (incx.eq.1.and.incy.eq.1)go to 20
!c
!c        code for unequal increments or equal increments
!c          not equal to 1
!c
      ix = 1
      iy = 1
      if (incx.lt.0)ix = (-n+1)*incx + 1
      if (incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
!c
!c        code for both increments equal to 1
!c
!c
!c        clean-up loop
!c
   20 m = mod(n,5)
      if ( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if ( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + &
        dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
!C    ------------------------------------------------------------------
!C
!C     **************************
!C     LINE SEARCH ROUTINE MCSRCH
!C     **************************
!C
      SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL,MAXFEV,INFO,NFEV,WA)
      
      USE constants , ONLY : dp 
      USE io,      ONLY :  stdout
     
      IMPLICIT NONE
      INTEGER N,MAXFEV,INFO,NFEV
      real(kind=dp) F,STP,FTOL,GTOL,XTOL,STPMIN,STPMAX
      real(kind=dp) X(N),G(N),S(N),WA(N)
      COMMON /LB3/GTOL,STPMIN,STPMAX
      SAVE
!C
!C                     SUBROUTINE MCSRCH
!C                
!C     A slight modification of the subroutine CSRCH of More' and Thuente.
!C     The changes are to allow reverse communication, and do not affect
!C     the performance of the routine. 
!C
!C     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
!C     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
!C
!C     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
!C     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
!C     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
!C     MINIMIZER OF THE MODIFIED FUNCTION
!C
!C          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
!C
!C     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
!C     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
!C     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
!C     CONTAINS A MINIMIZER OF F(X+STP*S).
!C
!C     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
!C     THE SUFFICIENT DECREASE CONDITION
!C
!C           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
!C
!C     AND THE CURVATURE CONDITION
!C
!C           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
!C
!C     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
!C     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
!C     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
!C     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
!C     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
!C     SATISFIES THE SUFFICIENT DECREASE CONDITION.
!C
!C     THE SUBROUTINE STATEMENT IS
!C
!C        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
!C     WHERE
!C
!C       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
!C         OF VARIABLES.
!C
!C       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!C         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
!C         X + STP*S.
!C
!C       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
!C         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
!C
!C       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
!C         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
!C         OF F AT X + STP*S.
!C
!C       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
!C         SEARCH DIRECTION.
!C
!C       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
!C         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
!C         STP CONTAINS THE FINAL ESTIMATE.
!C
!C       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
!C         communication implementation GTOL is defined in a COMMON
!C         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
!C         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
!C         SATISFIED.
!C
!C       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
!C         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!C         IS AT MOST XTOL.
!C
!C       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
!C         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
!C         communication implementatin they are defined in a COMMON
!C         statement).
!C
!C       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
!C         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
!C         MAXFEV BY THE END OF AN ITERATION.
!C
!C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!C
!C         INFO = 0  IMPROPER INPUT PARAMETERS.
!C
!C         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.
!C
!C         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
!C                   DIRECTIONAL DERIVATIVE CONDITION HOLD.
!C
!C         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
!C                   IS AT MOST XTOL.
!C
!C         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
!C
!C         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
!C
!C         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
!C
!C         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
!C                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
!C                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
!C                   TOLERANCES MAY BE TOO SMALL.
!C
!C       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
!C         CALLS TO FCN.
!C
!C       WA IS A WORK ARRAY OF LENGTH N.
!C
!C     SUBPROGRAMS CALLED
!C
!C       MCSTEP
!C
!C       FORTRAN-SUPPLIED...ABS,MAX,MIN
!C
!C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!C     JORGE J. MORE', DAVID J. THUENTE
!C
!C     **********
      INTEGER INFOC,J
      LOGICAL BRACKT,STAGE1
      real(kind=dp) DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM, &
             FINIT,FTEST1,FM,FX,FXM,FY,FYM,P5,P66,STX,STY, &
             STMIN,STMAX,WIDTH,WIDTH1,XTRAPF,ZERO
      DATA P5,P66,XTRAPF,ZERO /0.5_dp,0.66_dp,4.0_dp,0.0_dp/
      IF(INFO.EQ.-1) GO TO 45
      INFOC = 1
!C
!C     CHECK THE INPUT PARAMETERS FOR ERRORS.
!C
      IF (N .LE. 0 .OR. STP .LE. ZERO .OR. FTOL .LT. ZERO .OR. &
         GTOL .LT. ZERO .OR. XTOL .LT. ZERO .OR. STPMIN .LT. ZERO &
         .OR. STPMAX .LT. STPMIN .OR. MAXFEV .LE. 0) RETURN
!C
!C     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
!C     AND CHECK THAT S IS A DESCENT DIRECTION.
!C
      DGINIT = ZERO
      DO 10 J = 1, N
         DGINIT = DGINIT + G(J)*S(J)
   10    CONTINUE
      IF (DGINIT .GE. ZERO) then
         WRITE ( stdout , 15 )
   15    FORMAT(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
         RETURN
         ENDIF
!C
!C     INITIALIZE LOCAL VARIABLES.
!C
      BRACKT = .FALSE.
      STAGE1 = .TRUE.
      NFEV = 0
      FINIT = F
      DGTEST = FTOL*DGINIT
      WIDTH = STPMAX - STPMIN
      WIDTH1 = WIDTH/P5
      DO 20 J = 1, N
         WA(J) = X(J)
   20    CONTINUE
!C
!C     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
!C     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
!C     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
!C     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
!C     THE INTERVAL OF UNCERTAINTY.
!C     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
!C     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
!C
      STX = ZERO
      FX = FINIT
      DGX = DGINIT
      STY = ZERO
      FY = FINIT
      DGY = DGINIT
!C
!C     START OF ITERATION.
!C
   30 CONTINUE
!C
!C        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
!C        TO THE PRESENT INTERVAL OF UNCERTAINTY.
!C
         IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
         ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
            END IF
!C
!C        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
!C
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)
!C
!C        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
!C        STP BE THE LOWEST POINT OBTAINED SO FAR.
!C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX)) &
           .OR. NFEV .GE. MAXFEV-1 .OR. INFOC .EQ. 0 &
           .OR. (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX)) STP = STX
!C
!C        EVALUATE THE FUNCTION AND GRADIENT AT STP
!C        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
!C        We return to main program to obtain F and G.
!C
         DO 40 J = 1, N
            X(J) = WA(J) + STP*S(J)
   40       CONTINUE
         INFO=-1
         RETURN
!C
   45    INFO=0
         NFEV = NFEV + 1
         DG = ZERO
         DO 50 J = 1, N
            DG = DG + G(J)*S(J)
   50       CONTINUE
         FTEST1 = FINIT + STP*DGTEST
!C
!C        TEST FOR CONVERGENCE.
!C
         IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX)) .OR. INFOC .EQ. 0) INFO = 6
         IF (STP .EQ. STPMAX .AND. F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
         IF (STP .EQ. STPMIN .AND. (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
         IF (NFEV .GE. MAXFEV) INFO = 3
         IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
         IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1
!C
!C        CHECK FOR TERMINATION.
!C
         IF (INFO .NE. 0) RETURN
!C
!C        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
!C        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
!C
         IF (STAGE1 .AND. F .LE. FTEST1 .AND. DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.
!C
!C        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
!C        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
!C        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
!C        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
!C        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
!C
         IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
!C
!C           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
!C
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
!C
!C           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!C           AND TO COMPUTE THE NEW STEP.
!C
            CALL MCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,BRACKT,STMIN,STMAX,INFOC)
!C
!C           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
!C
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
         ELSE
!C
!C           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
!C           AND TO COMPUTE THE NEW STEP.
!C
            CALL MCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,BRACKT,STMIN,STMAX,INFOC)
            END IF
!C
!C        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
!C        INTERVAL OF UNCERTAINTY.
!C
         IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1) STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
            END IF
!C
!C        END OF ITERATION.
!C
         GO TO 30
!C
!C     LAST LINE OF SUBROUTINE MCSRCH.
!C
      END


      SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DPP,BRACKT,STPMIN,STPMAX,INFO)
      USE constants , ONLY : dp 
      IMPLICIT NONE 
      INTEGER INFO
      real(kind=dp) STX,FX,DX,STY,FY,DY,STP,FP,DPP,STPMIN,STPMAX
      LOGICAL BRACKT,BOUND
!C
!C     SUBROUTINE MCSTEP
!C
!C     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
!C     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
!C     A MINIMIZER OF THE FUNCTION.
!C
!C     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
!C     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
!C     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
!C     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
!C     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
!C     WITH ENDPOINTS STX AND STY.
!C
!C     THE SUBROUTINE STATEMENT IS
!C
!C       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DPP,BRACKT,
!C                        STPMIN,STPMAX,INFO)
!C
!C     WHERE
!C
!C       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
!C         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
!C         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
!C         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
!C         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
!C
!C       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
!C         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
!C         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
!C         UPDATED APPROPRIATELY.
!C
!C       STP, FP, AND DPP ARE VARIABLES WHICH SPECIFY THE STEP,
!C         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
!C         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
!C         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
!C
!C       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
!C         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
!C         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
!C         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
!C
!C       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
!C         AND UPPER BOUNDS FOR THE STEP.
!C
!C       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!C         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
!C         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
!C         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
!C
!C     SUBPROGRAMS CALLED
!C
!C       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT
!C
!C     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
!C     JORGE J. MORE', DAVID J. THUENTE
!C
      real(kind=dp) GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
      INFO = 0
!C
!C     CHECK THE INPUT PARAMETERS FOR ERRORS.
!C
      IF ((BRACKT .AND. (STP .LE. MIN(STX,STY) .OR. STP .GE. MAX(STX,STY))) .OR. DX*(STP-STX) .GE. 0.0 .OR. STPMAX .LT. STPMIN) RETURN
!C
!C     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
!C
      SGND = DPP*(DX/ABS(DX))
!C
!C     FIRST CASE. A HIGHER FUNCTION VALUE.
!C     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
!C     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
!C     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
!C
      IF (FP .GT. FX) THEN
         INFO = 1
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DPP
         S = MAX(ABS(THETA),ABS(DX),ABS(DPP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DPP/S))
         IF (STP .LT. STX) GAMMA = -GAMMA
         P = (GAMMA - DX) + THETA
         Q = ((GAMMA - DX) + GAMMA) + DPP
         R = P/Q
         STPC = STX + R*(STP - STX)
         STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
         IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
         ELSE
           STPF = STPC + (STPQ - STPC)/2
           END IF
         BRACKT = .TRUE.
!C
!C     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
!C     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
!C     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
!C     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
!C
      ELSE IF (SGND .LT. 0.0) THEN
         INFO = 2
         BOUND = .FALSE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DPP
         S = MAX(ABS(THETA),ABS(DX),ABS(DPP))
         GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DPP/S))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DPP) + THETA
         Q = ((GAMMA - DPP) + GAMMA) + DX
         R = P/Q
         STPC = STP + R*(STX - STP)
         STPQ = STP + (DPP/(DPP-DX))*(STX - STP)
         IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
         ELSE
            STPF = STPQ
            END IF
         BRACKT = .TRUE.
!C
!C     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
!C     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
!C     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
!C     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
!C     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
!C     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
!C     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
!C
      ELSE IF (ABS(DPP) .LT. ABS(DX)) THEN
         INFO = 3
         BOUND = .TRUE.
         THETA = 3*(FX - FP)/(STP - STX) + DX + DPP
         S = MAX(ABS(THETA),ABS(DX),ABS(DPP))
!C
!C        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
!C        TO INFINITY IN THE DIRECTION OF THE STEP.
!C
         GAMMA = S*SQRT(MAX(0.0_dp,(THETA/S)**2 - (DX/S)*(DPP/S)))
         IF (STP .GT. STX) GAMMA = -GAMMA
         P = (GAMMA - DPP) + THETA
         Q = (GAMMA + (DX - DPP)) + GAMMA
         R = P/Q
         IF (R .LT. 0.0 .AND. GAMMA .NE. 0.0) THEN
            STPC = STP + R*(STX - STP)
         ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
         ELSE
            STPC = STPMIN
            END IF
         STPQ = STP + (DPP/(DPP-DX))*(STX - STP)
         IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
         ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
               STPF = STPC
            ELSE
               STPF = STPQ
               END IF
            END IF
!C
!C     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
!C     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
!C     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
!C     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
!C
      ELSE
         INFO = 4
         BOUND = .FALSE.
         IF (BRACKT) THEN
            THETA = 3*(FP - FY)/(STY - STP) + DY + DPP
            S = MAX(ABS(THETA),ABS(DY),ABS(DPP))
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DPP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DPP) + THETA
            Q = ((GAMMA - DPP) + GAMMA) + DY
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
         ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
         ELSE
            STPF = STPMIN
            END IF
         END IF
!C
!C     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
!C     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
!C
      IF (FP .GT. FX) THEN
         STY = STP
         FY = FP
         DY = DPP
      ELSE
         IF (SGND .LT. 0.0) THEN
            STY = STX
            FY = FX
            DY = DX
            END IF
         STX = STP
         FX = FP
         DX = DPP
         END IF
!C
!C     COMPUTE THE NEW STEP AND SAFEGUARD IT.
!C
      STPF = MIN(STPMAX,STPF)
      STPF = MAX(STPMIN,STPF)
      STP = STPF
      IF (BRACKT .AND. BOUND) THEN
         IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
         ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
            END IF
         END IF
      RETURN
!C
!C     LAST LINE OF SUBROUTINE MCSTEP.
!C
      END

