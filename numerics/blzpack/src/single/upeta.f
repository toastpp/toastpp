      SUBROUTINE UPETA (J,AMAXN,BMAXN,BMINN,EPS1,ETAQ,ETAR) 
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    UPETA updates bounds for the partial reorthog. strategy          *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    J     (sii) : number of steps                                    *
C*    AMAXN (ari) : Euclidian norms of (ALPHA)                         *
C*    BMAXN (ari) : Euclidian norms of (BETA)                          *
C*    BMINN (sri) : minimum singular value of current (BETA)           *
C*    EPS1  (sri) : EPS*NVB*sqrt(N)                                    *
C*    ETAQ  (arb) : orthog. bounds among (Q) and Lanczos vectors       *
C*    ETAR  (arb) : orthog. bounds among (R) and Lanczos vectors       *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*    Recurrence implemented:                                          *
C*                                                                     *
C*    eta_(j+1,i) = { [||alpha_(i)||+||alpha_(j)||]*eta_(j,i) +        *
C*                    ||beta_(i+1)||*eta_(j  ,i+1) +                   *
C*                    ||beta_(  i)||*eta_(j  ,i-1) +                   *
C*                    ||beta_(  j)||*eta_(j-1,i  ) } / ||beta_(j+1)||  *
C*                                                                     *
C*    For details, see:                                                *
C*                                                                     *
C*  . `The Lanczos Algorithm with Partial Reorthogonalization', H. D.  *
C*    Simon, Mathematics of Computation, 42:115-142, 1984.             *
C*  . `A Shifted Block Lanczos Algorithm for Solving Sparse Symmetric  *
C*    Eigenvalue Problems', R. G. Grimes, J. G. Lewis and H. D. Simon, *
C*    SIAM J. Matrix Anal. Appl., 15:228-272, 1994.                    *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      REAL             TWO,ZERO
      PARAMETER        (TWO=2.0E0,ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          J
      REAL             AMAXN(*),BMAXN(*),BMINN,EPS1,ETAQ(*),ETAR(*)
C
C==== local variables ==================================================
C
      INTEGER          I
      REAL             TEMP
C
C**** executable statements ********************************************
C
      ETAQ(J) = ZERO
      ETAR(J) = EPS1
C
      IF      ( J .EQ. 2 ) THEN
C
              ETAQ(J-1) = ( (AMAXN(J)+AMAXN(J-1))*EPS1 + 
     &                      BMAXN(  J)*EPS1*TWO ) / BMINN
C
      ELSE IF ( J .GT. 2 ) THEN
C
              ETAQ(J-1) = ( (AMAXN(J)+AMAXN(J-1))*EPS1 + 
     &                      BMAXN(J-1)*ETAQ(J-2) +
     &                      BMAXN(  J)*EPS1*TWO ) / BMINN
C
              ETAQ(  1) = ( (AMAXN(1)+AMAXN(J))*ETAR(1) +
     &                      BMAXN(  2)*ETAR(  2) +
     &                      BMAXN(  J)*ETAQ(  1) ) / BMINN
C
              DO 10 I = 2,J-2
                 ETAQ(  I) = ( (AMAXN(I)+AMAXN(J))*ETAR(I) +
     &                         BMAXN(I+1)*ETAR(I+1) +
     &                         BMAXN(  I)*ETAR(I-1) + 
     &                         BMAXN(  J)*ETAQ(  I) ) / BMINN
   10         CONTINUE
C
      END IF
C
      DO 20 I = 1,J-1
         TEMP    = ETAQ(I)
         ETAQ(I) = ETAR(I)
         ETAR(I) = TEMP
   20 CONTINUE
C
      RETURN
C
C**** end of UPETA *****************************************************
C
      END
