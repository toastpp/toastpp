      DOUBLE PRECISION FUNCTION RESNRM (N,BETA,BIGNUM,S)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RESNRM computes the residual norm of a Ritz pair                 *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    N      (sii) : dimension of (BETA)                               *
C*    BETA   (ari) : matrix (BETA) in (R)=(Q)*(BETA) at j-th step      *
C*    BIGNUM (sri) : big number                                        *
C*    S      (ari) : bottom N elements of an eigenvector (s) of (TB)   *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,SQRT                                                         *
C*                                                                     *
C***********************************************************************
C     
C==== parameter ========================================================
C
      DOUBLE PRECISION ONE,ZERO
      PARAMETER        (ONE=1.0D0,ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          N
      DOUBLE PRECISION BETA(N,N),BIGNUM,S(N)
C
C==== local variables ==================================================
C
      INTEGER          I,J
      DOUBLE PRECISION SUM
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,SQRT
C
C**** executable statements ********************************************
C
C     resnrm_i = || (BETA)*(S)_i ||
C
      RESNRM = ZERO
C
C.... multiply (BETA) by (S) and store in RESNRM .......................
C
      DO 20 I = 1,N
         SUM = ZERO
         DO 10 J = I,N
            SUM = SUM + BETA(I,J)*S(J)
   10    CONTINUE
         RESNRM = RESNRM + SUM*SUM
   20 CONTINUE
C
C.... the residual norm is the square root .............................
C
      IF ( RESNRM*BIGNUM .GT. ONE ) THEN
         RESNRM = SQRT(RESNRM)
      ELSE
         RESNRM = ONE/BIGNUM
      END IF
C
      RETURN
C
C**** end of RESNRM ****************************************************
C
      END
