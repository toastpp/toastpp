      REAL             FUNCTION TBGMIN (N,REPS,THETA)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBGMIN looks for the minimum gap in (THETA)                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    N     (sii) : dimension of (THETA)                               *
C*    REPS  (sri) : sqrt(EPS)                                          *
C*    THETA (ari) : eigenvalues of the block tridiagonal matrix        *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,MIN                                                          *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          N
      REAL             REPS,THETA(N)
C
C==== local variables ==================================================
C
      INTEGER          I
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,MIN
C
C**** executable statements ********************************************
C
      TBGMIN = ABS(THETA(1)-THETA(2))
C
C.... gapmin is the minimum distance between two thetas ................
C
      DO 10 I = 2,N-1
         TBGMIN = MIN(TBGMIN,ABS(THETA(I)-THETA(I+1)))
   10 CONTINUE
C
      IF ( TBGMIN .LT. REPS ) TBGMIN = REPS
C
      RETURN 
C
C**** end of TBGMIN ****************************************************
C
      END
