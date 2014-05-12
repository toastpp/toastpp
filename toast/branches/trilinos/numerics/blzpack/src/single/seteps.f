      SUBROUTINE SETEPS (EPS)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SETEPS sets the roundoff unit                                    *
C*                                                                     *
C*    SETEPS is a modified version of the function `machar' available  *
C*           in Netlib and implemented by W. J. Cody.                  *
C*                                                                     *
C*  - Argument:                                                        *
C*                                                                     *
C*    EPS (sro) : roundoff unit                                        *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    REAL,INT                                                         *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      REAL             ONE,ZERO
      PARAMETER        (ONE=1.0E0,ZERO=0.0E0)
C
C==== argument =========================================================
C
      REAL             EPS
C
C==== local variables ==================================================
C
      INTEGER          I,IBETA,IRND,IT
      REAL             A,B,BETA,BETAIN,BETAM1,X
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        REAL,INT
C
C**** executable statements ********************************************
C
C.... determine IBETA and BETA .........................................
C
      A = ONE
   10 CONTINUE
      A = A + A
      X = A + ONE
      X = X - A
      X = X - ONE
      IF ( X .EQ. ZERO ) GO TO 10
C
      B = ONE
   20 CONTINUE
      B = B + B
      X = A + B
      X = X - A
      IF ( X .EQ. ZERO ) GO TO 20
C
      IBETA = INT((A+B)-A)
      BETA = REAL(IBETA)
C
C.... determine IT and IRND ............................................
C
      IT = 0
      B = ONE
   30 CONTINUE
      IT = IT + 1
      B = B*BETA
      X = B + ONE
      X = X - B
      X = X - ONE
      IF ( X .EQ. ZERO ) GO TO 30
C
      IRND = 0
      BETAM1 = BETA - ONE
      X = A + BETAM1
      X = X - A
      IF ( X .NE. ZERO ) IRND = 1
C
C.... determine EPS ....................................................
C
      A = ONE
      BETAIN = ONE/BETA
      DO 40 I = 1,IT+3
         A = A*BETAIN
   40 CONTINUE
C
   50 CONTINUE
      X = ONE + A
      X = X - ONE
      IF ( X .EQ. ZERO ) THEN
         A = A*BETA
         GO TO 50
      END IF
C
      EPS = A
C
      IF ( IBETA.NE.2 .AND. IRND.NE.0 ) THEN
         A = (A*(ONE+A))/(ONE+ONE)
         X = ONE + A
         X = X - ONE
         IF ( X .NE. ZERO ) EPS = A
      END IF
C
      RETURN
C
C**** end of SETEPS ****************************************************
C
      END
