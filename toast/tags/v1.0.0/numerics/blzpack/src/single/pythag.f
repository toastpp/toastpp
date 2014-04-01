      REAL             FUNCTION PYTHAG (A,B)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    PYTHAG finds SQRT(A**2+B**2) safely                              *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,MAX,MIN                                                      *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      REAL             FOUR,ONE,TWO,ZERO
      PARAMETER        (FOUR=4.0E0,ONE=1.0E0,TWO=2.0E0,ZERO=0.0E0)
C
C==== arguments ========================================================
C
      REAL             A,B
C
C==== local variables ==================================================
C
      REAL             P,R,S,T,U
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,MAX,MIN
C
C**** executable statements ********************************************
C
      P = MAX(ABS(A),ABS(B))
      IF (P .EQ. ZERO) GO TO 20
      R = (MIN(ABS(A),ABS(B))/P)**2
   10 CONTINUE
         T = FOUR + R
         IF (T .EQ. FOUR) GO TO 20
         S = R/T
         U = ONE + TWO*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
C
      RETURN
C
C**** END OF PYTHAG ****************************************************
C
      END
