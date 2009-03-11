      SUBROUTINE TBBRED (MB,N,A,D,E,Z,I0)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBBRED reduces a real symmetric band matrix to tridiagonal form  *
C*                                                                     *
C*    TBBRED is a modified version of the EISPACK subroutine BANDR     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    MB   (sii) : is the (half) band width of the matrix, defined as  *
C*                 the number of adjacent diagonals, including the     *
C*                 principal diagonal, required to specify the non     *
C*                 zero portion of the lower triangle of the matrix.   *
C*    N    (sii) : is the order of the matrix.                         *
C*    A    (ari) : contains the lower triangle of the symmetric band   *
C*                 input matrix stored as an N by MB array. Its lowest *
C*                 subdiagonal is stored in the last N+1-MB positions  *
C*                 of the first column, its next subdiagonal in the    *
C*                 the last N+2-MB positions of the second column,     *
C*                 further subdiagonals similarly, and finally its     *
C*                 principal diagonal in the N positions of the last   *
C*                 column. Other entries are arbitrary.                *
C*    D    (aro) : contains the diagonal elements of the               *
C*                 tridiagonal matrix.                                 *
C*    E    (aro) : contains the subdiagonal elements of the            *
C*                 tridiagonal matrix in its last N-1 positions.       *
C*    Z    (aro) : contains the orthogonal transformation matrix       *
C*                 produced in the reduction if I0 <= N, otherwise     *
C*                 Z is not referenced.                                *
C*    I0   (sii) : defines the part of (Z) to be computed, Z(I0:N,:).  *
C*                 If I0 > N then Z is not referenced.                 *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,SQRT                                                         *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      DOUBLE PRECISION HALF,ONE,TWO,ZERO
      PARAMETER        (HALF=0.50D0,ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          I0,MB,N
      DOUBLE PRECISION A(N,MB),D(N),E(N),Z(N,N)
C
C==== local variables ==================================================
C
      INTEGER          I1,I2,J,J1,J2,K,KR,L,MR,M1,N2,R,R1,UGL,MAXL,MAXR
      DOUBLE PRECISION B1,B2,C2,DMIN,DMINRT,G,F1,F2,S2,U
      LOGICAL          MATZ
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,SQRT
C
C**** executable statements ********************************************
C
      MATZ = I0.LE.N
      DMIN = TWO**(-64)
      DMINRT = TWO**(-32)
C     
      DO 30 J = 1, N
   30 D(J) = ONE
C
      IF (.NOT. MATZ) GO TO 60
C
      DO 50 J = I0, N
C
         DO 40 K = 1, N
   40    Z(J,K) = ZERO
C
         Z(J,J) = ONE
   50 CONTINUE
C
   60 M1 = MB - 1
      IF (M1 - 1) 900, 800, 70
   70 N2 = N - 2
C
      DO 700 K = 1, N2
C
         MAXR = MIN(M1,N-K)
C
         DO 600 R1 = 2, MAXR
            R = MAXR + 2 - R1
            KR = K + R
            MR = MB - R
            G = A(KR,MR)
            A(KR-1,1) = A(KR-1,MR+1)
            UGL = K
C
            DO 500 J = KR, N, M1
               J1 = J - 1
               J2 = J1 - 1
               IF (G .EQ. ZERO) GO TO 600
               B1 = A(J1,1) / G
               B2 = B1 * D(J1) / D(J)
               IF (ABS(B1) .GT. ONE) THEN
                  U = ONE / B1
                  S2 = U / (U + B2)
               ELSE 
                  S2 = ONE / (ONE + B1 * B2)
               ENDIF
C
               IF (S2 .GE. HALF) GO TO 450
               B1 = G / A(J1,1)
               B2 = B1 * D(J) / D(J1)
               C2 = ONE - S2
               D(J1) = C2 * D(J1)
               D(J) = C2 * D(J)
               F1 = TWO * A(J,M1)
               F2 = B1 * A(J1,MB)
               A(J,M1) = -B2 * (B1 * A(J,M1) - A(J,MB)) - F2 + A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J,MB) + F1) + A(J1,MB)
               A(J,MB) = B1 * (F2 - F1) + A(J,MB)
C
               DO 200 L = UGL, J2
                  I2 = MB - J + L
                  U = A(J1,I2+1) + B2 * A(J,I2)
                  A(J,I2) = -B1 * A(J1,I2+1) + A(J,I2)
                  A(J1,I2+1) = U
  200          CONTINUE
C
               UGL = J
               A(J1,1) = A(J1,1) + B2 * G
               IF (J .EQ. N) GO TO 350
               MAXL = MIN(M1,N-J1)
C
               DO 300 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = A(I1,I2) + B2 * A(I1,I2+1)
                  A(I1,I2+1) = -B1 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2) = U
  300          CONTINUE
C
               I1 = J + M1
               IF (I1 .GT. N) GO TO 350
               G = B2 * A(I1,1)
  350          IF (.NOT. MATZ) GO TO 500
C
               DO 400 L = I0, N
                  U = Z(L,J1) + B2 * Z(L,J)
                  Z(L,J) = -B1 * Z(L,J1) + Z(L,J)
                  Z(L,J1) = U
  400          CONTINUE
C
               GO TO 500
C
  450          U = D(J1)
               D(J1) = S2 * D(J)
               D(J) = S2 * U
               F1 = TWO * A(J,M1)
               F2 = B1 * A(J,MB)
               U = B1 * (F2 - F1) + A(J1,MB)
               A(J,M1) = B2 * (B1 * A(J,M1) - A(J1,MB)) + F2 - A(J,M1)
               A(J1,MB) = B2 * (B2 * A(J1,MB) + F1) + A(J,MB)
               A(J,MB) = U
C
               DO 460 L = UGL, J2
                  I2 = MB - J + L
                  U = B2 * A(J1,I2+1) + A(J,I2)
                  A(J,I2) = -A(J1,I2+1) + B1 * A(J,I2)
                  A(J1,I2+1) = U
  460          CONTINUE
C
               UGL = J
               A(J1,1) = B2 * A(J1,1) + G
               IF (J .EQ. N) GO TO 480
               MAXL = MIN(M1,N-J1)
C
               DO 470 L = 2, MAXL
                  I1 = J1 + L
                  I2 = MB - L
                  U = B2 * A(I1,I2) + A(I1,I2+1)
                  A(I1,I2+1) = -A(I1,I2) + B1 * A(I1,I2+1)
                  A(I1,I2) = U
  470          CONTINUE
C
               I1 = J + M1
               IF (I1 .GT. N) GO TO 480
               G = A(I1,1)
               A(I1,1) = B1 * A(I1,1)
  480          IF (.NOT. MATZ) GO TO 500
C
               DO 490 L = I0, N
                  U = B2 * Z(L,J1) + Z(L,J)
                  Z(L,J) = -Z(L,J1) + B1 * Z(L,J)
                  Z(L,J1) = U
  490          CONTINUE
C
  500       CONTINUE
C
  600    CONTINUE
C
         IF (MOD(K,64) .NE. 0) GO TO 700
C 
         DO 650 J = K, N
            IF (D(J) .GE. DMIN) GO TO 650
            MAXL = MAX(1,MB+1-J)
C
            DO 610 L = MAXL, M1
  610       A(J,L) = DMINRT * A(J,L)
C
            IF (J .EQ. N) GO TO 630
            MAXL = MIN(M1,N-J)
C
            DO 620 L = 1, MAXL
               I1 = J + L
               I2 = MB - L
               A(I1,I2) = DMINRT * A(I1,I2)
  620       CONTINUE
C
  630       IF (.NOT. MATZ) GO TO 645
C
            DO 640 L = I0, N
  640       Z(L,J) = DMINRT * Z(L,J)
C
  645       A(J,MB) = DMIN * A(J,MB)
            D(J) = D(J) / DMIN
  650    CONTINUE
C
  700 CONTINUE
C 
  800 DO 810 J = 2, N
  810 E(J) = SQRT(D(J))
C
      IF (.NOT. MATZ) GO TO 840
C
      DO 830 J = I0, N
C
         DO 820 K = 2, N
  820    Z(J,K) = E(K) * Z(J,K)
C
  830 CONTINUE
C
  840 U = ONE
C
      DO 850 J = 2, N
         A(J,M1) = U * E(J) * A(J,M1)
         U = E(J)
         A(J,MB) = D(J) * A(J,MB)
         D(J) = A(J,MB)
         E(J) = A(J,M1)
  850 CONTINUE
C
      D(1) = A(1,MB)
      E(1) = ZERO
      GO TO 1001
C
  900 DO 950 J = 1, N
         D(J) = A(J,MB)
         E(J) = ZERO
  950 CONTINUE
C
 1001 RETURN
C
C**** end of TBBRED ****************************************************
C
      END
