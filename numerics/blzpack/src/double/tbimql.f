      SUBROUTINE TBIMQL (N,D,E,Z,I0,INFO)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBIMQL finds the eigenpairs of a symmetric tridiagonal matrix    *
C*                                                                     *
C*    TBIMQL is a modified version of the EISPACK subroutine IMTQL2    *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    N    (sii) : is the order of the matrix.                         *
C*    D    (arb) : contains the diagonal elements of the input matrix. *
C*                 On output it contains the eigenvalues in ascending  *
C*                 order.                                              *
C*    E    (arb) : contains the subdiagonal elements of the input      *
C*                 matrix in its last N-1 positions. E(1) is not used. *
C*                 On output E has been destroyed.                     *
C*    Z    (arb) : contains the transformation matrix produced in the  *
C*                 reduction to tridiagonal form.                      *
C*    I0   (sii) : defines the part of (Z) to be computed, Z(I0:N,:).  *
C*                 If I0 > N then Z is not computed.                   *
C*    INFO (sio) : is set to zero for normal return or j if the j-th   *
C*                 eigenvalue has not been determined after 30 steps.  *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    PYTHAG                                                           *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,SIGN                                                         *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      DOUBLE PRECISION ONE,TWO,ZERO
      PARAMETER        (ONE=1.0D0,TWO=2.0D0,ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          I0,INFO,N
      DOUBLE PRECISION D(N),E(N),Z(N,N) 
C
C==== local variables ==================================================
C
      INTEGER          I,J,K,L,M,II,MML
      DOUBLE PRECISION B,C,F,G,P,R,S,TST1,TST2
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION PYTHAG        
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,SIGN
C
C**** executable statements ********************************************

      INFO = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      E(N) = ZERO
C
      DO 240 L = 1, N
         J = 0
C 
  105    DO 110 M = L, N
            IF (M .EQ. N) GO TO 120
            TST1 = ABS(D(M)) + ABS(D(M+1))
            TST2 = TST1 + ABS(E(M))
            IF (TST2 .EQ. TST1) GO TO 120
  110    CONTINUE
C
  120    P = D(L)
         IF (M .EQ. L) GO TO 240
         IF (J .EQ. 30) GO TO 1000
         J = J + 1
C
         G = (D(L+1) - P) / (TWO * E(L))
         R = PYTHAG(G,ONE)
         G = D(M) - P + E(L) / (G + SIGN(R,G))
         S = ONE
         C = ONE
         P = ZERO
         MML = M - L
C
         DO 200 II = 1, MML
            I = M - II
            F = S * E(I)
            B = C * E(I)
            R = PYTHAG(F,G)
            E(I+1) = R
            IF (R .EQ. ZERO) GO TO 210
            S = F / R
            C = G / R
            G = D(I+1) - P
            R = (D(I) - G) * S + TWO * C * B
            P = S * R
            D(I+1) = G + P
            G = C * R - B
C
            DO 180 K = I0, N
               F = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * F
               Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
C
  200    CONTINUE
C
         D(L) = D(L) - P
         E(L) = G
         E(M) = ZERO
         GO TO 105
C
  210    D(I+1) = D(I+1) - P
         E(M) = ZERO
         GO TO 105
  240 CONTINUE
C
      DO 300 II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 260 J = II, N
            IF (D(J) .GE. P) GO TO 260
            K = J
            P = D(J)
  260    CONTINUE
C
         IF (K .EQ. I) GO TO 300
         D(K) = D(I)
         D(I) = P
C
         DO 280 J = I0, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  280    CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C
 1000 INFO = L
 1001 RETURN
C
C**** END OF TBIMQL*****************************************************
C
      END
