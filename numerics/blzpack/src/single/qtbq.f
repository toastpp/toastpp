      SUBROUTINE QTBQ (JL,JT,NVB,NI,LFILE,BASIS,ORTH,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    QTBQ computes and prints the basis orthogonality level           *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL    (sii) : number of steps                                    *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    NVB   (sii) : number of vectors in a block                       *
C*    NI    (sii) : dimension of the vectors in (BASIS)                *
C*    LFILE (sii) : file unit for output                               *
C*    BASIS (ari) : Lanczos vectors array                              *
C*    ORTH  (arw) : stores (Q')*(B)*(Q)                                *
C*    GNRZD (sli) : problem type flag                                  *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    LZPRT5                                                           *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    SGEMM                                                            *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,LOG10                                                        *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      REAL             ONE,ZERO
      PARAMETER        (ONE=1.0E0,ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          JL,JT,LFILE,NI,NVB
      REAL             BASIS(*),ORTH(JT,JT)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          I,J,K,L
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,LOG10
C
C**** executable statements ********************************************
C
      IF ( GNRZD ) THEN   
C
C....... multiply (Q')*(B)*(Q) .........................................
C
         K = 1
         DO 20 I = 1,JT,NVB
            L = NI + 1
            DO 10 J = 1,JT,NVB
               CALL SGEMM ('T','N',NVB,NVB,NI,ONE,BASIS(K),NI*2,
     &                     BASIS(L),NI*2,ZERO,ORTH(I,J),JT)
               L = L + NI*NVB*2
   10       CONTINUE
            K = K + NI*NVB*2
   20    CONTINUE
C
      ELSE
C
C....... multiply (Q')*(Q) .............................................
C
         K = 1
         DO 40 I = 1,JT,NVB
            L = 1
            DO 30 J = 1,JT,NVB
               CALL SGEMM ('T','N',NVB,NVB,NI,ONE,BASIS(K),NI,
     &                     BASIS(L),NI,ZERO,ORTH(I,J),JT)
               L = L + NI*NVB
   30       CONTINUE
            K = K + NI*NVB
   40    CONTINUE
C
      END IF
C
C.... compute |log10(|ORTH|)| ..........................................
C
      DO 60 I = 1,JT
         DO 50 J = 1,JT
            IF ( ORTH(J,I) .NE. ZERO ) THEN
               ORTH(J,I) = ABS(LOG10(ABS(ORTH(J,I))))
            END IF
   50    CONTINUE
   60 CONTINUE
C
C.... print the orthogonality level ....................................
C
      CALL LZPRT5 (JL,JT,NVB,LFILE,ORTH)
C
      RETURN 
C
C**** end of QTBQ ******************************************************
C
      END
