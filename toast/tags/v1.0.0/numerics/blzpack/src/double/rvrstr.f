      SUBROUTINE RVRSTR (JT,NR,NVB,BIGNUM,REPS,ENDL,ENDR,RITZ,RNORM,
     &                   S,SSUM,SIGMA,IWORK,RWORK,GNRZD,SLICE)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RVRSTR identifies the quasi-converged Ritz vectors               *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    NR     (sio) : number of starting vectors available              *
C*    NVB    (sii) : number of vectors in a block                      *
C*    BIGNUM (sri) : big number                                        *
C*    REPS   (sri) : sqrt(EPS)                                         *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    RITZ   (ari) : Ritz values                                       *
C*    RNORM  (ari) : estimated residuals                               *
C*    S      (ari) : eigenvectors of the tridiagonal matrix            *
C*    SSUM   (aro) : sum of selected columns of (S)                    *
C*    SIGMA  (sri) : origin translation                                *
C*    IWORK  (arw) : work array (integer)                              *
C*    RWORK  (arw) : work array (real)                                 *
C*    GNRZD  (sli) : problem type flag                                 *
C*    SLICE  (sli) : spectrum slicing flag                             *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    DAXPY                                                            *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,MAX                                                          *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      DOUBLE PRECISION ONE,ZERO
      PARAMETER        (ONE=1.0D0,ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          IWORK(*),JT,NR,NVB
      DOUBLE PRECISION BIGNUM,ENDL,ENDR,REPS,RITZ(JT),RNORM(JT),
     &                 S(JT,JT),SIGMA,SSUM(JT,NVB),RWORK(*)
      LOGICAL          GNRZD,SLICE
C
C==== local variables ==================================================
C
      INTEGER          I,J,K,L
      DOUBLE PRECISION LOWER,SSCALE,SSUMJT,UPPER,WMIN
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,MAX,MIN
C
C**** executable statements ********************************************
C
      NR = 0
C
      DO 10 I = 1,NVB
         IWORK(I) = 0
  10  CONTINUE
C
      CALL SETTO0 (JT*NVB,SSUM,1)
C
C.... define bounds for the Ritz values ................................
C
      IF ( SLICE ) THEN
         LOWER = ENDL
         UPPER = ENDR 
      ELSE
         LOWER = -BIGNUM
         UPPER = +BIGNUM
      END IF
C
C.... criterion for identifying the vectors ............................
C
      IF ( GNRZD ) THEN
         DO 20 I = 1,JT
            RWORK(I) = ABS(RITZ(I)-SIGMA)
   20    CONTINUE
      ELSE
         DO 30 I = 1,JT
            RWORK(I) = ABS(RNORM(I))
   30    CONTINUE
      END IF
C
C.... sum up the quasi-converged vectors ...............................
C
      DO 60 I = 1,2
         DO 50 J = 1,NVB
            L = 0
            WMIN = BIGNUM
            DO 40 K = 1,JT 
               IF ( RWORK(K).LE.WMIN .AND. RNORM(K).LT.ZERO .AND.
     &              RITZ(K).GE.LOWER .AND. RITZ(K).LE.UPPER ) THEN
                    WMIN = RWORK(K)
                    L = K
               END IF
   40       CONTINUE
            IF ( L .NE. 0 ) THEN
               IF ( GNRZD ) THEN
                  SSCALE = ONE/MAX(REPS,RWORK(L))
               ELSE
                  SSCALE = MIN(RITZ(L),ONE/REPS)
               END IF
               CALL DAXPY (JT,SSCALE,S(1,L),1,SSUM(1,J),1)
               IWORK(J) = IWORK(J) + 1
               RNORM(L) = ZERO
               NR = MAX(J,NR)
            END IF
   50    CONTINUE
   60 CONTINUE
C
C.... if needed add perturbation to vectors ............................
C
      DO 70 I = 1,NVB
         IF ( IWORK(I) .EQ. 1 ) THEN
            SSUMJT = SSUM(JT,I)
            DO 80 J = 1,JT
               SSUM(J,I) = SSUM(J,I) + SSUMJT
   80       CONTINUE
         END IF
   70 CONTINUE
C
C.... mark the entries of RNORM to be taken by RVMNGR ..................
C
      DO 90 I = 1,NR
         RNORM(I) = ZERO
   90 CONTINUE
      DO 100 I = NR+1,JT
         RNORM(I) = -ONE
  100 CONTINUE
C
      RETURN 
C
C**** end of RVRSTR ****************************************************
C
      END
