      SUBROUTINE RVCOMP (JT,LBLAS,LDB,LNI,NI,NQ,NX,BASIS,S,X)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RVCOMP computes the Ritz vectors (X)=(X)+(BASIS)*(S)             *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    LBLAS (sii) : BLAS level setting                                 *
C*    LNB   (sii) : leading dimension of (BASIS)                       *
C*    LNI   (sii) : leading dimension of (X)                           *
C*    NI    (sii) : dimension of the vectors in (BASIS) and (X)        *
C*    NQ    (sii) : number of vectors in (BASIS)                       *
C*    NX    (sii) : number of vectors in (X)                           *
C*    BASIS (ari) : Lanczos vectors array                              *
C*    S     (ari) : Ritz coordinates                                   *
C*    X     (aro) : Ritz vectors                                       *
C*                                                                     *
C*  - BLAS kernels:                                                    *
C*                                                                     *
C*    DAXPY,DGEMM,DGEMV                                                *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C 
      DOUBLE PRECISION ONE
      PARAMETER        (ONE=1.0D0)
C
C==== arguments ========================================================
C
      INTEGER          JT,LBLAS,LDB,LNI,NI,NQ,NX
      DOUBLE PRECISION BASIS(LDB,NQ),S(JT,NX),X(LNI,NX)
C
C==== local variables ==================================================
C 
      INTEGER          I,J
C
C**** executable statements ********************************************
C
      IF      ( LBLAS .EQ. 1 ) THEN
C
C............ (X)=(X)+(BASIS)*(S) using BLAS 1 .........................
C
              DO 20 I = 1,NX
                 DO 10 J = 1,NQ
                    CALL DAXPY (NI,S(J,I),BASIS(1,J),1,X(1,I),1)
   10            CONTINUE
   20         CONTINUE
C
      ELSE IF ( LBLAS .EQ. 2 ) THEN
C
C............ (X)=(X)+(BASIS)*(S) using BLAS 2 .........................
C
              DO 30 I = 1,NX
                 CALL DGEMV ('N',NI,NQ,ONE,BASIS,LDB,
     &                       S(1,I),1,ONE,X(1,I),1)
   30         CONTINUE
C
      ELSE IF ( LBLAS .EQ. 3 ) THEN
C
C............ (X)=(X)+(BASIS)*(S) using BLAS 3 .........................
C
              CALL DGEMM ('N','N',NI,NX,NQ,ONE,BASIS,LDB,S,JT,ONE,X,LNI)
C
      END IF 
C
      RETURN
C
C**** end of RVCOMP ****************************************************
C
      END
