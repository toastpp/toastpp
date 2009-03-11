      SUBROUTINE RQBETA (R,Q,BETA,LBLAS,N,NVB)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RQBETA computes (R)=(R)-(Q)*(BETA')                              *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    R     (arb) : residual vectors                                   *
C*    Q     (ari) : Lanczos  vectors at previous step                  *
C*    BETA  (ari) : matrix (BETA) in (R)=(Q)*(BETA) at step JL-1       *
C*    LBLAS (sii) : BLAS level setting                                 *
C*    N     (sii) : dimension of the vectors in (Q) and (R)            *
C*    NVB   (sii) : number of vectors in a block                       *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SAXPY,SGEMV                                                      *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C 
      REAL             ONE
      PARAMETER        (ONE=1.0E0)
C 
C==== arguments ========================================================
C
      INTEGER          LBLAS,N,NVB
      REAL             BETA(NVB,NVB),Q(N,NVB),R(N,NVB)
C
C==== local variables ==================================================
C 
      INTEGER          I,J
C
C**** executable statements ********************************************
C
      IF      ( LBLAS .EQ. 1 ) THEN
C
C............ (R)=(R)-(Q)*(BETA') using BLAS 1 .........................
C
              DO 20 I = 1,NVB
                 DO 10 J = I,NVB
                    CALL SAXPY (N,-BETA(I,J),Q(1,J),1,R(1,I),1)
   10            CONTINUE
   20         CONTINUE
C
      ELSE IF ( LBLAS .GE. 2 ) THEN
C
C............ (R)=(R)-(Q)*(BETA') using BLAS 2 .........................
C
              J = NVB
C
              DO 30 I = 1,NVB
                 CALL SGEMV ('N',N,J,-ONE,Q(1,I),N,BETA(I,I),
     &                       NVB,ONE,R(1,I),1)
                 J = J - 1
   30         CONTINUE
C
      END IF 
C
      RETURN
C
C**** end of RQBETA ****************************************************
C
      END
