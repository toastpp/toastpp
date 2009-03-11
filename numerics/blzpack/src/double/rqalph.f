      SUBROUTINE RQALPH (N,R,LDR,NCR,Q,LDQ,NCQ,ALPHA,LBLAS)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RQALPH computes (R)=(R)-(Q)*(ALPHA)                              *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    N     (sii) : dimension of the vectors in (Q) and (R)            *
C*    R     (arb) : residual vectors                                   *
C*    LDR   (sii) : leading dimension of (R)                           *
C*    NCR   (sii) : number of columns in (R)                           *
C*    Q     (ari) : Lanczos  vectors                                   *
C*    LDQ   (sii) : leading dimension of (Q)                           *
C*    NCQ   (sii) : number of columns in (Q)                           *
C*    ALPHA (ari) : (Q')*(B)*(R)                                       *
C*    LBLAS (sii) : BLAS level setting                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
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
      INTEGER          LBLAS,LDQ,LDR,N,NCQ,NCR
      DOUBLE PRECISION ALPHA(NCQ,NCR),Q(LDQ,NCQ),R(LDR,NCR)
C
C==== local variables ==================================================
C 
      INTEGER          I,J
C
C**** executable statements ********************************************
C
      IF      ( LBLAS .EQ. 1 ) THEN
C
C............ (R)=(R)-(Q)*(ALPHA) using BLAS 1 .........................
C
              DO 20 I = 1,NCR
                 DO 10 J = 1,NCQ
                    CALL DAXPY (N,-ALPHA(J,I),Q(1,J),1,R(1,I),1)
   10            CONTINUE
   20         CONTINUE
C
      ELSE IF ( LBLAS .EQ. 2 ) THEN
C
C............ (R)=(R)-(Q)*(ALPHA) using BLAS 2 .........................
C
              DO 30 I = 1,NCR
                 CALL DGEMV ('N',N,NCQ,-ONE,Q,LDQ,ALPHA(1,I),
     &                       1,ONE,R(1,I),1)
   30         CONTINUE
C
      ELSE IF ( LBLAS .EQ. 3 ) THEN
C
C............ (R)=(R)-(Q)*(ALPHA) using BLAS 3 .........................
C
              CALL DGEMM ('N','N',N,NCR,NCQ,-ONE,Q,LDQ,
     &                    ALPHA,NCQ,ONE,R,LDR)
C
      END IF 
C
      RETURN
C
C**** end of RQALPH ****************************************************
C
      END
