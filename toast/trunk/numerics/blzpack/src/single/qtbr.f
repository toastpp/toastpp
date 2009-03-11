      SUBROUTINE QTBR (N,BQ,LDBQ,NCBQ,R,LDR,NCR,ALPHA,AMAX,
     &                 WORK,LBLAS,LCOMM,LRERR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    QTBR computes (ALPHA)=(Q')*(B)*(R)                               *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    N     (sii) : dimension of the vectors in (BQ) and (R)           *
C*    BQ    (ari) : (B)*(Q)                                            *
C*    LDBQ  (sii) : leading dimension of (BQ)                          *
C*    NCBQ  (sii) : number of columns in (BQ)                          *
C*    R     (ari) : residual vectors                                   *
C*    LDR   (sii) : leading dimension of (R)                           *
C*    NCR   (sii) : number of columns in (R)                           *
C*    ALPHA (aro) : (Q')*(B)*(R)                                       *
C*    AMAX  (sro) : largest absolute value in (ALPHA)                  *
C*    WORK  (arw) : workspace                                          *
C*    LBLAS (sii) : BLAS level setting                                 *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    PISRED                                                           *
C*                                                                     *
C*  - BLAS kernels:                                                    *
C*                                                                     *
C*    SDOT,SGEMM,SGEMV                                                 *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    ABS                                                              *
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
      INTEGER          LBLAS,LCOMM,LDBQ,LDR,LRERR,N,NCBQ,NCR 
      REAL             ALPHA(NCBQ,NCR),AMAX,BQ(LDBQ,NCBQ),
     &                 R(LDR,NCR),WORK(NCBQ,NCR)
C
C==== local variables ==================================================
C
      INTEGER          I,INFO,J
C
C==== BLAS kernel ======================================================
C
      REAL             SDOT
C
C==== intrinsic function ===============================================
C
      INTRINSIC        ABS
C
C**** executable statements ********************************************
C
      IF      ( LBLAS .EQ. 1 ) THEN
C
C............ (Q')*(B)*(R) using BLAS 1 ................................
C
              DO 20 I = 1,NCR
                 DO 10 J = 1,NCBQ
                    WORK(J,I) = SDOT(N,BQ(1,J),1,R(1,I),1)
   10            CONTINUE
   20         CONTINUE
C
      ELSE IF ( LBLAS .EQ. 2 ) THEN
C
C............ (Q')*(B)*(R) using BLAS 2 ................................
C
              DO 30 I = 1,NCR
                 CALL SGEMV ('T',N,NCBQ,ONE,BQ,LDBQ,R(1,I),1,
     &                       ZERO,WORK(1,I),1)
   30         CONTINUE
C
      ELSE IF ( LBLAS .EQ. 3 ) THEN
C
C............ (Q')*(B)*(R) using BLAS 3 ................................
C
              CALL SGEMM ('T','N',NCBQ,NCR,N,ONE,BQ,LDBQ,
     &                    R,LDR,ZERO,WORK,NCBQ)
C
      END IF
C
      CALL PISRED ('SUM',NCBQ*NCR,WORK,ALPHA,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
C.... largest entry ....................................................
C
      AMAX = ALPHA(1,1)
C
      DO 50 I = 1,NCR
         DO 40 J = 1,NCBQ
            IF ( ABS(ALPHA(J,I)) .GT. AMAX ) AMAX = ABS(ALPHA(J,I))
   40    CONTINUE
   50 CONTINUE
C
      RETURN
C
C**** end of QTBR ******************************************************
C
      END
