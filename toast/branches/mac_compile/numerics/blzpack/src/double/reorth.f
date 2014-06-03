      SUBROUTINE REORTH (LBLAS,LCOMM,LRERR,LRWRN,NI,NVB,NULLP,BQ,BR,Q,R,
     &                   AMAXN,BMAXN,BMINN,BETA,EPS,REPS,WORK,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    REORTH performs a local reorthogonalization                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LBLAS (aii) : BLAS level setting                                 *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    LRWRN (sio) : code for warning messages                          *
C*    NI    (sii) : dimension of the vectors in (Q) and (R)            *
C*    NVB   (sii) : number of vectors in a block                       *
C*    NULLP (sio) : number of null pivots in (BETA)                    *
C*    BQ    (arb) : (B)*(Q)                                            *
C*    BR    (arb) : (B)*(R)                                            *
C*    Q     (arb) : Lanczos  vectors at current step                   *
C*    R     (arb) : residual vectors at current step                   *
C*    AMAXN (sro) : maximum singular value of (ALPHA)                  *
C*    BMAXN (sro) : maximum singular value of (BETA)                   *
C*    BMINN (sro) : minimum singular value of (BETA)                   *
C*    BETA  (arb) : matrix (BETA) in (R)=(Q)*(BETA)                    *
C*    EPS   (sri) : roundoff unit                                      *
C*    REPS  (sri) : sqrt(EPS)                                          *
C*    WORK  (arw) : work array                                         *
C*    GNRZD (sli) : problem type flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    QTBR,RFACTR,RQALPH                                               *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LBLAS(*),LCOMM,LRERR,LRWRN,NI,NULLP,NVB
      DOUBLE PRECISION AMAXN,BETA(*),BMAXN,BMINN,BQ(*),BR(*),
     &                 EPS,Q(*),R(*),REPS,WORK(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          I,IWORK1,IWORK2,IWORK3
      DOUBLE PRECISION WMAX
C
C**** executable statements ********************************************
C
      IWORK1 = 1
      IWORK2 = IWORK1 + NVB*NVB
      IWORK3 = IWORK2 + NVB*NVB
C
      DO 10 I = 1,NVB
C
C....... reorthogonalization of (R) against (Q) ........................
C
         IF ( GNRZD ) THEN
            CALL QTBR   (NI,BQ,NI,NVB,R,NI,NVB,WORK(IWORK1),WMAX,
     &                   WORK(IWORK2),LBLAS(2),LCOMM,LRERR)
            CALL RQALPH (NI,BR,NI,NVB,BQ,NI,NVB,WORK,LBLAS(3))
            CALL RQALPH (NI, R,NI,NVB, Q,NI,NVB,WORK,LBLAS(3))
         ELSE
            CALL QTBR   (NI, Q,NI,NVB,R,NI,NVB,WORK(IWORK1),WMAX,
     &                   WORK(IWORK2),LBLAS(2),LCOMM,LRERR)
            CALL RQALPH (NI, R,NI,NVB, Q,NI,NVB,WORK,LBLAS(3))
         END IF
C
C....... (R)=(Q)*(BETA) updating .......................................
C
         CALL RFACTR (LBLAS(4),LCOMM,LRERR,LRWRN,NI,NVB,NULLP,BETA,
     &                BR,R,WORK(IWORK1),WORK(IWORK2),WORK(IWORK3),
     &                AMAXN,BMAXN,BMINN,EPS,GNRZD)
C
         IF ( NULLP.GT.0 .OR. WMAX.LE.REPS .OR.
     &        LRERR.GT.0 .OR. LRWRN.GT.0 ) RETURN   
C
   10 CONTINUE
C  
      RETURN 
C
C**** end of REORTH ****************************************************
C
      END
