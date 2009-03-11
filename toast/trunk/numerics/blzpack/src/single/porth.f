      SUBROUTINE PORTH (JL,NVB,NI,LBLAS,LCOMM,LRERR,NPORTH,NQMAX,BASIS,
     &                  BQ,BR,Q,R,EPS1,ETAQ,ETAR,ZETA,WORK,FHNDL,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    PORTH performs a partial reorthogonalization                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL     (sii) : number of steps                                   *
C*    NVB    (sii) : number of vectors in a block                      *
C*    NI     (sii) : dimension of the vectors in (Q) and (R)           *
C*    LBLAS  (aii) : BLAS level setting                                *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    LRERR  (sio) : code for error messages                           *
C*    NPORTH (sib) : number of partial reorthogonalizations performed  *
C*    NQMAX  (sii) : maximum number of vectors in (BASIS)              *
C*    BASIS  (ari) : Lanczos vectors array                             *
C*    BQ     (arb) : (B)*(Q)                                           *
C*    BR     (arb) : (B)*(R)                                           *
C*    Q      (arb) : Lanczos  vectors at current step                  *
C*    R      (arb) : residual vectors at current step                  *
C*    EPS1   (sri) : EPS*NVB*sqrt(N)                                   *
C*    ETAQ   (sri) : orthog. bounds among (Q) and Lanczos vectors      *
C*    ETAR   (sri) : orthog. bounds among (R) and Lanczos vectors      *
C*    ZETA   (arw) : work array                                        *
C*    WORK   (arw) : work array                                        *
C*    FHNDL  (aii) : file handle                                       *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZIOOP,QTBR,RQALPH                                               *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          FHNDL(*),JL,LBLAS(*),LCOMM,LRERR,
     &                 NI,NPORTH,NQMAX,NVB
      REAL             BASIS(*),BQ(NI,NVB),BR(NI,NVB),EPS1,ETAQ(*),
     &                 ETAR(*),Q(NI,NVB),R(NI,NVB),ZETA(*),WORK(*)
      LOGICAL GNRZD
C
C==== local variables ==================================================
C
      INTEGER          IBQ,IBQOUT,IQ,IQOUT,IZETA1,IZETA2,J,LDB
      REAL             DUMMY
C
C**** executable statements ********************************************
C
      NPORTH = NPORTH + 1
C
C.... set pointers .....................................................
C
      IQ = 1 
C
      IF ( GNRZD ) THEN
         LDB = NI*2 
         IBQ = IQ + NI
         IQOUT = IQ + NI*NVB*NQMAX*2
         IBQOUT = IQOUT + NI*NVB
      ELSE
         LDB = NI
         IBQ = IQ
         IQOUT = IQ + NI*NVB*NQMAX
         IBQOUT = IQOUT
      END IF
C
      IZETA1 = 1
      IZETA2 = IZETA1 + NVB*NVB
C
C.... loop on all Lanczos vectors ......................................
C
      DO 10 J = 1,JL-1
C
         ETAQ(J) = EPS1
         ETAR(J) = EPS1
C
C....... retrieve the corresponding old (B)*(Q) ........................
C
         IF ( J .GT. NQMAX ) THEN
            IQ = IQOUT
            IBQ = IBQOUT
            LDB = NI
            CALL LZIOOP (FHNDL,LCOMM,LRERR,J-NQMAX,NI*NVB,BASIS(IQ),
     &                   'Q ','GET')
            IF ( LRERR .NE. 0 ) RETURN
            IF ( GNRZD ) THEN
               CALL LZIOOP (FHNDL,LCOMM,LRERR,J-NQMAX,NI*NVB,BASIS(IBQ),
     &                      'BQ','GET')
               IF ( LRERR .NE. 0 ) RETURN
            END IF
	 END IF
C
C....... multiplication of old (B)*(Q) transpose by (Q) and (R) ........
C
         CALL QTBR (NI,BASIS(IBQ),LDB,NVB,Q,NI,NVB,ZETA(IZETA1),
     &              DUMMY,WORK,LBLAS(2),LCOMM,LRERR)
         CALL QTBR (NI,BASIS(IBQ),LDB,NVB,R,NI,NVB,ZETA(IZETA2),
     &              DUMMY,WORK,LBLAS(2),LCOMM,LRERR)
C
C....... orthogonalization of (Q) and (R) against old (Q) ..............
C
         CALL RQALPH (NI,Q,NI,NVB,BASIS(IQ),LDB,NVB,
     &                ZETA(IZETA1),LBLAS(3))
         CALL RQALPH (NI,R,NI,NVB,BASIS(IQ),LDB,NVB,
     &                ZETA(IZETA2),LBLAS(3))
C
C....... orthogonalization of (BQ) and (BR) against old (B)*(R) ........
C
         IF ( GNRZD ) THEN
            CALL RQALPH (NI,BQ,NI,NVB,BASIS(IBQ),LDB,NVB,
     &                   ZETA(IZETA1),LBLAS(3))
            CALL RQALPH (NI,BR,NI,NVB,BASIS(IBQ),LDB,NVB,
     &                   ZETA(IZETA2),LBLAS(3))
         END IF
C
C....... increment pointers ............................................
C
         IF ( J .LT. NQMAX ) THEN
            IF ( GNRZD ) THEN
               IQ = IQ + NI*NVB*2
               IBQ = IBQ + NI*NVB*2
            ELSE
               IQ = IQ + NI*NVB
               IBQ = IBQ + NI*NVB
            END IF
         END IF
C
   10 CONTINUE
C
      ETAQ(JL) = EPS1
      ETAR(JL) = EPS1
C
      RETURN 
C
C**** end of PORTH *****************************************************
C
      END
