      SUBROUTINE PURGE (IDXETA,IDXTAU,JL,JLMAX,JT,NVB,LTAU,LNI,NI,N,
     &                  FHNDL,LBLAS,LCOMM,LRERR,LRWRN,NPORTH,NSORTH,
     &                  NTEIG,NQMAX,NXMAX,NBXMAX,BQ,BR,BX,Q,R,X,
     &                  BASIS,BETAQ,BETAR,ANORM,BNORM,WORK,EPS,
     &                  EPS1,REPS,ETA,TAU,ABSETA,ABSTAU,
     &                  NULLDQ,NULLDR,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    PURGE performs reorthogonalizations                              *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    IDXETA (sio) : index of the maximum entry in (ETA)               *
C*    IDXTAU (sio) : index of the maximum entry in (TAU)               *
C*    JL     (sii) : number of steps                                   *
C*    JLMAX  (sii) : maximum number of steps                           *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    NVB    (sii) : number of vectors in a block                      *
C*    LTAU   (sii) : leading dimension of (TAU)                        *
C*    LNI    (sii) : leading dimension of (X)                          *
C*    NI     (sii) : dimension of the vectors in (Q), (R) and (X)      *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    FHNDL  (aii) : file handle                                       *
C*    LBLAS  (sii) : BLAS level setting                                *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    LRERR  (sio) : code for error messages                           *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NPORTH (sib) : number of partial reorthogonalizations performed  *
C*    NSORTH (sib) : number of selective orthogonalizations performed  *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    NQMAX  (sii) : maximum number of vectors in (BASIS)              *
C*    NXMAX  (aii) : maximum number of vectors in (X)                  *
C*    NBXMAX (aii) : maximum number of vectors in (BX)                 *
C*    BQ     (arb) : (B)*(Q)                                           *
C*    BR     (arb) : (B)*(R)                                           *
C*    BX     (arb) : (B)*(X)                                           *
C*    Q      (arb) : Lanczos  vectors at current step                  *
C*    R      (arb) : residual vectors at current step                  *
C*    X      (ari) : eigenvector approximations                        *
C*    BASIS  (ari) : Lanczos vectors array                             *
C*    BETAQ  (arb) : matrix (BETA) in (R)*(Q)*(BETA) at previous step  *
C*    BETAR  (arb) : matrix (BETA) in (R)*(Q)*(BETA) at current  step  *
C*    ANORM  (arb) : extreme singular values of (ALPHA)                *
C*    BNORM  (arb) : extreme singular values of (BETA)                 *
C*    WORK   (arw) : work array                                        *
C*    EPS    (sri) : roundoff unit                                     *
C*    EPS1   (sri) : EPS*NVB*sqrt(N)                                   *
C*    REPS   (sri) : sqrt(EPS)                                         *
C*    ETA    (arb) : orthog. bounds among (R) and Lanczos vectors      *
C*    TAU    (arb) : orthog. bounds among (R) and Ritz    vectors      *
C*    ABSETA (aro) : maximum entry in (ETA)                            *
C*    ABSTAU (aro) : maximum entry in (TAU)                            *
C*    NULLDQ (sio) : number of null pivots in (BETAQ)                  *
C*    NULLDR (sio) : number of null pivots in (BETAR)                  *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    PORTH,REORTH,RFACTR,SORTH                                        *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    ISAMAX                                                           *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    ABS                                                              *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C
      REAL             ZERO
      PARAMETER        (ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          FHNDL(*),IDXETA,IDXTAU,JL,JLMAX,JT,LBLAS(*),
     &                 LCOMM,LNI,LRERR,LRWRN,LTAU,N,NBXMAX,NI,NPORTH,
     &                 NQMAX,NSORTH,NTEIG,NULLDQ,NULLDR,NVB,NXMAX
      REAL             ABSETA,ABSTAU,ANORM(JLMAX,2),BASIS(*),BETAQ(*),
     &                 BETAR(*),BNORM(JLMAX,2),BQ(*),BR(*),BX(*),
     &                 ETA(JLMAX,2),EPS,EPS1,Q(*),R(*),REPS,
     &                 TAU(LTAU,2),WORK(*),X(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          IWRK1,IWRK2,IWRK3,JLM1
      LOGICAL          PORTON,SORTON
C
C==== BLAS kernel ======================================================
C
      INTEGER          ISAMAX
C
C==== intrinsic function ===============================================
C
      INTRINSIC        ABS
C
C**** executable statements ********************************************
C
      IWRK1 = 1
      IWRK2 = IWRK1 + NVB*NVB
      IWRK3 = IWRK2 + NVB*NVB
C
      PORTON = .FALSE.
      SORTON = .FALSE.
      ABSETA = ZERO
      ABSTAU = ZERO
C
C.... selective orthogonalization (Lanczos against Ritz) ...............
C
      IF ( NTEIG .NE. 0 ) THEN
C
C....... find the maximum ETA and check its magnitude ..................
C
         IDXTAU = ISAMAX(NTEIG,TAU(1,2),1)
         ABSTAU = ABS(TAU(IDXTAU,2)) 
C
         IF ( ABSTAU .GT. REPS ) THEN 
C
            SORTON = .TRUE.
C
            CALL SORTH (LBLAS(4),LCOMM,LRERR,LNI,NI,NSORTH,NTEIG,NVB,
     &                  NXMAX,NBXMAX,BR,BX,R,X,EPS1,REPS,TAU(1,1),
     &                  TAU(1,2),WORK(IWRK1),WORK(IWRK2),
     &                  FHNDL,GNRZD)
C
            IF ( LRERR.NE.0 ) RETURN
C
         ELSE
C
            ABSTAU = ZERO
            IDXTAU = 0
C
         END IF
C
      END IF
C
C.... partial reorthogonalization (Lanczos against Lanczos) ............
C
      IF ( JL .GT. 1 ) THEN
C
C....... find the maximum TAU and check its magnitude ..................
C
         IDXETA = ISAMAX(JL,ETA(1,2),1)
         ABSETA = ABS(ETA(IDXETA,2)) 
C
         IF ( ABSETA .GT. REPS ) THEN
C
            PORTON = .TRUE.
C
            CALL PORTH (JL,NVB,NI,LBLAS,LCOMM,LRERR,NPORTH,NQMAX,
     &                  BASIS,BQ,BR,Q,R,EPS1,ETA(1,1),ETA(1,2),
     &                  WORK(IWRK1),WORK(IWRK3),FHNDL,GNRZD)
C
            IF ( LRERR.NE.0 ) RETURN
C
         ELSE
C
            ABSETA = ZERO
            IDXETA = 0
C
         END IF
C
      END IF
C
C.... update factorizations ............................................
C
      IF      ( PORTON ) THEN
C
              JLM1 = JL - 1
C
C............ block with index j .......................................
C
              CALL RFACTR (LBLAS(4),LCOMM,LRERR,LRWRN,NI,NVB,
     &                     NULLDQ,BETAQ,BQ,Q,WORK(IWRK1),WORK(IWRK2),
     &                     WORK(IWRK3),ANORM(JLM1,2),BNORM(JLM1,2),
     &                     BNORM(JLM1,1),EPS,GNRZD)
C
              IF ( LRWRN.NE.0 ) RETURN
C
C............ block with index j+1 .....................................
C
              CALL RFACTR (LBLAS(4),LCOMM,LRERR,LRWRN,NI,NVB,
     &                     NULLDR,BETAR,BR,R,WORK(IWRK1),WORK(IWRK2),
     &                     WORK(IWRK3),ANORM(JL  ,2),BNORM(JL  ,2),
     &                     BNORM(JL  ,1),EPS,GNRZD)
C
              IF ( LRWRN.NE.0 ) RETURN
C
      ELSE IF ( SORTON ) THEN
C
C............ block with index j+1 .....................................
C
              CALL RFACTR (LBLAS(4),LCOMM,LRERR,LRWRN,NI,NVB,
     &                     NULLDR,BETAR,BR,R,WORK(IWRK1),WORK(IWRK2),
     &                     WORK(IWRK3),ANORM(JL  ,2),BNORM(JL  ,2),
     &                     BNORM(JL  ,1),EPS,GNRZD)
C
              IF ( LRWRN.NE.0 ) RETURN
C
      END IF
C
C.... local reorthogonalization (Lanczos against last Lanczos) .........
C
      IF ( PORTON .OR. SORTON .OR. ( N.GT.JT .AND. NVB.GT.1 ) ) THEN
C
         CALL REORTH (LBLAS,LCOMM,LRERR,LRWRN,NI,NVB,NULLDR,BQ,BR,Q,R,
     &                ANORM(JL,2),BNORM(JL,2),BNORM(JL,1),BETAR,
     &                EPS,REPS,WORK,GNRZD)
C
      END IF
C
      RETURN 
C
C**** end of PURGE *****************************************************
C
      END
