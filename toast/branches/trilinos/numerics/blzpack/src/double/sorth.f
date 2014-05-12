      SUBROUTINE SORTH (LBLAS,LCOMM,LRERR,LNI,NI,NSORTH,NTEIG,NVB,
     &                  NXMAX,NBXMAX,BR,BX,R,X,EPS1,REPS,TAUQ,
     &                  TAUR,ZETA,WORK,FHNDL,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SORTH performs a selective orthogonalization                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LBLAS  (sii) : BLAS level setting                                *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    LRERR  (sio) : code for error messages                           *
C*    LNI    (sii) : leading dimension of (X)                          *
C*    NI     (sii) : dimension of the vectors in (R) and (X)           *
C*    NSORTH (sib) : number of selective orthogonalizations performed  *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    NVB    (sii) : number of vectors in a block                      *
C*    NXMAX  (aii) : maximum number of vectors in (X)                  *
C*    NBXMAX (aii) : maximum number of vectors in (BX)                 *
C*    BR     (arb) : (B)*(R)                                           *
C*    BX     (arb) : (B)*(X)                                           *
C*    R      (arb) : residual vectors at current step                  *
C*    X      (arb) : eigenvector approximations                        *
C*    EPS1   (sri) : EPS*NVB*sqrt(N)                                   *
C*    REPS   (sri) : sqrt(EPS)                                         *
C*    TAUQ   (arb) : orthog. bounds among (Q) and Ritz vectors         *
C*    TAUR   (arb) : orthog. bounds among (R) and Ritz vectors         *
C*    ZETA   (arw) : work array                                        *
C*    WORK   (arw) : work array                                        *
C*    FHNDL  (aii) : file handle                                       *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZIOOP,QTBR,RQALPH                                               *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    ABS                                                              *
C*                                                                     *
C***********************************************************************
C
C==== parameter ========================================================
C 
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          FHNDL(*),LBLAS,LCOMM,LNI,LRERR,NBXMAX,NI,
     &                 NSORTH,NTEIG,NVB,NXMAX
      DOUBLE PRECISION BR(NI,NVB),BX(NI,*),EPS1,R(NI,NVB),REPS,
     &                 TAUQ(*),TAUR(*),X(LNI,*),ZETA(*),WORK(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          I,IBX,IX
      DOUBLE PRECISION ZMAX
C
C==== intrinsic function ===============================================
C
      INTRINSIC        ABS
C
C**** executable statements ********************************************
C
C.... loop on all eigenvectors .........................................
C
      DO 10 I = 1,NTEIG
C
         IF ( ABS(TAUR(I)) .LT. REPS ) GO TO 10        
C
         NSORTH = NSORTH + 1
C
C....... update (TAUQ) and (TAUR) ......................................
C
         IF ( TAUR(I) .LT. ZERO ) THEN
            TAUQ(I) = +EPS1
            TAUR(I) = +EPS1
         ELSE
            TAUQ(I) = -TAUQ(I)
            TAUR(I) = -TAUR(I)
         END IF
C
C....... retrieve the corresponding (B)*(X) ............................
C
         IF ( GNRZD ) THEN
            IF ( NBXMAX .GE. I ) THEN
               IBX = I
            ELSE
               IBX = NBXMAX + 1
               CALL LZIOOP (FHNDL,LCOMM,LRERR,I-NBXMAX,
     &                      NI,BX(1,IBX),'BX','GET')
               IF ( LRERR .NE. 0 ) RETURN
            END IF
         END IF
C
C....... retrieve the corresponding (X) ................................
C
         IF ( NXMAX .GT. 0 ) THEN
            IX = I
         ELSE
            IX = 1
            CALL LZIOOP (FHNDL,LCOMM,LRERR,I,NI,X(1,IX),'X ','GET')
            IF ( LRERR .NE. 0 ) RETURN
         END IF
C
C....... orthogonalization step ........................................
C
         IF ( GNRZD ) THEN 
            CALL QTBR   (NI,BX(1,IBX),NI,1,R,NI,NVB,ZETA,
     &                   ZMAX,WORK,LBLAS,LCOMM,LRERR)
            CALL RQALPH (NI,BR,NI,NVB,BX(1,IBX),NI,1,ZETA,LBLAS)
            CALL RQALPH (NI,R ,NI,NVB,X(1,IX)  ,NI,1,ZETA,LBLAS)
         ELSE
            CALL QTBR   (NI, X(1,IX) ,NI,1,R,NI,NVB,ZETA,
     &                   ZMAX,WORK,LBLAS,LCOMM,LRERR)
            CALL RQALPH (NI,R ,NI,NVB,X(1,IX)  ,NI,1,ZETA,LBLAS)
         END IF
C
   10 CONTINUE
C
      RETURN 
C
C**** end of SORTH *****************************************************
C
      END
