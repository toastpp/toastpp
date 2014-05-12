      SUBROUTINE SORTH0 (LBLAS,LCOMM,LRERR,LEIG,LNI,NI,NSORTH,NTEIG,
     &                   NVB,NXMAX,NBXMAX,BX,R,X,EIG,ENDL,ENDR,
     &                   EPS1,TAUQ,TAUR,ZETA,WORK,FHNDL,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SORTH0 purges the starting vectors and sets TAU                  *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LBLAS  (sii) : BLAS level setting                                *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    LRERR  (sio) : code for error messages                           *
C*    LEIG   (sii) : leading dimension of (EIG)                        *
C*    LNI    (sii) : leading dimension of (X)                          *
C*    NI     (sii) : dimension of the vectors in (R) and (X)           *
C*    NSORTH (sib) : number of selective orthogonalizations performed  *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    NVB    (sii) : number of vectors in a block                      *
C*    NXMAX  (aii) : maximum number of vectors in (X)                  *
C*    NBXMAX (aii) : maximum number of vectors in (BX)                 *
C*    BX     (arb) : (B)*(X)                                           *
C*    R      (arb) : residual vectors at current step                  *
C*    X      (arb) : eigenvector approximations                        *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    EPS1   (sri) : EPS*NVB*sqrt(N)                                   *
C*    TAUQ   (arb) : orthog. bounds among (Q) and Ritz vectors         *
C*    TAUR   (arb) : orthog. bounds among (R) and Ritz vectors         *
C*    ZETA   (arw) : workspace                                         *
C*    WORK   (arw) : workspace                                         *
C*    FHNDL  (aii) : file handle                                       *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZIOOP,QTBR,RQALPH                                               *
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
      INTEGER          FHNDL(*),LBLAS,LCOMM,LEIG,LNI,LRERR,NI,NSORTH,
     &                 NTEIG,NVB,NXMAX,NBXMAX
      DOUBLE PRECISION BX(NI,*),EIG(LEIG,2),ENDL,ENDR,EPS1,R(NI,NVB),
     &                 TAUQ(*),TAUR(*),X(LNI,*),WORK(*),ZETA(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          I,IBX,IX
      DOUBLE PRECISION ZMAX 
C
C**** executable statements ********************************************
C
C.... check all vectors in the neighbourhood ...........................
C
      DO 10 I = 1,NTEIG
C
         IF ( EIG(I,2).EQ.ZERO ) EIG(I,2) = EPS1
C
         IF ( EIG(I,1).GE.ENDL .AND. EIG(I,1).LE.ENDR ) THEN
C
              NSORTH = NSORTH + 1
C
C............ this vector is a candidate for orthogonalization .........
C
              TAUR(I) = -ONE
              TAUQ(I) = -ONE
C
C............ retrieve the corresponding (B)*(X) .......................
C
              IF ( GNRZD ) THEN
                 IF ( NBXMAX .GE. I ) THEN
                    IBX = I
                 ELSE
                    IBX = NBXMAX + 1
                    CALL LZIOOP (FHNDL,LCOMM,LRERR,I-NBXMAX,
     &                           NI,BX(1,IBX),'BX','GET')
                    IF ( LRERR .NE. 0 ) RETURN
                 END IF
              END IF
C
C............ retrieve the corresponding (X) ...........................
C
              IF ( NXMAX .GT. 0 ) THEN
                 IX = I
              ELSE
                 IX = 1
                 CALL LZIOOP (FHNDL,LCOMM,LRERR,I,NI,X(1,IX),'X ','GET')
                 IF ( LRERR .NE. 0 ) RETURN
              END IF
C
C............ orthogonalization step ...................................
C
              IF ( GNRZD ) THEN
                 CALL QTBR   (NI,BX(1,IBX),NI,1,R,NI,NVB,ZETA,
     &                        ZMAX,WORK,LBLAS,LCOMM,LRERR)
                 CALL RQALPH (NI,R,NI,NVB,X(1,IX),NI,1,ZETA,LBLAS)
              ELSE
                 CALL QTBR   (NI, X(1,IX) ,NI,1,R,NI,NVB,ZETA,
     &                        ZMAX,WORK,LBLAS,LCOMM,LRERR)
                 CALL RQALPH (NI,R,NI,NVB,X(1,IX),NI,1,ZETA,LBLAS)
              END IF
C
         ELSE
C
C............ this vector is beyond the bounds .........................
C
              TAUR(I) = EPS1
              TAUQ(I) = EPS1
C
         END IF 
C
   10 CONTINUE
C
      RETURN 
C
C**** end of SORTH0 ****************************************************
C
      END
