      SUBROUTINE LZSTP5 (JL,JT,JTMAX,JTMIN,NVB,N,LCOMM,NPE,LFILE,LPRNT,
     &                   LRMDE,LRWRN,NEWSIG,NDEIG,NSLS,NSRS,NSRLS,NSRRS,
     &                   NRITZ,NULLDQ,NULLDR,TIME,RITZ,BETAR,TB,S,THETA,
     &                   IWORK,RWORK,BIGNUM,EPS,REPS,SIGMA,THETA0,
     &                   THETAL,THETAR,THRSH,IDXETA,IDXTAU,ABSETA,
     &                   ABSTAU,EIGON,GNRZD,SLICE)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZSTP5 checks convergence                                        *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL     (sii) : number of steps                                   *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    JTMAX  (sii) : maximum dimension of the block tridiagonal matrix *
C*    JTMIN  (sii) : minimum dimension of the block tridiagonal matrix *
C*    NVB    (sii) : number of vectors in a block                      *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    NPE    (sii) : number of processes                               *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    LRMDE  (sii) : run mode                                          *
C*    LRWRN  (sii) : code for warning messages                         *
C*    NEWSIG (sio) : flag for a new starting point                     *
C*    NDEIG  (sii) : number of eigenpairs required in the run          *
C*    NSLS   (sio) : number of negative eigenvalues converged          *
C*    NSRS   (sio) : number of positive eigenvalues converged          *
C*    NSRLS  (sii) : number of eigenvalues required < SIGMA            *
C*    NSRRS  (sii) : number of eigenvalues required > SIGMA            *
C*    NRITZ  (sii) : number of eigenvalues of (TB) computed            *
C*    NULLDQ (sii) : number of zero diagonal entries in BETAQ          *
C*    NULLDR (sii) : number of zero diagonal entries in BETAR          *
C*    TIME   (arb) : time table                                        *
C*    RITZ   (aro) : Ritz values and estimated residuals               *
C*    BETAR  (ari) : matrix (BETA) in (R)*(Q)*(BETA) at current  step  *
C*    TB     (ari) : block tridiagonal matrix                          *
C*    S      (aro) : eigenvectors of the tridiagonal matrix            *
C*    THETA  (aro) : eigenvalues of the tridiagonal matrix             *
C*    IWORK  (arw) : integer workspace                                 *
C*    RWORK  (arw) : real workspace                                    *
C*    BIGNUM (sri) : big number                                        *
C*    EPS    (sri) : roundoff unit                                     *
C*    REPS   (sri) : sqrt(EPS)                                         *
C*    SIGMA  (sri) : origin translation                                *
C*    THETA0 (sri) : reference point for THETA                         *
C*    THETAL (sri) : inferior limit to converged eigenvalues           *
C*    THETAR (sri) : superior limit to converged eigenvalues           *
C*    THRSH  (sri) : threshold for convergence                         *
C*    IDXETA (sii) : index of the maximum entry in (ETA)               *
C*    IDXTAU (sii) : index of the maximum entry in (TAU)               *
C*    ABSETA (sri) : maximum entry in (ETA)                            *
C*    ABSTAU (sri) : maximum entry in (TAU)                            *
C*    EIGON  (slo) : eigenvectors computation flag                     *
C*    GNRZD  (sli) : problem type flag                                 *
C*    SLICE  (sli) : spectrum slicing flag                             *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZCHEK,LZPRT3,SITIME,SSTRSF,TBEIGP                               *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    DCOPY                                                            *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          IDXETA,IDXTAU,IWORK(*),JL,JT,JTMAX,JTMIN,
     &                 LCOMM,LFILE,LPRNT,LRMDE,LRWRN,N,NDEIG,
     &                 NEWSIG,NPE,NRITZ,NSLS,NSRS,NSRLS,
     &                 NSRRS,NULLDQ,NULLDR,NVB
      DOUBLE PRECISION ABSETA,ABSTAU,BETAR(*),BIGNUM,EPS,GRNRM,REPS,
     &                 RITZ(*),RWORK(*),S(*),SIGMA,TB(*),THETA(*),
     &                 THETA0,THETAL,THETAR,THRSH,TIME(*)
      LOGICAL          EIGON,GNRZD,SLICE
C
C==== local variables ==================================================
C
      DOUBLE PRECISION THRES,TIME0,TEIGTB
      LOGICAL          OUTER
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
C
C**** executable statements ********************************************
C
      NRITZ = 0
C
      IF ( JL.GT.0 .AND. ( IDXETA.GT.0 .OR. JT.EQ.JTMAX .OR.
     &                     NULLDR.GT.0 .OR. LPRNT.GT.10000 ) ) THEN
C
         NRITZ = JT
C
C....... solve the eigensystem (TB)*(S)=(S)*(THETA) ....................
C
         TIME0 = SITIME(ZERO)
C
         CALL TBEIGP (JT,JTMAX,NVB,LRWRN,NSLS,NSRS,NSRLS,NSRRS,NULLDR,
     &                BIGNUM,EPS,REPS,BETAR,TB,S,THETA,THETA0,THETAL,
     &                THETAR,THRES,THRSH,RWORK,IWORK,GNRZD,
     &                SLICE,OUTER,.FALSE.)
C
         CALL DCOPY  (JT*2,THETA,1,RITZ,1)
C
         TEIGTB = SITIME(TIME0)
C
         TIME( 6) = TIME( 6) + TEIGTB
         TIME(10) = TIME(10) + TEIGTB
C
C....... transform the spectrum ........................................
C
         IF ( GNRZD ) CALL SSTRSF (JT,0,RITZ,S,SIGMA,THETA0)
C
C....... check the numbers obtained so far .............................
C
         CALL LZCHEK (JT,JTMAX,JTMIN,NVB,N,LCOMM,NPE,LRMDE,LRWRN,NEWSIG,
     &                NDEIG,NSLS,NSRS,NSRLS,NSRRS,NULLDQ,NULLDR,BIGNUM,
     &                EPS,GRNRM,TIME,THETA,THETAL,THETAR,THRES,
     &                EIGON,GNRZD,OUTER,SLICE)
C
C....... print information .............................................
C
         CALL LZPRT3 (JL,JT,LFILE,LPRNT,IDXETA,IDXTAU,
     &                ABSETA,ABSTAU,THRES,RITZ)
C
      END IF

      RETURN 
C
C**** end of LZSTP5 ****************************************************
C
      END
