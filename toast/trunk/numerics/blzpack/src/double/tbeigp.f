      SUBROUTINE TBEIGP (JT,JTMAX,NVB,LRWRN,NSLS,NSRS,NSRLS,NSRRS,NULLD,
     &                   BIGNUM,EPS,REPS,BETA,TB,S,THETA,THETA0,THETAL,
     &                   THETAR,THRES,THRSH,RWORK,IWORK,
     &                   GNRZD,SLICE,OUTER,FULL)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBEIGP solves the reduced eigenproblem (TB)*(S)=(S)*(THETA)      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    JTMAX  (sii) : maximum dimension of the block tridiagonal matrix *
C*    NVB    (sii) : number of vectors in a block                      *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NSLS   (sio) : number of negative eigenvalues converged          *
C*    NSRS   (sio) : number of positive eigenvalues converged          *
C*    NSRLS  (sii) : number of negative eigenvalues required           *
C*    NSRRS  (sii) : number of positive eigenvalues required           *
C*    NULLD  (sii) : number of eigenvalues equal to zero               *
C*    BIGNUM (sri) : big number                                        *
C*    EPS    (sri) : roundoff unit                                     *
C*    REPS   (sri) : sqrt(EPS)                                         *
C*    BETA   (ari) : matrix (BETA) in (R)=(Q)*(BETA) at step JL        *
C*    TB     (ari) : block tridiagonal matrix                          *
C*    S      (aro) : eigenvectors of the block tridiagonal matrix      *
C*    THETA  (aro) : THETA(:,1) returns the eigenvalues of (TB)        *
C*                   THETA(:,2) returns the estimated residuals        *
C*    THETA0 (sri) : reference point for THETA                         *
C*    THETAL (sri) : inferior limit to converged eigenvalues           *
C*    THETAR (sri) : superior limit to converged eigenvalues           *
C*    THRES  (sro) : threshold for convergence                         *
C*    THRSH  (sri) : threshold for convergence (default)               *
C*    RWORK  (arw) : work array                                        *
C*    IWORK  (aiw) : work array                                        *
C*    GNRZD  (sli) : problem type flag                                 *
C*    SLICE  (sli) : spectrum slicing flag                             *
C*    OUTER  (slo) : convergence flag                                  *
C*    FULL   (sli) : eigenvalue computation flag                       *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    TBCONV,TBTRID                                                    *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MIN                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          IWORK(*),JT,JTMAX,LRWRN,NSLS,NSRS,
     &                 NSRLS,NSRRS,NULLD,NVB
      DOUBLE PRECISION BETA(NVB,NVB),BIGNUM,EPS,REPS,RWORK(*),
     &                 S(JT,JT),TB(JTMAX,NVB+1),THETA(JT,2),
     &                 THETA0,THETAL,THETAR,THRES,THRSH
      LOGICAL          FULL,GNRZD,OUTER,SLICE
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MIN
C
C**** executable statements ********************************************
C
      RWORK(1) = 0
      IWORK(1) = 0
C
      IF ( JT .EQ. 0 ) RETURN
C
C.... solution of the reduced eigenproblem .............................
C
      CALL TBTRID (JT,JTMAX,NVB,LRWRN,TB,S,THETA,RWORK,IWORK,FULL)
C
C.... look for converged eigenvalues ...................................
C
      CALL TBCONV (JT,NVB,NSLS,NSRS,NSRLS,NSRRS,NULLD,BIGNUM,EPS,REPS,
     &             BETA,THETA(1,2),S,THETA(1,1),THETA0,THETAL,THETAR,
     &             THRES,THRSH,GNRZD,SLICE,OUTER)
C
      RETURN 
C
C**** end of TBEIGP ****************************************************
C
      END
