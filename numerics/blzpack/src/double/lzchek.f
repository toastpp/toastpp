      SUBROUTINE LZCHEK (JT,JTMAX,JTMIN,NVB,N,LCOMM,NPE,LRMDE,LRWRN,
     &                   NEWSIG,NDEIG,NSLS,NSRS,NSRLS,NSRRS,NULLDQ,
     &                   NULLDR,BIGNUM,EPS,GRNRM,TIME,THETA,THETAL,
     &                   THETAR,THRES,EIGON,GNRZD,OUTER,SLICE)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZCHEK checks the status of the run                              *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    JTMAX  (sii) : maximum dimension of the block tridiagonal matrix *
C*    JTMIN  (sii) : minimum dimension of the block tridiagonal matrix *
C*    NVB    (sii) : number of vectors in a block                      *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    NPE    (sii) : number of processes                               *
C*    LRMDE  (sii) : run mode                                          *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NEWSIG (sio) : flag for a new starting point                     *
C*    NDEIG  (sii) : number of eigenvalues required in the run         *
C*    NSLS   (sii) : number of eigenvalues converged < than SIGMA      *
C*    NSRS   (sii) : number of eigenvalues converged > than SIGMA      *
C*    NSRLS  (sii) : number of eigenvalues required  < than SIGMA      *
C*    NSRRS  (sii) : number of eigenvalues required  > than SIGMA      *
C*    NULLDQ (sii) : number of zero diagonal entries in BETAQ          *
C*    NULLDR (sii) : number of zero diagonal entries in BETAR          *
C*    BIGNUM (sri) : big number                                        *
C*    EPS    (ari) : time table                                        *
C*    GRNRM  (srb) : global residual of unconverged eigenpairs         *
C*    TIME   (ari) : time table                                        *
C*    THETA  (ari) : eigenvalues of the block tridiagonal matrix       *
C*    THETAL (sri) : inferior limit to converged eigenvalues           *
C*    THETAR (sri) : superior limit to converged eigenvalues           *
C*    THRES  (sri) : threshold for convergence                         *
C*    EIGON  (slo) : eigenvectors computation flag                     *
C*    GNRZD  (sli) : problem type flag                                 *
C*    OUTER  (sli) : convergence flag                                  *
C*    SLICE  (sli) : spectrum slicing flag                             *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SETLRM,SSNEED,TBILLC                                             *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MAX                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          JT,JTMAX,JTMIN,LCOMM,LRMDE,LRWRN,N,NEWSIG,NDEIG,
     &                 NPE,NSLS,NSRS,NSRLS,NSRRS,NULLDQ,NULLDR,NVB
      DOUBLE PRECISION BIGNUM,EPS,GRNRM,THETA(JT,2),
     &                 THETAL,THETAR,THRES,TIME(*)
      LOGICAL          EIGON,GNRZD,OUTER,SLICE
C
C==== subprogram =======================================================
C
      LOGICAL          TBILLC
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MAX
C
C**** executable statements ********************************************
C
      EIGON = .FALSE.
C
      IF      ( JT .EQ. JTMAX ) THEN 
C
C............ activate the computation of eigenvectors .................
C
              NEWSIG = 2
              EIGON = .TRUE.
C
      ELSE IF ( NULLDR .NE. 0 ) THEN
C
C............ deal with an ill conditioned or invariant subspace .......
C
              IF ( (NULLDQ.GT.0) .OR.
     &             (NULLDR.EQ.NVB) .OR.
     &             ((JT+NVB).LT.N) ) THEN
                   CALL SETLRM (4,LRWRN)
                   EIGON = .TRUE.
              END IF
C
      ELSE IF ( SLICE ) THEN
C
C............ check for a new shift ....................................
C
              CALL SSNEED (JT,JTMIN,NVB,LCOMM,NPE,LRMDE,NEWSIG,NSLS,
     &                     NSRS,NSRLS,NSRRS,BIGNUM,EPS,GRNRM,THETA(1,2),
     &                     THETA(1,1),THETAL,THETAR,THRES,TIME)
C
              IF ( NEWSIG.EQ.3 .AND. LRMDE.GT.6 ) CALL SETLRM (2,LRWRN)
C
              EIGON = ( OUTER .AND. JT.GE.JTMIN ) .OR. NEWSIG.GT.1
C
      ELSE IF ( GNRZD ) THEN
C
              IF ( TBILLC(JT,EPS,THETA) ) THEN
C
C............... operator is assumed to be ill conditioned .............
C
                 CALL SETLRM (2,LRWRN)
                 EIGON = .TRUE.
C
              ELSE
C
C............... check convergence .....................................
C
                 EIGON = ( OUTER .AND. JT.GE.JTMIN ) .OR. 
     &                   ( NSLS.GE.NSRLS .AND. NSRS.GE.NSRRS )
C
              END IF
C
      ELSE
C
C............ standard problem .........................................
C
              EIGON = NDEIG.LE.(NSLS+NSRS)
C
      END IF
C
C.... check for a small problem and an empty interval ..................
C
      IF ( JT.GE.N .AND. (NSLS+NSRS).EQ.0 ) CALL SETLRM (1,LRWRN)
C
      RETURN
C
C**** end of LZCHEK ****************************************************
C
      END
