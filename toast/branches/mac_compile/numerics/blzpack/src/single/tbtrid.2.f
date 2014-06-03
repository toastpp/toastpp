      SUBROUTINE TBTRID (JT,JTMAX,NVB,LRWRN,TB,S,THETA,RWORK,IWORK,FULL)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    TBTRID solves the reduced eigenproblem (TB)*(S)=(S)*(THETA)      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    JTMAX (sii) : maximum dimension of the block tridiagonal matrix  *
C*    NVB   (sii) : number of vectors in a block                       *
C*    LRWRN (sio) : code for warning messages                          *
C*    TB    (ari) : block tridiagonal matrix                           *
C*    S     (aro) : eigenvectors of the block tridiagonal matrix       *
C*    THETA (aro) : eigenvalues  of the block tridiagonal matrix       *
C*    RWORK (arw) : real workspace                                     *
C*    IWORK (aiw) : integer workspace                                  *
C*    FULL  (sli) : eigenvector computation flag                       *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SETLRM,TBCOPY                                                    *
C*                                                                     *
C*  - LAPACK:                                                          *
C*                                                                     *
C*    SSBTRD,SSTEGR                                                    *
C*                                                                     *
C*  - BLAS kernels:                                                    *
C*                                                                     *
C*    SCOPY,SGEMV                                                      *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MIN                                                              *
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
      INTEGER          IWORK(*),JT,JTMAX,LRWRN,NFAIL,NVB
      REAL             RWORK(*),S(JT,JT),TB(JTMAX,NVB+1),THETA(JT)  
      LOGICAL          FULL
C
C==== local variables ==================================================
C
      INTEGER          HLFBND,INDD,INDE,INDIS,INDIW,INDQ,INDRW,
     &                 INDTB,INFO,J,LIWORK,LRWORK,M,N
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MIN
C
C**** executable statements ********************************************
C
      IF ( FULL ) CONTINUE
C
      N = JT
      NFAIL = 0
      RWORK(1) = 0
      IWORK(1) = 0
C
   10 CONTINUE
C
      LIWORK = 10*N
      LRWORK = 18*N
      HLFBND = MIN(NVB+1,N)
C
C.... set pointers .....................................................
C
C     real workspace: (2+N+MAX(HLFBND+1,18))*N 
C
      INDD  = 1
      INDE  = INDD + N
      INDQ  = INDE + N
      INDTB = INDQ + N*N
      INDRW = INDTB + HLFBND*N
C
C     integer workspace: 12*N 
C
      INDIS = 1
      INDIW = INDIS + 2*N
C
C.... reduction to tridiagonal form if necessary .......................
C
      IF ( NVB .GT. 1 ) THEN
         CALL TBCOPY (HLFBND,N,JTMAX,TB,RWORK(INDTB))
         CALL SSBTRD ('V','L',N,HLFBND-1,RWORK(INDTB),HLFBND,
     &                RWORK(INDD),RWORK(INDE),RWORK(INDQ),
     &                N,RWORK(INDRW),INFO)
      ELSE
         CALL SCOPY  (N  ,TB(1,1),1,RWORK(INDD),1)
         CALL SCOPY  (N-1,TB(1,2),1,RWORK(INDE),1)
      END IF
C
C.... solve the tridiagonal eigenproblem ...............................
C
      INDRW = INDTB
C
      CALL SSTEGR ('V','A',N,RWORK(INDD),RWORK(INDE),ZERO,ZERO,0,0,
     &             ZERO,M,THETA,S,N,IWORK(INDIS),RWORK(INDRW),
     &             LRWORK,IWORK(INDIW),LIWORK,INFO)
C
C.... apply orthogonal matrix used in reduction to tridiagonal form ....
C
      IF ( NVB .GT. 1 ) THEN
         DO 20 J = 1,M
            CALL SCOPY (N,S(1,J),1,RWORK(1),1)
            CALL SGEMV ('N',N,N,ONE,RWORK(INDQ),N,RWORK(1),1,
     &                  ZERO,S(1,J),1)
   20    CONTINUE
      END IF
C
C.... test the exit flag ...............................................
C
      IF ( INFO .GT. 0 ) THEN
         N = N - NVB
         NFAIL = NFAIL + 1
         CALL SETLRM (5,LRWRN)
         IF ( NFAIL .LE. 2 ) GO TO 10
      END IF
C
      RETURN 
C
C**** end of TBTRID ****************************************************
C
      END
