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
C*    IDENTY,SETLRM,TBBRED,TBCOPY,TBIMQL                               *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    DCOPY                                                            *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MIN                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          IWORK(*),JT,JTMAX,LRWRN,NVB
      DOUBLE PRECISION RWORK(*),S(JT,JT),TB(JTMAX,NVB+1),THETA(JT)  
      LOGICAL          FULL
C
C==== local variables ==================================================
C
      INTEGER          HLFBND,I0,INFO,N,NFAIL
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MIN
C
C**** executable statements ********************************************
C
      N = JT
      NFAIL = 0
      RWORK(1) = 0
      IWORK(1) = 0
C
   10 CONTINUE
C
C.... set variables ....................................................
C
      I0 = 1
C
      IF ( .NOT.FULL ) I0 = I0 + N - NVB
C
      HLFBND = MIN(NVB+1,N)
C
C.... initialization of (S) ............................................
C
      CALL IDENTY (N,S)
C
C.... reduction to tridiagonal form if necessary .......................
C 
      IF ( NVB .GT. 1 ) THEN
         CALL TBCOPY (HLFBND,N,JTMAX,TB,RWORK(1+N))
         CALL TBBRED (HLFBND,N,RWORK(1+N),THETA,RWORK,S,I0)
      ELSE
         CALL DCOPY  (N  ,TB(1,1),1,THETA(1),1)
         CALL DCOPY  (N-1,TB(1,2),1,RWORK(2),1)
      END IF
C
C.... solve the tridiagonal eigenproblem ...............................
C
      CALL TBIMQL (N,THETA,RWORK(1),S,I0,INFO)
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
