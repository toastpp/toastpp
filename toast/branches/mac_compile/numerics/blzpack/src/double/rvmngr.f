      SUBROUTINE RVMNGR (JL,JT,NVB,LEIG,LNI,NI,NQMAX,NXMAX,FHNDL,LCOMM,
     &                   LRERR,LRWRN,RITZ,S,SIGMA,THETA0,BETA,BASIS,Q,
     &                   NTEIG,EIG,X,WORK,GNRZD,PURFY)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RVMNGR computes Ritz vectors                                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL     (sii) : number of steps                                   *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    NVB    (sii) : number of vectors in a block                      *
C*    LEIG   (sii) : leading dimension of (EIG)                        *
C*    LNI    (sii) : leading dimension of (X)                          *
C*    NI     (sii) : dimension of the vectors in (Q) ane (X)           *
C*    NQMAX  (sii) : maximum number of vectors in (BASIS)              *
C*    NXMAX  (aii) : maximum number of vectors in (X)                  *
C*    FHNDL  (sii) : file handle                                       *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    LRERR  (sio) : code for error messages                           *
C*    LRWRN  (sio) : code for warning messages                         *
C*    RITZ   (ari) : Ritz values and estimated residuals               *
C*    S      (ari) : eigenvectors of the tridiagonal matrix            *
C*    SIGMA  (ari) : origin translation                                *
C*    THETA0 (sri) : reference point for THETA                         *
C*    BETA   (ari) : matrix (BETA) in (R)=(Q)*(BETA) at step JL        *
C*    BASIS  (ari) : Lanczos vectors array                             *
C*    Q      (ari) : Lanczos vectors at step JL+1                      *
C*    NTEIG  (sib) : number of computed eigenpairs                     *
C*    EIG    (arb) : eigenvalue approximations and estimated residuals *
C*    X      (arb) : eigenvector approximations                        *
C*    WORK   (arw) : work array                                        *
C*    GNRZD  (sli) : problem type flag                                 *
C*    PURFY  (sli) : eigenvectors purification flag                    *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZIOOP,RVCOMP,SETLRM,SETTO0                                      *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    DCOPY                                                            *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    MAX,MIN                                                          *
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
      INTEGER          FHNDL(*),JL,JT,LCOMM,LEIG,LNI,LRERR,LRWRN,NI,
     &                 NQMAX,NTEIG,NVB,NXMAX
      DOUBLE PRECISION BASIS(*),BETA(NVB,NVB),EIG(LEIG,2),Q(NI,NVB),
     &                 RITZ(JT,2),SIGMA,S(JT,*),THETA0,WORK(JT,*),
     &                 X(LNI,*)
      LOGICAL          GNRZD,PURFY
C
C==== local variables ==================================================
C
      INTEGER          I,INDEX,J,K,LBLAS,LDB,NCEIG,NEWR0,
     &                 NQ,NRITZ,NRMAX,NX
      DOUBLE PRECISION DELTA,SUM
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        MAX,MIN
C
C**** executable statements ********************************************
C
      NEWR0 = 0
      NRITZ = 0
      IF ( GNRZD ) THEN 
         LDB = NI*2
      ELSE
         LDB = NI
      END IF
      IF ( NXMAX .GT. 0 ) THEN 
         NCEIG = NTEIG
      ELSE
         NCEIG = 0
      END IF
C
C.... maximum number of Ritz vectors that can be computed at once ......
C
C     number of columns of S to be stored in WORK
C
      NRMAX = MIN(JT,MAX(1,NXMAX))
C
C.... loop on the JT Ritz values .......................................
C
      DO 50 INDEX = 1,JT
C
         IF      ( RITZ(INDEX,2) .EQ. ZERO ) THEN
C
C............... computation of restarting vectors .....................
C
                 NCEIG = 0
                 NEWR0 = NEWR0 + 1 
                 CALL DCOPY (JT,S(1,INDEX),1,WORK(1,NEWR0),1)
C
         ELSE IF ( RITZ(INDEX,2) .GT. ZERO ) THEN
C
C............... computation of Ritz vectors ...........................
C
   	         IF ( NTEIG .EQ. LEIG ) THEN
                    CALL SETLRM (9,LRWRN)
                    RETURN
   	         ELSE 
                    NRITZ = NRITZ + 1
                    NTEIG = NTEIG + 1
                    EIG(NTEIG,1) = RITZ(INDEX,1)
                    EIG(NTEIG,2) = RITZ(INDEX,2)
                    CALL DCOPY (JT,S(1,INDEX),1,WORK(1,NRITZ),1)
         	 END IF
C
         END IF 
C
C....... compute the Ritz vectors : (X)=(Q)*(s) ........................
C
         IF ( ( (NRITZ.NE.0  ) .AND. (NRITZ.EQ.NRMAX) ) .OR.
     &        ( (NRITZ.NE.0  ) .AND. (NTEIG.EQ.LEIG ) ) .OR.
     &        ( (NRITZ.NE.0  ) .AND. (INDEX.EQ.JT   ) ) .OR.
     &        ( (NEWR0.NE.0  ) .AND. (INDEX.EQ.JT   ) ) .OR.  
     &          (NEWR0.EQ.NVB) ) THEN
C
            CALL SETTO0 (LNI*MAX(NEWR0,NRITZ),X(1,NCEIG+1),1)
C
            NQ = MIN(JL,NQMAX)*NVB
            NX = MAX(NRITZ,NEWR0)
            LBLAS = MIN(3,NX)
C
C.......... sum up the the j-th (Q) stored in BASIS ....................
C
            IF ( NQ .GT. 0 ) CALL RVCOMP (JT,LBLAS,LDB,LNI,NI,NQ,NX,
     &                                    BASIS,WORK,X(1,NCEIG+1))
C
C.......... sum up the the j-th (Q) retrieved from secondary storage ...
C
            IF ( JL .GT. NQMAX) THEN
               J = 1 + NVB*NQMAX
               K = 1 + NI*NVB*NQMAX
               DO 10 I = 1,JL-NQMAX
                  CALL LZIOOP (FHNDL,LCOMM,LRERR,I,NI*NVB,
     &                         BASIS(K),'Q ','GET')
                  IF ( LRERR .NE. 0 ) RETURN
                  CALL RVCOMP (JT,LBLAS,NI,LNI,NI,NVB,NX,BASIS(K),
     &                         WORK(J,1),X(1,NCEIG+1))
                  J = J + NVB
   10          CONTINUE
            END IF
C
C.......... if required refine the Ritz vectors ........................
C
            IF ( PURFY ) THEN
               DO 40 I = 1,NX 
                  IF ( THETA0 .EQ. ZERO ) THEN
                     DELTA = RITZ(I,1) - SIGMA
                  ELSE
                     DELTA = (RITZ(I,1)-SIGMA)/SIGMA
                  END IF
                  DO 30 J = 1,NVB
                     SUM = ZERO
                     DO 20 K = J,NVB
                        SUM = SUM + BETA(J,K)*WORK(JT-NVB+K,I)
   20                CONTINUE
                     WORK(J,I) = SUM*DELTA
   30             CONTINUE
   40          CONTINUE
               CALL RVCOMP (JT,LBLAS,LDB,LNI,NI,NVB,NX,
     &                      Q,WORK,X(1,NCEIG+1))
            END IF
C
C.......... if required store (X) in a file ............................
C
            IF ( NXMAX.EQ.0 .AND. NRITZ.GT.0 ) THEN
               CALL LZIOOP (FHNDL,LCOMM,LRERR,NTEIG,NI,
     &                      X(1,NCEIG+1),'X ','PUT')
               IF ( LRERR .NE. 0 ) RETURN
            END IF
C
            NCEIG = MIN(NTEIG,NXMAX)
            NEWR0 = 0           
            NRITZ = 0           
C
         END IF
C
   50 CONTINUE
C
      RETURN 
C
C**** end of RVMNGR ****************************************************
C
      END
