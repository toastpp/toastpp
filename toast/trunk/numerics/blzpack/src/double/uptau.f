      SUBROUTINE UPTAU (LEIG,NTEIG,NVB,ALPHA,BMAXN,BMINN,DELTA,EIG,
     &                  ENDL,ENDR,EPS,EPS1,SIGMA,TAUQ,TAUR,
     &                  THETA0,WORK,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    UPTAU updates bounds for the selective orthog. strategy          *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LEIG   (sii) : leading dimension of (EIG)                        *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    NVB    (sii) : number of vectors in a block                      *
C*    ALPHA  (ari) : (Q')*(B)*(R) at current step                      *
C*    BMAXN  (sri) : maximum singular value of current (BETA)          *
C*    BMINN  (sri) : minimum singular value of current (BETA)          *
C*    DELTA  (arw) : work array                                        *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    EPS    (sri) : roundoff unit                                     *
C*    EPS1   (sri) : EPS*NVB*sqrt(N)                                   *
C*    SIGMA  (sri) : origin translation                                *
C*    TAUQ   (arb) : orthog. bounds among (Q) and Ritz vectors         *
C*    TAUR   (arb) : orthog. bounds among (R) and Ritz vectors         *
C*    THETA0 (sri) : reference point for THETA                         *
C*    WORK   (arw) : work array                                        *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    NORM2A                                                           *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*    Recurrence implemented:                                          *
C*                                                                     *
C*    tau_(j+1,i) = { ||I/[eig_(i)-sigma]-alpha_(j)||*tau_(j,i) +      *
C*                    ||beta_(j+1)||*tau_(j-1,i) +                     *
C*                    ||error_(i)|| } / ||beta_(j+1)||                 *
C*                                                                     *
C*    For details, see:                                                *
C*                                                                     *
C*  . `The Lanczos Algorithm with Selective Orthogonalization', B. N.  *
C*    Parlett and D. S. Scott, Mathematics of Computation, 33:217-238, *
C*    1984.                                                            *
C*  . `A Shifted Block Lanczos Algorithm for Solving Sparse Symmetric  *
C*    Eigenvalue Problems', R. G. Grimes, J. G. Lewis and H. D. Simon, *
C*    SIAM J. Matrix Anal. Appl., 15:228-272, 1994.                    *
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
      INTEGER          LEIG,NTEIG,NVB
      DOUBLE PRECISION ALPHA(NVB,NVB),BMAXN,BMINN,DELTA(NVB,NVB),
     &                 EIG(LEIG,2),ENDL,ENDR,EPS,EPS1,SIGMA,
     &                 TAUQ(*),TAUR(*),THETA0,WORK(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C
      INTEGER          I,J
      DOUBLE PRECISION SVMAX,SVMIN,TEMPR,THETA
C
C**** executable statements ********************************************
C
C.... update bounds for existing Ritz pairs ............................
C
      DO 10 I = 1,NTEIG
C
         IF      ( EIG(I,1).LT.ENDL .OR. EIG(I,1).GT.ENDR ) THEN
C
C............... eigenvalue out of bounds ..............................
C
                 TAUR(I) = EPS1 
                 TAUQ(I) = EPS1
C
         ELSE IF ( TAUR(I) .GT. ZERO ) THEN        
C
C............... compute the norm of [theta_(i)*(I)-(alpha_(j))] .......
C
                 IF ( GNRZD .AND. EIG(I,1).NE.SIGMA ) THEN
                    IF ( THETA0 .EQ. ZERO ) THEN
                       THETA = ONE/(EIG(I,1)-SIGMA)
                    ELSE
                       THETA = EIG(I,1)/(EIG(I,1)-SIGMA)
                    END IF
                 ELSE
                    THETA = EIG(I,1)
                 END IF
C
                 CALL DCOPY (NVB*NVB,ALPHA,1,DELTA,1)
C
                 DO 20 J = 1,NVB
                    DELTA(J,J) = DELTA(J,J) - THETA
   20            CONTINUE
C
                 CALL NORM2A (NVB,DELTA,WORK,EPS,SVMAX,SVMIN)
C
C............... update TAUR and TAUQ ..................................
C
                 TEMPR   = TAUR(I)
                 TAUR(I) = SVMAX*TAUR(I)+BMAXN*TAUQ(I)+EIG(I,2)+EPS1
                 TAUR(I) = TAUR(I)/BMINN
                 TAUQ(I) = TEMPR
C
         END IF 
C
   10 CONTINUE
C
      RETURN 
C
C**** end of UPTAU *****************************************************
C
      END
