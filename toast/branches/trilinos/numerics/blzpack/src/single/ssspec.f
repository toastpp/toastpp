      SUBROUTINE SSSPEC (JT,LFILE,LPRNT,N,NEWSIG,NNEIG,NONEWS,NSLS,NSRS,
     &                   NSRLS,NSRRS,ENDL,ENDR,RADIUS,RITZ,RITZL,RITZR,
     &                   RNORM,SFARL,SFARR,SIGMA,SIGMAL,SIGMAR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSSPEC examines the Ritz values spectrum                         *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JT     (sii) : dimension of the basis                            *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    NEWSIG (sib) : flag for a new starting point                     *
C*    NNEIG  (sii) : number of eigenvalues less than SIGMA             *
C*    NONEWS (sii) : number of runs without any converged eingenvalue  *
C*    NSLS   (sii) : number of eigenvalues converged < SIGMA           *
C*    NSRS   (sii) : number of eigenvalues converged > SIGMA           *
C*    NSRLS  (sii) : number of eigenvalues required  < SIGMA           *
C*    NSRRS  (sii) : number of eigenvalues required  > SIGMA           *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    RADIUS (srb) : radius of convergence                             *
C*    RITZ   (ari) : Ritz values                                       *
C*    RITZL  (sro) : farthest Ritz value < SIGMA                       *
C*    RITZR  (sro) : farthest Ritz value > SIGMA                       *
C*    RNORM  (ari) : estimated residuals                               *
C*    SFARL  (sii) : farthest SIGMA to the left  of ORIGIN             *
C*    SFARR  (sii) : farthest SIGMA to the right of ORIGIN             *
C*    SIGMA  (sri) : origin translation                                *
C*    SIGMAL (sro) : origin translation  < SIGMA                       *
C*    SIGMAR (sro) : origin translation  > SIGMA                       *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZPRT7,SIBTST,SSSIGL,SSSIGR                                      *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    ABS,LOG10,MAX,MIN,SQRT                                           *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*     SIGMAL           RITZL               RITZR          SIGMAR      *
C*                                                                     *
C*     left                       sigma                     right      *
C*     sigma  1  2     3 4 5  6  7  |  8  9 10 11     12 13 sigma      *
C*    ---^----.--.-----.-x-x--x--x--|--x--x--x-x------.--.----^--->    *
C*                       |          |          |                       *
C*                       |   RADL   |   RADR   |                       *
C*                 (-)  -+----------+----------+-  (+)                 *
C*                        converged eigenvalues                        *
C*                                                                     *
C*                       RADIUS = max(RADL,RADR)                       *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      REAL             TEN,ZERO
      PARAMETER        (TEN=10.0E0,ZERO=0.0E0)
      REAL             TGAP
      PARAMETER        (TGAP=0.00010E0)
C
C==== arguments ========================================================
C
      INTEGER          JT,LFILE,LPRNT,N,NEWSIG,NNEIG,NONEWS,
     &                 NSLS,NSRS,NSRLS,NSRRS
      REAL             ENDL,ENDR,RADIUS,RITZ(JT),RITZL,RITZR,RNORM(JT),
     &                 SFARL,SFARR,SIGMA,SIGMAL,SIGMAR
C
C==== local variables ==================================================
C
      INTEGER          I,ILFAR,ILNEAR,IRFAR,IRNEAR,NRVL,NRVR
      REAL             RADL,RADR,RCWL,RCWR
      LOGICAL          KEEPLS,KEEPRS
C
C==== subprogram =======================================================
C
      LOGICAL          SIBTST
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        ABS,LOG10,MAX,MIN
C
C**** executable statements ********************************************
C
      NRVL = 0
      NRVR = 0
      RADL = ZERO
      RADR = ZERO
      RCWL = ZERO
      RCWR = ZERO
C
      RITZL  = SIGMA
      RITZR  = SIGMA
      SIGMAL = SIGMA
      SIGMAR = SIGMA
C
      KEEPLS = .FALSE.
      KEEPRS = .FALSE.
C
C.... closest and farthest eigenvalue to the left ......................
C
      ILFAR  = 0
      ILNEAR = 0
C
      DO 10 I = JT,1,-1
	 IF ( RITZ(I) .LT. SIGMA ) THEN
            NRVL = NRVL + 1
            RCWL = RCWL + LOG10(ABS(RITZ(I)))
	    IF ( RNORM(I) .GT. ZERO ) THEN
               IF      ( ILNEAR .EQ. 0 ) THEN
                       ILNEAR = I
               ELSE IF ( NRVL .LE. NSRLS ) THEN
                       ILFAR  = I
               END IF
            END IF
         ELSE
            RADR = ABS(RITZ(I)-SIGMA)
         END IF 
   10 CONTINUE
C
      IF ( ILFAR .EQ. 0 ) ILFAR = ILNEAR
C
C.... closest and farthest eigenvalue to the right .....................
C
      IRFAR  = 0
      IRNEAR = 0
C
      DO 20 I = 1,JT,+1
	 IF ( RITZ(I) .GT. SIGMA ) THEN
            NRVR = NRVR + 1
            RCWR = RCWR + LOG10(ABS(RITZ(I)))
	    IF ( RNORM(I) .GT. ZERO ) THEN
               IF      ( IRNEAR .EQ. 0 ) THEN
                       IRNEAR = I
               ELSE IF ( NRVR .LE. NSRRS ) THEN
                       IRFAR  = I
               END IF
            END IF
         ELSE
            RADL = ABS(RITZ(I)-SIGMA)
         END IF 
   20 CONTINUE
C
      IF ( IRFAR .EQ. 0 ) IRFAR = IRNEAR
C
C.... `convergence radius' .............................................
C
      IF ( ILFAR .GT. 0 ) THEN
         RADL  = ABS(RITZ(ILFAR )-SIGMA)
         RCWL  = RCWL/NRVL
         RITZL = RITZ(ILFAR)
      END IF
C
      IF ( IRFAR .GT. 0 ) THEN
         RADR  = ABS(RITZ(IRFAR )-SIGMA)
         RCWR  = RCWR/NRVR
         RITZR = RITZ(IRFAR)
      END IF
C
      RADIUS = MAX(RADIUS,RADL,RADR)
C  
      RCWL = TEN**RCWL
      RCWR = TEN**RCWR
C  
C.... define a sigma to the left .......................................
C
      IF      ( NNEIG .EQ. 0 ) THEN
C
C............ lower limit ..............................................
C  
	      SIGMAL = SIGMA
              RITZL  = SIGMA
C
      ELSE IF ( ILNEAR .EQ. 0 ) THEN
C  
C............ no converged eigenvalue to the left ......................
C
              SIGMAL = SIGMA - RADIUS*TEN**NONEWS
C
              IF ( ABS(SIGMA-SIGMAL).LE.ABS(SIGMA)*TGAP ) THEN
                 KEEPLS = NEWSIG.LE.2 .AND. NSRLS.GT.0
                 SIGMAL = SIGMA 
              END IF
C
      ELSE IF ( NSRLS .GT. 0 ) THEN
C  
C............ choose a shift based on Ritz values ......................
C
              SIGMAL = SIGMA
C
              CALL SSSIGL (JT,NRVL,NSLS,NSRLS,RITZ,SIGMAL,TGAP)
C
              IF ( SIGMAL .EQ. SIGMA ) THEN
                 IF ( RADIUS .GE. RCWL ) THEN
                    SIGMAL = SIGMA - SQRT(RADIUS*RCWL)*2
                 ELSE
                    SIGMAL = SIGMA - RADIUS*2
                 END IF
              END IF
C
              KEEPLS = NEWSIG.LE.2 .AND.
     &                 SIGMA.GT.SFARL .AND. SIGMAL.LE.SFARL
C
      END IF
C  
C.... define a sigma to the right ......................................
C
      IF      ( NNEIG .EQ. N ) THEN
C
C............ upper limit ..............................................
C  
              RITZR  = SIGMA
              SIGMAR = SIGMA
C
      ELSE IF ( IRNEAR.EQ.0 .AND. RADIUS*TGAP.GE.ABS(SIGMA) ) THEN
C  
C............ no converged eigenvalue to the right .....................
C
              SIGMAR = SIGMA + RADIUS*TEN**NONEWS
C
              IF ( ABS(SIGMA-SIGMAR).LE.ABS(SIGMA)*TGAP ) THEN
                 KEEPRS = NEWSIG.LE.2 .AND. NSRRS.GT.0
                 SIGMAR = SIGMA 
              END IF
C
      ELSE IF ( NSRRS .GT. 0 ) THEN
C  
C............ choose a shift based on Ritz values ......................
C
              SIGMAR = SIGMA
C
              CALL SSSIGR (JT,NRVL,NSRS,NSRRS,RITZ,SIGMAR,TGAP)
C
              IF ( SIGMAR .EQ. SIGMA ) THEN
                 IF ( RADIUS .GE. RCWR ) THEN
                    SIGMAR = SIGMA + SQRT(RADIUS*RCWR)*2
                 ELSE
                    SIGMAR = SIGMA + RADIUS*2
                 END IF
              END IF
C
              KEEPRS = NEWSIG.LE.2. AND.
     &                 SIGMA.LT.SFARR .AND. SIGMAR.GE.SFARR
C
      END IF 
C
C.... check bounds .....................................................
C
      SIGMAL = MAX(SIGMAL,ENDL)
      SIGMAR = MIN(SIGMAR,ENDR)
C
      IF ( KEEPLS .OR. KEEPRS ) NEWSIG = 1
C
      IF ( SIBTST(6,LPRNT) ) CALL LZPRT7 (LFILE,RITZL,RITZR,
     &                                    SIGMAL,SIGMAR)
C
      RETURN 
C
C**** end of SSSPEC ****************************************************
C
      END
