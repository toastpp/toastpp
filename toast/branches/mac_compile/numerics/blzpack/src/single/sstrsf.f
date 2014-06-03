      SUBROUTINE SSTRSF (NVAL,NVEC,RITZ,S,SIGMA,THETA0)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSTRSF transforms the eigenvalue spectrum from [THETA] to [RITZ] *
C*                                                                     *
C*           THETA0.eq.0: RITZ = SIGMA + 1/RITZ                        *
C*           THETA0.ne.0: RITZ = SIGMA*RITZ/(RITZ-THETA0)              *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    NVAL   (sii) : number of eigenvalues                             *
C*    NVEC   (sii) : number of eigenvectors                            *
C*    RITZ   (arb) : Ritz values and estimated residuals               *
C*    S      (arb) : eigenvectors of the tridiagonal matrix            *
C*    SIGMA  (sri) : origin translation                                *
C*    THETA0 (sri) : reference point for THETA                         *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    EIGSRT                                                           *
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
      INTEGER          NVAL,NVEC
      REAL             RITZ(NVAL,2),SIGMA,S(NVEC,NVEC),THETA0
C
C==== local variables ==================================================
C
      INTEGER          I
C
C**** executable statements ********************************************
C
C.... compute the reciprocals ..........................................
C
      IF ( THETA0 .EQ. ZERO ) THEN
         DO 10 I = 1,NVAL
            IF ( RITZ(I,1).NE.ZERO ) RITZ(I,1) = 
     &                               SIGMA+ONE/RITZ(I,1)
   10    CONTINUE
      ELSE
         DO 20 I = 1,NVAL
            IF ( RITZ(I,1).NE.ONE  ) RITZ(I,1) =
     &                               SIGMA*RITZ(I,1)/(RITZ(I,1)-THETA0)
   20    CONTINUE
      END IF
C
C.... sort the eigenpairs in ascending order ...........................
C
      CALL EIGSRT (NVAL,NVEC,NVAL,RITZ,S)
C
      RETURN 
C
C**** end of SSTRSF ****************************************************
C
      END
