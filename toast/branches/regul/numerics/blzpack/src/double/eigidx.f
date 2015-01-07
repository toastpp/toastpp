      SUBROUTINE EIGIDX (INDXP,INDXS,NNEIG,NPEIG,NREIG,NTEIG,
     &                   BIGNUM,EIG,EIGL,EIGR,ORIGIN)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    EIGIDX finds the closest or farthest converged eigenvalue        *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    INDXP  (sio) : index for the eigenpairs (printing)               *
C*    INDXS  (sio) : index for the eigenpairs (algebraic)              *
C*    NNEIG  (sii) : number of eigenvalues less than ORIGIN            *
C*    NPEIG  (sii) : number of eigenpairs to be printed                *
C*    NREIG  (sii) : number of required eigenpairs                     *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    BIGNUM (sri) : big number                                        *
C*    EIG    (ari) : eigenvalues                                       *
C*    EIGL   (sri) : inferior bound for eigenvalues                    *
C*    EIGR   (sri) : superior bound for eigenvalues                    *
C*    ORIGIN (sri) : starting point (first SIGMA)                      *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    ABS                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          INDXP,INDXS,NNEIG,NPEIG,NREIG,NTEIG
      DOUBLE PRECISION BIGNUM,EIG(*),EIGL,EIGR,ORIGIN
C
C==== local variables ==================================================
C
      INTEGER          I,J,K
      DOUBLE PRECISION DELTA,BOUND 
C
C==== intrinsic function ===============================================
C
      INTRINSIC        ABS 
C
C**** executable statements ********************************************
C
      NPEIG = 0
      INDXP = NTEIG
      INDXS = NNEIG + 1
C
C.... determine INDXP and INDXS ........................................
C
      DO 20 I = 1,NREIG 
         K = 0
         BOUND = BIGNUM
         DO 10 J = 1,NTEIG 
            DELTA = ABS(ORIGIN-EIG(J))
	    IF ( ( DELTA  .LE. BOUND ) .AND.
     &           ( EIG(J) .GT. EIGL  ) .AND.
     &           ( EIG(J) .LT. EIGR  ) ) THEN
                 BOUND = DELTA    
                 K = J
            END IF
   10    CONTINUE
         IF ( K .NE. 0 ) THEN
            IF ( EIG(K) .LT. ORIGIN ) INDXS = INDXS - 1
            IF ( K .LT. INDXP ) INDXP = K
            NPEIG  = NPEIG + 1
            EIG(K) = BIGNUM
         END IF
   20 CONTINUE
C
      RETURN 
C
C**** end of EIGIDX ****************************************************
C
      END
