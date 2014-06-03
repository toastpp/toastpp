      SUBROUTINE RFACTR (LBLAS,LCOMM,LRERR,LRWRN,NI,NVB,NULLP,
     &                   BETA,BR,R,DELTA,ENORM,WORK,AMAXN,
     &                   BMAXN,BMINN,EPS,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RFACTR performs the factorization (R)=(Q)*(BETA)                 *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LBLAS (sii) : BLAS level setting                                 *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    LRWRN (sio) : code for warning messages                          *
C*    NI    (sii) : dimension of the vectors in (Q) and (R)            *
C*    NVB   (sii) : number of vectors in a block                       *
C*    NULLP (sio) : number of null pivots in (BETA)                    *
C*    BETA  (arb) : matrix (BETA) in (R)=(Q)*(BETA)                    *
C*    BR    (arb) : (B)*(R)                                            *
C*    R     (arb) : residual vectors transformed to Lanczos vectors    *
C*    DELTA (arw) : matrix (BETA) in each iteration                    *
C*    ENORM (arw) : diagonals of (DELTA)                               *
C*    WORK  (arw) : workspace                                          *
C*    AMAXN (sri) : minimum singular value of (ALPHA)                  *
C*    BMAXN (sro) : maximum singular value of (BETA)                   *
C*    BMINN (sro) : minimum singular value of (BETA)                   *
C*    EPS   (sri) : roundoff unit                                      *
C*    GNRZD (sli) : problem type flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    MGSCHM,NORM2A,SETLRM,SETTO0,UPBETA                               *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C 
      REAL             THRQT,HNDRD,ONE,ZERO
      PARAMETER        (THRQT=0.750E0,HNDRD=100.0E0,ONE=1.0E0,
     &                  ZERO=0.0E0)
C
C     GSTOL as suggested by Daniel, Cragg, Kaufman and Stewart
C
      REAL             GSTOL
      PARAMETER        (GSTOL=0.70710680E0)
C
C==== arguments ========================================================
C
      INTEGER          LBLAS,LCOMM,LRERR,LRWRN,NI,NULLP,NVB
      REAL             AMAXN,BETA(NVB,NVB),BMAXN,BMINN,BR(NI,NVB),
     &                 DELTA(NVB,NVB),ENORM(NVB),EPS,
     &                 R(NI,NVB),WORK(*)
      LOGICAL          GNRZD
C
C==== local variables ==================================================
C 
      INTEGER          I,J,K
      REAL             BETMAX,BETMIN
      LOGICAL          LEAVE,NPIVT,RANKD
C
C**** executable statements ********************************************
C
      NULLP = 0
      RANKD = .FALSE.
      LEAVE = .FALSE.
C
C.... factorization repeated up to NVB*2 times, if necessary ...........
C
      DO 50 I = 1,NVB*2
C
C....... modified Gram-Schimdt factorization ...........................
C
         CALL MGSCHM (LBLAS,LCOMM,LRERR,NI,NVB,BR,R,
     &                DELTA,WORK,GNRZD,NPIVT)
C
C....... updating: (BETA) := (DELTA)*(BETA) ............................
C
         CALL UPBETA (NVB,BETA,DELTA)
C
C....... check for negative and small pivots ...........................
C
         IF      ( NPIVT ) THEN
                 CALL SETLRM (15,LRWRN)
                 RETURN
         ELSE IF ( I .EQ. 1 ) THEN
                 DO 10 J = 1,NVB
                    ENORM(J) = DELTA(J,J)
   10            CONTINUE
         ELSE
                 BETMIN = ONE
                 BETMAX = ZERO
                 LEAVE = .TRUE.
                 DO 20 J = 1,NVB
	            IF ( DELTA(J,J).LT.ENORM(J)*GSTOL ) LEAVE = .FALSE.
                    IF ( BETA(J,J).GT.BETMAX ) BETMAX = BETA(J,J)
                    IF ( BETA(J,J).LT.BETMIN ) BETMIN = BETA(J,J)
                    ENORM(J) = DELTA(J,J)
   20            CONTINUE
                 RANKD = BETMIN.LE.BETMAX*EPS
         END IF
C
C....... check the conditioning of (BETA) ..............................
C
         IF ( LEAVE .OR. RANKD .OR. I.EQ.NVB*2 .OR. NVB.EQ.1 ) THEN
C
 	    CALL NORM2A (NVB,BETA,WORK,EPS,BMAXN,BMINN) 
C
            IF ( BMINN.LE.BMAXN*EPS*HNDRD .OR. 
     &           BMAXN.LE.AMAXN*EPS**THRQT ) THEN
C
C............... (BETA) is ill-conditioned, set some entries to zero ...
C
                 DO 40 J = 1,NVB
                    IF ( BETA(J,J) .LE. EPS*HNDRD ) THEN
                       IF ( GNRZD ) CALL SETTO0 (NI,BR(1,J),1)
	                            CALL SETTO0 (NI, R(1,J),1)
                       NULLP = NULLP + 1
                       DO 30 K = 1,NVB
                          BETA(J,K)  = ZERO 
   30                  CONTINUE
	            END IF
   40            CONTINUE
C
            END IF
C
            RETURN    
C
         END IF
C
   50 CONTINUE
C
      RETURN
C
C**** end of RFACTR ****************************************************
C
      END
