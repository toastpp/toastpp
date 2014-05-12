      SUBROUTINE SSCLSD (LRMDE,NNSPNT,NNTRTL,NNTRTR,NREIG,NREIGL,
     &                   NREIGR,NSINT,NSLXI,NSRXI,NTEIG,EIG,
     &                   EIGL,EIGR,ORIGIN,TRUSTL,TRUSTR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSCLSD checks whether all required eigenvalues were found        *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LRMDE  (sib) : run mode                                          *
C*    NNSPNT (sii) : number of eigenvalues less than ORIGIN            *
C*    NNTRTL (sii) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sii) : number of eigenvalues less than TRUSTR            *
C*    NREIG  (sii) : number of required eigenvalues                    *
C*    NREIGL (sii) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sii) : number of required eigenvalues greater than EIGL  *
C*    NSINT  (sii) : number of subintervals                            *
C*    NSLXI  (sio) : number of eigenvalues less    than EIGR           *
C*    NSRXI  (sio) : number of eigenvalues greater than EIGL           *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    EIGL   (sri) : inferior bound for eigenvalues                    *
C*    EIGR   (sri) : superior bound for eigenvalues                    *
C*    ORIGIN (sri) : starting-point (first SIGMA)                      *
C*    TRUSTL (sri) : inferior trust bound                              *
C*    TRUSTR (sri) : superior trust bound                              *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    NEIGAB                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LRMDE,NNSPNT,NNTRTL,NNTRTR,NREIG,NREIGL,
     &                 NREIGR,NSINT,NSLXI,NSRXI,NTEIG
      DOUBLE PRECISION EIG(*),EIGL,EIGR,ORIGIN,TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      DOUBLE PRECISION RANGE
C
C==== subprogram =======================================================
C
      INTEGER          NEIGAB
C
C**** executable statements ********************************************
C
C.... check whether all required eigenvalues were found ................
C
      IF      ( NSINT .EQ. 0 ) THEN
C
C............ no interval left (exit) ..................................
C
              TRUSTL = EIGL
              TRUSTR = EIGR
              LRMDE = 0
C
      ELSE IF ( ORIGIN.NE.EIGL .AND. ORIGIN.NE.EIGR ) THEN
C
C............ run around a reference value .............................
C
	      IF ( NTEIG .EQ. (NREIGL+NREIGR) ) THEN
C
C............... check based on the number of eigenvalues converged ....
C
                 NSLXI = NEIGAB(NTEIG,TRUSTL,ORIGIN,EIG)
                 NSRXI = NEIGAB(NTEIG,ORIGIN,TRUSTR,EIG)
C
	         IF ( NSLXI.GE.NREIGL .AND. NSRXI.GE.NREIGR ) LRMDE = 0
C
	      ELSE
C
C............... check based on the trust regions ......................
C
	         RANGE = MIN(ORIGIN-TRUSTL,TRUSTR-ORIGIN)
                 NSLXI = NEIGAB(NTEIG,ORIGIN-RANGE,ORIGIN,EIG)
                 NSRXI = NEIGAB(NTEIG,ORIGIN+RANGE,ORIGIN,EIG)
C
	         IF ( NSLXI+NSRXI .GE. NREIG )                LRMDE = 0
C
	      END IF 
C
      ELSE IF ( NSINT .EQ. 1 ) THEN
C
C............ one interval left ........................................
C
              IF      ( NREIGR .EQ. 0 ) THEN
                      NSRXI = NEIGAB(NTEIG,TRUSTL,ORIGIN,EIG)
	              IF ( NSRXI.GE.NREIGL .AND. 
     &                     NSRXI.EQ.NNSPNT-NNTRTL )           LRMDE = 0
              ELSE IF ( NREIGL .EQ. 0 ) THEN
                      NSLXI = NEIGAB(NTEIG,ORIGIN,TRUSTR,EIG)
	              IF ( NSLXI.GE.NREIGR .AND.
     &                     NSLXI.EQ.NNTRTR-NNSPNT )           LRMDE = 0
              END IF
C
      END IF
C
      IF ( LRMDE .EQ. 0 ) NSINT = 0
C
      RETURN 
C
C**** end of SSCLSD ****************************************************
C
      END
