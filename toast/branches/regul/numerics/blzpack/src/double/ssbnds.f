      SUBROUTINE SSBNDS (LRWRN,N,NNEIG,NDEIG,NFARL,NFARR,NNSPNT,NNTRTL,
     &                   NNTRTR,NREIGL,NREIGR,NSIGMA,NSRLS,NSRRS,NTEIG,
     &                   RSINT,ISINT,INDSI,NSINT,NSIMAX,EIG,EIGL,EIGR,
     &                   ENDL,ENDR,SFARL,SFARR,SIGMA,ORIGIN,THETA0,
     &                   THETAL,THETAR,TRUSTL,TRUSTR,TIME,BIGNUM)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSBNDS initializes the computational subintervals                *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LRWRN  (sio) : code for warning messages                         *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    NNEIG  (sii) : number of eigenvalues less than SIGMA             *
C*    NDEIG  (sio) : number of eigenvalues required in the run         *
C*    NFARL  (sio) : number of eigenvalues less than SFARL             *
C*    NFARR  (sio) : number of eigenvalues less than SFARR             *
C*    NNSPNT (sio) : number of eigenvalues less than ORIGIN            *
C*    NNTRTL (sio) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sio) : number of eigenvalues less than TRUSTR            *
C*    NREIGL (sii) : number of eigenvalues required less    than EIGR  *
C*    NREIGR (sii) : number of eigenvalues required greater than EIGL  *
C*    NSIGMA (sii) : number of origin translations                     *
C*    NSRLS  (sio) : number of eigenvalues required less    than SIGMA *
C*    NSRRS  (sio) : number of eigenvalues required greater than SIGMA *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    RSINT  (arb) : lower and upper limits of each subinterval        *
C*    ISINT  (arb) : inertias of the lower and upper limits            *
C*    INDSI  (sii) : index of the subinterval to be closed             *
C*    NSINT  (sii) : number of subintervals                            *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    EIGL   (sri) : inferior bound for eigenvalues                    *
C*    EIGR   (sri) : superior bound for eigenvalues                    *
C*    ENDL   (sro) : inferior bound for Ritz values                    *
C*    ENDR   (sro) : superior bound for Ritz values                    *
C*    SFARL  (sro) : farthest SIGMA to the left  of ORIGIN             *
C*    SFARR  (sro) : farthest SIGMA to the right of ORIGIN             *
C*    SIGMA  (sri) : origin translation                                *
C*    ORIGIN (sro) : starting point (first SIGMA)                      *
C*    THETA0 (sri) : reference point for THETA                         *
C*    THETAL (sro) : inferior limit to converged eigenvalues           *
C*    THETAR (sro) : superior limit to converged eigenvalues           *
C*    TRUSTL (sro) : inferior trust bound                              *
C*    TRUSTR (sro) : superior trust bound                              *
C*    TIME   (arb) : time table                                        *
C*    BIGNUM (sri) : big number                                        *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    NEIGAB,SETLRM,SETTO0,SITIME,SSMOVB                               *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    SIGN                                                             *
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
      INTEGER          ISINT(2,*),INDSI,LRWRN,N,NFARL,NFARR,NDEIG,NNEIG,
     &                 NNSPNT,NNTRTL,NNTRTR,NREIGL,NREIGR,NSIGMA,NSIMAX,
     &                 NSINT,NSRLS,NSRRS,NTEIG
      DOUBLE PRECISION BIGNUM,EIG(*),EIGL,EIGR,ENDL,ENDR,ORIGIN,
     &                 SFARL,SFARR,SIGMA,RSINT(6,*),THETA0,
     &                 THETAL,THETAR,TIME(*),TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      INTEGER          NELSI,NERSI,NSLXI,NSRXI
      DOUBLE PRECISION A,B
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
      INTEGER          NEIGAB
C
C==== intrinsic function ===============================================
C
      INTRINSIC        SIGN
C
C**** executable statements ********************************************
C
C.... increment factorization counter ..................................
C
      NSIGMA = NSIGMA + 1 
C
      TIME(12) = TIME(11) 
      TIME( 8) = SITIME(TIME(11))
      TIME( 3) = TIME(3) + TIME(8)
C
      IF ( NSIGMA.EQ.1 ) THEN
C
C....... set basic bounds and number of required eigenpairs ............
C
         INDSI  = 1
         NSINT  = 2
         NFARL  = NNEIG
         NFARR  = NNEIG
         NNSPNT = NNEIG
         NNTRTL = NNEIG
         NNTRTR = NNEIG
         NREIGL = MIN(NREIGL,  NNEIG)
         NREIGR = MIN(NREIGR,N-NNEIG)
         ORIGIN = SIGMA
         TRUSTL = SIGMA
         TRUSTR = SIGMA
         SFARL  = SIGMA
         SFARR  = SIGMA
C
         IF ( (NREIGL.EQ.0) .AND. (NREIGR.EQ.0) ) CALL SETLRM (1,LRWRN)
C
         CALL SETTO0 (NSIMAX*6,RSINT,1)
C
         IF      ( ORIGIN.NE.EIGL .AND. ORIGIN.NE.EIGR ) THEN
C
C............... run around ORIGIN .....................................
C
                 ISINT(1,1) = 0
                 ISINT(2,1) = N
                 RSINT(1,1) = -BIGNUM
                 RSINT(2,1) = -BIGNUM
                 RSINT(3,1) = -BIGNUM
                 RSINT(4,1) = +BIGNUM
                 RSINT(5,1) = +BIGNUM
                 RSINT(6,1) = +BIGNUM
C
         ELSE IF ( NREIGR .EQ. 0 ) THEN
C
C............... run from the right to the left ........................
C
                 ISINT(1,1) = 0
                 ISINT(2,1) = NNEIG
                 RSINT(1,1) = EIGL
                 RSINT(2,1) = EIGL
                 RSINT(3,1) = EIGL
                 RSINT(4,1) = SIGMA
                 RSINT(5,1) = SIGMA
                 RSINT(6,1) = SIGMA
C
         ELSE IF ( NREIGL .EQ. 0 ) THEN
C
C............... run from the left to the right ........................
C
                 ISINT(1,1) = NNEIG
                 ISINT(2,1) = N
                 RSINT(1,1) = SIGMA
                 RSINT(2,1) = SIGMA
                 RSINT(3,1) = SIGMA
                 RSINT(4,1) = EIGR
                 RSINT(5,1) = EIGR
                 RSINT(6,1) = EIGR
C
         END IF
C
      END IF
C
C.... update bounds ....................................................
C
      CALL SSMOVB(-1,NSINT,INDSI+1,ISINT,RSINT)
C
      ISINT(2,INDSI  ) = NNEIG 
      ISINT(1,INDSI+1) = NNEIG 
C
      RSINT(4,INDSI  ) = SIGMA 
      RSINT(5,INDSI  ) = SIGMA 
      RSINT(6,INDSI  ) = SIGMA 
      RSINT(1,INDSI+1) = SIGMA 
      RSINT(2,INDSI+1) = SIGMA 
      RSINT(3,INDSI+1) = SIGMA 
C
C.... move the inferior bound if possible ..............................
C
      IF ( (SIGMA.LE.ORIGIN) .AND. 
     &     ((NREIGL.LE.(NNSPNT-NNEIG)).OR.(SIGMA.LE.EIGL)) ) THEN
           ISINT(1,    1) = NNEIG 
           RSINT(1,    1) = SIGMA 
           RSINT(2,    1) = SIGMA 
           RSINT(3,    1) = SIGMA 
           NFARL = NNEIG
           SFARL = SIGMA
           EIGL  = SIGMA
      END IF
C
C.... move the superior bound if possible ..............................
C
      IF ( (SIGMA.GE.ORIGIN) .AND. 
     &     ((NREIGR.LE.(NNEIG-NNSPNT)).OR.(SIGMA.GE.EIGR)) ) THEN
           ISINT(2,NSINT) = NNEIG 
           RSINT(6,NSINT) = SIGMA 
           RSINT(5,NSINT) = SIGMA 
           RSINT(4,NSINT) = SIGMA 
           NFARR = NNEIG
           SFARR = SIGMA
           EIGR  = SIGMA
      END IF
C
C.... number of required eigenvalues to the left  of SIGMA .............
C
      A = RSINT(1,INDSI  )
      B = RSINT(6,INDSI  )
      NSLXI = NEIGAB(NTEIG,A,B,EIG)
C
      IF      ( A.NE.EIGL .OR. ( A.EQ.EIGL .AND. A.EQ.SFARL ) ) THEN
              NELSI = ISINT(2,INDSI  ) - ISINT(1,INDSI  )
              NSRLS = NELSI - NSLXI
      ELSE IF ( B.EQ.SIGMA ) THEN
              NSRLS = NREIGL - (NSLXI+(NNSPNT-NNEIG))
      ELSE
              NSRLS = 0 
      END IF
C
C.... number of required eigenvalues to the right of SIGMA .............
C
      A = RSINT(1,INDSI+1)
      B = RSINT(6,INDSI+1)
      NSRXI = NEIGAB(NTEIG,A,B,EIG)
C
      IF      ( B.NE.EIGR .OR. ( B.EQ.EIGR .AND. B.EQ.SFARR ) ) THEN
              NERSI = ISINT(2,INDSI+1) - ISINT(1,INDSI+1)
              NSRRS = NERSI - NSRXI
      ELSE IF ( A.EQ.SIGMA ) THEN
              NSRRS = NREIGR - (NSRXI+(NNEIG-NNSPNT))
      ELSE
              NSRRS = 0 
      END IF
C
C.... active interval ..................................................
C
      ENDL = RSINT(1,INDSI  )
      ENDR = RSINT(6,INDSI+1)
C
C.... check the number of required eigenvalues .........................
C
      IF ( NSRLS .LE. 0 ) THEN
         ENDL  = SIGMA
         NSRLS = 0
      END IF
      IF ( NSRRS .LE. 0 ) THEN
         ENDR  = SIGMA
         NSRRS = 0
      END IF
C
      NDEIG = NSRLS + NSRRS
C
C.... define limits for accepting thetas ...............................
C
      IF ( SIGMA .EQ. ENDL ) THEN
         IF ( THETA0 .EQ. ZERO ) THEN
            THETAL = -BIGNUM
         ELSE
            THETAL = -BIGNUM*SIGN(ONE,SIGMA)
         END IF
      ELSE
         IF ( THETA0 .EQ. ZERO ) THEN
            THETAL = ONE/(ENDL-SIGMA)
         ELSE
            THETAL = ENDL/(ENDL-SIGMA)
         END IF
      END IF
C
      IF ( SIGMA .EQ. ENDR ) THEN
         IF ( THETA0 .EQ. ZERO ) THEN
            THETAR = +BIGNUM
         ELSE
            THETAR = +BIGNUM*SIGN(ONE,SIGMA)
         END IF
      ELSE
         IF ( THETA0 .EQ. ZERO ) THEN
            THETAR = ONE/(ENDR-SIGMA)
         ELSE
            THETAR = ENDR/(ENDR-SIGMA)
         END IF
      END IF
C
C.... check for an empty computational interval ........................
C
      IF ( SIGMA.LE.EIGL .AND. SIGMA.LE.EIGR .AND.
     &     SIGMA.LT.ORIGIN .AND. NNEIG.EQ.NNSPNT ) CALL SETLRM (1,LRWRN)
C
      IF ( SIGMA.GE.EIGL .AND. SIGMA.GE.EIGR .AND.
     &     SIGMA.GT.ORIGIN .AND. NNEIG.EQ.NNSPNT ) CALL SETLRM (1,LRWRN)
C
      RETURN 
C
C**** end of SSBNDS ****************************************************
C
      END
