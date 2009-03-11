      SUBROUTINE SSORGN (LRMDE,LRWRN,NSIMAX,NESIKL,NESIKR,NEWSIG,
     &                   NFARL,NFARR,NDEIG,NNTRTL,NNTRTR,NSLS,NSRS,
     &                   NSRLS,NSRRS,NSSIKL,NSSIKR,BIGNUM,EIGL,EIGR,
     &                   ORIGIN,RADIUS,SFARL,SFARR,SIGMA,INDSI,
     &                   NSINT,RSINT,TRUSTL,TRUSTR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSORGN sets the origin translation                               *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LRMDE  (sib) : run mode                                          *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    NESIKL (sii) : number of eigenvalues in the left  subinterval    *
C*    NESIKR (sii) : number of eigenvalues in the right subinterval    *
C*    NEWSIG (sii) : flag for a new starting point                     *
C*    NFARL  (sii) : number of eigenvalues less than SFARL             *
C*    NFARR  (sii) : number of eigenvalues less than SFARR             *
C*    NDEIG  (sio) : number of eigenvalues required in the run         *
C*    NNTRTL (sii) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sii) : number of eigenvalues less than TRUSTR            *
C*    NSLS   (sii) : number of eigenvalues converged < than SIGMA      *
C*    NSRS   (sii) : number of eigenvalues converged > than SIGMA      *
C*    NSRLS  (sib) : number of eigenvalues required  < than SIGMA      *
C*    NSRRS  (sib) : number of eigenvalues required  > than SIGMA      *
C*    NSSIKL (sii) : number of eigenvalues in the left  subinterval    *
C*    NSSIKR (sii) : number of eigenvalues in the right subinterval    *
C*    BIGNUM (sri) : big number                                        *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    EIGL   (sri) : inferior bound for eigenvalues                    *
C*    EIGR   (sri) : superior bound for eigenvalues                    *
C*    ORIGIN (sri) : starting-point (first SIGMA)                      *
C*    RADIUS (sri) : radius of convergence                             *
C*    SFARL  (sii) : farthest SIGMA to the left  of ORIGIN             *
C*    SFARR  (sii) : farthest SIGMA to the right of ORIGIN             *
C*    SIGMA  (sro) : origin translation                                *
C*    INDSI  (sii) : index of the subinterval to be closed             *
C*    NSINT  (sib) : number of subintervals                            *
C*    RSINT  (ari) : lower and upper limits of each subinterval        *
C*    TRUSTL (srb) : inferior trust bound                              *
C*    TRUSTR (srb) : superior trust bound                              *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SETLRM,SSRSTR,SSSPLT                                             *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    MAX,MIN                                                          *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          INDSI,LRMDE,LRWRN,NSIMAX,NESIKL,NESIKR,NEWSIG,
     &                 NFARL,NFARR,NDEIG,NNTRTL,NNTRTR,NSINT,NSLS,
     &                 NSRS,NSRLS,NSRRS,NSSIKL,NSSIKR
      DOUBLE PRECISION BIGNUM,EIGL,EIGR,ORIGIN,RADIUS,SFARL,SFARR,
     &                 SIGMA,RSINT(6,NSIMAX),TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      INTEGER          I
      DOUBLE PRECISION LOWER,UPPER
      LOGICAL          EXTNDL,EXTNDR

C==== subprogram =======================================================
C
      LOGICAL          SSRSTR
C
C==== intrinsic functions ==============================================
C
      INTRINSIC        MAX,MIN
C
C**** executable statements ********************************************
C
      IF ( LRMDE .EQ. 0 ) THEN
C
C....... work done .....................................................
C
         TRUSTL = MIN(TRUSTL,SFARL)
         TRUSTR = MAX(TRUSTR,SFARR)
         NNTRTL = MIN(NNTRTL,NFARL)
         NNTRTR = MAX(NNTRTR,NFARR)
         NEWSIG = 0
         NDEIG  = 0
         NSINT  = 0
C
         RETURN
C
      ELSE
C
C....... make sure SIGMA is a bound of a valid interval ................
C
         DO 10 I = 1,NSINT
            IF ( SIGMA.EQ.RSINT(1,I) .OR. SIGMA.EQ.RSINT(6,I) ) GO TO 20
   10    CONTINUE
C
C        SIGMA is not a bound: a new sigma is mandatory 
C
         NEWSIG = 4
C
      END IF
C
   20 CONTINUE
C
      IF      ( SSRSTR(LRMDE,NESIKL,NESIKR,NEWSIG,NSLS,NSRS,
     &                 NSRLS,NSRRS,NSSIKL,NSSIKR,
     &                 ORIGIN,SIGMA) ) THEN
C
C............ reinitialize with the same sigma .........................
C
              IF ( NSRLS .GT. 0 ) THEN
                 NSRLS = MAX(0,NSRLS-NSLS)
              ELSE
                 NSRLS = MAX(0,NESIKL-NSSIKL)
              END IF
              IF ( NSRRS .GT. 0 ) THEN
                 NSRRS = MAX(0,NSRRS-NSRS)
              ELSE
                 NSRRS = MAX(0,NESIKR-NSSIKR)
              END IF
              NDEIG  = NSRLS + NSRRS
              LRMDE  = 5 
              NEWSIG = 1
C
      ELSE IF ( NSINT .LT. NSIMAX ) THEN
C
C............ extend the bounds or split a subinterval .................
C
              NDEIG  = 0
              LRMDE  = 6
              NEWSIG = 2
C
              LOWER  = RSINT(1,INDSI)
              UPPER  = RSINT(6,INDSI)
              EXTNDL = SIGMA.GT.EIGL .AND. LOWER.LT.ORIGIN .AND.
     &                 SFARL.NE.EIGL
              EXTNDR = SIGMA.LT.EIGR .AND. UPPER.GT.ORIGIN .AND.
     &                 SFARR.NE.EIGR
C
              IF      ( LOWER .EQ. -BIGNUM ) THEN
C
C.................... infinite bound to the left .......................
C
                      IF ( RSINT(6,INDSI) .EQ. RSINT(5,INDSI) ) THEN
                         SIGMA = RSINT(5,INDSI) - RADIUS
                      ELSE
                         SIGMA = RSINT(4,INDSI)
                      END IF
C
              ELSE IF ( UPPER .EQ. +BIGNUM ) THEN
C
C.................... infinite bound to the right ......................
C
                      IF ( RSINT(1,INDSI) .EQ. RSINT(2,INDSI) ) THEN
                         SIGMA = RSINT(2,INDSI) + RADIUS
                      ELSE
                         SIGMA = RSINT(3,INDSI)
                      END IF
C
              ELSE IF ( NSINT.EQ.1 .AND. (EXTNDL.OR.EXTNDR) ) THEN
C
C.................... extend the bound .................................
C
                      IF      ( LOWER .LT. ORIGIN ) THEN
C
C............................ bound extended to the left ...............
C
                              SIGMA = RSINT(4,INDSI)
C
                              IF ( SIGMA .EQ. UPPER ) THEN
                                 SIGMA = SIGMA - RADIUS
                              END IF
C
                      ELSE IF ( UPPER .GT. ORIGIN ) THEN
C
C............................ bound extended to the right ..............
C
                              SIGMA = RSINT(3,INDSI)
C
                              IF ( SIGMA .EQ. LOWER ) THEN
                                 SIGMA = SIGMA + RADIUS
                              END IF
C
                      END IF
C
                      IF ( SIGMA .LE. EIGL ) LRMDE = 7
                      IF ( SIGMA .GE. EIGR ) LRMDE = 8
C
              ELSE
C
C.................... split a subinterval ..............................
C
                      CALL SSSPLT (INDSI,SIGMA,RSINT)
C
                      LRMDE = 4
C
              END IF
C
              NSINT = NSINT + 1
C
      ELSE
C
C............ algorithm unable to deal with subinterval ................
C
              CALL SETLRM (8,LRWRN)
              NEWSIG = 0
C
      END IF
C
      RETURN
C
C**** end of SSORGN ****************************************************
C
      END
