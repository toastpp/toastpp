      SUBROUTINE SSICHK (LFILE,LPRNT,LRERR,LRMDE,LRWRN,NNEIG,NDEIG,
     &                   NEWSIG,NFARL,NFARR,NNSPNT,NNTRTL,NNTRTR,NREIG,
     &                   NREIGL,NREIGR,NRUNMX,NSFAIL,NSIGMA,NSLOG,NSRLS,
     &                   NSRRS,NTEIG,RSINT,ISINT,INDSI,NSINT,NSIMAX,EIG,
     &                   EIGL,EIGR,ENDL,ENDR,ORIGIN,RADIUS,SFARL,SFARR,
     &                   SIGMA,TRUSTL,TRUSTR,SSLOG,TIME,BIGNUM)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSICHK examines the inertia information                          *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    LRERR  (sio) : code for error messages                           *
C*    LRMDE  (sib) : run mode                                          *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NNEIG  (sii) : number of eigenvalues less than SIGMA             *
C*    NDEIG  (sii) : number of eigenvalues required in the run         *
C*    NEWSIG (sib) : flag for a new starting point                     *
C*    NFARL  (sii) : number of eigenvalues less than SFARL             *
C*    NFARR  (sii) : number of eigenvalues less than SFARR             *
C*    NNSPNT (sii) : number of eigenvalues less than ORIGIN            *
C*    NNTRTL (sii) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sii) : number of eigenvalues less than TRUSTR            *
C*    NREIG  (sii) : number of required eigenvalues                    *
C*    NREIGL (sii) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sii) : number of required eigenvalues greater than EIGL  *
C*    NRUNMX (sii) : maximum number of runs                            *
C*    NSFAIL (sib) : number of factorizations failed                   *
C*    NSIGMA (sii) : number of origin translations                     *
C*    NSLOG  (sib) : number of subintervals recorded in SSLOG          *
C*    NSRLS  (sii) : number of eigenvalues required less    than SIGMA *
C*    NSRRS  (sii) : number of eigenvalues required greater than SIGMA *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    RSINT  (ari) : lower and upper limits of each subinterval        *
C*    ISINT  (ari) : inertia of the lower and upper limits             *
C*    INDSI  (sib) : index of the subinterval to be closed             *
C*    NSINT  (sib) : number of subintervals                            *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    EIGL   (sri) : inferior bound for eigenvalues                    *
C*    EIGR   (sri) : superior bound for eigenvalues                    *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    ORIGIN (sri) : starting-point (first SIGMA)                      *
C*    RADIUS (sri) : radius of convergence                             *
C*    SFARL  (sii) : farthest SIGMA to the left  of ORIGIN             *
C*    SFARR  (sii) : farthest SIGMA to the right of ORIGIN             *
C*    SIGMA  (srb) : origin translation                                *
C*    TRUSTL (sri) : inferior trust bound                              *
C*    TRUSTR (sri) : superior trust bound                              *
C*    SSLOG  (arb) : spectrum slicing history                          *
C*    TIME   (ari) : time table                                        *
C*    BIGNUM (sri) : big number                                        *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZPRT2,SETLRM,SETSSL,SSBACK,SSCHEK,SSORGN                        *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MAX                                                              *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C  
      INTEGER          ISINT(2,*),INDSI,LFILE,LPRNT,LRERR,LRMDE,LRWRN,
     &                 NEWSIG,NFARL,NFARR,NDEIG,NNEIG,NNSPNT,NNTRTL,
     &                 NNTRTR,NREIG,NREIGL,NREIGR,NRUNMX,NSFAIL,
     &                 NSIGMA,NSIMAX,NSINT,NSLOG,
     &                 NSRLS,NSRRS,NTEIG
      DOUBLE PRECISION BIGNUM,EIG(*),EIGL,EIGR,ENDL,ENDR,ORIGIN,
     &                 RADIUS,RSINT(6,*),SFARL,SFARR,SIGMA,
     &                 SSLOG(8,*),TIME(*),TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      INTEGER          NESIKL,NESIKR,NSLS,NSRS,NSSIKL,NSSIKR
      DOUBLE PRECISION SIGMA0 
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MAX
C
C**** executable statements ********************************************
C
      NSLS = 0
      NSRS = 0
      SIGMA0 = SIGMA
C
C.... print the subintervals ...........................................
C
      CALL LZPRT2 (LFILE,LPRNT,NNEIG,NDEIG,NSIGMA,RSINT,
     &             ISINT,NSINT,ENDL,ENDR,SIGMA)
C
C.... check the subintervals ...........................................
C
      CALL SSCHEK (LRERR,LRMDE,NESIKL,NESIKR,NNSPNT,NNTRTL,NNTRTR,
     &             NREIG,NREIGL,NREIGR,NSSIKL,NSSIKR,NTEIG,EIG,
     &             EIGL,EIGR,ORIGIN,SIGMA,INDSI,NSINT,RSINT,
     &             ISINT,TRUSTL,TRUSTR)
C
C.... check for the number of required eigenvalues .....................
C
      IF ( NDEIG .EQ. 0 ) THEN
C
C....... define a new origin ...........................................
C
         CALL SSORGN (LRMDE,LRWRN,NSIMAX,NESIKL,NESIKR,NEWSIG,
     &                NFARL,NFARR,NDEIG,NNTRTL,NNTRTR,NSLS,NSRS,
     &                NSRLS,NSRRS,NSSIKL,NSSIKR,BIGNUM,EIGL,EIGR,
     &                ORIGIN,RADIUS,SFARL,SFARR,SIGMA,INDSI,
     &                NSINT,RSINT,TRUSTL,TRUSTR)
C
         IF ( NSINT .NE. 0 ) THEN
            NSFAIL = NSFAIL + 1
         ELSE
            NSFAIL = 0
         END IF
C
      ELSE
C
C....... check whether SIGMA was applied too far .......................
C
         CALL SSBACK (LFILE,LPRNT,LRWRN,NSIMAX,NNEIG,NNSPNT,NNTRTL,
     &                NNTRTR,NREIGL,NREIGR,NSFAIL,NSIGMA,NSINT,
     &                ORIGIN,SIGMA,TRUSTL,TRUSTR)
C
      END IF
C
C.... store the information ............................................
C
      NSLOG = NSLOG + 1
C
      IF ( NSLOG .GT. NRUNMX ) THEN
         CALL SETLRM (7,LRWRN)
      ELSE
         CALL SETSSL (0,0,NNEIG,0,ENDL,ENDR,SIGMA0,NSLOG,SSLOG,TIME)
      END IF
C
C.... check for failure ................................................
C
      IF ( NSFAIL .GE. 3 ) CALL SETLRM (8,LRWRN)
C
      RETURN 
C
C**** end of SSICHK ****************************************************
C
      END
