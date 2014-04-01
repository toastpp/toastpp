      SUBROUTINE BLZEIG (ISTOR,RSTOR,LFLAG,SIGMA,NNEIG,EIG,
     &                   X,U,AGAIN,ENDON,RBXON,STRON)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZEIG deals with the converged eigenvalues and eigenvectors     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aib) : array for integer variables                        *
C*    RSTOR (arb) : array for real variables                           *
C*    LFLAG (sio) : reverse communication flag                         *
C*    SIGMA (srb) : origin translation                                 *
C*    NNEIG (sri) : number of eigenvalues less than SIGMA              *
C*    EIG   (arb) : eigenvalue approximations and estimated residuals  *
C*    X     (arb) : Ritz vectors                                       *
C*    U     (ari) : Lanczos vectors at step JL+1                       *
C*    AGAIN (slo) : loop control flag                                  *
C*    ENDON (slo) : finalization flag                                  *
C*    RBXON (slo) : (B)*(X) computation flag                           *
C*    STRON (slo) : run starting flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    RVMNGR,SETLRM,SETSSL,SITIME,SSBEXT,                              *
C*    SSCHEK,SSORGN,SSSPEC,SSTRSF,TBEIGP                               *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    MAX                                                              *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO=0.0D0)
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*),LFLAG,NNEIG
      DOUBLE PRECISION EIG(*),RSTOR(*),SIGMA,U(*),X(*)
      LOGICAL          AGAIN,ENDON,RBXON,STRON
C
C==== local variables ==================================================
C
      INTEGER          BASIS,BETAR,IWORK,NCEIG,NESIKL,NESIKR,
     &                 NSLS,NSRS,NSSIKL,NSSIKR,RITZ,RNORM,
     &                 RSINT,RWORK,SSLOG,S,TB,THETA,TIME
      DOUBLE PRECISION RITZL,RITZR,SIGMAL,SIGMAR,THRES,T0,T1,T2
      LOGICAL          GNRZD,OUTER,PURFY,SLICE
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
C
C==== intrinsic function ===============================================
C
      INTRINSIC        MAX
C
C**** executable statements ********************************************
C
      TIME  = ISTOR(ITIME)
      RSINT = ISTOR(IRSINT)
      SSLOG = ISTOR(ISSLOG)
      RITZ  = ISTOR(IRITZ)
      TB    = ISTOR(ITB)
      BETAR = ISTOR(IBETAR)
      THETA = ISTOR(ITHETA)
      S     = ISTOR(IS)
      BASIS = ISTOR(IBASIS)
      RWORK = ISTOR(IRWORK)
      IWORK = ISTOR(IIWORK)
C
      RNORM = RITZ + ISTOR(JT)
C
      ISTOR(NBX) = ISTOR(NTEIG)
C
      GNRZD = ISTOR(LOPTS  ).GT.0
      SLICE = ISTOR(LOPTS+1).GT.0
      PURFY = ISTOR(LOPTS+2).GT.0
C
C.... solve the reduced eigenproblem ...................................
C
      T0 = SITIME(ZERO)
C
      CALL TBEIGP (ISTOR(JT)    ,ISTOR(JTMAX) ,ISTOR(NVB)   ,
     &             ISTOR(LRWRN) ,NSLS         ,NSRS         ,
     &             ISTOR(NSRLS) ,ISTOR(NSRRS) ,ISTOR(NULLDR),
     &             RSTOR(BIGNUM),RSTOR(EPS)   ,RSTOR(REPS)  ,
     &             RSTOR(BETAR) ,RSTOR(TB)    ,RSTOR(S)     ,
     &             RSTOR(THETA) ,RSTOR(THETA0),RSTOR(THETAL),
     &             RSTOR(THETAR),THRES        ,RSTOR(THRSH) ,
     &             RSTOR(RWORK) ,ISTOR(IWORK) ,GNRZD        ,
     &             SLICE        ,OUTER        ,.TRUE.       )
C
      IF ( GNRZD ) CALL SSTRSF (ISTOR(JT)    ,ISTOR(JT)    ,
     &                          RSTOR(THETA) ,RSTOR(S)     ,
     &                          SIGMA        ,RSTOR(THETA0))
C
      T1 = SITIME(ZERO)
C
C.... compute the eigenvectors .........................................
C
      CALL RVMNGR (ISTOR(JL)    ,ISTOR(JT)    ,ISTOR(NVB)   ,
     &             ISTOR(LEIG)  ,ISTOR(LNI)   ,ISTOR(NI)    ,
     &             ISTOR(NQMAX) ,ISTOR(NXMAX) ,ISTOR(FHNDL) ,
     &             ISTOR(LCOMM) ,ISTOR(LRERR) ,ISTOR(LRWRN) ,
     &             RSTOR(THETA) ,RSTOR(S)     ,SIGMA        ,
     &             RSTOR(THETA0),RSTOR(BETAR) ,RSTOR(BASIS) ,
     &             U            ,ISTOR(NTEIG) ,EIG          ,
     &             X            ,RSTOR(RWORK) ,GNRZD        ,
     &             PURFY        )
C
      T2 = SITIME(ZERO)
C
      RSTOR(TIME+5) = RSTOR(TIME+5) + (T1-T0)
      RSTOR(TIME+6) = RSTOR(TIME+6) + (T2-T1)
C
      IF      ( ISTOR(LRERR).NE.0 .OR. ISTOR(LRWRN).NE.0 ) THEN
C
C............ prepare to exit ..........................................
C
              ISTOR(NEWSIG) = 0
              ISTOR(LRMDE) = 0
              ISTOR(NDEIG) = 0
C
      ELSE IF ( SLICE ) THEN
C
C............ store the information ....................................
C
              NCEIG = NSLS + NSRS
C
              CALL SETSSL (ISTOR(JL)    ,ISTOR(JT)    ,NNEIG        ,
     &                     NCEIG        ,RSTOR(ENDL)  ,RSTOR(ENDR)  ,
     &                     SIGMA        ,ISTOR(NSLOG) ,RSTOR(SSLOG) ,
     &                     RSTOR(TIME)  )
C
C............ check the spectrum .......................................
C
              CALL SSSPEC (ISTOR(JT)    ,ISTOR(LFILE) ,ISTOR(LPRNT) ,
     &                     ISTOR(N)     ,ISTOR(NEWSIG),NNEIG        ,
     &                     ISTOR(NONEWS),NSLS         ,NSRS         ,
     &                     ISTOR(NSRLS) ,ISTOR(NSRRS) ,RSTOR(ENDL)  ,
     &                     RSTOR(ENDR)  ,RSTOR(RADIUS),RSTOR(RITZ)  ,
     &                     RITZL        ,RITZR        ,RSTOR(RNORM) ,
     &                     RSTOR(SFARL) ,RSTOR(SFARR) ,SIGMA        ,
     &                     SIGMAL       ,SIGMAR       )
C
C............ extend the boundaries ....................................
C
              CALL SSBEXT (ISTOR(NSINT) ,RITZL        ,RITZR        ,
     &                     SIGMA        ,SIGMAL       ,SIGMAR       ,
     &                     RSTOR(RSINT) )
C
C............ check the subintervals ...................................
C
              CALL SSCHEK (ISTOR(LRERR) ,ISTOR(LRMDE) ,NESIKL       ,
     &                     NESIKR       ,ISTOR(NNSPNT),ISTOR(NNTRTL),
     &                     ISTOR(NNTRTR),ISTOR(NREIG) ,ISTOR(NREIGL),
     &                     ISTOR(NREIGR),NSSIKL       ,NSSIKR       ,
     &                     ISTOR(NTEIG) ,EIG          ,RSTOR(EIGL)  ,
     &                     RSTOR(EIGR)  ,RSTOR(ORIGIN),SIGMA        ,
     &                     ISTOR(INDSI) ,ISTOR(NSINT) ,RSTOR(RSINT) ,
     &                     ISTOR(ISINT) ,RSTOR(TRUSTL),RSTOR(TRUSTR))
C
C............ define a new origin ......................................
C
              CALL SSORGN (ISTOR(LRMDE) ,ISTOR(LRWRN) ,ISTOR(NSIMAX),
     &                     NESIKL       ,NESIKR       ,ISTOR(NEWSIG),
     &                     ISTOR(NFARL) ,ISTOR(NFARR) ,ISTOR(NDEIG) ,
     &                     ISTOR(NNTRTL),ISTOR(NNTRTR),NSLS         ,
     &                     NSRS         ,ISTOR(NSRLS) ,ISTOR(NSRRS) ,
     &                     NSSIKL       ,NSSIKR       ,RSTOR(BIGNUM),
     &                     RSTOR(EIGL)  ,RSTOR(EIGR)  ,RSTOR(ORIGIN),
     &                     RSTOR(RADIUS),RSTOR(SFARL) ,RSTOR(SFARR) ,
     &                     SIGMA        ,ISTOR(INDSI) ,ISTOR(NSINT) ,
     &                     RSTOR(RSINT) ,RSTOR(TRUSTL),RSTOR(TRUSTR))
C
              IF ( ISTOR(NEWSIG).EQ.1 ) ISTOR(NSLOG) = ISTOR(NSLOG) + 1
C
              LFLAG = 3
C
      ELSE IF ( GNRZD ) THEN
C
C............ generalized eigenproblem .................................
C
              IF ( OUTER ) THEN
                 ISTOR(NSRLS) = 0
                 ISTOR(NSRRS) = 0
                 ISTOR(NDEIG) = 0
              ELSE
                 ISTOR(NSRLS) = MAX(0,ISTOR(NSRLS)-NSLS)
                 ISTOR(NSRRS) = MAX(0,ISTOR(NSRRS)-NSRS)
                 ISTOR(NDEIG) = ISTOR(NSRLS) + ISTOR(NSRRS)
              END IF
C
      ELSE
C
C............ standard eigenproblem ....................................
C
              ISTOR(NDEIG) = MAX(0,ISTOR(NREIG)-ISTOR(NTEIG))
C
      END IF
C
      RSTOR(TIME+11) = ZERO 
C
      IF ( ISTOR(JT) .GE. ISTOR(N) ) ISTOR(NDEIG) = 0
C
      IF ( MAX(ISTOR(NRUN),ISTOR(NSLOG)) .GE. ISTOR(NRUNMX) ) THEN
         CALL SETLRM (7,ISTOR(LRWRN))
         RETURN
      END IF
C
C.... set flags accordingly ............................................
C
      STRON = ISTOR(NDEIG).NE.0 .AND. .NOT.GNRZD
      RBXON = ( ISTOR(NDEIG).NE.0 .AND. GNRZD ) .OR. PURFY
      ENDON = ( ISTOR(NSINT).EQ.0 .AND. SLICE .AND. .NOT.RBXON ) .OR. 
     &        ( ISTOR(NDEIG).EQ.0 .AND. .NOT.SLICE .AND. .NOT.PURFY )
C
      STRON = ISTOR(LRERR).EQ.0 .AND. ISTOR(LRWRN).EQ.0 .AND. STRON
      RBXON = ISTOR(LRERR).EQ.0 .AND. ISTOR(LRWRN).EQ.0 .AND. RBXON
      ENDON = ISTOR(LRERR).NE.0 .OR.  ISTOR(LRWRN).NE.0 .OR.  ENDON
      AGAIN = ENDON .OR. RBXON .OR. STRON
C
      RETURN 
C
C**** end of BLZEIG ****************************************************
C
      END
