      SUBROUTINE BLZKLL (ISTOR,RSTOR,SIGMA,EIG,X,U)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZKLL deals with a premature finalization                       *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aib) : array for integer variables                        *
C*    RSTOR (arb) : array for real variables                           *
C*    SIGMA (sri) : origin translation                                 *
C*    EIG   (arb) : eigenvalue approximations and estimated residuals  *
C*    X     (arb) : Ritz vectors                                       *
C*    U     (ari) : Lanczos vectors at step JL+1                       *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZPRT4,RVMNGR,SETLRM,SSTRSF,TBEIGP                               *
C*                                                                     *
C*  - BLAS kernels:                                                    *
C*                                                                     *
C*    SCOPY,ISAMAX                                                     *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
      REAL             ZERO
      PARAMETER        (ZERO=0.0E0)
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*)
      REAL             EIG(*),RSTOR(*),SIGMA,U(*),X(*)
C
C==== local variables ==================================================
C
      REAL             THRES
      INTEGER          BASIS,BETAR,IDXRES,IWORK,NSLS,NSRS,
     &                 RITZ,RNORM,RWORK,S,TB,THETA
      LOGICAL          GNRZD,OUTER
C
C==== BLAS kernel ======================================================
C
      INTEGER          ISAMAX
C
C**** executable statements ********************************************
C
      RITZ  = ISTOR(IRITZ)
      TB    = ISTOR(ITB)
      BETAR = ISTOR(IBETAR)
      RWORK = ISTOR(IRWORK)
      THETA = ISTOR(ITHETA)
      S     = ISTOR(IS)
      BASIS = ISTOR(IBASIS)
      IWORK = ISTOR(IIWORK)
C
      RNORM = RITZ + ISTOR(JT)
C
      GNRZD = ISTOR(LOPTS).GT.0
C
C.... solve the reduced eigenproblem ...................................
C
      CALL TBEIGP (ISTOR(JT)    ,ISTOR(JTMAX) ,ISTOR(NVB)   , 
     &             ISTOR(LRWRN) ,NSLS         ,NSRS         ,
     &             ISTOR(NSRLS) ,ISTOR(NSRRS) ,ISTOR(NULLDR),
     &             RSTOR(BIGNUM),RSTOR(EPS)   ,RSTOR(REPS)  ,
     &             RSTOR(BETAR) ,RSTOR(TB)    ,RSTOR(S)     ,
     &             RSTOR(THETA) ,RSTOR(THETA0),RSTOR(THETAL),
     &             RSTOR(THETAR),THRES        ,RSTOR(THRSH) ,
     &             RSTOR(RWORK) ,ISTOR(IWORK) ,GNRZD        ,
     &             .FALSE.      ,OUTER        ,.TRUE.       )
C
      IF ( GNRZD ) CALL SSTRSF (ISTOR(JT)    ,ISTOR(JT)    ,
     &                          RSTOR(THETA) ,RSTOR(S)     ,
     &                          SIGMA        ,RSTOR(THETA0))
C
      IDXRES = ISAMAX(ISTOR(JT),RSTOR(RNORM),1)
C
      IF ( RSTOR(RNORM+IDXRES-1) .EQ. ZERO ) CALL SCOPY (ISTOR(JT)*2 ,
     &                                                   RSTOR(THETA),1,
     &                                                   RSTOR(RITZ) ,1)
C
C.... compute the eigenvectors .........................................
C
      CALL RVMNGR (ISTOR(JL)    ,ISTOR(JT)    ,ISTOR(NVB)   ,
     &             ISTOR(LEIG)  ,ISTOR(LNI)   ,ISTOR(NI)    ,
     &             ISTOR(NQMAX) ,ISTOR(NXMAX) ,ISTOR(FHNDL) ,
     &             ISTOR(LCOMM) ,ISTOR(LRERR) ,ISTOR(LRWRN) ,
     &             RSTOR(RITZ)  ,RSTOR(S)     ,SIGMA        ,
     &             ZERO         ,RSTOR(BETAR) ,RSTOR(BASIS) ,
     &             U            ,ISTOR(NTEIG) ,EIG          ,
     &             X            ,RSTOR(RWORK) ,GNRZD        ,
     &             .FALSE.      )
C
C.... print Ritz values ................................................
C
      CALL LZPRT4 (0            ,ISTOR(JT)    ,ISTOR(LFILE) ,
     &             ISTOR(LPRNT) ,0            ,0            ,
     &             0            ,0            ,RSTOR(RITZ)  ,
     &             .FALSE.      )
C
      CALL SETLRM (17,ISTOR(LRWRN))
C
      RETURN 
C
C**** end of BLZKLL ****************************************************
C
      END
