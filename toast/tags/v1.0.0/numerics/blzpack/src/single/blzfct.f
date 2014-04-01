      SUBROUTINE BLZFCT (ISTOR,RSTOR,SIGMA,NNEIG,EIG,
     &                   AGAIN,ENDON,RBXON,STRON)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZFCT updates all computational subintervals                    *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aib) : array for integer variables                        *
C*    RSTOR (aib) : array for real variables                           *
C*    SIGMA (sri) : origin translation                                 *
C*    NNEIG (sri) : number of eigenvalues less than SIGMA              *
C*    EIG   (ari) : eigenvalue approximations and estimated residuals  *
C*    AGAIN (slo) : loop control flag                                  *
C*    ENDON (slo) : finalization flag                                  *
C*    RBXON (slo) : (B)*(X) computation flag                           *
C*    STRON (slo) : run starting flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SETLRM,SSBNDS,SSICHK                                             *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'

C==== arguments ========================================================
C
      INTEGER          ISTOR(*),NNEIG 
      REAL             EIG(*),RSTOR(*),SIGMA
      LOGICAL          AGAIN,ENDON,RBXON,STRON
C
C==== local variables ==================================================
C
      REAL             DIFF
      INTEGER          RSINT,SSLOG,TIME
C
C**** executable statements ********************************************
C
      DIFF  = SIGMA
      TIME  = ISTOR(ITIME)
      RSINT = ISTOR(IRSINT)
      SSLOG = ISTOR(ISSLOG)
C
      IF ( NNEIG .LT. 0 ) THEN
C
C....... NNEIG can not be negative .....................................
C
         CALL SETLRM (24,ISTOR(LRERR))
C
      ELSE
C
C....... set bounds and number of required eigenpairs ..................
C
         CALL SSBNDS (ISTOR(LRWRN) ,ISTOR(N)     ,NNEIG        ,
     &                ISTOR(NDEIG) ,ISTOR(NFARL) ,ISTOR(NFARR) ,
     &                ISTOR(NNSPNT),ISTOR(NNTRTL),ISTOR(NNTRTR),
     &                ISTOR(NREIGL),ISTOR(NREIGR),ISTOR(NSIGMA),
     &                ISTOR(NSRLS) ,ISTOR(NSRRS) ,ISTOR(NTEIG) ,
     &                RSTOR(RSINT) ,ISTOR(ISINT) ,ISTOR(INDSI) ,
     &                ISTOR(NSINT) ,ISTOR(NSIMAX),EIG          ,
     &                RSTOR(EIGL)  ,RSTOR(EIGR)  ,RSTOR(ENDL)  ,
     &                RSTOR(ENDR)  ,RSTOR(SFARL) ,RSTOR(SFARR) ,
     &                SIGMA        ,RSTOR(ORIGIN),RSTOR(THETA0),
     &                RSTOR(THETAL),RSTOR(THETAR),RSTOR(TRUSTL),
     &                RSTOR(TRUSTR),RSTOR(TIME)  ,RSTOR(BIGNUM))
C
C....... check whether the shift applied is acceptable .................
C
         CALL SSICHK (ISTOR(LFILE) ,ISTOR(LPRNT) ,ISTOR(LRERR) ,
     &                ISTOR(LRMDE) ,ISTOR(LRWRN) ,NNEIG        ,
     &                ISTOR(NDEIG) ,ISTOR(NEWSIG),ISTOR(NFARL) ,
     &                ISTOR(NFARR) ,ISTOR(NNSPNT),ISTOR(NNTRTL),
     &                ISTOR(NNTRTR),ISTOR(NREIG) ,ISTOR(NREIGL),
     &                ISTOR(NREIGR),ISTOR(NRUNMX),ISTOR(NSFAIL),
     &                ISTOR(NSIGMA),ISTOR(NSLOG) ,ISTOR(NSRLS) ,
     &                ISTOR(NSRRS) ,ISTOR(NTEIG) ,RSTOR(RSINT) ,
     &                ISTOR(ISINT) ,ISTOR(INDSI) ,ISTOR(NSINT) ,
     &                ISTOR(NSIMAX),EIG          ,RSTOR(EIGL)  ,
     &                RSTOR(EIGR)  ,RSTOR(ENDL)  ,RSTOR(ENDR)  ,
     &                RSTOR(ORIGIN),RSTOR(RADIUS),RSTOR(SFARL) ,
     &                RSTOR(SFARR) ,SIGMA        ,RSTOR(TRUSTL),
     &                RSTOR(TRUSTR),RSTOR(SSLOG) ,RSTOR(TIME)  ,
     &                RSTOR(BIGNUM))
C
      END IF
C
      DIFF = DIFF - SIGMA
C
C.... set flags accordingly ............................................
C
      ENDON = ISTOR(LRMDE).EQ.0
      RBXON = ISTOR(NDEIG).GT.0 .AND. DIFF.EQ.0 .AND.
     &        ISTOR(NTEIG).GT.0 .AND. ISTOR(NTEIG).NE.ISTOR(NBX)
      STRON = ISTOR(LRMDE).NE.0 .AND. DIFF.EQ.0 .AND. .NOT.RBXON
C
      RBXON = ISTOR(LRERR).EQ.0 .AND. ISTOR(LRWRN).EQ.0 .AND. RBXON
      STRON = ISTOR(LRERR).EQ.0 .AND. ISTOR(LRWRN).EQ.0 .AND. STRON
      ENDON = ISTOR(LRERR).NE.0 .OR.  ISTOR(LRWRN).NE.0 .OR.  ENDON
C
      AGAIN = ENDON .OR. RBXON .OR. STRON
C
      RETURN 
C
C**** end of BLZFCT ****************************************************
C
      END
