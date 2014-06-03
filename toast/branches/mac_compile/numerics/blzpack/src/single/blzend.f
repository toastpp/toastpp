      SUBROUTINE BLZEND (ISTOR,RSTOR,LFLAG,EIG,X,ABORT)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZEND finishes the Lanczos algorithm run                        *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aib) : array for integer variables                        *
C*    RSTOR (aib) : array for real variables                           *
C*    LFLAG (sio) : reverse communication flag                         *
C*    EIG   (ari) : eigenvalue approximations and estimated residuals  *
C*    X     (ari) : eigenvector approximations                         *
C*    ABORT (sli) : lock control flag                                  *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    EIGCHK,EIGPRT,LZIOOP,LZPRT6,SETLRM                               *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*),LFLAG
      REAL             EIG(*),RSTOR(*),X(*)
      LOGICAL          ABORT
C
C==== local variables ==================================================
C
      INTEGER          IEIG,INDXP,INDXS,IX,NFEIG,NPEIG,RWORK,SSLOG,TIME
      REAL             DUMMY(1)
      LOGICAL          GNRZD,SLICE
C
C**** executable statements ********************************************
C
      TIME  = ISTOR(ITIME)
      SSLOG = ISTOR(ISSLOG)
      RWORK = ISTOR(IRWORK)
C
      DUMMY(1) = 0
      GNRZD = ISTOR(LOPTS  ).GT.0
      SLICE = ISTOR(LOPTS+1).GT.0
C
      IF ( ABORT ) CALL SETLRM (1,ISTOR(LRERR))
C
C.... number of eigenpairs found .......................................
C
      IEIG = 1 + ISTOR(NEPIN)
      IX = 1 + ISTOR(NEPIN)*ISTOR(LNI)
C
      NFEIG = ISTOR(NTEIG) - ISTOR(NEPIN)
C
C.... close temporary files ............................................
C
      CALL LZIOOP (ISTOR(FHNDL),ISTOR(LCOMM),ISTOR(LRERR),
     &             0,1,DUMMY,'BQ','DEL')
      CALL LZIOOP (ISTOR(FHNDL),ISTOR(LCOMM),ISTOR(LRERR),
     &             0,1,DUMMY,'BX','DEL')
      CALL LZIOOP (ISTOR(FHNDL),ISTOR(LCOMM),ISTOR(LRERR),
     &             0,1,DUMMY,'Q ','DEL')
C
C.... sort the eigenpairs ..............................................
C
      IF ( ISTOR(NXMAX) .GT. 0 ) 
     &   CALL EIGCHK (INDXP        ,INDXS        ,ISTOR(LEIG)  ,
     &                ISTOR(LNI)   ,ISTOR(LRWRN) ,ISTOR(NNSPNT),
     &                NPEIG        ,ISTOR(NREIG) ,NFEIG        ,
     &                RSTOR(BIGNUM),EIG(IEIG)    ,X(IX)        ,
     &                RSTOR(EIGL)  ,RSTOR(EIGR)  ,RSTOR(ORIGIN),
     &                RSTOR(TRUSTL),RSTOR(TRUSTR),RSTOR(RWORK) ,
     &                GNRZD        ,SLICE        )
C
C.... print time table and exit messages ...............................
C
      CALL LZPRT6 (ISTOR(LFILE) ,ISTOR(LPRNT) ,ISTOR(LRERR) ,
     &             ISTOR(LRWRN) ,ISTOR(NMOPA) ,ISTOR(NMOPB) ,
     &             ISTOR(NPORTH),ISTOR(NSORTH),ISTOR(NRUN)  ,
     &             ISTOR(NSIGMA),ISTOR(NSLOG) ,NFEIG        ,
     &             RSTOR(SSLOG) ,RSTOR(TIME)  ,SLICE        )
C
C.... print the eigenpairs .............................................
C
      IF ( ISTOR(NXMAX) .GT. 0 ) 
     &   CALL EIGPRT (INDXP        ,INDXS        ,NPEIG        ,
     &                ISTOR(NPE)   ,ISTOR(LFILE) ,ISTOR(LPRNT) ,
     &                ISTOR(LEIG)  ,ISTOR(LNI)   ,ISTOR(NI)    ,
     &                EIG(IEIG)    ,X(IX)        )
C
C.... set output flag ..................................................
C
      LFLAG = -ISTOR(LRERR) 
C
      RETURN 
C
C**** end of BLZEND ****************************************************
C
      END
