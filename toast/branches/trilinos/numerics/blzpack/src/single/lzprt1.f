      SUBROUTINE LZPRT1 (LFILE,LPRNT,LRMDE,NRUN)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT1 prints the run mode                                       *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE (sii) : file unit for output                               *
C*    LPRNT (sii) : level of printing                                  *
C*    LRMDE (sii) : run mode                                           *
C*    NRUN  (sii) : number of runs                                     *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   LFILE,LPRNT,LRMDE,NRUN
C
C==== local variables ==================================================
C
      CHARACTER REASON(9)*36
C
C==== subprogram =======================================================
C
      LOGICAL   SIBTST
C
C==== initialization of variables ======================================
C
      DATA      REASON /'(origin in EIGL=EIGR)               ',
     &                  '(origin in EIGL, moving to EIGR)    ',
     &                  '(origin in EIGR, moving to EIGL)    ',
     &                  '(splitting a subinterval)           ',
     &                  '(restart with the previous origin)  ',
     &                  '(restart with a new origin)         ',
     &                  '(restart with EIGL)                 ',
     &                  '(restart with EIGR)                 ',
     &                  '                                    '/
C
C**** executable statements ********************************************
C
      IF ( SIBTST(4,LPRNT) ) WRITE (LFILE,1000) NRUN,REASON(LRMDE)
C
      RETURN
C
 1000 FORMAT (/,71('-'),/,'run:',I3,2X,A36,/,71('-'))
C
C**** end of LZPRT1 ****************************************************
C
      END
