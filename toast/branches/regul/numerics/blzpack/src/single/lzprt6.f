      SUBROUTINE LZPRT6 (LFILE,LPRNT,LRERR,LRWRN,NMOPA,NMOPB,
     &                   NPORTH,NSORTH,NRUN,NSIGMA,NSLOG,
     &                   NFEIG,SSLOG,TIME,SLICE)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT6 prints statistics, time table and exit messages           *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    LRERR  (sii) : code for error messages                           *
C*    LRWRN  (sii) : code for warning messages                         *
C*    NMOPA  (sii) : number of op(A)*vector performed                  *
C*    NMOPB  (sii) : number of op(B)*vector performed                  *
C*    NPORTH (sii) : number of partial reorthogonalizations performed  *
C*    NSORTH (sii) : number of selective orthogonalizations performed  *
C*    NRUN   (sii) : number of runs                                    *
C*    NSIGMA (sii) : number of origin translations                     *
C*    NSLOG  (sii) : number of subintervals recorded in SSLOG          *
C*    NFEIG  (sii) : number of computed eigenpairs                     *
C*    SSLOG  (arb) : spectrum slicing history                          *
C*    TIME   (ari) : time table                                        *
C*    SLICE  (sli) : spectrum slicing flag                             *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    LZERRS,LZSTTS,LZTIME,LZWRNS,SIBTST,SSHIST                        *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LFILE,LPRNT,LRERR,LRWRN,NFEIG,NMOPA,NMOPB,
     &                 NPORTH,NRUN,NSIGMA,NSLOG,NSORTH
      REAL             SSLOG(8,*),TIME(*)
      LOGICAL          SLICE
C
C==== subprogram =======================================================
C
      LOGICAL          SIBTST
C
C**** executable statements ********************************************
C
      IF ( SIBTST(3,LPRNT) ) THEN 
         WRITE (LFILE,1000)
         IF ( LRERR.EQ.0 ) THEN
            IF ( SLICE ) CALL SSHIST (LFILE,NSLOG,SSLOG)
            CALL LZSTTS (LFILE,NMOPA,NMOPB,NPORTH,NRUN,
     &                   NSIGMA,NSORTH,NFEIG)
            CALL LZTIME (LFILE,TIME)
         END IF
      END IF
C
      IF ( SIBTST(1,LPRNT) ) THEN
         CALL LZWRNS (LFILE,LRWRN)
         CALL LZERRS (LFILE,LRERR)
         WRITE (LFILE,1001)
      END IF
C
      RETURN
C
 1000 FORMAT (/,71('*'))
 1001 FORMAT (/,'BLZPACK exit ',58('*'))
C
C**** end of LZPRT6 ****************************************************
C
      END
