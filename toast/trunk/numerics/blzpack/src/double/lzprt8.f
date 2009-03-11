      SUBROUTINE LZPRT8 (LFILE,NSIGMA,SIGMA,SGNEW)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT8 prints warning about a discarded factorization            *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE  (sii) : file unit for output                              *
C*    NSIGMA (sii) : number of origin translations                     *
C*    SIGMA  (srb) : current origin translation                        *
C*    SGNEW  (srb) : new origin translation                            *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C  
      INTEGER          LFILE,NSIGMA
      DOUBLE PRECISION SGNEW,SIGMA
C
C**** executable statements ********************************************
C
      WRITE (LFILE,1000) NSIGMA,SIGMA,SGNEW
C
      RETURN 
C
 1000 FORMAT (/,'* Warning: factorization',I4,' will be discarded',
     &        /,11X,'SIGMA =',1P,E12.4,' was applied too far away',
     &        /,11X,'SIGMA =',1P,E12.4,' will be used instead')
C
C**** end of LZPRT8 ****************************************************
C
      END
