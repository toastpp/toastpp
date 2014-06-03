      SUBROUTINE LZSTTS (LFILE,NMOPA,NMOPB,NPORTH,
     &                   NRUN,NSIGMA,NSORTH,NTEIG) 
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZSTTS prints statistics for the Lanczos algorithm               *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE  (sii) : file unit for output                              *
C*    NMOPA  (sii) : number of op(A)*vector performed                  *
C*    NMOPB  (sii) : number of op(B)*vector performed                  *
C*    NPORTH (sii) : number of partial reorthogonalizations performed  *
C*    NRUN   (sii) : number of runs                                    *
C*    NSIGMA (sii) : number of origin translations                     *
C*    NSORTH (sii) : number of selective orthogonalizations performed  *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER LFILE,NMOPA,NMOPB,NPORTH,NRUN,NSIGMA,NSORTH,NTEIG
C
C**** executable statements ********************************************
C
      WRITE  (LFILE,1000) NRUN,NSIGMA,NTEIG,NPORTH,NSORTH,NMOPA,NMOPB
C
      RETURN 
C
 1000 FORMAT (/,'statistics',/,44('-'),
     &        /,'number of runs .................. =',I09,  
     &        /,'number of origin translations ... =',I09,  
     &        /,'number of converged eigenpairs .. =',I09,  
     &        /,'number of products for p.o. ..... =',I09,  
     &        /,'number of products for s.o. ..... =',I09,  
     &        /,'number of op(A)*vector .......... =',I09,  
     &        /,'number of op(B)*vector .......... =',I09,  
     &        /,44('-'))  
C
C**** end of LZSTTS ****************************************************
C
      END
