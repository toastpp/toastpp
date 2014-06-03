      SUBROUTINE LZPRT2 (LFILE,LPRNT,NNEIG,NDEIG,NSIGMA,RSINT,
     &                   ISINT,NSINT,ENDL,ENDR,SIGMA)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT2 prints the computational subintervals                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    NNEIG  (sii) : number of eigenvalues less than SIGMA             *
C*    NDEIG  (sii) : number of eigenvalues required in the run         *
C*    NSIGMA (sii) : number of origin translations                     *
C*    RSINT  (ari) : lower and upper limits of each subinterval        *
C*    ISINT  (ari) : inertia of the lower and upper limits             *
C*    NSINT  (sii) : number of subintervals                            *
C*    ENDL   (sri) : inferior bound for Ritz values                    *
C*    ENDR   (sri) : superior bound for Ritz values                    *
C*    SIGMA  (sri) : origin translation                                *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C  
      INTEGER          ISINT(2,NSINT),LFILE,LPRNT,NDEIG,
     &                 NNEIG,NSIGMA,NSINT
      DOUBLE PRECISION ENDL,ENDR,RSINT(6,NSINT),SIGMA
C
C==== local variable ===================================================
C
      INTEGER          I
C
C==== subprogram =======================================================
C
      LOGICAL          SIBTST
C
C**** executable statements ********************************************
C
      IF ( SIBTST(4,LPRNT) ) WRITE (LFILE,1000) NSIGMA,SIGMA,NDEIG,
     &                                          NNEIG,ENDL,ENDR
C
      IF ( SIBTST(6,LPRNT) ) WRITE (LFILE,1001) (I,ISINT(1,I),
     &                                             RSINT(1,I),
     &                                             RSINT(2,I),
     &                                             RSINT(3,I),
     &                                             ISINT(2,I),
     &                                             RSINT(6,I),
     &                                             RSINT(5,I),
     &                                             RSINT(4,I),
     &                                             I=1,NSINT)
C
      RETURN 
C
 1000 FORMAT (/,71('='),
     &        /,'factorization:',S,I3,
     &          5X,'SIGMA:',1P,SP,E12.4,
     &          5X,'eigenvalues required:',S,I5,
     &        /,71('='),/
     &        /,'eigenvalues smaller than SIGMA:',S,I8,
     &        /,'lower boundary for eigenvalues:',1P,SP,E12.4,
     &        /,'upper boundary for eigenvalues:',1P,SP,E12.4)
 1001 FORMAT (/,'interval',2X,'bound',3X,'nneig',5X,'xi',10X,
     &          'ritz',8X,'new xi',
     &        /,(S,I5,5X,'lower',S,I8,2X,1P,SP,3E12.4,/,10X,
     &                   'upper',S,I8,2X,1P,SP,3E12.4))
C
C**** end of LZPRT2 ****************************************************
C
      END
