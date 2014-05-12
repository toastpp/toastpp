      SUBROUTINE LZPRT3 (JL,JT,LFILE,LPRNT,IDXETA,IDXTAU,
     &                   ABSETA,ABSTAU,THRES,RITZ)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT3 prints the Ritz values at every step                      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL     (sii) : number of steps                                   *
C*    JT     (sii) : dimension of the block tridiagonal matrix         *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    IDXETA (sii) : index of the maximum entry in (ETA)               *
C*    IDXTAU (sii) : index of the maximum entry in (TAU)               *
C*    ABSETA (sri) : maximum entry in (ETA)                            *
C*    ABSTAU (sri) : maximum entry in (TAU)                            *
C*    THRES  (sri) : threshold for convergence                         *
C*    RITZ   (ari) : Ritz values and estimated residuals               *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          IDXETA,IDXTAU,JL,JT,LFILE,LPRNT 
      DOUBLE PRECISION ABSETA,ABSTAU,RITZ(JT,2),
     &                 THRES
C
C==== local variables ==================================================
C
      INTEGER          I
C
C==== subprograms ======================================================
C
      LOGICAL          SIBTST
C
C**** executable statements ********************************************
C
      IF ( SIBTST(6,LPRNT) ) THEN
         WRITE (LFILE,1000) JL,JT
         IF ( IDXTAU.GT.0 ) WRITE (LFILE,1001) IDXTAU,ABSTAU
         IF ( IDXETA.GT.0 ) WRITE (LFILE,1002) IDXETA,ABSETA
         WRITE (LFILE,1003) THRES,(RITZ(I,1),RITZ(I,2),I=1,JT)
      END IF 
C
      RETURN
C
 1000 FORMAT (20('.'),'  step:',I4,',  basis size:',I4,2X,20('.'))
 1001 FORMAT ('selective orthogonalization on:',I4,
     &        1P,SP,' (',E9.2,') ')
 1002 FORMAT ('partial reorthogonalization on:',I4,
     &        1P,SP,' (',E9.2,') ')
 1003 FORMAT ('eigenvalue approximations and estimated residuals ',
     &        1P,'(tolerance:',E9.2,')',/,
     &        (1P,SP,3(E11.4,' (',E9.2,') ')))
C
C**** end of LZPRT3 ****************************************************
C
      END
