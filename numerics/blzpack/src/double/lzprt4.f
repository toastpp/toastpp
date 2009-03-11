      SUBROUTINE LZPRT4 (JL,JT,LFILE,LPRNT,NSLS,NSRS,
     &                   NSRLS,NSRRS,RITZ,SLICE)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZPRT4 prints the Ritz values at the end of the run              *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL    (sii) : number of steps                                    *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    LFILE (sii) : file unit for output                               *
C*    LPRNT (sii) : level of printing                                  *
C*    NSLS  (sii) : number of eigenvalues converged < than SIGMA       *
C*    NSRS  (sii) : number of eigenvalues converged > than SIGMA       *
C*    NSRLS (sii) : number of eigenvalues required  < than SIGMA       *
C*    NSRRS (sii) : number of eigenvalues required  > than SIGMA       *
C*    RITZ  (ari) : Ritz values and estimated residuals                *
C*    SLICE (sli) : spectrum slicing flag                              *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          JL,JT,LFILE,LPRNT,NSLS,NSRS,NSRLS,NSRRS
      DOUBLE PRECISION RITZ(JT,2)
      LOGICAL          SLICE
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
      IF ( SIBTST(4,LPRNT) ) THEN
         IF ( SLICE ) WRITE (LFILE,1000) NSLS,NSRLS,NSRS,NSRRS
         IF ( JL.GT.0 ) WRITE (LFILE,1001) JL,NSLS+NSRS
         IF ( JT.GT.0 ) WRITE (LFILE,1002) (RITZ(I,1),RITZ(I,2),I=1,JT)
      END IF
C
      RETURN
C
 1000 FORMAT (/,'eigenvalues < 0:',I5,4X,
     &          'eigenvalues required < 0:',I5,
     &        /,'eigenvalues > 0:',I5,4X,
     &          'eigenvalues required > 0:',I5)
 1001 FORMAT (/,'steps performed     :',I5,
     &        /,'eigenvalues accepted:',I5)
 1002 FORMAT (/,'eigenvalue approximations and estimated residuals:',
     &        /,(1P,SP,3(E11.4,' (',E9.2,') ')))
C
C**** end of LZPRT4 ****************************************************
C
      END
