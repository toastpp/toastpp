      SUBROUTINE SSBACK (LFILE,LPRNT,LRWRN,NSIMAX,NNEIG,NNSPNT,NNTRTL,
     &                   NNTRTR,NREIGL,NREIGR,NSFAIL,NSIGMA,NSINT,
     &                   ORIGIN,SIGMA,TRUSTL,TRUSTR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSBACK checks whether SIGMA was applied too far                  *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE  (sii) : file unit for output                              *
C*    LPRNT  (sii) : level of printing                                 *
C*    LRWRN  (sio) : code for warning messages                         *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    NNEIG  (sii) : number of eigenvalues less than SIGMA             *
C*    NNSPNT (sii) : number of eigenvalues less than ORIGIN            *
C*    NNTRTL (sii) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sii) : number of eigenvalues less than TRUSTR            *
C*    NREIGL (sii) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sii) : number of required eigenvalues greater than EIGL  *
C*    NSFAIL (sib) : number of factorizations failed                   *
C*    NSIGMA (sii) : number of origin translations                     *
C*    NSINT  (sib) : number of subintervals                            *
C*    ORIGIN (sri) : starting-point (first SIGMA)                      *
C*    SIGMA  (srb) : origin translation                                *
C*    TRUSTL (sri) : inferior trust bound                              *
C*    TRUSTR (sri) : superior trust bound                              *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SETLRM,SIBTST                                                    *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C  
      INTEGER          LFILE,LPRNT,LRWRN,NSIMAX,NNEIG,NNSPNT,NNTRTL,
     &                 NNTRTR,NREIGL,NREIGR,NSFAIL,NSIGMA,NSINT
      DOUBLE PRECISION ORIGIN,SIGMA,TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      DOUBLE PRECISION SGNEW
      LOGICAL          MVBACK
C
C==== subprogram =======================================================
C
      LOGICAL          SIBTST
C
C**** executable statements ********************************************
C
      MVBACK = .FALSE.
C
C.... check the number of required eigenpairs ..........................
C
      IF      ((SIGMA.LT.ORIGIN).AND.(NNSPNT-NNEIG.GT.NREIGL+30)) THEN
C
C............ move SIGMA back to the right .............................
C
              SGNEW = TRUSTL - 30*(TRUSTL-SIGMA)/(NNTRTL-NNEIG)
              MVBACK = .TRUE.    
C
      ELSE IF ((SIGMA.GT.ORIGIN).AND.(NNEIG-NNSPNT.GT.NREIGR+30)) THEN
C
C............ move SIGMA back to the left ..............................
C
              SGNEW = TRUSTR + 30*(SIGMA-TRUSTR)/(NNEIG-NNTRTR)
              MVBACK = .TRUE.    
C
      END IF
C
C.... discard SIGMA if necessary .......................................
C
      IF ( MVBACK ) THEN
C
         NSINT = NSINT + 1
         NSFAIL = NSFAIL + 1
C
         IF ( NSINT .GT. NSIMAX ) CALL SETLRM (10,LRWRN)
         IF ( SIBTST(4,LPRNT) ) CALL LZPRT8 (LFILE,NSIGMA,SIGMA,SGNEW)
C
         SIGMA = SGNEW
C
      ELSE
C
         NSFAIL = 0
C
      END IF
C
      RETURN 
C
C**** end of SSBACK ****************************************************
C
      END
