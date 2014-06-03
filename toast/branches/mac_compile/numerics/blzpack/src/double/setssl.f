      SUBROUTINE SETSSL (JL,JT,NNEIG,NCEIG,ENDL,ENDR,
     &                   SIGMA,NSLOG,SSLOG,TIME)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SETSSL sets the entries of SSLOG                                 *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JL    (sii) : number of steps                                    *
C*    JT    (sii) : dimension of the block tridiagonal matrix          *
C*    NNEIG (sii) : number of eigenvalues less than SIGMA              *
C*    NCEIG (sii) : number of computed eigenpairs                      *
C*    ENDL  (sri) : inferior bound for Ritz values                     *
C*    ENDR  (sri) : superior bound for Ritz values                     *
C*    SIGMA (sri) : origin translation                                 *
C*    NSLOG (sii) : number of subintervals recorded in SSLOG           *
C*    SSLOG (aro) : spectrum slicing history                           *
C*    TIME  (ari) : time table                                         *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SITIME                                                           *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    DBLE                                                             *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C  
      INTEGER          JL,JT,NCEIG,NNEIG,NSLOG
      DOUBLE PRECISION ENDL,ENDR,SIGMA,SSLOG(8,*),TIME(*)
C
C==== subprogram =======================================================
C
      DOUBLE PRECISION SITIME
C
C==== intrinsic function ===============================================
C
      INTRINSIC        DBLE
C
C**** executable statements ********************************************
C
      SSLOG(1,NSLOG) = SIGMA
      SSLOG(2,NSLOG) = ENDL
      SSLOG(3,NSLOG) = ENDR
      SSLOG(4,NSLOG) = SITIME(TIME(12))
      SSLOG(5,NSLOG) = DBLE(NNEIG)
      SSLOG(6,NSLOG) = DBLE(NCEIG)
      SSLOG(7,NSLOG) = DBLE(JL)
      SSLOG(8,NSLOG) = DBLE(JT)
C
      RETURN
C
C**** end of SETSSL ****************************************************
C
      END
