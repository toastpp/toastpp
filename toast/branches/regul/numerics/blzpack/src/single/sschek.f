      SUBROUTINE SSCHEK (LRERR,LRMDE,NESIKL,NESIKR,NNSPNT,NNTRTL,NNTRTR,
     &                   NREIG,NREIGL,NREIGR,NSSIKL,NSSIKR,NTEIG,EIG,
     &                   EIGL,EIGR,ORIGIN,SIGMA,INDSI,NSINT,RSINT,
     &                   ISINT,TRUSTL,TRUSTR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSCHEK checks convergence in all subintervals                    *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LRERR  (sio) : code for error messages                           *
C*    LRMDE  (sio) : run mode                                          *
C*    NESIKL (sio) : number of eigenvalues in the left  subinterval    *
C*    NESIKR (sio) : number of eigenvalues in the right subinterval    *
C*    NNSPNT (sii) : number of eigenvalues less than ORIGIN            *
C*    NNTRTL (sib) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sib) : number of eigenvalues less than TRUSTR            *
C*    NREIG  (sii) : number of required eigenvalues                    *
C*    NREIGL (sii) : number of required eigenvalues less    than EIGR  *
C*    NREIGR (sii) : number of required eigenvalues greater than EIGL  *
C*    NSSIKL (sio) : number of eigenvalues in the left  subinterval    *
C*    NSSIKR (sio) : number of eigenvalues in the right subinterval    *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    EIGL   (sri) : inferior bound for eigenvalues                    *
C*    EIGR   (sri) : superior bound for eigenvalues                    *
C*    ORIGIN (sri) : starting point (first SIGMA)                      *
C*    SIGMA  (sri) : origin translation                                *
C*    INDSI  (sio) : index of the subinterval to be closed             *
C*    NSINT  (sib) : number of subintervals                            *
C*    RSINT  (arb) : lower and upper limits of each subinterval        *
C*    ISINT  (ari) : inertia of the lower and upper limits             *
C*    TRUSTL (srb) : inferior trust bound                              *
C*    TRUSTR (srb) : superior trust bound                              *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    SSCLSD,SSSLCT,SSSTAT                                             *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          INDSI,ISINT(2,*),LRERR,LRMDE,NESIKL,NESIKR,
     &                 NNSPNT,NNTRTL,NNTRTR,NREIG,NREIGL,NREIGR,
     &                 NSINT,NSSIKL,NSSIKR,NTEIG
      REAL             EIG(*),EIGL,EIGR,ORIGIN,SIGMA,
     &                 RSINT(6,*),TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      INTEGER          KEPTLS,KEPTRS,NSLXI,NSRXI
C
C**** executable statements ********************************************
C
C.... check all subintervals ...........................................
C
      CALL SSSTAT (KEPTLS,KEPTRS,LRERR,NESIKL,NESIKR,NNTRTL,
     &             NNTRTR,NSSIKL,NSSIKR,NTEIG,NSINT,RSINT,
     &             ISINT,EIG,ORIGIN,SIGMA,TRUSTL,TRUSTR)
C
C.... check whether all required eigenvalues were found ................
C
      CALL SSCLSD (LRMDE,NNSPNT,NNTRTL,NNTRTR,NREIG,NREIGL,
     &             NREIGR,NSINT,NSLXI,NSRXI,NTEIG,EIG,
     &             EIGL,EIGR,ORIGIN,TRUSTL,TRUSTR)
C
C.... select a subinterval .............................................
C
      CALL SSSLCT (KEPTLS,KEPTRS,LRMDE,NREIGL,NREIGR,NSLXI,
     &             NSRXI,NTEIG,INDSI,NSINT,RSINT,ISINT,
     &             EIG,ORIGIN,SIGMA)   
C
      RETURN 
C
C**** end of SSCHEK ****************************************************
C
      END
