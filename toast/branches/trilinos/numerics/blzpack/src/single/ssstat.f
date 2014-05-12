      SUBROUTINE SSSTAT (KEPTLS,KEPTRS,LRERR,NESIKL,NESIKR,NNTRTL,
     &                   NNTRTR,NSSIKL,NSSIKR,NTEIG,NSINT,RSINT,
     &                   ISINT,EIG,ORIGIN,SIGMA,TRUSTL,TRUSTR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSSTAT checks the status of each subinterval                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    KEPTLS (sio) : index of the left  subinterval                    *
C*    KEPTRS (sio) : index of the right subinterval                    *
C*    LRERR  (sii) : code for error messages                           *
C*    NESIKL (sio) : number of eigenvalues in the left  subinterval    *
C*    NESIKR (sio) : number of eigenvalues in the right subinterval    *
C*    NNTRTL (sib) : number of eigenvalues less than TRUSTL            *
C*    NNTRTR (sib) : number of eigenvalues less than TRUSTR            *
C*    NSSIKL (sio) : number of eigenvalues in the left  subinterval    *
C*    NSSIKR (sio) : number of eigenvalues in the right subinterval    *
C*    NTEIG  (sii) : number of computed eigenpairs                     *
C*    NSINT  (sib) : number of subintervals                            *
C*    RSINT  (ari) : lower and upper limits of each subinterval        *
C*    ISINT  (ari) : inertia of the lower and upper limits             *
C*    EIG    (ari) : eigenvalue approximations and estimated residuals *
C*    ORIGIN (sri) : starting point (first SIGMA)                      *
C*    SIGMA  (sri) : origin translation                                *
C*    TRUSTL (srb) : inferior trust bound                              *
C*    TRUSTR (srb) : superior trust bound                              *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    NEIGAB,SETLRM,SSMOVB                                             *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          ISINT(2,*),KEPTLS,KEPTRS,LRERR,NESIKL,NESIKR,
     &                 NNTRTL,NNTRTR,NSINT,NSSIKL,NSSIKR,NTEIG
      REAL             EIG(*),ORIGIN,RSINT(6,*),SIGMA,TRUSTL,TRUSTR
C
C==== local variables ==================================================
C
      INTEGER          I,J,NERNG,NRSINT,NSRNG,NTSINV
C
C==== subprogram =======================================================
C
      INTEGER          NEIGAB
C
C**** executable statements ********************************************
C
      NTSINV = 0
      KEPTRS = 0 
      KEPTLS = 0 
      NESIKL = 0
      NESIKR = 0
      NSSIKL = 0
      NSSIKR = 0
C
      NRSINT = NSINT
C
C.... check the status of each subinterval .............................
C
      DO 10 I = 1,NSINT
C
         NERNG = ISINT(2,NTSINV+1) - ISINT(1,NTSINV+1)
         NSRNG = NEIGAB(NTEIG,RSINT(1,NTSINV+1),RSINT(6,NTSINV+1),EIG)
C
         IF      ( NSRNG .GT. NERNG ) THEN
C
C............... something went wrong ..................................
C
                 CALL SETLRM (27,LRERR)
                 RETURN
C
         ELSE IF ( NSRNG .LT. NERNG ) THEN
C
C............... keep  this subinterval ................................
C
                 NTSINV = NTSINV + 1
C
                 IF ( SIGMA .EQ. RSINT(6,NTSINV) ) THEN
                    KEPTLS = NTSINV 
                    NESIKL = NERNG
                    NSSIKL = NSRNG
                 END IF
C
                 IF ( SIGMA .EQ. RSINT(1,NTSINV) ) THEN
                    KEPTRS = NTSINV 
                    NESIKR = NERNG
                    NSSIKR = NSRNG
                 END IF
C
         ELSE 
C
C............... update bounds .........................................
C
                 IF ( RSINT(1,NTSINV+1) .LT. ORIGIN ) THEN
                    TRUSTL = RSINT(1,NTSINV+1)
                    NNTRTL = ISINT(1,NTSINV+1)
                 END IF
                 IF ( RSINT(6,NTSINV+1) .GT. ORIGIN ) THEN
                    TRUSTR = RSINT(6,NTSINV+1)
                    NNTRTR = ISINT(2,NTSINV+1)
                 END IF
C
C............... close this subinterval ................................
C
                 CALL SSMOVB (+1,NTSINV+1,NRSINT-1,ISINT,RSINT)
C
	         NRSINT = NRSINT - 1
C
         END IF
C
   10 CONTINUE
C
      NSINT = NTSINV 
C
C.... update bounds ....................................................
C
      DO 20 I = 1,NSINT
         J = NSINT - I + 1
         IF ( RSINT(1,I) .LT. ORIGIN ) THEN
            TRUSTL = RSINT(6,I)
            NNTRTL = ISINT(2,I)
         END IF
         IF ( RSINT(6,J) .GT. ORIGIN ) THEN
            TRUSTR = RSINT(1,J)
            NNTRTR = ISINT(1,J)
         END IF
   20 CONTINUE
C
      RETURN 
C
C**** end of SSSTAT ****************************************************
C
      END
