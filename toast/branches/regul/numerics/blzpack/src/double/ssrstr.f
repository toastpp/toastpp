      LOGICAL FUNCTION SSRSTR (LRMDE,NESIKL,NESIKR,NEWSIG,NSLS,NSRS,
     &                         NSRLS,NSRRS,NSSIKL,NSSIKR,
     &                         ORIGIN,SIGMA)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSRSTR checks for a restart with SIGMA (origin translation)      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LRMDE  (sii) : run mode                                          *
C*    NESIKL (sii) : number of eigenvalues in the left  subinterval    *
C*    NESIKR (sii) : number of eigenvalues in the right subinterval    *
C*    NEWSIG (sii) : flag for a new starting point                     *
C*    NSLS   (sii) : number of eigenvalues converged < than SIGMA      *
C*    NSRS   (sii) : number of eigenvalues converged > than SIGMA      *
C*    NSRLS  (sib) : number of eigenvalues required  < than SIGMA      *
C*    NSRRS  (sib) : number of eigenvalues required  > than SIGMA      *
C*    NSSIKL (sii) : number of eigenvalues in the left  subinterval    *
C*    NSSIKR (sii) : number of eigenvalues in the right subinterval    *
C*    ORIGIN (sri) : starting-point (first SIGMA)                      *
C*    SIGMA  (sri) : origin translation                                *
C*                                                                     *
C***********************************************************************
C*                                                                     *
C*    condition for SSRSTR = .TRUE. :                                  *
C*                                                                     *
C*    eigenvalues_missing_in_gap < eigenvalues_computed_in_gap * scale *
C*                                                                     *
C*    "scale" : function of a new factorization cost + a new run cost  *
C*    "scale" : set to 2 in the present configuration                  *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LRMDE,NESIKL,NESIKR,NEWSIG,NSLS,NSRS,
     &                 NSRLS,NSRRS,NSSIKL,NSSIKR
      DOUBLE PRECISION ORIGIN,SIGMA
C
C==== local variable ===================================================
C
      LOGICAL          HOLDON,HOLDLS,HOLDRS
C
C**** executable statements ********************************************
C
      HOLDON = .FALSE.
C
      IF      ( LRMDE .EQ. 4 ) THEN
C
C............ splitting a subinterval ..................................
C
              HOLDLS = NESIKL.GT.NSSIKL .AND. (NESIKL-NSSIKL).LT.NSLS
     &                 .AND. NSRLS.GT.NSLS
              HOLDRS = NESIKR.GT.NSSIKR .AND. (NESIKR-NSSIKR).LT.NSRS
     &                 .AND. NSRRS.GT.NSRS
C
              HOLDON = NEWSIG.LE.2 .AND. (HOLDLS.OR.HOLDRS)
C
      ELSE IF ( (SIGMA.LT.ORIGIN) .AND. (NESIKR.GT.NSSIKR) ) THEN
C
C............ SIGMA is a  lower bound for the subinterval ..............
C
              HOLDON = NEWSIG.LE.2 .AND. (NESIKR-NSSIKR).LT.NSRS*2
C
      ELSE IF ( (SIGMA.GT.ORIGIN) .AND. (NESIKL.GT.NSSIKL) ) THEN
C
C............ SIGMA is an upper bound for the subinterval ..............
C
              HOLDON = NEWSIG.LE.2 .AND. (NESIKL-NSSIKL).LT.NSLS*2
C
      END IF
C
C.... condition for keeping SIGMA ......................................
C
      SSRSTR = HOLDON .OR. NEWSIG.EQ.1
C
      RETURN 
C
C**** end of SSRSTR ****************************************************
C
      END
