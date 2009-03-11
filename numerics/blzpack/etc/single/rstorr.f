      REAL             FUNCTION RSTORR (RSTOR,VNAME)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    RSTORR is an auxiliar function for retrieving internal variables *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    RSTOR (ari) : array for real variables                           *
C*    VNAME (sci) : variable name (case sensitive)                     *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
C==== arguments ========================================================
C
      REAL             RSTOR(*)
      CHARACTER        VNAME*(*)
C
C==== local variable ===================================================
C
      INTEGER          RSTART
C
C**** executable statements ********************************************
C
      RSTART = RINIT - 1
C
      IF      ( VNAME .EQ. 'BIGNUM' ) THEN
              RSTORR = RSTOR(RSTART+BIGNUM)
      ELSE IF ( VNAME .EQ. 'EIGL'   ) THEN
              RSTORR = RSTOR(RSTART+EIGL)
      ELSE IF ( VNAME .EQ. 'EIGR'   ) THEN
              RSTORR = RSTOR(RSTART+EIGR)
      ELSE IF ( VNAME .EQ. 'ENDL'   ) THEN
              RSTORR = RSTOR(RSTART+ENDL)
      ELSE IF ( VNAME .EQ. 'ENDR'   ) THEN
              RSTORR = RSTOR(RSTART+ENDR)
      ELSE IF ( VNAME .EQ. 'EPS'    ) THEN
              RSTORR = RSTOR(RSTART+EPS)
      ELSE IF ( VNAME .EQ. 'EPS1'   ) THEN
              RSTORR = RSTOR(RSTART+EPS1)
      ELSE IF ( VNAME .EQ. 'GRNRM'  ) THEN
              RSTORR = RSTOR(RSTART+GRNRM)
      ELSE IF ( VNAME .EQ. 'ORIGIN' ) THEN
              RSTORR = RSTOR(RSTART+ORIGIN)
      ELSE IF ( VNAME .EQ. 'RADIUS' ) THEN
              RSTORR = RSTOR(RSTART+RADIUS)
      ELSE IF ( VNAME .EQ. 'REPS'   ) THEN
              RSTORR = RSTOR(RSTART+REPS)
      ELSE IF ( VNAME .EQ. 'SFARL'  ) THEN
              RSTORR = RSTOR(RSTART+SFARL)
      ELSE IF ( VNAME .EQ. 'SFARR'  ) THEN
              RSTORR = RSTOR(RSTART+SFARR)
      ELSE IF ( VNAME .EQ. 'THETA0' ) THEN
              RSTORR = RSTOR(RSTART+THETA0)
      ELSE IF ( VNAME .EQ. 'THETAL' ) THEN
              RSTORR = RSTOR(RSTART+THETAL)
      ELSE IF ( VNAME .EQ. 'THETAR' ) THEN
              RSTORR = RSTOR(RSTART+THETAR)
      ELSE IF ( VNAME .EQ. 'THRSH'  ) THEN
              RSTORR = RSTOR(RSTART+THRSH)
      ELSE IF ( VNAME .EQ. 'TRUSTL' ) THEN
              RSTORR = RSTOR(RSTART+TRUSTL)
      ELSE IF ( VNAME .EQ. 'TRUSTR' ) THEN
              RSTORR = RSTOR(RSTART+TRUSTR)
      ELSE
              WRITE (*,'(2A)') '* RSTORR, variable not defined: ',VNAME
              RSTORR = 0
      END IF
C
      RETURN
C
C**** end of RSTORR ****************************************************
C
      END
