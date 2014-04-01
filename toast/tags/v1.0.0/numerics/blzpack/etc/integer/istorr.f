      INTEGER FUNCTION ISTORR (ISTOR,VNAME)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    ISTORR is an auxiliar function for retrieving internal variables *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aii) : array for integer variables                        *
C*    VNAME (sci) : variable name (case sensitive)                     *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE   'addresses.h'
C
C==== arguments ========================================================
C
      INTEGER   ISTOR(*)
      CHARACTER VNAME*(*)
C
C==== local variable ===================================================
C
      INTEGER   ISTART,RSTART
C
C**** executable statements ********************************************
C
      ISTART = IINIT - 1
      RSTART = RINIT - 1
C
      IF      ( VNAME .EQ. 'INDSI'  ) THEN
              ISTORR = ISTOR(ISTART+INDSI)
      ELSE IF ( VNAME .EQ. 'JL'     ) THEN
              ISTORR = ISTOR(ISTART+JL)
      ELSE IF ( VNAME .EQ. 'JLMAX'  ) THEN
              ISTORR = ISTOR(ISTART+JLMAX)
      ELSE IF ( VNAME .EQ. 'JT'     ) THEN
              ISTORR = ISTOR(ISTART+JT)
      ELSE IF ( VNAME .EQ. 'JTMAX'  ) THEN
              ISTORR = ISTOR(ISTART+JTMAX)
      ELSE IF ( VNAME .EQ. 'JTMIN'  ) THEN
              ISTORR = ISTOR(ISTART+JTMIN)
      ELSE IF ( VNAME .EQ. 'LCOMM'  ) THEN
              ISTORR = ISTOR(ISTART+LCOMM)
      ELSE IF ( VNAME .EQ. 'LEIG'   ) THEN
              ISTORR = ISTOR(ISTART+LEIG)
      ELSE IF ( VNAME .EQ. 'LFILE'  ) THEN
              ISTORR = ISTOR(ISTART+LFILE)
      ELSE IF ( VNAME .EQ. 'LNI'    ) THEN
              ISTORR = ISTOR(ISTART+LNI)
      ELSE IF ( VNAME .EQ. 'LPRNT'  ) THEN
              ISTORR = ISTOR(ISTART+LPRNT)
      ELSE IF ( VNAME .EQ. 'LRERR'  ) THEN
              ISTORR = ISTOR(ISTART+LRERR)
      ELSE IF ( VNAME .EQ. 'LRMDE'  ) THEN
              ISTORR = ISTOR(ISTART+LRMDE)
      ELSE IF ( VNAME .EQ. 'LRWRN'  ) THEN
              ISTORR = ISTOR(ISTART+LRWRN)
      ELSE IF ( VNAME .EQ. 'LTAU'   ) THEN
              ISTORR = ISTOR(ISTART+LTAU)
      ELSE IF ( VNAME .EQ. 'MYPE'   ) THEN
              ISTORR = ISTOR(ISTART+MYPE)
      ELSE IF ( VNAME .EQ. 'N'      ) THEN
              ISTORR = ISTOR(ISTART+N)
      ELSE IF ( VNAME .EQ. 'NBX'    ) THEN
              ISTORR = ISTOR(ISTART+NBX)
      ELSE IF ( VNAME .EQ. 'NBXMAX' ) THEN
              ISTORR = ISTOR(ISTART+NBXMAX)
      ELSE IF ( VNAME .EQ. 'NDEIG'  ) THEN
              ISTORR = ISTOR(ISTART+NDEIG)
      ELSE IF ( VNAME .EQ. 'NEPIN'  ) THEN
              ISTORR = ISTOR(ISTART+NEPIN)
      ELSE IF ( VNAME .EQ. 'NEWSIG' ) THEN
              ISTORR = ISTOR(ISTART+NEWSIG)
      ELSE IF ( VNAME .EQ. 'NFARL'  ) THEN
              ISTORR = ISTOR(ISTART+NFARL)
      ELSE IF ( VNAME .EQ. 'NFARR'  ) THEN
              ISTORR = ISTOR(ISTART+NFARR)
      ELSE IF ( VNAME .EQ. 'NI'     ) THEN
              ISTORR = ISTOR(ISTART+NI)
      ELSE IF ( VNAME .EQ. 'NMOPA'  ) THEN
              ISTORR = ISTOR(ISTART+NMOPA)
      ELSE IF ( VNAME .EQ. 'NMOPB'  ) THEN
              ISTORR = ISTOR(ISTART+NMOPB)
      ELSE IF ( VNAME .EQ. 'NNSPNT' ) THEN
              ISTORR = ISTOR(ISTART+NNSPNT)
      ELSE IF ( VNAME .EQ. 'NNTRTL' ) THEN
              ISTORR = ISTOR(ISTART+NNTRTL)
      ELSE IF ( VNAME .EQ. 'NNTRTR' ) THEN
              ISTORR = ISTOR(ISTART+NNTRTR)
      ELSE IF ( VNAME .EQ. 'NONEWS' ) THEN
              ISTORR = ISTOR(ISTART+NONEWS)
      ELSE IF ( VNAME .EQ. 'NPE'    ) THEN
              ISTORR = ISTOR(ISTART+NPE)
      ELSE IF ( VNAME .EQ. 'NPORTH' ) THEN
              ISTORR = ISTOR(ISTART+NPORTH)
      ELSE IF ( VNAME .EQ. 'NQMAX'  ) THEN
              ISTORR = ISTOR(ISTART+NQMAX)
      ELSE IF ( VNAME .EQ. 'NREIG'  ) THEN
              ISTORR = ISTOR(ISTART+NREIG)
      ELSE IF ( VNAME .EQ. 'NREIGL' ) THEN
              ISTORR = ISTOR(ISTART+NREIGL)
      ELSE IF ( VNAME .EQ. 'NREIGR' ) THEN
              ISTORR = ISTOR(ISTART+NREIGR)
      ELSE IF ( VNAME .EQ. 'NRITZ'  ) THEN
              ISTORR = ISTOR(ISTART+NRITZ)
      ELSE IF ( VNAME .EQ. 'NRUN'   ) THEN
              ISTORR = ISTOR(ISTART+NRUN)
      ELSE IF ( VNAME .EQ. 'NRUNMX' ) THEN
              ISTORR = ISTOR(ISTART+NRUNMX)
      ELSE IF ( VNAME .EQ. 'NSFAIL' ) THEN
              ISTORR = ISTOR(ISTART+NSFAIL)
      ELSE IF ( VNAME .EQ. 'NSIGMA' ) THEN
              ISTORR = ISTOR(ISTART+NSIGMA)
      ELSE IF ( VNAME .EQ. 'NSIMAX' ) THEN
              ISTORR = ISTOR(ISTART+NSIMAX)
      ELSE IF ( VNAME .EQ. 'NSINT'  ) THEN
              ISTORR = ISTOR(ISTART+NSINT)
      ELSE IF ( VNAME .EQ. 'NSORTH' ) THEN
              ISTORR = ISTOR(ISTART+NSORTH)
      ELSE IF ( VNAME .EQ. 'NSRLS'  ) THEN
              ISTORR = ISTOR(ISTART+NSRLS)
      ELSE IF ( VNAME .EQ. 'NSRRS'  ) THEN
              ISTORR = ISTOR(ISTART+NSRRS)
      ELSE IF ( VNAME .EQ. 'NSVIN'  ) THEN
              ISTORR = ISTOR(ISTART+NSVIN)
      ELSE IF ( VNAME .EQ. 'NTEIG'  ) THEN
              ISTORR = ISTOR(ISTART+NTEIG)
      ELSE IF ( VNAME .EQ. 'NULLDQ' ) THEN
              ISTORR = ISTOR(ISTART+NULLDQ)
      ELSE IF ( VNAME .EQ. 'NULLDR' ) THEN
              ISTORR = ISTOR(ISTART+NULLDR)
      ELSE IF ( VNAME .EQ. 'NVB'    ) THEN
              ISTORR = ISTOR(ISTART+NVB)
      ELSE IF ( VNAME .EQ. 'NWBSY'  ) THEN
              ISTORR = ISTOR(ISTART+NWBSY)
      ELSE IF ( VNAME .EQ. 'NWMAX'  ) THEN
              ISTORR = ISTOR(ISTART+NWMAX)
      ELSE IF ( VNAME .EQ. 'NWMIN'  ) THEN
              ISTORR = ISTOR(ISTART+NWMIN)
      ELSE IF ( VNAME .EQ. 'NXMAX'  ) THEN
              ISTORR = ISTOR(ISTART+NXMAX)
      ELSE IF ( VNAME .EQ. 'ITIME'  ) THEN
              ISTORR = RSTART + ISTOR(ISTART+ITIME)
      ELSE IF ( VNAME .EQ. 'IRSINT' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IRSINT)
      ELSE IF ( VNAME .EQ. 'ISSLOG' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+ISSLOG)
      ELSE IF ( VNAME .EQ. 'IRITZ'  ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IRITZ)
      ELSE IF ( VNAME .EQ. 'ITB'    ) THEN
              ISTORR = RSTART + ISTOR(ISTART+ITB)
      ELSE IF ( VNAME .EQ. 'IALPHA' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IALPHA)
      ELSE IF ( VNAME .EQ. 'IBETAQ' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IBETAQ)
      ELSE IF ( VNAME .EQ. 'IBETAR' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IBETAR)
      ELSE IF ( VNAME .EQ. 'IANORM' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IANORM)
      ELSE IF ( VNAME .EQ. 'IBNORM' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IBNORM)
      ELSE IF ( VNAME .EQ. 'IETA'   ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IETA)
      ELSE IF ( VNAME .EQ. 'ITAU'   ) THEN
              ISTORR = RSTART + ISTOR(ISTART+ITAU)
      ELSE IF ( VNAME .EQ. 'IR'     ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IR)
      ELSE IF ( VNAME .EQ. 'ITHETA' ) THEN
              ISTORR = RSTART + ISTOR(ITHETA)
      ELSE IF ( VNAME .EQ. 'IS'     ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IS)
      ELSE IF ( VNAME .EQ. 'IBASIS' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IBASIS)
      ELSE IF ( VNAME .EQ. 'IBX'    ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IBX)
      ELSE IF ( VNAME .EQ. 'IRWORK' ) THEN
              ISTORR = RSTART + ISTOR(ISTART+IRWORK)
      ELSE IF ( VNAME .EQ. 'IIWORK' ) THEN
              ISTORR = ISTOR(ISTART+IIWORK)
      ELSE IF ( VNAME .EQ. 'INDR'   ) THEN
              ISTORR = ISTOR(ISTART+INDR)
      ELSE IF ( VNAME .EQ. 'LOPTS'  ) THEN
              ISTORR = ISTOR(ISTART+LOPTS)
      ELSE IF ( VNAME .EQ. 'LBLAS'  ) THEN
              ISTORR = ISTOR(ISTART+LBLAS)
      ELSE IF ( VNAME .EQ. 'FHNDL'  ) THEN
              ISTORR = ISTOR(ISTART+FHNDL)
      ELSE IF ( VNAME .EQ. 'ISINT'  ) THEN
              ISTORR = ISTOR(ISTART+ISINT)
      ELSE
              WRITE (*,'(2A)') '* ISTORR, variable not defined: ',VNAME
              ISTORR = 0
      END IF
C
      RETURN
C
C**** end of ISTORR ****************************************************
C
      END
