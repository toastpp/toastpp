      SUBROUTINE LZMMRY (JLMAX,JLSET,JTMAX,NVB,NVSET,LTAU,NI,N,LCOMM,
     &                   NPE,LRERR,NQMAX,NBXMAX,NRUNMX,NSIMAX,LISTOR,
     &                   LRSTOR,IBUSY,RBUSY,GNRZD)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZMMRY checks the workspace available                            *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    JLMAX  (sio) : maximum number of steps                           *
C*    JLSET  (sii) : maximum number of steps set by the user           *
C*    JTMAX  (sio) : maximum dimension of the block tridiagonal matrix *
C*    NVB    (sib) : number of vectors in a block                      *
C*    NVSET  (sii) : number of vectors in a block set by the user      *
C*    LTAU   (sii) : leading dimension of (TAU)                        *
C*    NI     (sii) : dimension of the vectors in (U), (V) and (X)      *
C*    N      (sii) : dimension of the eigenvalue problem               *
C*    LCOMM  (sii) : communicator for the parallel version             *
C*    NPE    (sii) : number of processes                               *
C*    LRERR  (sio) : code for error messages                           *
C*    NQMAX  (sio) : maximum number of vectors in (BASIS)              *
C*    NBXMAX (sio) : maximum number of vectors in (BX)                 *
C*    NRUNMX (sii) : maximum number of runs                            *
C*    NSIMAX (sii) : maximum number of subintervals                    *
C*    LISTOR (sib) : workspace available in ISTOR                      *
C*    LRSTOR (sib) : workspace available in RSTOR                      *
C*    IBUSY  (sib) : number of positions already used in ISTOR         *
C*    RBUSY  (sib) : number of positions already used in RSTOR         *
C*    GNRZD  (sli) : problem type flag                                 *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    PIIRED,SETLRM                                                    *
C*                                                                     *
C*  - Intrinsic Functions:                                             *
C*                                                                     *
C*    MAX,MIN,MOD                                                      *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   IBUSY,JLMAX,JLSET,JTMAX,LCOMM,LISTOR,LRERR,LRSTOR,LTAU,
     &          N,NBXMAX,NI,NPE,NQMAX,NRUNMX,NSIMAX,NVB,NVSET,RBUSY
      LOGICAL   GNRZD
C
C==== local variables ==================================================
C
      INTEGER   FACTOR,INFO,MEMPTY,MWORK,TEMP
C
C==== intrinsic functions ==============================================
C
      INTRINSIC MAX,MIN,MOD
C
C**** executable statements ********************************************
C
      IF ( GNRZD ) THEN
         FACTOR = 2
      ELSE
         FACTOR = 1
      END IF
C
      NQMAX  = 0
      NBXMAX = 0
      MEMPTY = LRSTOR
C
C.... minimum workspace ................................................
C
      CALL PIIRED ('MIN',1,LISTOR,TEMP,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
      IF ( TEMP .LT. 0 ) LISTOR = -1
      CALL PIIRED ('MIN',1,LRSTOR,TEMP,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
      IF ( TEMP .LT. 0 ) LRSTOR = -1
C
C.... set limits .......................................................
C
      IF ( NVSET .GT. 0 ) NVB   = NVSET
      IF ( JLSET .GT. 0 ) JTMAX = JLSET*NVB
      IF ( JTMAX .GT. N ) JTMAX = N
C
C.... check the space available in the integer work array ..............
C
      MWORK = IBUSY    
      MWORK = MWORK + NSIMAX*2
      MWORK = MWORK + JTMAX*12
C
      IF ( MWORK.GT.LISTOR .AND. LISTOR.GT.0 ) THEN
         CALL SETLRM (22,LRERR)
         RETURN
      END IF
C
      IBUSY = MWORK
C
C.... check the space available in the real work array .................
C
   10 CONTINUE
C
      IF ( MOD(JTMAX,NVB) .NE. 0 ) THEN
         JLMAX = JTMAX/NVB + 1
      ELSE
         JLMAX = JTMAX/NVB
      END IF
C
      JTMAX = JLMAX*NVB
C
      TEMP = JLMAX
      CALL PIIRED ('MIN',1,TEMP,JLMAX,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
      TEMP = JTMAX
      CALL PIIRED ('MIN',1,TEMP,JTMAX,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
      MWORK = RBUSY + 15
      MWORK = MWORK + 6*NSIMAX
      MWORK = MWORK + 8*NRUNMX
      MWORK = MWORK + JTMAX*2      
      MWORK = MWORK + JTMAX*(NVB+1)
      MWORK = MWORK + NVB*NVB
      MWORK = MWORK + NVB*NVB
      MWORK = MWORK + NVB*NVB
      MWORK = MWORK + JLMAX*2
      MWORK = MWORK + JLMAX*2
      MWORK = MWORK + JLMAX*2
      MWORK = MWORK + LTAU*2
      MWORK = MWORK + NI*NVB*4
      MWORK = MWORK + JTMAX*2
      MWORK = MWORK + JTMAX*JTMAX
      MWORK = MWORK + NI*NVB*FACTOR
      MWORK = MWORK + NI*(FACTOR-1)
      MWORK = MWORK + (2+JTMAX+MAX(NVB+2,18))*JTMAX
C
      IF ( MWORK .GT. LRSTOR ) THEN
C
C....... reduce the block size .........................................
C
         IF      ( LRSTOR.LE.0 ) THEN
                 RBUSY = MWORK
                 GO TO 30
         ELSE IF ( NVB.EQ.1 .OR. NVSET.GT.0 ) THEN
                 CALL SETLRM (23,LRERR)
                 GO TO 20
         END IF
C
         NVB = NVB - 1
C
      ELSE  
C
C....... space for Lanczos vectors .....................................
C
         MEMPTY = MAX(0,MEMPTY-MWORK)
C
         NQMAX  = MIN(JLMAX,MEMPTY/(NI*NVB*FACTOR))
C
C....... space for Ritz vectors ........................................
C
         MEMPTY = MAX(0,MEMPTY-NQMAX*NI*NVB*FACTOR)
C
         NBXMAX = MIN(N*(FACTOR-1),LTAU,MEMPTY/NI)
C
      END IF
C
      IF ( MEMPTY .EQ. LRSTOR ) GO TO 10
C
      IF ( NBXMAX .EQ. LTAU-1 ) NBXMAX = LTAU
      IF ( NQMAX .EQ. JLMAX-1 ) NQMAX = JLMAX
C
   20 CONTINUE
C
C.... cross check error ................................................
C
      CALL PIIRED ('MAX',1,LRERR,TEMP,LCOMM,INFO)
      IF ( TEMP .NE. 0 ) CALL SETLRM (23,LRERR)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
C.... find the minima ..................................................
C
      TEMP = NVB
      CALL PIIRED ('MIN',1,TEMP,NVB,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
      TEMP = NQMAX
      CALL PIIRED ('MIN',1,TEMP,NQMAX,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
      TEMP = NBXMAX
      CALL PIIRED ('MIN',1,TEMP,NBXMAX,LCOMM,INFO)
      IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
C
   30 CONTINUE
C
C.... the current code does not support parallel IO operations .........
C
      IF ( NPE .GT. 1 ) THEN
         IF      ( LRSTOR.LE.0 ) THEN
                 RBUSY = RBUSY + NI*NVB*(JLMAX-1)*FACTOR
                 RBUSY = RBUSY + NI*(LTAU-1)*(FACTOR-1)
         ELSE IF ( JLMAX.GT.NQMAX .OR. (GNRZD.AND.LTAU.GT.NBXMAX) ) THEN
                 CALL SETLRM (23,LRERR)
         END IF
      END IF
C
      RETURN 
C
C**** end of LZMMRY ****************************************************
C
      END
