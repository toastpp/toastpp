      SUBROUTINE BLZRBX (ISTOR,RSTOR,LFLAG,NVOPU,X,U,V,
     &                   AGAIN,ENDON,RBXON,STRON)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    BLZRBX deals with the product (B)*(X)                            *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    ISTOR (aib) : array for integer variables                        *
C*    RSTOR (aib) : array for real variables                           *
C*    LFLAG (sib) : reverse communication flag                         *
C*    NVOPU (sib) : number of vectors for reverse communication        *
C*    X     (ari) : eigenvector approximations                         *
C*    U     (arb) : array for reverse communication, U(LNI,NVB)        *
C*    V     (arb) : array for reverse communication, V(LNI,NVB)        *
C*    AGAIN (slo) : loop control flag                                  *
C*    ENDON (slo) : finalization flag                                  *
C*    RBXON (slo) : (B)*(X) computation flag                           *
C*    STRON (slo) : run starting flag                                  *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    LZNRMX                                                           *
C*                                                                     *
C***********************************************************************
C
C==== parameters =======================================================
C
      INCLUDE          'addresses.h'
C
C==== arguments ========================================================
C
      INTEGER          ISTOR(*),LFLAG,NVOPU
      REAL             RSTOR(*),U(*),V(*),X(*)
      LOGICAL          AGAIN,ENDON,RBXON,STRON
C
C==== local variables ==================================================
C
      INTEGER          BX,TIME
      LOGICAL          PURFY,SLICE
C
C**** executable statements ********************************************
C
      BX    = ISTOR(IBX)
      TIME  = ISTOR(ITIME)
      SLICE = ISTOR(LOPTS+1).GT.0
      PURFY = ISTOR(LOPTS+2).GT.0
C
C.... normalize/select vectors for multiplication ......................
C
      CALL LZNRMX (ISTOR(NVB)   ,ISTOR(LNI)   ,ISTOR(NI)    ,
     &             ISTOR(FHNDL) ,ISTOR(LCOMM) ,ISTOR(LRERR) ,
     &             ISTOR(LRWRN) ,ISTOR(NMOPB) ,ISTOR(NBX)   ,
     &             ISTOR(NTEIG) ,NVOPU        ,ISTOR(NXMAX) ,
     &             ISTOR(NBXMAX),RSTOR(BX)    ,X            ,
     &             U            ,V            ,RSTOR(TIME)  ,
     &             PURFY        )
C
      IF      ( NVOPU .GT. 0 ) THEN
              LFLAG = 2
      ELSE IF ( ISTOR(NEWSIG).EQ.0 ) THEN
              ENDON = .TRUE.
      ELSE IF ( SLICE ) THEN
              LFLAG = 3
      END IF
C
C.... set flags accordingly .. .........................................
C
      STRON = ISTOR(NDEIG).NE.0 .AND. ISTOR(NTEIG).EQ.ISTOR(NBX) 
      RBXON = ISTOR(LRERR).EQ.0 .AND. ISTOR(NTEIG).NE.ISTOR(NBX)
      ENDON = ENDON .OR. (.NOT.RBXON .AND. .NOT.SLICE .AND. .NOT.STRON)
C
      STRON = ISTOR(LRERR).EQ.0 .AND. ISTOR(LRWRN).EQ.0 .AND. STRON
      RBXON = ISTOR(LRERR).EQ.0 .AND. ISTOR(LRWRN).EQ.0 .AND. RBXON
      ENDON = ISTOR(LRERR).NE.0 .OR.  ISTOR(LRWRN).NE.0 .OR.  ENDON
C
      AGAIN = ENDON .OR. STRON
C
      RETURN 
C
C**** end of BLZRBX ****************************************************
C
      END
