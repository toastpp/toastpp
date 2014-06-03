      PROGRAM DRVSP3
C     ==============
C
C***********************************************************************
C*                                                                     *
C*    DRVSP3 is a driver model for the standard eigenvalue problem     *
C*                                                                     *
C*           (A)*(x)-eig*(x)=(0)                                       *
C*                                                                     *
C*    where (A) is a symmetric matrix of dimension N defined as        *
C*                                                                     *
C*           (A) = (C)^T*(C)                                           *
C*                                                                     *
C*    where (C) is an M-by-N matrix. Therefore, the pairs eig,(x)      *
C*    lead to the singular value decomposition of (C),                 *
C*                                                                     *
C*           (C) = (U)*(S)*(V)^T                                       *
C*                                                                     *
C*           eig = s^2                                                 *
C*                                                                     *
C*           (x) = (v)                                                 *
C*                                                                     *
C*    Data for DRVSP3 are read from the standard input (see below).    *
C*    The entries of (C) are read from the file MATRXC (coordinate     *
C*    format).                                                         *
C*                                                                     *
C*    Examples: files `drvsp3.dat' and `C1.dat'                        *
C*                                                                     *
C***********************************************************************
C
C.... parameters .......................................................
C
      INTEGER          LEIG
      PARAMETER        (LEIG   =     30)
      INTEGER          LISTOR
      PARAMETER        (LISTOR =  10000)
      INTEGER          LN
      PARAMETER        (LN     =    100)
      INTEGER          LRSTOR
      PARAMETER        (LRSTOR =  10000)
      INTEGER          MAXC
      PARAMETER        (MAXC   =  10000)
      INTEGER          NCUV
      PARAMETER        (NCUV   =      3)
      REAL             ZERO
      PARAMETER        (ZERO   =  0.0E0)
C
C.... work variables ...................................................
C
      INTEGER          I,J,K
      REAL             AJK
C
C.... BLZDRS variables .................................................
C
      INTEGER          ISTOR(LISTOR),LFLAG,NNEIG,NVOPU
      REAL             EIG(LEIG,2),RSTOR(LRSTOR),SIGMA,U(LN,NCUV),
     &                 V(LN,NCUV),X(LN,LEIG),Z(LN,NCUV)
C
C.... matrix (C) .......................................................
C
      INTEGER          IRN(MAXC),JCN(MAXC),M,N,NE
      REAL             C(MAXC)
      CHARACTER        MATRXC*16
C
C=======================================================================
C
C.... read data ........................................................
C
      READ (*,ERR=1,FMT=*) ISTOR( 3)
      READ (*,ERR=1,FMT=*) ISTOR( 5)
      READ (*,ERR=1,FMT=*) ISTOR( 6)
      READ (*,ERR=1,FMT=*) ISTOR( 7)
      READ (*,ERR=1,FMT=*) ISTOR( 8)
      READ (*,ERR=1,FMT=*) ISTOR( 9)
      READ (*,ERR=1,FMT=*) ISTOR(10)
      READ (*,ERR=1,FMT=*) ISTOR(11)
      READ (*,ERR=1,FMT=*) ISTOR(12)
      READ (*,ERR=1,FMT=*) ISTOR(13)
      READ (*,ERR=1,FMT=*) RSTOR( 1)
      READ (*,ERR=1,FMT=*) RSTOR( 2)
      READ (*,ERR=1,FMT=*) RSTOR( 3)
      READ (*,ERR=1,FMT=*) MATRXC
C
C.... read (C) .........................................................
C
      M = 0
      N = 0
      NE = 0
      OPEN (UNIT=10,ERR=2,STATUS='OLD',FILE=MATRXC)
      DO 10 I = 1,LN*LN
         READ (UNIT=10,ERR=2,END=20,FMT=*) J,K,AJK
         NE = NE + 1
         IF ( NE .GT. MAXC ) STOP '* Error: NE > MAXC *'
         M = MAX(J,M)
         N = MAX(K,N)
         IRN(NE) = J
         JCN(NE) = K
         C(NE) = AJK
   10 CONTINUE
   20 CLOSE (UNIT=10,ERR=2)
C
C.... check the dimensions .............................................
C
      IF ( N .GT. LN ) STOP '* Error: N > LN *'
      IF ( ISTOR(5) .GT. NCUV ) STOP '* Error: NVB > NCUV *'
C
      ISTOR( 1) = N
      ISTOR( 2) = LN
      ISTOR( 4) = LEIG
      ISTOR(14) = 0     
      ISTOR(15) = LISTOR    
      RSTOR( 4) = LRSTOR    
C
C.... reverse communication strategy ...................................
C
      LFLAG = 0
C
   30 CONTINUE
C
C     ***********************************************************
      CALL BLZDRS (ISTOR,RSTOR,SIGMA,NNEIG,U,V,LFLAG,NVOPU,EIG,X)
C     ***********************************************************
C
      IF      ( LFLAG .LT. 0 ) THEN
C
C............ early  finalization ......................................
C
              WRITE (*,'(/A)') 'execution finished: abnormal exit' 
C
      ELSE IF ( LFLAG .EQ. 0 ) THEN
C
C............ normal finalization ......................................
C
              WRITE (*,'(/A)') 'execution finished: standard exit'
C
      ELSE IF ( LFLAG .EQ. 1 ) THEN
C
C............ given (U), compute (V) = (C)^T*(C)*(U) ...................
C
              DO 80 I = 1,NVOPU
C
C............... (Z) = (C)*(U) .........................................
C
                 DO 40 J = 1,M
                    Z(J,I) = ZERO
   40            CONTINUE
C
                 DO 50 J = 1,NE
                    Z(IRN(J),I) = Z(IRN(J),I) + C(J)*U(JCN(J),I)
   50            CONTINUE 
C
C............... (V) = (C)^T*(Z) .......................................
C
                 DO 60 J = 1,N
                    V(J,I) = ZERO
   60            CONTINUE
C
                 DO 70 J = 1,NE
                    V(JCN(J),I) = V(JCN(J),I) + C(J)*Z(IRN(J),I)
   70            CONTINUE
C
   80         CONTINUE
C
              GO TO 30
C
      ELSE 
C
C............ other flags should not be used here ......................
C
              STOP '* Error: LFLAG does not apply in this case *' 
C
      END IF
C
      STOP
    1 STOP '* IO error: standard input *' 
    2 STOP '* IO error: file MATRXC *'
C
C**** end of DRVSP3 ****************************************************
C
      END
