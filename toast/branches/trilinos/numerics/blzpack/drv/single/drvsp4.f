      PROGRAM DRVSP4
C     ==============
C
C***********************************************************************
C*                                                                     *
C*    DRVSP4 is a driver model for the standard eigenvalue problem     *
C*                                                                     *
C*           (A)*(x)-eig*(x)=(0)                                       *
C*                                                                     *
C*    where (A) is a symmetric matrix of dimension M+N defined as      *
C*                                                                     *
C*                 |  (0)  (C) |                                       *
C*           (A) = |           |                                       *
C*                 | (C)^T (0) |                                       *
C*                                                                     *
C*    where (C) is an M-by-N matrix. Therefore, the pairs eig,(x)      *
C*    lead to the singular value decomposition of (C),                 *
C*                                                                     *
C*           (C) = (U)*(S)*(V)^T                                       *
C*                                                                     *
C*           eig = +/- s                                               *
C*                                                                     *
C*                 | (u) |                                             *
C*           (x) = |     |*(1/sqrt(2))                                 *
C*                 | (v) |                                             *
C*                                                                     *
C*    Data for DRVSP4 are read from the standard input (see below).    *
C*    The entries of (C) are read from the file MATRXC (coordinate     *
C*    format).                                                         *
C*                                                                     *
C*    Examples: files `drvsp4.dat' and `C2.dat'                        *
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
      LOGICAL          ODD
C
C.... BLZDRS variables .................................................
C
      INTEGER          ISTOR(LISTOR),LFLAG,NNEIG,NRUN,NVOPU
      REAL             EIG(LEIG,2),RSTOR(LRSTOR),SIGMA,
     &                 U(LN,NCUV),V(LN,NCUV),X(LN,LEIG)
C
C.... matrix (C) .......................................................
C
      INTEGER          IRN(MAXC),JCN(MAXC),M,MPN,N,NE
      REAL             C(MAXC)
      CHARACTER        MATRXC*16
C
C.... external function ................................................
C
      REAL             SIRAND
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
      MPN = M + N
C
C.... check the dimensions .............................................
C
      IF ( MPN .GT. LN ) STOP '* Error: MPN > LN *'
      IF ( ISTOR(5) .GT. NCUV ) STOP '* Error: NVB > NCUV *'
      IF ( ISTOR(5) .NE. ISTOR(7) ) STOP '* Set ISTOR(7)=ISTOR(5) *'
C
      ISTOR( 1) = MPN
      ISTOR( 2) = LN
      ISTOR( 4) = LEIG
      ISTOR(14) = 0     
      ISTOR(15) = LISTOR    
      RSTOR( 4) = LRSTOR    
C
C.... reverse communication strategy ...................................
C
      NRUN = 0
      LFLAG = 0
      ODD = M.LE.N
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
C............ given (U), compute (V)=(A)*(U) ...........................
C
              DO 70 I = 1,NVOPU
                 DO 40 J = 1,MPN
                    V(J,I) = ZERO
   40            CONTINUE
                 IF ( ODD ) THEN
                    DO 50 J = 1,NE
                       V(JCN(J)+M,I) = V(JCN(J)+M,I)+C(J)*U(IRN(J)  ,I)
   50               CONTINUE
                 ELSE
                    DO 60 J = 1,NE
                       V(IRN(J)  ,I) = V(IRN(J)  ,I)+C(J)*U(JCN(J)+M,I)
   60               CONTINUE
                 END IF
   70         CONTINUE
C
              ODD = .NOT.ODD
C
              GO TO 30
C
      ELSE IF ( LFLAG .EQ. 4 ) THEN
C
C............ starting vectors (SIRAND is borrowed from BLZPACK) .......
C
              NRUN = NRUN + 1
              V(1,1) = SIRAND(NRUN*2)
C
              DO 110 I = 1,NVOPU
                 DO 80 J = 1,MPN
                    V(J,I) = ZERO
   80            CONTINUE
                 IF ( ODD ) THEN
                    DO 90 J = 1,M
                       V(J  ,I) = SIRAND(0)
   90               CONTINUE
                 ELSE
                    DO 100 J = 1,N
                       V(J+M,I) = SIRAND(0)
  100               CONTINUE
                 END IF
  110         CONTINUE
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
C**** end of DRVSP4 ****************************************************
C
      END
