      PROGRAM DRVSP1
C     ==============
C
C***********************************************************************
C*                                                                     *
C*    DRVSP1 is a driver model for the standard eigenvalue problem     *
C*                                                                     *
C*           (A)*(x)-eig*(x)=(0)                                       *
C*                                                                     *
C*    where (A) is a symmetric matrix of dimension N.                  *
C*                                                                     *
C*    Data for DRVSP1 are read from the standard input (see below).    *
C*    The upper triangle of (A) is read from the file MATRXA (in       *
C*    coordinate format).                                              *
C*                                                                     *
C*    Examples: files `drvsp1.dat' and `A.dat'                         *
C*    BLAS kernel used: SSYMM                                          *
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
      INTEGER          NCUV
      PARAMETER        (NCUV   =      3)
      REAL             ONE
      PARAMETER        (ONE    =  1.0E0)
      REAL             ZERO
      PARAMETER        (ZERO   =  0.0E0)
C
C.... work variables ...................................................
C
      INTEGER          I,J,K
C
C.... BLZDRS variables .................................................
C
      INTEGER          ISTOR(LISTOR),LFLAG,NNEIG,NVOPU
      REAL             EIG(LEIG,2),RSTOR(LRSTOR),SIGMA,
     &                 U(LN,NCUV),V(LN,NCUV),X(LN,LEIG)
C
C.... matrix (A) .......................................................
C
      INTEGER          N
      REAL             A(LN,LN)
      CHARACTER        MATRXA*16
C
C=======================================================================
C
C.... initialize (A) ...................................................
C
      DO 10 I = 1,LN
         DO 20 J = 1,LN
            A(J,I) = ZERO
   20    CONTINUE
   10 CONTINUE
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
      READ (*,ERR=1,FMT=*) MATRXA
C
C.... read (A) .........................................................
C
      N = 0
      OPEN (UNIT=10,ERR=2,STATUS='OLD',FILE=MATRXA)
      DO 30 I = 1,LN*LN
         READ (UNIT=10,ERR=2,END=40,FMT=*) J,K,A(J,K)
         N = MAX(J,K,N)
   30 CONTINUE 
   40 CLOSE (UNIT=10,ERR=2)
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
   50 CONTINUE
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
              CALL SSYMM ('L','U',N,NVOPU,ONE,A,LN,U,LN,ZERO,V,LN)
C
              GO TO 50
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
    2 STOP '* IO error: file MATRXA *'
C
C**** end of DRVSP1 ****************************************************
C
      END
