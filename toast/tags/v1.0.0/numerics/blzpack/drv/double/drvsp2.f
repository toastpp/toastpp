      PROGRAM DRVSP2
C     ==============
C
C***********************************************************************
C*                                                                     *
C*    DRVSP2 is a driver model for the standard eigenvalue problem     *
C*                                                                     *
C*           (A)*(x)-eig*(x)=(0)                                       *
C*                                                                     *
C*    where (A) is a symmetric matrix of dimension N. DRVSP2 uses      *
C*    eigenvalues/eigenvectors of (A) as input.                        *
C*                                                                     *
C*    Data for DRVSP2 are read from the standard input (see below).    *
C*    The upper triangle of (A) is read from the file MATRXA (in       *
C*    coordinate format) and the eigenvalues/eigenvectors of (A)       *
C*    are read from the file EIGENA.                                   *
C*                                                                     *
C*    Examples: files `drvsp2.dat', `A.dat' and `A.eig.dat'            *
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
      INTEGER          MAXA
      PARAMETER        (MAXA   =  10000)
      INTEGER          NCUV
      PARAMETER        (NCUV   =      3)
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO   =  0.0D0)
C
C.... work variables ...................................................
C
      INTEGER          I,J,K,NGEIG
      DOUBLE PRECISION AJK
C
C.... BLZDRD variables .................................................
C
      INTEGER          ISTOR(LISTOR),LFLAG,NNEIG,NVOPU
      DOUBLE PRECISION EIG(LEIG,2),RSTOR(LRSTOR),SIGMA,
     &                 U(LN,NCUV),V(LN,NCUV),X(LN,LEIG)
C
C.... matrix (A) .......................................................
C
      INTEGER          IRN(MAXA),JCN(MAXA),N,NE
      DOUBLE PRECISION A(MAXA)
      CHARACTER        EIGENA*16,MATRXA*16
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
      READ (*,ERR=1,FMT=*) MATRXA
      READ (*,ERR=1,FMT=*) EIGENA
C
C.... read (A) .........................................................
C
      N = 0
      NE = 0
      OPEN (UNIT=10,ERR=2,STATUS='OLD',FILE=MATRXA)
      DO 10 I = 1,LN*LN
         READ (UNIT=10,ERR=2,END=20,FMT=*) J,K,AJK
         NE = NE + 1
         IF ( NE .GT. MAXA ) STOP '* Error: NE > MAXA *'
         N = MAX(J,K,N)
         IRN(NE) = J
         JCN(NE) = K
         A(NE) = AJK
   10 CONTINUE
   20 CLOSE (UNIT=10,ERR=2)
C
C.... read eigenvalues and eigenvectors of (A) .........................
C
      NGEIG = ISTOR(8)
C
      OPEN (UNIT=10,ERR=3,STATUS='OLD',FILE=EIGENA)
      DO 30 I = 1,NGEIG
         READ (UNIT=10,ERR=3,FMT=*) EIG(I,1),EIG(I,2)
         READ (UNIT=10,ERR=3,FMT=*) (X(J,I),J=1,N)
   30 CONTINUE
      CLOSE (UNIT=10,ERR=3)
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
   40 CONTINUE
C
C     ***********************************************************
      CALL BLZDRD (ISTOR,RSTOR,SIGMA,NNEIG,U,V,LFLAG,NVOPU,EIG,X)
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
                 DO 50 J = 1,N
                    V(J,I) = ZERO
   50            CONTINUE
                 DO 60 J = 1,NE
                    IF ( IRN(J) .EQ. JCN(J) ) THEN
                       V(IRN(J),I) = V(IRN(J),I) + A(J)*U(JCN(J),I)
                    ELSE
                       V(IRN(J),I) = V(IRN(J),I) + A(J)*U(JCN(J),I)
                       V(JCN(J),I) = V(JCN(J),I) + A(J)*U(IRN(J),I)
                    END IF
   60            CONTINUE
   70         CONTINUE
C
              GO TO 40
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
    3 STOP '* IO error: file EIGENA *'
C
C**** end of DRVSP2 ****************************************************
C
      END
