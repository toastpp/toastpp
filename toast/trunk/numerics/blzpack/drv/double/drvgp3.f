      PROGRAM DRVGP3
C     ==============
C
C***********************************************************************
C*                                                                     *
C*    DRVGP3 is a driver model for the generalized eigenvalue problem  *
C*                                                                     *
C*           (A)*(x)-eig*(B)*(x)=(0)                                   *
C*                                                                     *
C*    where (A) is a symmetric matrix of dimension N and (B) has the   *
C*    same pattern of (A).                                             *
C*                                                                     *
C*    Data for DRVGP3 are read from the standard input (see below).    *
C*    The upper triangle of (A) is read from the file MATRXA (in       *
C*    coordinate format) and the upper triangle of (B) is read         *
C*    from the file MATRXB (in coordinate format).                     *
C*                                                                     *
C*    Examples: files `drvgp3.dat', `A.dat' and `B2.dat'               *
C*    Factorization and solver used: MA47AD, MA47BD, MA47CD, MA47ID    *
C*    BLAS kernels used: DAXPY,DCOPY                                   *
C*                                                                     *
C***********************************************************************
C
C.... parameters .......................................................
C
      INTEGER          LEIG
      PARAMETER        (LEIG   =     40)
      INTEGER          LISTOR
      PARAMETER        (LISTOR =  10000)
      INTEGER          LN
      PARAMETER        (LN     =    110)
      INTEGER          LRSTOR
      PARAMETER        (LRSTOR =  10000)
      INTEGER          MAXA
      PARAMETER        (MAXA   =  10000)
      INTEGER          MAXIW1
      PARAMETER        (MAXIW1 =  10000)
      INTEGER          MAXNE
      PARAMETER        (MAXNE  =   1100)
      INTEGER          NCUV
      PARAMETER        (NCUV   =      3)
      DOUBLE PRECISION ZERO
      PARAMETER        (ZERO   =  0.0D0)
C
C.... work variables ...................................................
C
      INTEGER          I,J,K,L
      DOUBLE PRECISION AJK,BJK 
C
C.... BLZDRD variables .................................................
C
      INTEGER          ISTOR(LISTOR),LFLAG,NNEIG,NVOPU
      DOUBLE PRECISION EIG(LEIG,2),RSTOR(LRSTOR),SIGMA,
     &                 U(LN,NCUV),V(LN,NCUV),X(LN,LEIG)
C
C.... matrices (A) and (B) .............................................
C
      INTEGER          ICNTL(7),INFO(24),IRN(MAXNE),IW1(MAXIW1),
     &                 IW2(LN*2+2),JCN(MAXNE),KEEP(MAXNE+LN*5+2),N,NE
      DOUBLE PRECISION A(MAXA),B(MAXA),C(MAXA),CNTL(2),RINFO(4),W(LN)
      CHARACTER        MATRXA*16,MATRXB*16
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
      READ (*,ERR=1,FMT=*) MATRXB
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
         IF ( NE .GT. MAXNE ) STOP '* Error: NE > MAXNE *'
         N = MAX(J,K,N)
         IRN(NE) = J
         JCN(NE) = K
         A(NE) = AJK
         B(I) = ZERO
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
C.... read (B), same pattern of (A) ....................................
C
      L = 1
      OPEN (UNIT=10,ERR=3,STATUS='OLD',FILE=MATRXB)
   30 CONTINUE
      READ (UNIT=10,ERR=3,END=50,FMT=*) J,K,BJK
      DO 40 I = L,NE
         IF ( J.EQ.IRN(I) .AND. K.EQ.JCN(I) ) THEN
            B(I) = BJK
            L = I + 1
            GO TO 30
         END IF
   40 CONTINUE
      STOP '* Error: Patterns of (A) and (B) differ *'
   50 CLOSE (UNIT=10,ERR=3)
C
C.... set default parameters for MA47 and analyse sparsity pattern .....
C
      CALL MA47ID (CNTL,ICNTL)
C
      CALL MA47AD (N,NE,IRN,JCN,IW1,MAXIW1,KEEP,ICNTL,RINFO,INFO)
C
      IF ( INFO(1).NE.0 ) STOP '* Error: MA47AD, INFO(1)>0 *'
C
C.... reverse communication strategy ...................................
C
      LFLAG = 0
C
   60 CONTINUE
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
C............ given (U), solve (C)*(V)=(U) .............................
C
              DO 70 I = 1,NVOPU
                 CALL DCOPY  (N,U(1,I),1,V(1,I),1)
                 CALL MA47CD (N,C,MAXA,IW1,MAXIW1,W,V(1,I),IW2,ICNTL)
   70         CONTINUE
C
              GO TO 60
C
      ELSE IF ( LFLAG .EQ. 2 ) THEN
C
C............ given (U), compute (V)=(B)*(U) ...........................
C
              DO 100 I = 1,NVOPU
                 DO 80 J = 1,N
                    V(J,I) = ZERO
   80            CONTINUE
                 DO 90 J = 1,NE
                    IF ( IRN(J) .EQ. JCN(J) ) THEN
                       V(IRN(J),I) = V(IRN(J),I) + B(J)*U(JCN(J),I)
                    ELSE
                       V(IRN(J),I) = V(IRN(J),I) + B(J)*U(JCN(J),I)
                       V(JCN(J),I) = V(JCN(J),I) + B(J)*U(IRN(J),I)
                    END IF
   90            CONTINUE
  100         CONTINUE
C
              GO TO 60
C
      ELSE IF ( LFLAG .EQ. 3 ) THEN
C
C............ given SIGMA, form (C)=(A)-SIGMA*(B) ......................
C
              CALL DCOPY (NE,A,1,C,1)
              CALL DAXPY (NE,-SIGMA,B,1,C,1)
C
C............ factor (C)=(L)*(D)*(L') ..................................
C
              CALL MA47BD (N,NE,JCN,C,MAXA,IW1,MAXIW1,KEEP,
     &                     CNTL,ICNTL,IW2,RINFO,INFO)
C
              IF ( INFO( 1).NE.0 ) STOP '* Error: MA47BD, INFO(1)>0 *'
              IF ( INFO(24).GT.0 ) STOP '* Error: MA47BD, INFO(24)>0 *'
C
              NNEIG = INFO(23)
C
              GO TO 60
C
      END IF
C
      STOP
    1 STOP '* IO error: standard input *'
    2 STOP '* IO error: file MATRXA *'
    3 STOP '* IO error: file MATRXB *'
C
C**** end of DRVGP3 ****************************************************
C
      END
