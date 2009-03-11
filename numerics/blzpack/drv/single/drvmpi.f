      PROGRAM DRVMPI
C     ==============
C
C***********************************************************************
C*                                                                     *
C*    DRVMPI is a driver model for the standard eigenvalue problem     *
C*                                                                     *
C*           (A)*(x)-eig*(x)=(0)                                       *
C*                                                                     *
C*    where (A) is a symmetric matrix of dimension N. DRVMPI is the    *
C*    parallel, MPI based, version of DRVSP1.                          *
C*                                                                     *
C*    Data for DRVMPI are read from `drvsp1.dat' (see below). The      *
C*    upper triangle of (A) is read from the file MATRXA (in           *
C*    coordinate format).                                              *
C*                                                                     *
C*    Examples: files `drvsp1.dat' and `A.dat'                         *
C*    BLAS kernel used: SGEMM                                          *
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
      INTEGER          NNZMAX
      PARAMETER        (NNZMAX =   4000)
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
      REAL             EIG(LEIG,2),RSTOR(LRSTOR),SIGMA,T(LN*NCUV),
     &                 U(LN*NCUV),V(LN*NCUV),X(LN*LEIG)
C
C.... matrix (A) .......................................................
C
      INTEGER          ICOL(NNZMAX),ICOLL,ICOLR,IROW(NNZMAX),
     &                 N,NCOLL,NCOLP,NNZ,RMNDR
      REAL             A(LN,LN),S(NNZMAX)
      CHARACTER        MATRXA*16
C
C.... MPI ..............................................................
C
      INTEGER          II,INFO,JJ,MYPE,NI,NJ,NPE,ROOT,STATUS(10)
C
      INCLUDE          'mpif.h'
C
C=======================================================================
C
      CALL MPI_INIT      (INFO)
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,NPE ,INFO)
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,MYPE,INFO)
C
      IF ( MYPE .EQ. 0 ) THEN
C 
C....... read data .....................................................
C
         OPEN  (UNIT=10,ERR=1,STATUS='OLD',FILE='drvsp1.dat')
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR( 3)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR( 5)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR( 6)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR( 7)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR( 8)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR( 9)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR(10)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR(11)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR(12)
         READ  (UNIT=10,ERR=1,FMT=*) ISTOR(13)
         READ  (UNIT=10,ERR=1,FMT=*) RSTOR( 1)
         READ  (UNIT=10,ERR=1,FMT=*) RSTOR( 2)
         READ  (UNIT=10,ERR=1,FMT=*) RSTOR( 3)
         READ  (UNIT=10,ERR=1,FMT=*) MATRXA
         CLOSE (UNIT=10,ERR=1)
C
C....... read (A) ......................................................
C
         N = 0
         NNZ = 0
         OPEN  (UNIT=11,ERR=2,STATUS='OLD',FILE=MATRXA)
         DO 10 I = 1,LN*LN
            READ (UNIT=11,ERR=2,END=20,FMT=*) IROW(I),ICOL(I),S(I)
            N = MAX(N,IROW(I),ICOL(I))
            NNZ = NNZ + 1
   10    CONTINUE
   20    CLOSE (UNIT=11,ERR=2)
C
         NCOLL = N/NPE
         ICOLL = N + 1
         RMNDR = MOD(N,NPE)
C
      END IF
C
C.... broadcast ISTOR and RSTOR ........................................
C
      CALL MPI_BCAST (ISTOR,13,MPI_INTEGER         ,0,
     &                MPI_COMM_WORLD,INFO)
      CALL MPI_BCAST (RSTOR,3 ,MPI_REAL            ,0,
     &                MPI_COMM_WORLD,INFO)
C
C.... make sure each process contains the right piece of (A) ...........
C
      DO 60 I = 1,NPE
C
         IF      ( MYPE .EQ. 0 ) THEN
C
C                PE 0 distributes (IROW,ICOL,S) among the processes
C
                 NCOLP = NCOLL
                 IF ( RMNDR .GT. NPE-I ) NCOLP = NCOLP + 1
                 ICOLR = ICOLL - 1
                 ICOLL = ICOLR - NCOLP + 1
                 DO 40 J = 1,LN
                    DO 30 K = 1,LN
                       A(K,J) = ZERO
   30               CONTINUE
   40            CONTINUE
                 DO 50 J = 1,NNZ
                    IF ( ICOLL.LE.IROW(J) .AND. IROW(J).LE.ICOLR )
     &                 A(ICOL(J),IROW(J)-ICOLL+1) = S(J)
                    IF ( ICOLL.LE.ICOL(J) .AND. ICOL(J).LE.ICOLR )
     &                 A(IROW(J),ICOL(J)-ICOLL+1) = S(J)
   50            CONTINUE
                 IF ( I .EQ. NPE ) THEN
                    II = 1
                    NI = NCOLP
                 ELSE
                    CALL MPI_SEND (ICOLL,1       ,MPI_INTEGER         ,
     &                             NPE-I,I,MPI_COMM_WORLD,INFO)
                    CALL MPI_SEND (NCOLP,1       ,MPI_INTEGER         ,
     &                             NPE-I,I,MPI_COMM_WORLD,INFO)
                    CALL MPI_SEND (A    ,LN*NCOLP,MPI_REAL            ,
     &                             NPE-I,I,MPI_COMM_WORLD,INFO)
                 END IF
C
         ELSE IF ( MYPE .EQ. NPE-I ) THEN
C
C                PE NPE-I receives data from PE 0
C
                 CALL MPI_RECV (II     ,1    ,MPI_INTEGER         ,0,I,
     &                          MPI_COMM_WORLD,STATUS,INFO)
                 CALL MPI_RECV (NI     ,1    ,MPI_INTEGER         ,0,I,
     &                          MPI_COMM_WORLD,STATUS,INFO)
                 CALL MPI_RECV (A      ,LN*NI,MPI_REAL            ,0,I,
     &                          MPI_COMM_WORLD,STATUS,INFO)
C
         END IF
C
   60 CONTINUE
C
      ISTOR( 1) = NI
      ISTOR( 2) = NI
      ISTOR( 4) = LEIG
      ISTOR(14) = MPI_COMM_WORLD
      ISTOR(15) = LISTOR
      RSTOR( 4) = LRSTOR
C
C.... check the dimensions of (U) and (V) ..............................
C
      IF ( ISTOR(1) .GT. LN ) STOP '* Error: ISTOR(1) > LN *'
      IF ( ISTOR(5) .GT. NCUV ) STOP '* Error: ISTOR(5) > NCUV *'
C
C.... reverse communication strategy ...................................
C
      LFLAG = 0
C
   70 CONTINUE
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
              DO 80 I = 1,NPE
C
                 ROOT = I - 1
C
C............... each process deals with a pice of (A), (U) and (V) ....
C
                 IF ( NPE .EQ. 1 ) THEN
C
                    JJ = 1
                    NJ = NI
C
                 ELSE 
C
C.................. process I tells the others what it needs ...........
C
C                   NI rows of (U) and (V): A(II:II+NI-1,:) 
C
                    JJ = II
                    NJ = NI
                    CALL MPI_BCAST (JJ,1,MPI_INTEGER,ROOT,
     &                              MPI_COMM_WORLD,INFO)
                    CALL MPI_BCAST (NJ,1,MPI_INTEGER,ROOT,
     &                              MPI_COMM_WORLD,INFO)
C
                 END IF
C
C............... each process computes (T)=A(JJ:JJ+NJ-1,:)*U ...........
C
                 CALL SGEMM ('N','N',NJ,NVOPU,NI,ONE,A(JJ,1),LN,
     &                        U,NI,ZERO,T,NJ)
C
C............... process I gets the sum of (T) and puts it in (V) ......
C
                 CALL MPI_REDUCE (T,V,NJ*NVOPU,MPI_REAL            ,
     &                            MPI_SUM,ROOT,MPI_COMM_WORLD,INFO)
C
   80         CONTINUE

              GO TO 70
C
      ELSE 
C
C............ other flags should not be used here ......................
C
              STOP '* Error: LFLAG does not apply in this case *' 
C
      END IF
C
      CALL MPI_BARRIER  (MPI_COMM_WORLD,INFO)
      CALL MPI_FINALIZE (INFO)
C
      STOP
    1 STOP '* IO error: file drvsp1.dat *'
    2 STOP '* IO error: file MATRXA *'
C
C**** end of DRVMPI ****************************************************
C
      END
