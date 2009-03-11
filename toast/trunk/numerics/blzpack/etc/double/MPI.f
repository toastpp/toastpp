C***********************************************************************
C*                                                                     *
C*    This file contains:                                              *
C*                                                                     *
C*    PIDGTR : interface for a parallel gather operation               *
C*    PIDRED : interface for a parallel reduce operation               *
C*                                                                     *
C***********************************************************************
C
      SUBROUTINE PIDGTR (N,X,Y,LCOMM,INFO)
C
      INTEGER          INFO,LCOMM,N
      DOUBLE PRECISION X(*),Y(*)
C
C     PIDGTR is an interface for a parallel gather operation 
C     ====== 
C
      INTEGER I,IINFO,J,K,PID,NPE,STATUS(10)
C
      include 'mpif.h'
C
      INFO = 0
C
      J = N + 1
C
      CALL MPI_COMM_SIZE (LCOMM,NPE,IINFO)
      CALL MPI_COMM_RANK (LCOMM,PID,IINFO)
C
      IF ( IINFO.NE.0 ) INFO = 1
C
      IF ( PID .EQ. 0 ) THEN
         DO 10 I = 1,N
            Y(I) = X(I)
   10    CONTINUE
      END IF
C
      IF ( NPE .EQ. 1 ) RETURN
C
      DO 20 I = 1,NPE-1
         IF      ( PID.EQ.0 ) THEN
                 CALL MPI_RECV (K   ,1,MPI_INTEGER         ,I,1,
     &                          LCOMM,STATUS,IINFO)
                 CALL MPI_RECV (Y(J),K,MPI_DOUBLE_PRECISION,I,2,
     &                          LCOMM,STATUS,IINFO)
                 J = J + K
         ELSE IF ( PID.EQ.I ) THEN
                 CALL MPI_SEND (N,1,MPI_INTEGER         ,0,1,
     &                          LCOMM,IINFO)
                 CALL MPI_SEND (X,N,MPI_DOUBLE_PRECISION,0,2,
     &                          LCOMM,IINFO)
         END IF
         IF ( IINFO.NE.0 ) INFO = 2 
   20 CONTINUE
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PIDRED (FUNC,N,X,Y,LCOMM,INFO)
C
      INTEGER          INFO,LCOMM,N
      DOUBLE PRECISION X(N),Y(N)
      CHARACTER        FUNC*3
C
C     PIDRED is an interface for a parallel reduce operation
C     ====== 
C
C            if FUNC = 'MAX', (Y) = max(X)
C            if FUNC = 'MIN', (Y) = min(X)
C            if FUNC = 'SUM', (Y) = sum(X)
C
      INTEGER IINFO
C
      include 'mpif.h'
C
      INFO = 0
C
      IF      ( FUNC .EQ. 'MAX' ) THEN
              CALL MPI_ALLREDUCE (X,Y,N,MPI_DOUBLE_PRECISION,MPI_MAX,
     &                            LCOMM,IINFO)
      ELSE IF ( FUNC .EQ. 'MIN' ) THEN
              CALL MPI_ALLREDUCE (X,Y,N,MPI_DOUBLE_PRECISION,MPI_MIN,
     &                            LCOMM,IINFO)
      ELSE IF ( FUNC .EQ. 'SUM' ) THEN
              CALL MPI_ALLREDUCE (X,Y,N,MPI_DOUBLE_PRECISION,MPI_SUM,
     &                            LCOMM,IINFO)
      ELSE
              INFO = 1
      END IF
C
      IF ( IINFO.NE.0 ) INFO = 2 
C
      RETURN
      END
