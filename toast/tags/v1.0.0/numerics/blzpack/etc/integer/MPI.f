C***********************************************************************
C*                                                                     *
C*    This file contains:                                              *
C*                                                                     *
C*    PIINPE : interface for NPE                                       *
C*    PIIPID : interface for PID                                       *
C*    PIIRED : interface for a parallel reduce operation               *
C*                                                                     *
C***********************************************************************
C
      SUBROUTINE PIINPE (NPE,LCOMM,INFO)
C
      INTEGER INFO,LCOMM,NPE
C
C     PIINPE is an interface for obtaining the number of processes
C     ======                        
C
      INTEGER IINFO
C
      INCLUDE 'mpif.h'
C
      INFO = 0
C
      CALL MPI_COMM_SIZE (LCOMM,NPE,IINFO)
C
      IF ( IINFO.NE.0 ) INFO = 1
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PIIPID (PID,LCOMM,INFO)
C
      INTEGER INFO,LCOMM,PID
C
C     PIIPID is an interface for obtaining the calling process ID
C     ====== 
C
      INTEGER IINFO
C
      INCLUDE 'mpif.h'
C
      INFO = 0
C
      CALL MPI_COMM_RANK (LCOMM,PID,IINFO)
C
      IF ( IINFO.NE.0 ) INFO = 1 
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PIIRED (FUNC,N,X,Y,LCOMM,INFO)
C
      INTEGER   INFO,LCOMM,N
      INTEGER   X(N),Y(N)
      CHARACTER FUNC*3
C
C     PIIRED is an interface for a parallel reduce operation
C     ====== 
C
C            if FUNC = 'MAX', (Y) = max(X)
C            if FUNC = 'MIN', (Y) = min(X)
C            if FUNC = 'SUM', (Y) = sum(X)
C
      INTEGER IINFO
C
      INCLUDE 'mpif.h'
C
      INFO = 0
C
      IF      ( FUNC .EQ. 'MAX' ) THEN
              CALL MPI_ALLREDUCE (X,Y,N,MPI_INTEGER,MPI_MAX,
     &                            LCOMM,IINFO)
      ELSE IF ( FUNC .EQ. 'MIN' ) THEN
              CALL MPI_ALLREDUCE (X,Y,N,MPI_INTEGER,MPI_MIN,
     &                            LCOMM,IINFO)
      ELSE IF ( FUNC .EQ. 'SUM' ) THEN
              CALL MPI_ALLREDUCE (X,Y,N,MPI_INTEGER,MPI_SUM,
     &                            LCOMM,IINFO)
      ELSE
              INFO = 1
      END IF
C
      IF ( IINFO.NE.0 ) INFO = 2
C
      RETURN
      END
