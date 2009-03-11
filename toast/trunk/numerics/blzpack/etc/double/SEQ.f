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
      INTEGER I
C
      INFO = 0
C
      DO 10 I = 1,N
         Y(I) = X(I)
   10 CONTINUE
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
C     PIDMAX is an interface for a parallel reduce operation
C     ====== 
C
      INTEGER I
C
      INFO = 0
C
      DO 10 I = 1,N
         Y(I) = X(I)
   10 CONTINUE
C
      RETURN
      END
