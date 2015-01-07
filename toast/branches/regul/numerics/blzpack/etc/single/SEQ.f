C***********************************************************************
C*                                                                     *
C*    This file contains:                                              *
C*                                                                     *
C*    PISGTR : interface for a parallel gather operation               *
C*    PISRED : interface for a parallel reduce operation               *
C*                                                                     *
C***********************************************************************
C
      SUBROUTINE PISGTR (N,X,Y,LCOMM,INFO)
C
      INTEGER          INFO,LCOMM,N
      REAL             X(*),Y(*)
C
C     PISGTR is an interface for a parallel gather operation 
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
      SUBROUTINE PISRED (FUNC,N,X,Y,LCOMM,INFO)
C
      INTEGER          INFO,LCOMM,N
      REAL             X(N),Y(N)
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
