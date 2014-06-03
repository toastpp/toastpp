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
      INTEGER   INFO,LCOMM,N
      REAL      X(*),Y(*)
C
C     PISGTR is an interface for a parallel gather operation 
C     ======
C
      INTEGER   I,IINFO,J,K,PID,NPE,TID0,TIDI
      CHARACTER GROUP*13
C
      INCLUDE 'fpvm3.h'
C
      INFO = 0
C
      J = N + 1
C
      CALL PIIGRP (GROUP,LCOMM)
      CALL PIIPID (PID,LCOMM,IINFO)
      CALL PIINPE (NPE,LCOMM,IINFO)
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
         IF      ( PID .EQ. 0 ) THEN
                 CALL PIITID       (I,TIDI,LCOMM)
                 CALL PVMFRECV     (TIDI,I,IINFO)
                 CALL PVMFUNPACK   (INTEGER8,K,1,1,IINFO)
                 CALL PVMFUNPACK   (REAL8,Y(J),K,1,IINFO)
                 J = J + K
         ELSE IF ( PID .EQ. I ) THEN
                 CALL PIITID       (0,TID0,LCOMM)
                 CALL PVMFINITSEND (PVMDATARAW,IINFO)
                 CALL PVMFPACK     (INTEGER8,N,1,1,IINFO)
                 CALL PVMFPACK     (REAL8,X,N,1,IINFO)
                 CALL PVMFSEND     (TID0,I,IINFO)
         END IF
   20 CONTINUE
C
      CALL PVMFBARRIER (GROUP,NPE,IINFO)
C
      IF ( IINFO .LT. 0 ) INFO = 2
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
C     PISRED is an interface for a parallel reduce operation
C     ====== 
C
C            if FUNC = 'MAX', (Y) = max(X)
C            if FUNC = 'MIN', (Y) = min(X)
C            if FUNC = 'SUM', (Y) = sum(X)
C
      INTEGER   I,IINFO,NPE,PID
      CHARACTER GROUP*13
C
      INCLUDE 'fpvm3.h'
C
      EXTERNAL PVMMAX
      EXTERNAL PVMMIN
      EXTERNAL PVMSUM
C
      INFO = 0
C
      CALL PIIGRP (GROUP,LCOMM)
      CALL PIIPID (PID,LCOMM,IINFO)
      CALL PIINPE (NPE,LCOMM,IINFO)
C
      DO 10 I = 1,N
         Y(I) = X(I)
   10 CONTINUE
C
      IF ( NPE .EQ. 1 ) RETURN
C
      IF      ( FUNC .EQ. 'MAX' ) THEN
              CALL PVMFREDUCE  (PVMMAX,Y,N,REAL8,1,GROUP,0,IINFO)
      ELSE IF ( FUNC .EQ. 'MIN' ) THEN
              CALL PVMFREDUCE  (PVMMIN,Y,N,REAL8,1,GROUP,0,IINFO)
      ELSE IF ( FUNC .EQ. 'SUM' ) THEN
              CALL PVMFREDUCE  (PVMSUM,Y,N,REAL8,1,GROUP,0,IINFO)
      ELSE
              INFO = 1
      END IF
C
      CALL PVMFBARRIER (GROUP,NPE,IINFO)
C
      IF ( IINFO .LT. 0 ) INFO = 2 
C
      IF ( PID .EQ. 0 ) THEN
         CALL PVMFINITSEND (PVMDATARAW,IINFO)
         CALL PVMFPACK     (REAL8,Y,N,1,IINFO)
         CALL PVMFBCAST    (GROUP,8,IINFO)
      ELSE
         CALL PVMFRECV     (-1,8,IINFO)
         CALL PVMFUNPACK   (REAL8,Y,N,1,IINFO)
      END IF
C
      IF ( IINFO .LT. 0 ) INFO = 3 
C
      RETURN
      END
