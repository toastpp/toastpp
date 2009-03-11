C***********************************************************************
C*                                                                     *
C*    This file contains:                                              *
C*                                                                     *
C*    PIIGRP : defines GROUP                                           *
C*    PIINPE : interface for NPE                                       *
C*    PIIPID : interface for PID                                       *
C*    PIIRED : interface for a parallel reduce operation               *
C*    PIITID : interface for TID                                       *
C*                                                                     *
C***********************************************************************
C
      SUBROUTINE PIIGRP (GROUP,LCOMM)
C
      INTEGER   LCOMM
      CHARACTER GROUP*13
C
C     PIIGRP defines GROUP
C     ======        
C
      INCLUDE 'fpvm3.h'
C
      WRITE (GROUP,'(A12,I1)') 'blzpack.pvm.',LCOMM
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PIINPE (NPE,LCOMM,INFO)
C
      INTEGER INFO,LCOMM,NPE
C
C     PIIPID is an interface for obtaining the number of processes
C     ======                        
C
      INTEGER   IINFO
      CHARACTER GROUP*13
C
      INCLUDE 'fpvm3.h'
C
      INFO = 0
C
      CALL PIIGRP    (GROUP,LCOMM)
      CALL PVMFGSIZE (GROUP,NPE)
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
      INTEGER   TID
      CHARACTER GROUP*13
C
      INCLUDE 'fpvm3.h'
C
      INFO = 0
C
      CALL PIIGRP      (GROUP,LCOMM)
      CALL PVMFMYTID   (TID)
      CALL PVMFGETINST (GROUP,TID,PID)
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
      CALL PVMFBARRIER (GROUP,NPE,IINFO)
C
      IF      ( FUNC .EQ. 'MAX' ) THEN
              CALL PVMFREDUCE  (PVMMAX,Y,N,INTEGER4,1,GROUP,0,IINFO)
      ELSE IF ( FUNC .EQ. 'MIN' ) THEN
              CALL PVMFREDUCE  (PVMMIN,Y,N,INTEGER4,1,GROUP,0,IINFO)
      ELSE IF ( FUNC .EQ. 'SUM' ) THEN
              CALL PVMFREDUCE  (PVMSUM,Y,N,INTEGER4,1,GROUP,0,IINFO)
      ELSE
              INFO = 1
      END IF
C
      IF ( PID .EQ. 0 ) THEN
         CALL PVMFINITSEND (PVMDATARAW,IINFO)
         CALL PVMFPACK     (INTEGER4,Y,N,1,IINFO)
         CALL PVMFBCAST    (GROUP,4,IINFO)
      ELSE
         CALL PVMFRECV     (-1,4,IINFO)
         CALL PVMFUNPACK   (INTEGER4,Y,N,1,IINFO)
      END IF
C
      IF ( IINFO .LT. 0 ) INFO = 2
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PIITID (PID,TID,LCOMM)
C
      INTEGER LCOMM,PID,TID
C
C     PIITID is an interface for obtaining the TID from PID
C     ====== (used in send/receive operations)
C
      CHARACTER GROUP*13
C
      INCLUDE 'fpvm3.h'
C
      CALL PIIGRP     (GROUP,LCOMM)
      CALL PVMFGETTID (GROUP,PID,TID)
C
      RETURN
      END
