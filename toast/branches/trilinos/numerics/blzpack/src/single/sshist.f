      SUBROUTINE SSHIST (LFILE,NSLOG,SSLOG)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    SSHIST prints the spectrum slicing history                       *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE (sii) : file unit for output                               *
C*    NSLOG (sii) : number of subintervals recorded in SSLOG           *
C*    SSLOG (ari) : spectrum slicing history                           *
C*                                                                     *
C*  - Intrinsic Function:                                              *
C*                                                                     *
C*    NINT                                                             *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LFILE,NSLOG
      REAL             SSLOG(8,*)
C
C==== local variables ==================================================
C
      INTEGER          J
C
C==== intrinsic function ===============================================
C
      INTRINSIC        NINT
C
C**** executable statements ********************************************
C
      IF ( NSLOG .EQ. 0 ) RETURN
C
      WRITE (LFILE,1000) 
      WRITE (LFILE,1001) 
      WRITE (LFILE,1002) 
      WRITE (LFILE,1001) 
      WRITE (LFILE,1003) (SSLOG(1,J),NINT(SSLOG(5,J)),SSLOG(2,J),
     &                    SSLOG(3,J),NINT(SSLOG(6,J)),NINT(SSLOG(7,J)),
     &                    NINT(SSLOG(8,J)),SSLOG(4,J),J=1,NSLOG)
      WRITE (LFILE,1001) 
C
      RETURN 
C
 1000 FORMAT (/,'spectrum slicing history')
 1001 FORMAT (71('-'))
 1002 FORMAT (2X,'sigma',7X,'nneig',4X,'lower',7X,'upper',
     &        5X,'nritz',3X,'jl',3X,'jt',4X,'time')
 1003 FORMAT (1P,SP,E11.4,S,I8,1X,SP,2E12.4,S,I6,2I5,2X,1P,E9.2)
C
C**** end of SSHIST ****************************************************
C
      END
