      SUBROUTINE SETLRM (INDEX,LRMSG)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C
C*    SETLRM sets error and warning messages                           *
C*                                                                     *
C*  - Arguments:                                                       *
C
C*    INDEX (sii) : message index                                      *
C*    LRMSG (sib) : message control flag                               *
C*                                                                     *
C*  - Subprogram:                                                      *
C
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER INDEX,LRMSG
C
C==== function =========================================================
C
      LOGICAL SIBTST
C
C**** executable statements ********************************************
C
      IF ( .NOT.SIBTST(INDEX,LRMSG) ) LRMSG = LRMSG + 2**INDEX
C
      RETURN 
C
C**** end of SETLRM ****************************************************
C
      END
