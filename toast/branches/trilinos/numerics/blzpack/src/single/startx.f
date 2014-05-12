      SUBROUTINE STARTX (LCOMM,LRERR,LNI,NI,NX,X)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    STARTX checks the starting eigenvectors defined by the user      *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    LNI   (sii) : leading dimension of (X)                           *
C*    NI    (sii) : dimension of the vectors in (X)                    *
C*    NX    (sii) : number of vectors in (X)                           *
C*    X     (ari) : starting eigenvectors                              *
C*                                                                     *
C*  - Subprograms:                                                     *
C*                                                                     *
C*    PISRED,PIIRED,SETLRM                                             *
C*                                                                     *
C*  - BLAS kernel:                                                     *
C*                                                                     *
C*    SASUM                                                            *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          LCOMM,LNI,LRERR,NI,NX
      REAL             X(LNI,NX)
C
C==== local variables ==================================================
C
      REAL             GSUM,SUMX
      INTEGER          GERR,I,INFO
C
C==== BLAS kernel ======================================================
C
      REAL             SASUM
C
C**** executable statements ********************************************
C
C.... loop on the columns of X .........................................
C
      DO 10 I = 1,NX
C
C....... sum up all entries of X .......................................
C
         SUMX = SASUM(NI,X(1,I),1)
C
C....... check for NaNs ................................................
C
         IF ( SUMX .NE. SUMX ) CALL SETLRM (26,LRERR)
         CALL PIIRED ('MAX',1,LRERR,GERR,LCOMM,INFO)
         IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
         IF ( GERR .NE. 0 ) CALL SETLRM (26,LRERR)
C
C....... check for a zero X ............................................
C
         IF ( LRERR .EQ. 0 ) THEN
            CALL PISRED ('MAX',1,SUMX,GSUM,LCOMM,INFO)
            IF ( INFO .NE. 0 ) CALL SETLRM (32,LRERR)
            IF ( GSUM .EQ. 0 ) CALL SETLRM (26,LRERR)
         END IF
C
   10 CONTINUE
C
      RETURN
C
C**** end of STARTX ****************************************************
C
      END
