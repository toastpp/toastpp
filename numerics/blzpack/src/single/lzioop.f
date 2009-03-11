      SUBROUTINE LZIOOP (FHNDL,LCOMM,LRERR,IY,N,Y,ARRAY,TASK)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZIOOP performs I/O operations                                   *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    FHNDL (sii) : IO handler                                         *
C*    LCOMM (sii) : communicator for the parallel version              *
C*    LRERR (sio) : code for error messages                            *
C*    IY    (sii) : index of (Y)                                       *
C*    N     (sii) : dimension of (Y)                                   *
C*    Y     (arb) : array copied from/to the file                      *
C*    ARRAY (sci) : formal name of (Y)                                 *
C*                  - `BQ' for (B) * Lanczos vectors                   *
C*                  - `BX' for (B) * Ritz vectors                      *
C*                  - `Q ' for Lanczos vectors                         *
C*                  - `X ' for Ritz vectors                            *
C*    TASK  (sci) : task to be performed                               *
C*                  - if `DEL', delete the file                        *
C*                  - if `OPN', open   the file                        *
C*                  - if `GET', (Y) is copied from the file            *
C*                  - if `PUT', (Y) is copied to   the file            *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SETLRM                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER          FHNDL(*),LCOMM,LRERR,IY,N
      REAL             Y(N)
      CHARACTER        ARRAY*2,TASK*3 
C
C==== local variables ==================================================
C
      INTEGER          ID,IOERR,J
      LOGICAL          ONLINE
      CHARACTER        LABEL*16 
C
C**** executable statements ********************************************
C
      IF ( LRERR .NE. 0 ) RETURN
C
C.... define units and file names ......................................
C
      IF      ( ARRAY .EQ. 'BQ' ) THEN
              ID = FHNDL(1)
              IOERR = 26
      ELSE IF ( ARRAY .EQ. 'BX' ) THEN
              ID = FHNDL(2)
              IOERR = 27
      ELSE IF ( ARRAY .EQ. 'Q ' ) THEN
              ID = FHNDL(3)
              IOERR = 28
      ELSE IF ( ARRAY .EQ. 'X ' ) THEN
              ID = FHNDL(4)
              IOERR = 29
      END IF     
C
      WRITE (LABEL,'(A8,I2.2,A1,A2)') 'blzpack.',ID,'.',ARRAY
C
C.... inquire file .....................................................
C
      IF ( TASK.EQ.'OPN' .OR. TASK.EQ.'DEL' .OR. IY.EQ.1 ) THEN
C
         INQUIRE  (UNIT   = ID           ,
     &             ERR    = 10           ,
     &             OPENED = ONLINE       )
C
         IF ( ONLINE .AND. TASK.EQ.'DEL' ) THEN
C
C...........close file .................................................
C
            CLOSE (UNIT   = ID           ,
     &             ERR    = 10           ,
     &             STATUS = 'DELETE'     )
C
            RETURN
C
         END IF
C
      END IF
C
C.... open/read/write ..................................................
C
      IF      ( (ARRAY.EQ.'BQ') .OR. (ARRAY.EQ.'Q ') ) THEN
C
C............ operations related with BQ and Q .........................
C
              IF ( .NOT.ONLINE .AND. TASK.EQ.'OPN' .OR. IY.EQ.1 ) THEN
C
C............... open  file ............................................
C
                 OPEN  (UNIT   = ID           ,
     &                  ERR    = 10           ,
     &                  FILE   = LABEL        ,
     &                  STATUS = 'UNKNOWN'    ,
     &                  ACCESS = 'SEQUENTIAL' ,
     &                  FORM   = 'UNFORMATTED')
C
              END IF
C
              IF      ( TASK .EQ. 'GET' ) THEN
C
C.................... read  data .......................................
C
	              IF ( IY .EQ. 1 ) REWIND (UNIT=ID,ERR=10)  
C
                      READ  (UNIT=ID,ERR=10) (Y(J),J=1,N)
C
              ELSE IF ( TASK .EQ. 'PUT' ) THEN
C
C.................... write data .......................................
C
	              IF ( IY .EQ. 1 ) REWIND (UNIT=ID,ERR=10)  
C
                      WRITE (UNIT=ID,ERR=10) (Y(J),J=1,N)
C
              END IF
C
      ELSE IF ( (ARRAY.EQ.'BX') .OR. (ARRAY.EQ.'X ') ) THEN
C
C............ operations related with BX and X .........................
C
              IF ( .NOT.ONLINE .AND. TASK.EQ.'OPN' .OR. IY.EQ.1 ) THEN
C
C............... open  file ............................................
C
                 OPEN  (UNIT   = ID           ,
     &                  ERR    = 10           ,
     &                  RECL   = N*8          ,
     &                  FILE   = LABEL        ,
     &                  ACCESS = 'DIRECT'     ,
     &                  STATUS = 'UNKNOWN'    ,
     &                  FORM   = 'UNFORMATTED')
C
              END IF
C
              IF      ( TASK .EQ. 'GET' ) THEN
C
C.................... read  data .......................................
C
  	              READ  (REC=IY,UNIT=ID,ERR=10) (Y(J),J=1,N)
C
              ELSE IF ( TASK .EQ. 'PUT' ) THEN
C
C.................... write data .......................................
C
                      WRITE (REC=IY,UNIT=ID,ERR=10) (Y(J),J=1,N)
C
              END IF
C
      END IF
C
      RETURN 
C
C.... IO error .........................................................
C
   10 CONTINUE
C
      CALL SETLRM (IOERR,LRERR)
C
      RETURN 
C
C**** end of LZIOOP ****************************************************
C
      END
