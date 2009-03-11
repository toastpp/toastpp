      SUBROUTINE LZWRNS (LFILE,LRWRN)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZWRNS prints warning messages                                   *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE (sii) : file unit for output                               *
C*    LRWRN (sii) : code for warning messages                          *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   LFILE,LRWRN
C
C==== local variables ==================================================
C
      INTEGER   IBIT          
      CHARACTER WM(17)*50
C
C==== subprogram =======================================================
C
      LOGICAL   SIBTST
C
C=======================================================================
C
      DATA WM(01) /'no eigenvalue in computational interval           '/
      DATA WM(02) /'SIGMA is too close to an eigenvalue               '/
      DATA WM(03) /'invariant or ill-conditioned starting block       '/
      DATA WM(04) /'invariant or ill-conditioned subspace             '/
      DATA WM(05) /'reduced eigenproblem can not be solved            '/
      DATA WM(06) /'three runs without any new information            '/
      DATA WM(07) /'maximum number of runs has been reached           '/
      DATA WM(08) /'algorithm unable to deal with subinterval         '/
      DATA WM(09) /'not enough space to continue eigensolution        '/
      DATA WM(10) /'execution interrupted due to a big cluster        '/
      DATA WM(11) /'no eigenpair has been computed                    '/
      DATA WM(12) /'no eigenpair has been checked                     '/
      DATA WM(13) /'less than NREIG eigenpairs have been accepted     '/
      DATA WM(14) /'less than NREIG eigenpairs have been checked      '/
      DATA WM(15) /'matrix (B) seems to be indefinite                 '/
      DATA WM(16) /'orthogonality check not performed                 '/
      DATA WM(17) /'LFLAG=5 (early termination)                       '/
C
C**** executable statements ********************************************
C
      IF ( LRWRN .GT. 0 ) THEN
         WRITE (LFILE,1000)
         DO 10 IBIT = 1,30
            IF ( SIBTST(IBIT,LRWRN) ) WRITE (LFILE,1001) IBIT,WM(IBIT)
   10    CONTINUE
      END IF
C
      RETURN 
C
 1000 FORMAT ()
 1001 FORMAT ('* Warning (',I2,'): ',A50)
C
C**** end of LZWRNS ****************************************************
C
      END
