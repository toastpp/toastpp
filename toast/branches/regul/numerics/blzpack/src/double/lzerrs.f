      SUBROUTINE LZERRS (LFILE,LRERR)
C
C***********************************************************************
C*                                                                     *
C*  - Purpose:                                                         *
C*                                                                     *
C*    LZERRS prints error messages                                     *
C*                                                                     *
C*  - Arguments:                                                       *
C*                                                                     *
C*    LFILE (sii) : file unit for output                               *
C*    LRERR (sii) : code for error messages                            *
C*                                                                     *
C*  - Subprogram:                                                      *
C*                                                                     *
C*    SIBTST                                                           *
C*                                                                     *
C***********************************************************************
C
C==== arguments ========================================================
C
      INTEGER   LFILE,LRERR
C
C==== local variables ==================================================
C
      INTEGER   IBIT          
      CHARACTER EM(32)*50
C
C==== subprogram =======================================================
C
      LOGICAL   SIBTST
C
C=======================================================================
C
      DATA EM( 1) /'illegal data, LFLAG                               '/
      DATA EM( 2) /'illegal data, dimension of (U), (V), (X)          '/
      DATA EM( 3) /'illegal data, leading dimension of (U), (V), (X)  '/
      DATA EM( 4) /'illegal data, leading dimension of (EIG)          '/
      DATA EM( 5) /'illegal data, number of required eigenpairs       '/
      DATA EM( 6) /'illegal data, Lanczos algorithm block size        '/
      DATA EM( 7) /'illegal data, maximum number of steps             '/
      DATA EM( 8) /'illegal data, number of starting vectors          '/
      DATA EM( 9) /'illegal data, number of eigenpairs provided       '/
      DATA EM(10) /'illegal data, problem type flag                   '/
      DATA EM(11) /'illegal data, spectrum slicing flag               '/
      DATA EM(12) /'illegal data, eigenvectors purification flag      '/
      DATA EM(13) /'illegal data, level of output                     '/
      DATA EM(14) /'illegal data, output file unit                    '/
      DATA EM(15) /'illegal data, LCOMM (MPI or PVM)                  '/
      DATA EM(16) /'illegal data, dimension of ISTOR                  '/
      DATA EM(17) /'illegal data, convergence threshold               '/
      DATA EM(18) /'illegal data, dimension of RSTOR                  '/
      DATA EM(19) /'illegal data on at least one PE                   '/
      DATA EM(20) /'ISTOR(3:14) must be equal on all PEs              '/
      DATA EM(21) /'RSTOR(1:3) must be equal on all PEs               '/
      DATA EM(22) /'not enough space in ISTOR to start eigensolution  '/
      DATA EM(23) /'not enough space in RSTOR to start eigensolution  '/
      DATA EM(24) /'illegal data, number of negative eigenvalues      '/
      DATA EM(25) /'illegal data, entries of V                        '/
      DATA EM(26) /'illegal data, entries of X                        '/
      DATA EM(27) /'failure in computational subinterval              '/
      DATA EM(28) /'file I/O error, blzpack.__.BQ                     '/
      DATA EM(29) /'file I/O error, blzpack.__.BX                     '/
      DATA EM(30) /'file I/O error, blzpack.__.Q                      '/
      DATA EM(31) /'file I/O error, blzpack.__.X                      '/
      DATA EM(32) /'parallel interface error                          '/
C
C**** executable statements ********************************************
C
      IF ( LRERR .GT. 0 ) THEN
         WRITE (LFILE,1000)
         DO 10 IBIT = 1,32
            IF ( SIBTST(IBIT,LRERR) ) WRITE (LFILE,1001) IBIT,EM(IBIT)
   10    CONTINUE
      END IF
C
      RETURN 
C
 1000 FORMAT ()
 1001 FORMAT ('* Error (',I2,'): ',A50)
C
C**** end of LZERRS ****************************************************
C
      END
