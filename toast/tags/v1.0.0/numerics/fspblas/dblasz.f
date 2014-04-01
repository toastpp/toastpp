      SUBROUTINE   DAXPYIZ   ( NZ, A, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  DAXPYIZ -- INDEXED DOUBLE PRECISION ELEMENTARY          ====
C     ====            VECTOR OPERATION                              ====
C     ====  MS: Modified version for zero-based index lists for use ====
C     ====  with C calling functions                                ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DAXPYIZ ADDS A DOUBLE PRECISION SCALAR MULTIPLE OF 
C             A DOUBLE PRECISION SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         TO  
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF  Y  WHOSE INDICES ARE LISTED IN  INDX 
C         ARE REFERENCED OR MODIFIED.  THE VALUES IN  INDX  MUST BE 
C         DISTINCT TO ALLOW CONSISTENT VECTOR OR PARALLEL EXECUTION.
C
C         ALTHOUGH DISTINCT INDICES WILL ALLOW VECTOR OR PARALLEL
C         EXECUTION, MOST COMPILERS FOR HIGH-PERFORMANCE MACHINES WILL
C         BE UNABLE TO GENERATE BEST POSSIBLE CODE WITHOUT SOME 
C         MODIFICATION, SUCH AS COMPILER DIRECTIVES, TO THIS CODE.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         A       DOUBLE      SCALAR MULTIPLIER OF  X.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  IT IS ASSUMED THAT
C                             THE ELEMENTS IN  INDX  ARE DISTINCT.
C                             *** INDICES ARE ZERO-BASED ***
C
C     UPDATED ...
C
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ON OUTPUT
C                             ONLY THE ELEMENTS CORRESPONDING TO THE
C                             INDICES IN  INDX  HAVE BEEN MODIFIED.
C
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C         
      DOUBLE PRECISION    Y (*), X (*), A
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      IF  ( NZ .LE. 0 )  RETURN
C
      IF  ( A .EQ. 0.0D0 )  RETURN
C
      DO 10 I = 1, NZ
          Y(INDX(I)+1)  = Y(INDX(I)+1) + A * X(I)
   10 CONTINUE
C
      RETURN
      END
      DOUBLE PRECISION FUNCTION   DDOTIZ   ( NZ, X, INDX, Y )
C
C     ==================================================================
C     ==================================================================
C     ====  DDOTIZ -- DOUBLE PRECISION INDEXED DOT PRODUCT          ====
C     ====  MS: Modified version for zero-based index lists for use ====
C     ====  with C calling functions                                ====
C     ==================================================================
C     ==================================================================
C
C     PURPOSE
C     -------
C
C         DDOTIZ COMPUTES THE VECTOR INNER PRODUCT OF 
C             A DOUBLE PRECISION SPARSE VECTOR  X
C             STORED IN COMPRESSED FORM  (X,INDX) 
C         WITH 
C             A DOUBLE PRECISION VECTOR  Y  IN FULL STORAGE FORM.
C
C         ONLY THE ELEMENTS OF Y WHOSE INDICES ARE LISTED IN INDX 
C         ARE REFERENCED.
C
C     ARGUMENTS
C     ---------
C
C     INPUT ...
C
C         NZ      INTEGER     NUMBER OF ELEMENTS IN THE COMPRESSED FORM.
C         X       DOUBLE      ARRAY CONTAINING THE VALUES OF THE 
C                             COMPRESSED FORM.
C         INDX    INTEGER     ARRAY CONTAINING THE INDICES OF THE
C                             COMPRESSED FORM.  
C                             *** INDICES ARE ZERO-BASED ***
C         Y       DOUBLE      ARRAY, ON INPUT, WHICH CONTAINS THE 
C                             VECTOR  Y  IN FULL STORAGE FORM.  ONLY
C                             THE ELEMENTS  CORRESPONDING TO THE
C                             INDICES IN  INDX  WILL BE ACCESSED.
C
C     OUTPUT ...
C
C         DDOTIZ  DOUBLE      DOUBLE PRECISION FUNCTION VALUE EQUAL TO
C                             THE VECTOR INNER PRODUCT.  
C                             IF  NZ .LE. 0  DDOTI IS SET TO ZERO.
C
C     SPARSE BASIC LINEAR ALGEBRA SUBPROGRAM
C
C     FORTRAN VERSION WRITTEN OCTOBER 1984
C     ROGER G GRIMES, BOEING COMPUTER SERVICES
C
C     ==================================================================
C
C     -------------
C     ... ARGUMENTS
C     -------------
C
      INTEGER             NZ, INDX (*)
C
      DOUBLE PRECISION    X (*), Y (*)
C
C     -------------------
C     ... LOCAL VARIABLES
C     -------------------
C
      INTEGER             I
C
C     ==================================================================
C
      DDOTIZ = 0.0D0
      IF  ( NZ .LE. 0 )  RETURN
C
      DO 10 I = 1, NZ
          DDOTIZ = DDOTIZ  +  X(I) * Y(INDX(I)+1)
   10 CONTINUE
C
      RETURN
      END
