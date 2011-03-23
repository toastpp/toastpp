/*----------------------------------------------------------------------
|     Read a Harwell-Boeing matrix.
|     Use readmtc first time to determine sizes of arrays.
|     Read in values on the second call.
|---------------------------------------------------------------------*/

    nrhs = 0;
    tmp0 = 0;

    if ((fp=fopen(fname,"r"))==NULL) {
        fprintf(STDERR," file %s not found\n",fname);
        exit(0);
    }
    fclose(fp);

    READMTC(&tmp0,&tmp0,&tmp0,fname,fullmat.a,fullmat.ja,fullmat.ia,
	    rhs,&nrhs,rhstyp,&n,&nc,&nz,title,key,type,
	    &nrhsix,rhsptr,rhsind,rhsval,&ierr,fnamelen,2,72,8,3);
    // if a right hand side is given, then use these
    if (nrhs>0)
       mynrhs=nrhs;

    if (ierr) {
        fprintf(STDERR," ierr = %d\n",ierr);
        fprintf(STDERR," error in reading the matrix, stop.\n");
	switch(ierr) {
	case 1:
	  fprintf(STDERR,"too many columns\n");
	  break;  
	case 2:
	  fprintf(STDERR,"too many nonzeros\n");
	  break;  
	case 3:
	  fprintf(STDERR,"too many columns and nonzeros\n");
	  break;  
	case 4:
	  fprintf(STDERR,"right hand side has incompatible type\n");
	  break;  
	case 5:
	  fprintf(STDERR,"too many right hand side entries\n");
	  break;  
	case 6:
	  fprintf(STDERR,"wrong type (real/complex)\n");
	  break;  
	}
        exit(ierr);
    }
    printf("Matrix: %s: size (%d,%d), nnz=%d(%4.1lfav.)\n", fname, n,nc,
	   nz,((double)nz)/n);

    if (fname[fnamelen-1]=='5')
      fo = fopen("out_mc64","aw");
    else
      fo = fopen("out_normal","aw");

    fprintf(fo,"%s|%7.1e|%7.1e|",foname,DROP_TOL,CONDEST);

    m=1;
    if (nrhs>0) {
       printf ("Number of right hand sides supplied: %d \n", nrhs) ;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') {
	 m++;
	 printf ("Initial solution(s) offered\n") ;
       }
       else
	 printf ("\n") ;

       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	 m++;
	 printf ("Exact solution(s) provided\n") ;
       }
       else
	 printf ("\n") ;
    }
    else 
      printf("\n\n\n");
    printf("\n");

    rhsptr=NULL;
    rhsind=NULL;
    rhsval=NULL;
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
       rhsptr=(integer *)  MALLOC((size_t)(nrhs+1)*sizeof(integer),"main:rhsptr");
       rhsind=(integer *)  MALLOC((size_t)nrhsix*sizeof(integer),  "main:rhsind");
       rhsval=(FLOAT *)MALLOC((size_t)nrhsix*sizeof(FLOAT),"main:rhsval");
       // no dense right hand side
       m--;
       m*=n*MAX(mynrhs,nrhs);
       // in any case we need at least one vector for the r.h.s.
       m+=n;
    }
    else
       m*=n*MAX(mynrhs,nrhs);
    fullmat.ia=(integer *)  MALLOC((size_t)(n+1)*sizeof(integer),"main:fullmat.ia");
    fullmat.ja=(integer *)  MALLOC((size_t)nz*sizeof(integer),   "main:fullmat.ja");
    fullmat.a =(FLOAT *)MALLOC((size_t)nz*sizeof(FLOAT), "main:fullmat.a");
    fullmat.nr=n;
    fullmat.nc=n;
    rhs       =(FLOAT *)MALLOC((size_t)m*sizeof(FLOAT),  "main:rhs");
    // advance pointer to reserve space when uncompressing the right hand side
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       rhs+=n;
    
    tmp = 3;
    tmp2 = n;
    tmp3 = nz;

    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       m-=n;
    READMTC(&tmp2,&tmp3,&tmp,fname,fullmat.a,fullmat.ja,fullmat.ia,
	    rhs,&m,rhstyp,&n,&nc,&nz,title,key,type,
	    &nrhsix,rhsptr,rhsind,rhsval,&ierr,fnamelen,2,72,8,3);
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       m+=n;



    if (ierr) {
        fprintf(STDERR," ierr = %d\n",ierr);
        fprintf(STDERR," error in reading the matrix, stop.\n");
	fprintf(fo,"out of memory\n");fclose(fo); 
        exit(ierr);
    }
    
/*----------------------------------------------------------------------
|     Convert the matrix from csc to csr format.  First, allocate
|     space in (fullmat).  After conversion, free space occupied by
|     initial format.
|---------------------------------------------------------------------*/

    A.nr=A.nc=n;
    if (transpose) {
       A.ia=fullmat.ia; 
       A.ja=fullmat.ja;
       A.a =fullmat.a;
    }
    else {
       A.ia=(integer *)  MALLOC((size_t)(n+1)*sizeof(integer),"main:A.ia");
       A.ja=(integer *)  MALLOC((size_t)nz*sizeof(integer),   "main:A.ja");
       A.a =(FLOAT *)MALLOC((size_t)nz*sizeof(FLOAT), "main:A.a");
       
       tmp = 1;
       CSRCSC(&n,&tmp,&tmp,fullmat.a,fullmat.ja,fullmat.ia,A.a,A.ja,A.ia);
       
       free(fullmat.a);
       free(fullmat.ja);
       free(fullmat.ia);
    }

