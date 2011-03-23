/*----------------------------------------------------------------------
|     Read a Harwell-Boeing matrix.
|     Use readmtc first time to determine sizes of arrays.
|     Read in values on the second call.
|---------------------------------------------------------------------*/

    nrhs = 0;
    tmp0 = 0;

    if ((fp=fopen(pmatrixname,"r"))==NULL) {
        fprintf(STDERR," file %s not found\n",pmatrixname);
        exit(0);
    }
    fclose(fp);

    READMTC(&tmp0,&tmp0,&tmp0,pmatrixname,C.a,C.ja,C.ia,
	    rhs,&nrhs,rhstyp,&n,&nc,&nz,title,key,type,
	    &nrhsix,rhsptr,rhsind,rhsval,&ierr,strlen(pmatrixname),2,72,8,3);

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
    printf("preconditioning Matrix: %s: size (%d,%d), nnz=%d(%4.1lfav.)\n", pmatrixname, n,nc,
	   nz,((double)nz)/n);

    m=1;
    if (nrhs>0) {
       if (rhstyp[1]=='G' || rhstyp[1]=='g') {
	 m++;
       }

       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	 m++;
       }
    }

    rhsptr=NULL;
    rhsind=NULL;
    rhsval=NULL;
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
       rhsptr=(integer *)  MALLOC((size_t)(nrhs+1)*sizeof(integer),"main:rhsptr");
       rhsind=(integer *)  MALLOC((size_t)nrhsix*sizeof(integer),  "main:rhsind");
       rhsval=(FLOAT *)MALLOC((size_t)nrhsix*sizeof(FLOAT),"main:rhsval");
       // no dense right hand side
       m--;
       m*=n*nrhs;
       // in any case we need at least one vector for the r.h.s.
       m+=n;
    }
    else
       m*=n*nrhs;
    C.ia=(integer *)  MALLOC((size_t)(n+1)*sizeof(integer),"main:C.ia");
    C.ja=(integer *)  MALLOC((size_t)nz*sizeof(integer),   "main:C.ja");
    C.a =(FLOAT *)MALLOC((size_t)nz*sizeof(FLOAT), "main:C.a");
    C.nr=n;
    C.nc=n;
    rhs       =(FLOAT *) MALLOC((size_t)m*sizeof(FLOAT),  "main:rhs");
    // advance pointer to reserve space when uncompressing the right hand side
    i=0;
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       i=n;
    
    tmp = 3;
    tmp2 = n;
    tmp3 = nz;

    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       m-=n;
    READMTC(&tmp2,&tmp3,&tmp,pmatrixname,C.a,C.ja,C.ia,
	    rhs+i,&m,rhstyp,&n,&nc,&nz,title,key,type,
	    &nrhsix,rhsptr,rhsind,rhsval,&ierr,strlen(pmatrixname),2,72,8,3);
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m'))
       m+=n;



    if (ierr) {
        fprintf(STDERR," ierr = %d\n",ierr);
        fprintf(STDERR," error in reading the matrix, stop.\n");
        exit(ierr);
    }

    free(rhs);
    if (rhsptr!=NULL) free(rhsptr);
    if (rhsind!=NULL) free(rhsind);
    if (rhsval!=NULL) free(rhsval);
    m=nrhs=tmp=tmp2=tmp3=nz=0;
