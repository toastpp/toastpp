#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
SP=(SAMGlevelmat *)&PRE;
#else
SP=(CAMGlevelmat *)&PRE;
#endif
    printf("   final elbow space factor=%8.2f\n",param.elbow+0.005);
    printf("   final condest on level 1=%8.2f\n",param.condest+0.005);
    printf("ILUPACK,   multilevel structure\n");

    if (param.mixedprecision)
       fprintf(fo,"%3d|",SP->nlev);
    else
       fprintf(fo,"%3d|",PRE.nlev);

    if (param.mixedprecision)
       sn=SP;
    else
       next=&PRE;
    nnzU=0;
    tmp=0;
    tmp0=0;
    tmp2=0;

    for (i=1; i<=((param.mixedprecision)?SP->nlev:PRE.nlev); i++) {
        // fill-in LU
        printf("level %3d, block size %7d\n",i,
	       ((param.mixedprecision)?sn->LU.nr:next->LU.nr)); fflush(stdout);
        printf("           E_L=%8.1le, E_S=%8.1le\n",
	       ((param.mixedprecision)?sn->errorL:next->errorL),
	       ((param.mixedprecision)?sn->errorS:next->errorS));
	fflush(stdout);

	if (!(param.flags&DISCARD_MATRIX)) {
	   if (i<((param.mixedprecision)?SP->nlev:PRE.nlev)) {
	      printf("  system size=%6d, fill-in=%8.1f(av.)\n",
		     ((param.mixedprecision)?sn->n:next->n),
		     ((param.mixedprecision)?
		      (sn->A.ia[sn->n]-1)/((double)sn->n):
		      (next->A.ia[next->n]-1)/((double)next->n))
		     );
	      fflush(stdout);
	      if (i>1) 
		 tmp0+=((param.mixedprecision)?
			sn->A.ia[sn->n]-1:
		        next->A.ia[next->n]-1);
	   }
	   else {
	      if (((param.mixedprecision)?sn->LU.ja:next->LU.ja)!=NULL) {
	         printf("  system size=%6d, fill-in=%8.1f(av.)\n",
			((param.mixedprecision)?sn->n:next->n),
			((param.mixedprecision)?
			 (sn->A.ia[sn->n]-1)/((double)sn->n):
			 (next->A.ia[next->n]-1)/((double)next->n)));
		 fflush(stdout);
		 tmp0+=((param.mixedprecision)?
		        sn->A.ia[sn->n]-1:next->A.ia[next->n]-1);
	     }
	   }
	}
	l=nnzU;
	if (i<((param.mixedprecision)?SP->nlev:PRE.nlev) || ((param.mixedprecision)?sn->LU.ja:next->LU.ja)!=NULL) {
	   nnzU+=((param.mixedprecision)?
		  sn->LU.ja[sn->LU.nr-1]-sn->LU.ja[0]+sn->nB:
		  next->LU.ja[next->LU.nr-1]-next->LU.ja[0]+next->nB);
	}
	if (i==((param.mixedprecision)?SP->nlev:PRE.nlev)) {
	   if (((param.mixedprecision)?sn->LU.ja:next->LU.ja)==NULL) {
	      printf("switched to full matrix processing\n");fflush(STDOUT);
	      tmp=-1;
	      j=((param.mixedprecision)?sn->LU.nr:next->LU.nr);
	      nnzU+=(j*(j-1))/2;
	   }
	}
	printf("  local fill-in %7d(%5.1lfav)\n",
	       ((param.mixedprecision)?
		nnzU-l+sn->LU.nr:nnzU-l+next->LU.nr),
	       ((param.mixedprecision)?
		(1.0*(nnzU-l+sn->LU.nr))/sn->LU.nr:
		(1.0*(nnzU-l+next->LU.nr))/next->LU.nr));

	if (i<((param.mixedprecision)?SP->nlev:PRE.nlev)) {
 	   if (param.flags&COARSE_REDUCE) {
	      // fill-in F
	      nnzU+=((param.mixedprecision)?
		     sn->F.ia[sn->F.nr]-1:
		     next->F.ia[next->F.nr]-1);
	      printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,
		     ((param.mixedprecision)?sn->LU.nr:next->LU.nr),
		     ((param.mixedprecision)?sn->F.nc:next->F.nc));
	      printf("  local fill-in F %7d(%5.1lfav pr)\n",
		     ((param.mixedprecision)?
		      sn->F.ia[sn->F.nr]-1:
		      next->F.ia[next->F.nr]-1),
		     ((param.mixedprecision)?
		      (1.0*(sn->F.ia[sn->F.nr]-1))/sn->LU.nr:
		      (1.0*(next->F.ia[next->F.nr]-1))/next->LU.nr));
	   } // end if
	   else {
	      printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,
		     ((param.mixedprecision)?sn->LU.nr:next->LU.nr),
		     ((param.mixedprecision)?sn->F.nc:next->F.nc));
	      fflush(stdout);
	   } // end else
	}
	if (param.mixedprecision)
	   sn=sn->next;
	else
	   next=next->next;
    }
    printf("\ntotal fill-in sum%8d(%5.1lfav)\n",
	   nnzU+n+tmp0,(1.0*(nnzU+n+tmp0))/n);
    printf("fill-in factor:      %5.1lf\n",(1.0*(nnzU+tmp0+n))/nz);
    
    if (tmp) {
       // nnzU-j*(j+1)/2+n-j memory for sparse data structures
       //                    indices (weight 1/3) and values (weight 2/3)
       // j*(j+1)/2          memory for dense data, no indices (weight 2/3)
       printf("memory usage factor: %5.1lf\n",(1.0*(tmp0+nnzU-(j*(j+1))/2+n-j))/nz
	                                     +(j*(j+1))/(3.0*nz));
       fprintf(fo,"%5.1f|",(1.0*(tmp0+nnzU-(j*(j+1))/2+n-j))/nz
                          +(j*(j+1))/(3.0*nz));
    }
    else {
       printf("memory usage factor: %5.1lf\n",(1.0*(nnzU+tmp0+n))/nz);
       fprintf(fo,"%5.1f|",(1.0*(nnzU+tmp0+n))/nz);
    }
    printf("total time: %8.1le [sec]\n",  (double)secnds); 
    printf("            %8.1le [sec]\n\n",(double)ILUPACK_secnds[7]); 
    fprintf(fo,"%7.1le|",(double)secnds);

#ifdef USE_MPILUC
    printf("refined timings for M-ILUPACK multilevel factorization\n"); 
#else
    printf("refined timings for   ILUPACK multilevel factorization\n"); 
#endif
    printf("initial preprocessing:       %8.1le [sec]\n",ILUPACK_secnds[0]); 
    printf("reorderings remaining levels:%8.1le [sec]\n",ILUPACK_secnds[1]); 
#ifdef USE_MPILUC
    printf("MPILDLC(sum over all levels):%8.1le [sec]\n",ILUPACK_secnds[2]); 
#else
    printf("PILDLC (sum over all levels):%8.1le [sec]\n",ILUPACK_secnds[2]); 
#endif
    printf("ILDLC (if used):             %8.1le [sec]\n",ILUPACK_secnds[3]); 
#ifdef USE_LAPACK_DRIVER
    printf("PPTRF (if used):             %8.1le [sec]\n",ILUPACK_secnds[4]); 
#else
    printf("LDLP (if used):              %8.1le [sec]\n",ILUPACK_secnds[4]); 
#endif
    printf("remaining parts:             %8.1le [sec]\n\n",MAX(0.0,(double)secnds
	                                                   -ILUPACK_secnds[0]
							   -ILUPACK_secnds[1]
							   -ILUPACK_secnds[2]
							   -ILUPACK_secnds[3]
							   -ILUPACK_secnds[4]));

    fflush(STDOUT);

    /*
    next=&PRE;
    for (i=1; i<=PRE.nlev; i++) {
      printf("%d,%d\n",next->n,next->nB);
      for (j=0; j<next->n; j++)
	printf("%8d",next->p[j]);
      printf("\n");
      for (j=0; j<next->n; j++)
	printf("%8d",next->invq[j]);
      printf("\n");
        next=next->next;
    }
    fflush(STDOUT);

    printf("multilevel structure\n");
    fflush(STDOUT);
    printf("number of levels: %d\n",PRE.nlev);
    fflush(STDOUT);
    next=&PRE;
    for (i=1; i<=PRE.nlev; i++) {
      printf("total size %d\n",next->n);
      fflush(STDOUT);
      printf("leading block %d\n",next->nB);
      fflush(STDOUT);

      printf("row permutation\n");
      for (j=0; j<next->n; j++)
	printf("%4d",next->p[j]);
      printf("\n");
      fflush(STDOUT);
      printf("inverse column permutation\n");
      for (j=0; j<next->n; j++)
	printf("%4d",next->invq[j]);
      printf("\n");
      fflush(STDOUT);

      if (PRE.nlev==1) {
	printf("row scaling\n");
	for (j=0; j<next->n; j++)
	  printf("%12.4e",next->rowscal[j]);
	printf("\n");
	fflush(STDOUT);
	printf("column scaling\n");
	for (j=0; j<next->n; j++)
	  printf("%12.4e",next->colscal[j]);
	printf("\n");
	fflush(STDOUT);
      }

      printf("(2,1) block E (%d,%d)\n",next->E.nr,next->E.nc);
      fflush(STDOUT);
      for (k=0; k<next->E.nr; k++) {
	printf("%3d: ",k+1);
	for (j=next->E.ia[k]-1; j<next->E.ia[k+1]-1;j++)
	  fprintf(STDOUT,"%12d",next->E.ja[j]);
	printf("\n");
	fflush(STDOUT);
	printf("     ");
	for (j=next->E.ia[k]-1; j<next->E.ia[k+1]-1;j++)
	  fprintf(STDOUT,"%12.4e",next->E.a[j]);
	printf("\n");
	fflush(STDOUT);
      }
      printf("(1,2) block F (%d,%d)\n",next->F.nr,next->F.nc);
      fflush(STDOUT);
      for (k=0; k<next->F.nr; k++) {
	printf("%3d: ",k+1);
	for (j=next->F.ia[k]-1; j<next->F.ia[k+1]-1;j++)
	  fprintf(STDOUT,"%12d",next->F.ja[j]);
	printf("\n");
	fflush(STDOUT);
	printf("     ");
	for (j=next->F.ia[k]-1; j<next->F.ia[k+1]-1;j++)
	  fprintf(STDOUT,"%12.4e",next->F.a[j]);
	printf("\n");
	fflush(STDOUT);
      }

      printf("ILU...\n");
      printf("Diagonal part\n");
      for (k=0; k<next->LU.nr; k++) {
	fprintf(STDOUT,"%12.4e",next->LU.a[k]);
      }
      printf("\n");
      printf("L part\n");
      for (k=0; k<next->LU.nr; k++) {
	printf("col %3d: ",k+1);
	for (j=next->LU.ja[k]-1; j<next->LU.ia[k]-1;j++)
	  fprintf(STDOUT,"%12d",next->LU.ja[j]);
	printf("\n");
	fflush(STDOUT);
	printf("         ");
	for (j=next->LU.ja[k]-1; j<next->LU.ia[k]-1;j++)
	  fprintf(STDOUT,"%12.4e",next->LU.a[j]);
	printf("\n");
	fflush(STDOUT);
      }
      printf("U part\n");
      for (k=0; k<next->LU.nr; k++) {
	printf("row %3d: ",k+1);
	for (j=next->LU.ia[k]-1; j<next->LU.ja[k+1]-1;j++)
	  fprintf(STDOUT,"%12d",next->LU.ja[j]);
	printf("\n");
	fflush(STDOUT);
	printf("         ");
	for (j=next->LU.ia[k]-1; j<next->LU.ja[k+1]-1;j++)
	  fprintf(STDOUT,"%12.4e",next->LU.a[j]);
	printf("\n");
	fflush(STDOUT);
      }

      next=next->next;
    }
    */
