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
    ierr=0;
    tmp0=0;
    tmp2=0;
    tmp3=0;
    
    if (param.ind!=NULL && (param.flags&SADDLE_POINT)) {
       param.nibuff=MAX(2*A.nc,param.nibuff);
       param.ibuff =(integer *)REALLOC(param.ibuff,param.nibuff*sizeof(integer),"mainsym:param.ibuff");
       ptr=param.ibuff;
       k=0;
       for (j=0; j<A.nc; j++) {
	   ptr[j]=param.ind[j];
	   if (ptr[j]<0) k++;
       }
       printf("initial saddle point block structure (%7d,%7d)\n",A.nc-k,k);
    }

    for (i=1; i<=((param.mixedprecision)?SP->nlev:PRE.nlev); i++) {
        // fill-in LU
        printf("level %3d, block size %7d\n",i,
	       ((param.mixedprecision)?sn->LU.nr:next->LU.nr)); fflush(stdout);
        printf("           E_L=%8.1le, E_S=%8.1le\n",
	       ((param.mixedprecision)?sn->errorL:next->errorL),
	       ((param.mixedprecision)?sn->errorS:next->errorS));
	fflush(stdout);
	if (param.ind!=NULL && (param.flags&SADDLE_POINT)) {
	   // permute index vector
	   for (j=0; j<((param.mixedprecision)?sn->n:next->n); j++) 
	       ptr[j+((param.mixedprecision)?sn->n:next->n)]=ptr[((param.mixedprecision)?sn->p[j]:next->p[j])-1];
	   for (j=0; j<((param.mixedprecision)?sn->n:next->n); j++) 
	       ptr[j]=ptr[j+((param.mixedprecision)?sn->n:next->n)];
	   k=0;
	   for (j=0; j<((param.mixedprecision)?sn->nB:next->nB); j++) 
	       if (ptr[j]<0) k++;
	   printf("           saddle point block structure (%6d,%6d),",next->nB-k,k);
	   k=0;
	   for (; j<((param.mixedprecision)?sn->n:next->n); j++) 
	       if (ptr[j]<0) k++;
	   printf("(%6d,%6d)\n",((param.mixedprecision)?sn->n-sn->nB:next->n-next->nB)-k,k);
	   ptr+=((param.mixedprecision)?sn->nB:next->nB);
	}
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
	   nnzU+=(param.mixedprecision)?
	         sn->LU.nnz:next->LU.nnz;

	   k=0;
	   for (m=0; m<((param.mixedprecision)?sn->LU.nr:next->LU.nr); m++) {
	       
	       if (((param.mixedprecision)?sn->isblock:next->isblock)) {
		  if (((param.mixedprecision)?sn->blocksize[m]:next->blocksize[m])>0)
		     k+=1;   
	       }
	       else {
	          if (((param.mixedprecision)?sn->LU.ja[sn->LU.nr+1+m]:next->LU.ja[next->LU.nr+1+m])>0) 
		    k+=2;
	       }
	   } // end for m
	   ierr+=k;
	   tmp3+=((param.mixedprecision)?sn->LU.nr:next->LU.nr);
	   if (((param.mixedprecision)?sn->isblock:next->isblock)) {
	      printf("           number of diagonal blocks %6d\n",k); 
	   }
	   else {
	      printf("           2x2 pivots %5.1f%%\n",(100.0*k)/((param.mixedprecision)?sn->LU.nr:next->LU.nr)); 
	   }
	   fflush(stdout);
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
	   }
	   else {
	      printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,
		     ((param.mixedprecision)?sn->LU.nr:next->LU.nr),
		     ((param.mixedprecision)?sn->F.nc:next->F.nc));
	      fflush(stdout);
	   }
	}
	if (param.mixedprecision)
	   sn=sn->next;
	else
	   next=next->next;
    }
    printf("\ntotal fill-in sum%8d(%5.1lfav)\n",
	   nnzU+tmp0,(1.0*(nnzU+tmp0))/n);
    printf("fill-in factor:      %5.1lf\n",(1.0*nnzU+tmp0)/nz);
    if (((param.mixedprecision)?SP->isblock:PRE.isblock)) {
       printf("total number of diagonal blocks %6d\n",ierr); 
    }
    else
       printf("total number of sparse  2x2 pivots %5.1f%%\n",(100.0*ierr)/tmp3); 
    fflush(stdout);
    ierr=tmp3=0;

    if (tmp) {
       // nnzU-j*(j+1)/2+n-j memory for sparse data structures
       //                    indices (weight 1/3) and values (weight 2/3)
       // j*(j+1)/2          memory for dense data, no indices (weight 2/3)
       printf("memory usage factor: %5.1lf\n",(1.0*(tmp0+nnzU-(j*(j+1))/2-j))/nz
	      +(j*(j+1))/(3.0*nz));
       fprintf(fo,"%5.1f|",(1.0*(tmp0+nnzU-(j*(j+1))/2-j))/nz
	       +(j*(j+1))/(3.0*nz));
    }
    else {
       printf("memory usage factor: %5.1lf\n",(1.0*(nnzU+tmp0))/nz);
       fprintf(fo,"%5.1f|",(1.0*(nnzU+tmp0))/nz);
    }
    printf("total time: %8.1le [sec]\n",  (double)secnds); 
    printf("            %8.1le [sec]\n\n",(double)ILUPACK_secnds[7]); 
    fprintf(fo,"%7.1le|",(double)secnds);

    printf("refined timings for   ILUPACK multilevel factorization\n"); 
    printf("initial preprocessing:         %8.1le [sec]\n",ILUPACK_secnds[0]); 
    printf("reorderings remaining levels:  %8.1le [sec]\n",ILUPACK_secnds[1]); 
    printf("SYMPILUC (sum over all levels):%8.1le [sec]\n",ILUPACK_secnds[2]); 
    printf("SYMILUC (if used):             %8.1le [sec]\n",ILUPACK_secnds[3]); 
    printf("SPTRF, LAPACK (if used):       %8.1le [sec]\n",ILUPACK_secnds[4]); 
    printf("remaining parts:               %8.1le [sec]\n\n",MAX(0.0,(double)secnds
	                                                   -ILUPACK_secnds[0]
							   -ILUPACK_secnds[1]
							   -ILUPACK_secnds[2]
							   -ILUPACK_secnds[3]
							   -ILUPACK_secnds[4]));

    fflush(STDOUT);

    /*
    next=&PRE;
    for (i=1; i<=nlev; i++) {
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
    printf("number of levels: %d\n",nlev);
    fflush(STDOUT);
    next=&PRE;
    for (i=1; i<=nlev; i++) {
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

      if (nlev==1) {
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

#ifdef PRINT_INFO
    next=&PRE;
    for (i=1; i<=PRE.nlev; i++) {
        // fill-in LU
        printf("level %3d, block size %7d\n",i,next->LU.nr); fflush(stdout);
	if (i<PRE.nlev || next->LU.ja!=NULL) {
	   printf("U-factor");
	   printf("\n");fflush(stdout);
	   for (l=0; l<next->LU.nr; ) {
	       if (next->LU.ja[next->LU.nr+1+l]==0){
		  for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		      printf("%8d",next->LU.ja[j-1]);
		  }
		  printf("\n");fflush(stdout);
		  for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		      printf("%8.1le",next->LU.a[j-1]);
		  }
		  l++;
	       }
	       else {
		  for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		      printf("%8d",next->LU.ja[j-1]);
		  }
		  printf("\n");fflush(stdout);
		  for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		      printf("%8.1le",
			     next->LU.a[next->LU.ja[l]+2*(j-next->LU.ja[l])-1]);
		  }
		  printf("\n");fflush(stdout);
		  for (j=next->LU.ja[l];j<next->LU.ja[l+1]; j++) {
		      printf("%8.1le",
			     next->LU.a[next->LU.ja[l]+2*(j-next->LU.ja[l])]);
		  }
		  l+=2;
	       }
	       printf("\n");fflush(stdout);
	   }

	   printf("Block diagonal factor\n");
	   for (l=0; l<next->LU.nr;) {
	       if (next->LU.ja[next->LU.nr+1+l]==0){
		  printf("%8.1le",next->LU.a[l]);
		  l++;
	       }
	       else {
		  printf("%8.1le%8.1le",next->LU.a[l],
			 next->LU.a[next->LU.nr+1+l]);
		  l+=2;
	       }
	   }
	   printf("\n");fflush(stdout);
	   for (l=0; l<next->LU.nr; ) {
	       if (next->LU.ja[next->LU.nr+1+l]==0) {
		 printf("        ");
		 l++;
	       }
	       else {
		 printf("%8.1le%8.1le",next->LU.a[next->LU.nr+1+l],
			next->LU.a[l+1]);
		 l+=2;
	       }
	   }
	   printf("\n");fflush(stdout);
	   
	}
	if (i==PRE.nlev) {
	   if (next->LU.ja==NULL) {
	      printf("switched to full matrix processing\n");fflush(STDOUT);
	   }
	}

	if (i<PRE.nlev) {
	   // fill-in F
	   nnzU+=next->F.ia[next->F.nr]-1;
	   printf("level %3d->%3d, block size (%7d,%7d)\n",i,i+1,next->LU.nr,
		  next->F.nc);
	   printf("  local fill-in F %7d(%5.1lfav pr)\n",
		  next->F.ia[next->F.nr]-1,
		  (1.0*(next->F.ia[next->F.nr]-1))/next->LU.nr);
	}
	next=next->next;
    }
#endif
