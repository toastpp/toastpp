    // -------   compute final residual   ------
    MYSMATVEC(A,sol+A.nr*l,param.dbuff);

    for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
        param.dbuff[i]-=rhs[i+A.nr*l];
#else
	param.dbuff[i].r-=rhs[i+A.nr*l].r;
	param.dbuff[i].i-=rhs[i+A.nr*l].i;
#endif
    } // end for i

    i=1;
    val=NRM(&n, param.dbuff, &i);
    // -----------------------------------------
    printf("current: %8.1le\n",val);


    // release part of rhs that may store the uncompressed rhs
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
       rhs+=n;
       rhs+=A.nr*l;
    }


    if (nrhs==0) {
       for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
	               param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR))
  	      sol[i+A.nr*l]-=mytestvector[i];
	   else
	      sol[i+A.nr*l]-=1.0+i*l;
#else
	   if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
	               param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR)) {
  	      sol[i+A.nr*l].r-=mytestvector[i].r;
  	      sol[i+A.nr*l].i-=mytestvector[i].i;
	   }
	   else {
	      sol[i+A.nr*l].r-=1.0+i*(l-1);
	      sol[i+A.nr*l].i-=0.1*(A.nr-i+l-1);
	   }
#endif
       }
       i=1;
       val=NRM(&n,sol+A.nr*l,&i);
       for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
	               param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR))
  	      sol[i+A.nr*l]+=mytestvector[i];
	   else
	      sol[i+A.nr*l]+=1.0+i*l;
#else
	   if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
	               param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR)) {
  	      sol[i+A.nr*l].r+=mytestvector[i].r;
  	      sol[i+A.nr*l].i+=mytestvector[i].i;
	   }
	   else {
	      sol[i+A.nr*l].r+=1.0+i*(l-1);
	      sol[i+A.nr*l].i+=0.1*(A.nr-i+l-1);
	   }
#endif
       }
       i=1;
       vb=NRM(&n,sol+A.nr*l,&i);
       // if (nrhs!=0 && ((rhstyp[2]=='X' || rhstyp[2]=='x') || (rhstyp[0]!='M' && rhstyp[0]!='m')))
       printf("rel. error in the solution: %8.1le\n\n",val/vb);
       //else printf("\n");
    }
    else {
       // exact solution provided
       m=1;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') m++;
       
       if (rhstyp[0]=='M' || rhstyp[0]=='m')
	  rhs+=n+(m-1)*n*nrhs;
       else
	  rhs+=n*m*nrhs;
       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	      sol[i+A.nr*l]-=rhs[i+A.nr*l];
#else
	      sol[i+A.nr*l].r-=rhs[i+A.nr*l].r;
	      sol[i+A.nr*l].i-=rhs[i+A.nr*l].i;
#endif
	  } // end for i
       } // end if

       i=1;
       if (rhstyp[2]=='X' || rhstyp[2]=='x') {
	  printf("rel. error in the solution: %8.1le\n\n",NRM(&n,sol+A.nr*l,&i)/NRM(&n,rhs+A.nr*l,&i));
	  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	      sol[i+A.nr*l]+=rhs[i+A.nr*l];
#else
	      sol[i+A.nr*l].r+=rhs[i+A.nr*l].r;
	      sol[i+A.nr*l].i+=rhs[i+A.nr*l].i;
#endif
	  } // end for i
       }
       else printf("\n");

       // re-adjust rhs
       if (nrhs!=0) {
	  m=1;
	  if (rhstyp[1]=='G' || rhstyp[1]=='g') m++;
	  if (rhstyp[0]=='M' || rhstyp[0]=='m')
	     rhs-=n+(m-1)*n*nrhs;
	  else
	     rhs-=n*m*nrhs;
       }
    }
