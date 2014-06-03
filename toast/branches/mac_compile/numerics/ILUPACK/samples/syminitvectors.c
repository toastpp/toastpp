    // seek for an exact solution
    if (nrhs!=0 && (rhstyp[2]=='X' || rhstyp[2]=='x')) {
       j=nrhs*n;
       if (rhstyp[1]=='G' || rhstyp[1]=='g') 
	 j*=2;
       if (rhstyp[0]=='M' || rhstyp[0]=='m')
	 j-=nrhs*n;
       for (i=0; i<n; i++,j++)
	   sol[i+A.nr*l]=rhs[j+A.nr*l];
    }


    // release part of rhs that may store the uncompressed rhs
    if (nrhs!=0 && (rhstyp[0]=='M' || rhstyp[0]=='m')) {
       rhs-=n;
       rhs-=A.nr*l;
    }

    // set up right hand side
    if (nrhs==0) {
       // prescribe artificial solution
       for (i=0; i<n; i++)  {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
  	   // for testing try test vector as solution if l=0
	   // In this case the preconditioned method must stop
	   // after a few steps.
	   if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
	               param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR)) {
  	      sol[i+A.nr*l]=mytestvector[i];
	   }
	   else
	      sol[i+A.nr*l]=1.0+i*l;
#else
	   if (l==0 && param.flags&DIAGONAL_COMPENSATION &&
	               param.flags&(STATIC_TESTVECTOR|DYNAMIC_TESTVECTOR))
	      sol[i+A.nr*l]=mytestvector[i+2*n];
	   else {
	      sol[i+A.nr*l].r=1.0+i*(l-1);
	      sol[i+A.nr*l].i=0.1*(A.nr-i+l-1);
	   } 
#endif
       } // end for i

       // construct associated right hand side rhs=A*sol
       MYSMATVEC(A,sol+A.nr*l,rhs+A.nr*l);

    } // end if (nrhs==0)
    else {
       if (rhstyp[0]=='M' || rhstyp[0]=='m') {
	  for (i=0; i<n; i++) {
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_ 
	      rhs[i+A.nr*l]=0;
#else
	      rhs[i+A.nr*l].r=rhs[i+A.nr*l].i=0;
#endif
	  }// end for i

	  // uncompress rhs
	  for (i=0; i<rhsptr[l+1]-rhsptr[l]; i++) {
	      j=rhsind[i]-1;
	      rhs[j+A.nr*l]=rhsval[i];
	  } // end for i
       } // end if
    } // end if-else


    // initial solution
    if (nrhs!=0 && (rhstyp[1]=='G' || rhstyp[1]=='g')) {
       j=nrhs*n;
       if (rhstyp[0]=='M' || rhstyp[0]=='m')
	  j=n;
       for (i=0; i<n; i++,j++)
	   sol[i+A.nr*l]=rhs[j+A.nr*l];
    }
    else
       for (i=0; i<n; i++)
#if defined _DOUBLE_REAL_ || defined _SINGLE_REAL_
	   sol[i+A.nr*l]=0;
#else
           sol[i+A.nr*l].r=sol[i+A.nr*l].i=0;
#endif

    // ------   compute initial residual   -----
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
