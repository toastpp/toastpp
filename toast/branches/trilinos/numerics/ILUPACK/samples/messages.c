    if (param.matching) {
       if (!strncmp("yes",param.FCpart,3)) {
          if (strncmp("none",param.typetv,4)) {
	     // 1. maximum weight matching
 	     // 2. fine/coarse grid partitioning
	     // 3. ... based on static or dynamic test vector
	     if (!strncmp("amd",param.ordering,3))
	        fprintf(fo,"mw/ad/FCv|");
	     else if (!strncmp("metisn",param.ordering,6))
	        fprintf(fo,"mw/mn/FCv|");
	     else if (!strncmp("metise",param.ordering,6))
	        fprintf(fo,"mw/me/FCv|");
	     else if (!strncmp("rcm",param.ordering,3))
	        fprintf(fo,"mw/rc/FCv|");
	     else if (!strncmp("amf",param.ordering,3))
	        fprintf(fo,"mw/af/FCv|");
	     else if (!strncmp("mmd",param.ordering,3))
	        fprintf(fo,"mw/md/FCv|");
	     else if (!strncmp("indset",param.ordering,6))
	        fprintf(fo,"mw/is/FCv|");
	     else if (!strncmp("pq",param.ordering,2))
	        fprintf(fo,"  /pq/   |");
	     else 
	        fprintf(fo,"  /  /   |");
	  }
	  // 3. no test vector
	  else {
	     // 1. maximum weight matching
 	     // 2. fine/coarse grid partitioning
	     // 3. no test vector
	     if (!strncmp("amd",param.ordering,3))
	        fprintf(fo,"mw/ad/FC |");
	     else if (!strncmp("metisn",param.ordering,6))
	        fprintf(fo,"mw/mn/FC |");
	     else if (!strncmp("metise",param.ordering,6))
	        fprintf(fo,"mw/me/FC |");
	     else if (!strncmp("rcm",param.ordering,3))
	        fprintf(fo,"mw/rc/FC |");
	     else if (!strncmp("amf",param.ordering,3))
	        fprintf(fo,"mw/af/FC |");
	     else if (!strncmp("mmd",param.ordering,3))
	        fprintf(fo,"mw/md/FC |");
	     else if (!strncmp("indset",param.ordering,6))
	        fprintf(fo,"mw/is/FC |");
	     else if (!strncmp("pq",param.ordering,2))
	        fprintf(fo,"  /pq/   |");
	     else 
	        fprintf(fo,"  /  /   |");
	  }
       }	  
       // 2. no FCpart
       else {
	  // 1. maximum weight matching
	  // 2. no fine/coarse grid partitioning
	  // 3. test vector is not used in this case, anyway
	  if (!strncmp("amd",param.ordering,3))
	     fprintf(fo,"mw/ad/   |");
	  else if (!strncmp("metisn",param.ordering,6))
	     fprintf(fo,"mw/mn/   |");
	  else if (!strncmp("metise",param.ordering,6))
	     fprintf(fo,"mw/me/   |");
	  else if (!strncmp("rcm",param.ordering,3))
	     fprintf(fo,"mw/rc/   |");
	  else if (!strncmp("amf",param.ordering,3))
	     fprintf(fo,"mw/af/   |");
	  else if (!strncmp("mmd",param.ordering,3))
	     fprintf(fo,"mw/md/   |");
	  else if (!strncmp("indset",param.ordering,6))
	     fprintf(fo,"mw/is/   |");
	  else if (!strncmp("pq",param.ordering,2))
	     fprintf(fo,"  /pq/   |");
	  else 
	     fprintf(fo,"  /  /   |");
       }
    }
    // 1. no matching
    else {
       if (!strncmp("yes",param.FCpart,3)) {
          if (strncmp("none",param.typetv,4)) {
	     // 1. no matching
 	     // 2. fine/coarse grid partitioning
	     // 3. ... based on static or dynamic test vector
	     if (!strncmp("amd",param.ordering,3))
	        fprintf(fo,"  /ad/FCv|");
	     else if (!strncmp("metisn",param.ordering,6))
	        fprintf(fo,"  /mn/FCv|");
	     else if (!strncmp("metise",param.ordering,6))
	        fprintf(fo,"  /me/FCv|");
	     else if (!strncmp("rcm",param.ordering,3))
	        fprintf(fo,"  /rc/FCv|");
	     else if (!strncmp("amf",param.ordering,3))
	        fprintf(fo,"  /af/FCv|");
	     else if (!strncmp("mmd",param.ordering,3))
	        fprintf(fo,"  /md/FCv|");
	     else if (!strncmp("indset",param.ordering,6))
	        fprintf(fo,"  /is/FCv|");
	     else if (!strncmp("pq",param.ordering,2))
	        fprintf(fo,"  /pq/   |");
	     else 
	        fprintf(fo,"  /  /   |");
	  }
	  // 3. no test vector
	  else {
	     // 1. no matching
 	     // 2. fine/coarse grid partitioning
	     // 3. no test vector
	     if (!strncmp("amd",param.ordering,3))
	        fprintf(fo,"  /ad/FC |");
	     else if (!strncmp("metisn",param.ordering,6))
	        fprintf(fo,"  /mn/FC |");
	     else if (!strncmp("metise",param.ordering,6))
	        fprintf(fo,"  /me/FC |");
	     else if (!strncmp("rcm",param.ordering,3))
	        fprintf(fo,"  /rc/FC |");
	     else if (!strncmp("amf",param.ordering,3))
	        fprintf(fo,"  /af/FC |");
	     else if (!strncmp("mmd",param.ordering,3))
	        fprintf(fo,"  /md/FC |");
	     else if (!strncmp("indset",param.ordering,6))
	        fprintf(fo,"  /is/FC |");
	     else if (!strncmp("pq",param.ordering,2))
	        fprintf(fo,"  /pq/   |");
	     else 
	        fprintf(fo,"  /  /   |");
	  }
       }	  
       // 2. no FCpart
       else {
	  // 1. no matching
	  // 2. no fine/coarse grid partitioning
	  // 3. test vector is not used in this case, anyway
	  if (!strncmp("amd",param.ordering,3))
	     fprintf(fo,"  /ad/   |");
	  else if (!strncmp("metisn",param.ordering,6))
	     fprintf(fo,"  /mn/   |");
	  else if (!strncmp("metise",param.ordering,6))
	     fprintf(fo,"  /me/   |");
	  else if (!strncmp("rcm",param.ordering,3))
	     fprintf(fo,"  /rc/   |");
	  else if (!strncmp("amf",param.ordering,3))
	     fprintf(fo,"  /af/   |");
	  else if (!strncmp("mmd",param.ordering,3))
	     fprintf(fo,"  /md/   |");
	  else if (!strncmp("indset",param.ordering,6))
	     fprintf(fo,"  /is/   |");
	  else if (!strncmp("pq",param.ordering,2))
	     fprintf(fo,"  /pq/   |");
	  else 
	     fprintf(fo,"  /  /   |");
       }
    }


    printf("ILUPACK PARAMETERS:\n");
    printf("   droptol=%8.1le\n",          param.droptol);
    printf("   condest=%8.1le\n",          param.condest);
    printf("   elbow space factor=%5.1f\n",param.elbow);
    if (!strncmp("ilu",param.typecoarse,3))
       printf("   simple ILU-type coarse grid system\n");
    else if (!strncmp("amg",param.typecoarse,3))
       printf("   medium AMG-type coarse grid system\n");

    if (!strncmp("static",param.typetv,6))
       printf("   diagonal compensation using a prescribed static test vector\n");
    else if (!strncmp("none",param.typetv,4))
       printf("   no diagonal compensation\n");
    else
       printf("   diagonal compensation using a dynamically generated test vector\n");

    printf("   type of preconditioner: ");
    if (!strncmp("ilu",param.amg,3))
       printf("multilevel ILU\n");
    else if (!strncmp("amli",param.amg,4)) {
       printf("AMLI-like\n");
       printf("   number of recursive calls: %d\n", param.ncoarse);
    }
    else if (!strncmp("mg",param.amg,2)) {
       printf("multigrid\n");
       printf("   number of pre-smoothing steps:  %d\n", param.npresmoothing);
       printf("   number of recursive calls:      %d\n", param.ncoarse);
       printf("   number of post-smoothing steps: %d\n", param.npostsmoothing);
       printf("   type of pre-smoother:  ");
       if (!strncmp("gsf",param.presmoother,3))
	  printf("Gauss-Seidel forward\n");
       else if (!strncmp("gsb",param.presmoother,3))
	  printf("Gauss-Seidel backward\n");
       else if (!strncmp("j",param.presmoother,1))
	  printf("(damped) Jacobi\n");
       else if (!strncmp("ilu",param.presmoother,3))
	  printf("ILU\n");
       else
	  printf("custom\n");
       printf("   type of post-smoother: ");
       if (!strncmp("gsf",param.postsmoother,3))
	  printf("Gauss-Seidel forward\n");
       else if (!strncmp("gsb",param.postsmoother,3))
	  printf("Gauss-Seidel backward\n");
       else if (!strncmp("j",param.postsmoother,1))
	  printf("(damped) Jacobi\n");
       else if (!strncmp("ilu",param.postsmoother,3))
	  printf("ILU\n");
       else
	  printf("custom\n");
    }



    if (param.matching)
       printf("   maximum weight matching prior to reorderings\n");
    else
       printf("   NO maximum weight matching\n");

#if !defined _SINGLE_REAL_ && !defined _SINGLE_COMPLEX_
    if (param.mixedprecision!=0)
       printf("   use mixed precision\n");
    else
       printf("   NO mixed precision\n");
#endif

    if (!strncmp("amd",param.ordering,3))
       printf("   reorder systems based on AMD\n");
    else if (!strncmp("metisn",param.ordering,6))
       printf("   reorder systems based on METIS nested dissection by nodes\n");
    else if (!strncmp("metise",param.ordering,6))
       printf("   reorder systems based on METIS nested dissection by edges\n");
    else if (!strncmp("rcm",param.ordering,3))
       printf("   reorder systems based on RCM\n");
    else if (!strncmp("amf",param.ordering,3))
       printf("   reorder systems based on HALOAMD\n");
    else if (!strncmp("mmd",param.ordering,3))
       printf("   reorder systems based on MMD\n");
    else if (!strncmp("indset",param.ordering,6))
       printf("   reorder systems based on independent sets\n");
    else if (!strncmp("pq",param.ordering,2))
       printf("   reorder systems based on ddPQ\n");
    else
       printf("   custom reordering strategy\n");

    if (!strncmp("yes",param.FCpart,3))
       if (strncmp("none",param.typetv,4))
	  printf("   a priori fine/coarse grid partitioning using test vector\n");
       else
	  printf("   a priori fine/coarse grid partitioning\n");
    else
       printf("   NO a priori fine/coarse grid partitioning\n");

    if (!strncmp("ilu",param.typecoarse,3))
       fprintf(fo,"ILUPACK(S)|");
    else if (!strncmp("amg",param.typecoarse,3))
       fprintf(fo,"ILUPACK(M)|");
    else
       fprintf(fo,"ILUPACK   |");

