/***************************************************************************
 * genmat_angint_sdm_3D_PN.cc         Simon Arridge            29.07.10    *
 *                                                                         *
 ***************************************************************************/

#ifdef __BORLANDC__
#include <strstrea.h>
#include <conio.h>
#include <process.h>

typedef unsigned pid_t;
#else
#include <sstream>
#include <unistd.h>

#endif
#include <time.h>
#include <stream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include <toast.h>
#include <rte3D.h>
#include <sphints.h>

#define USE_INTONSPHERE
    using namespace toast;

const double Y00N  = sqrt(1/M_PI)/2;
const double Y1m1N = sqrt(3/(2*M_PI))/2;
const double Y10N  = sqrt(3/M_PI)/2;
const double Y1p1N = -Y1m1N;

const double Y22[3][3] = 
  {{ Y1m1N*CGmep(1,(-1)), Y1m1N*CGmep(1,0), Y1m1N*CGmep(1,1)},
   { Y10N*CGzp(1,(-1)),   Y10N*CGzp(1,0), Y10N*CGzp(1,1)},
   { Y1p1N*CGpep(1,-1), Y1p1N*CGpep(1,0), Y1p1N*CGpep(1,1)}
  };         // coefficients of Y2 harmonics as products of Y1
 
void genmat_angint_sdm_3D_PN(RCompRowMatrix& Aintscsc,  RCompRowMatrix& Aintscss,  RCompRowMatrix& Aintscc,  RCompRowMatrix& Aintssss,  RCompRowMatrix& Aintssc,  RCompRowMatrix& Aintcc, RCompRowMatrix& Anvec_sc, RCompRowMatrix& Anvec_ss, RCompRowMatrix& Anvec_c, CCompRowMatrix& RY, const int sporder)
{
  // this function returns (L+1)^2 X (L+1) matrices of angular shape integrals.

  const int nso = sporder+1;
  const int nsp = nso*nso;
  int i,j;
  int nnz ; // number of nonzeros
  int rp = 0, ci = 0;
  int *rowptr ;
  int *colidx ;
  toast::complex *val;
  double *rval;
  CCompRowMatrix tmp(nsp,nsp);
  toast::complex val00;

  Aintscsc.New(nsp,nsp);// integrals of products of spherical harmonics with xx

  Aintscss.New(nsp,nsp);// integrals of products of spherical harmonics with xy

  Aintscc.New(nsp,nsp);// integrals of products of spherical harmonics with xz

  Aintssss.New(nsp,nsp);// integrals of products of spherical harmonics with yy

  Aintssc.New(nsp,nsp);// integrals of products of spherical harmonics with yz

  Aintcc.New(nsp,nsp);  // integrals of products of spherical harmonics with zz

  Anvec_sc.New(nsp,nsp);  // this one is just Aintsc

  Anvec_ss.New(nsp,nsp);  // and this one is just Aintss

  Anvec_c.New(nsp,nsp);  // and this one is just Aintc

  //  RDenseMatrix Y22mat(3,3);
  //  Y22mat = Y22;
  cerr << "Y22 coefficients :\n" ;
  for (i = 0; i < 3; i++ ){
    for(j = 0; j < 3; j++)
      cerr << Y22[i][j] << ' ';
    cerr << endl;
  }
  /*--------- Identity Matrix  -----*/
  CCompRowMatrix idnspC(nsp,nsp);
    nnz = nsp;
    rowptr = new int [nsp+1];
    colidx = new int [nnz];
    val = new toast::complex [nnz];
    for (i = 0; i < nsp; i++){
 //      cerr << "i " << i;
      rowptr[i] = rp++;
      val[ci] = 1.0;
      colidx[ci++] = i;
 //      cerr << " " << rowptr[i] << " " << colidx[i] << " " << val[i] <<endl;
    }
    rowptr[nsp] = rp;
    idnspC.Initialise(rowptr,colidx,val);
    //    cerr << "\nAint struture \n" << Aint << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;  
/*---------   Rotation matrix from complex to real spherical harmonics ------*/
/*- this is same as in genmat_angint_3D_PN so does not need recalculating  --*/
    RY.New(nsp,nsp);

    /* set sparsity structure for RY */
    nnz = 2*(nsp - sporder)+sporder; // number of nonzeros
    rowptr = new int [nsp+1];
    colidx = new int [nnz];
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    for (int l = 0; l < nso ; l++) {
	for (int m = -l; m < 0; m++) {
	  rowptr[lmYindex(l,m)] = rp;
	  colidx[ci++] = lmYindex(l,(-m));
	  colidx[ci++] = lmYindex(l,(+m));
	  rp += 2;
	}
	rowptr[lmYindex(l,0)] = rp++;
	colidx[ci++] = lmYindex(l,0);
	
	for (int m = 1; m <= l; m++) {
          int sgn =  pow(-1,m);
	  rowptr[lmYindex(l,m)] = rp;
	  colidx[ci++] = lmYindex(l,(-m));
	  colidx[ci++] = lmYindex(l,(+m));
	  rp += 2;
	}
    }
    rowptr[nsp] = rp;
    RY.Initialise(rowptr,colidx);
    //    cerr << "rotation structure \n" << RY << endl;

    delete [] colidx;
    delete [] rowptr;
    delete [] val;
    /* set the values */
    for (int l = 0 ; l < nso ; l++) {
	for (int m = -l; m < 0; m++) {
	  RY(lmYindex(l,m),lmYindex(l,m)) = irt2;
	  RY(lmYindex(l,m),lmYindex(l,-m)) = pow(-1,m)* irt2;
	}
	RY(lmYindex(l,0),lmYindex(l,0)) = 1;
	for (int m = 1; m <= l; m++) {
          int sgn =  pow(-1,m);
	  RY(lmYindex(l,m),lmYindex(l,-m)) = (sgn < 0 ? -iirt2 : iirt2);
	  RY(lmYindex(l,m),lmYindex(l,m)) = -iirt2;
	}
    }
    cerr << "rotation \n" << Sparse2Dense(RY) << endl;
    CCompRowMatrix RYT = transp(RY);
    for (i = 0; i < nsp; i++)
      for(j = 0; j < nsp; j++)
	if(RYT.Exists(i,j))
	  RYT(i,j) = conj(RYT(i,j));
    cerr << "conjg transp rotation \n" << Sparse2Dense(RYT) << endl;
    
    /*---------------------------- Y2m2YY integral -------------------------*/
    cout << "\nBuilding Y2m2YY integral \n";
    CCompRowMatrix Y2m2YYC(nsp,nsp);
   /* set sparsity structure for Y22YYC */
    nnz = 3*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    const double iY2m2N = (SQR(Y1m1N)/Y22[0][0]);

    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      cerr << "row " << k << ":(" << l << ',' << m << ")" << endl;
     
      if(l>1) { // lowering by 2
	if(lmYindex((l-2),(m-2)) >= 0 ) {
	  cerr <<"\t(" <<l <<',' <<m <<")->(" <<l-2 <<',' <<m-2 << ") val = " ;
	  val[ci] = CGmem(l,m)*CGmem((l-1),(m-1))*iY2m2N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l-2),(m-2));
	  rp++;ci++;
	}
      }
      if(lmYindex((l),(m-2))>= 0) {
   	  val00 = toast::complex(0,0);
	  if(lmYindex((l-1),(m-1)) >= 0)
	    val00 += CGmem(l,m)*CGmep((l-1),(m-1));
	  if(lmYindex((l+1),(m-1)) >= 0)
	    val00 += CGmep(l,m)*CGmem((l+1),(m-1));
	  cerr <<"\t(" <<l << ',' <<m <<")->(" <<l << ',' <<m-2 <<") val = " ;
	  val[ci] = val00*iY2m2N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l),(m-2));
	  rp++;ci++;
	}

      if(l <nso-2) {        // raise by 2
	  cerr << "\t("<<l << ','<<m <<")->(" <<l+2 << ',' <<m-2 <<") val = " ;
	  val[ci] = CGmep(l,m)*CGmep((l+1),(m-1))*iY2m2N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l+2),(m-2));
	  rp++;ci++;
      }
      //     for(int vv = 0; vv < ci; vv++)
      //	cerr << " " << val[vv];
    }
    rowptr[nsp] = rp;
    Y2m2YYC.Initialise(rowptr,colidx,val);
    cerr << "\nY2m2YY struture \n" << Y2m2YYC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;
   /*---------------------------- Y2m1YY integral -------------------------*/
    cout << "\nBuilding Y2m1YY integral \n";
    CCompRowMatrix Y2m1YYC(nsp,nsp);
   /* set sparsity structure for Y2m1YYC */
    nnz = 3*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    const double iY2m1N = ((Y1m1N*Y10N)/Y22[0][1]);

    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      cerr << "row " << k << ":(" << l << ',' << m << ")" << endl;
     
      if(l>1) { // lowering by 2
	if(lmYindex((l-2),(m-1)) >= 0 ) {
	  cerr << "\t(" << l << ',' << m << ")->(" << l-2 << ',' << m-1 << ") val = " ;
	  val[ci] = CGmem(l,m)*CGzm((l-1),(m-1))*iY2m1N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l-2),(m-1));
	  rp++;ci++;
	}
      }
      if(lmYindex((l),(m-1))>= 0) {
   	  val00 = toast::complex(0,0);
	  if(lmYindex((l-1),(m-1)) >= 0)
	    val00 += CGmem(l,m)*CGzp((l-1),(m-1));
	  if(lmYindex((l+1),(m-1)) >= 0)
	    val00 += CGmep(l,m)*CGzm((l+1),(m-1));
	  //	  cerr << "\t\t" << iY00N << "\t" << CGpem(l,m)*CGmep((l-1),(m+1))*iY20N << "\t" << CGpep(l,m)*CGmem((l+1),(m+1))*iY20N<< "\t" << CGmep(l,m)*CGpem((l+1),(m-1))*iY20N<< "\t" << CGmem(l,m)*CGpep((l-1),(m-1))*iY20N << endl;
	  cerr << "\t(" << l << ',' << m << ")->(" << l << ',' << m-1 <<") val = " ;
	  val[ci] = val00*iY2m1N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l),(m-1));
	  rp++;ci++;
	}

      if(l <nso-2) {        // raise by 2

      	  cerr << "\t(" << l << ',' << m << ")->(" << l+2 << ',' << m-1 <<") val = " ;
	  val[ci] = CGmep(l,m)*CGzp((l+1),(m-1))*iY2m1N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l+2),(m-1));
	  rp++;ci++;

      }
      //     for(int vv = 0; vv < ci; vv++)
      //	cerr << " " << val[vv];
    }
    rowptr[nsp] = rp;
    Y2m1YYC.Initialise(rowptr,colidx,val);
    cerr << "\nY2m1YY struture \n" << Y2m1YYC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;
   /*---------------------------- Y20YY integral -------------------------*/
    cout << "\nBuilding Y20YY integral \n";
    CCompRowMatrix Y20YYC(nsp,nsp);
   /* set sparsity structure for Y20YYC */
    nnz = 3*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    const double iY20N = ((Y1p1N*Y1m1N)/Y22[0][2]);
    const double iY00N = -((Y1p1N*Y00N*sqrt(2.0/3.0))/Y22[0][2]);

    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      cerr << "row " << k << ":(" << l << ',' << m << ")" << endl;
     
      if(l>1) { // lowering by 2
	if(lmYindex((l-2),(m)) >= 0 ) {
	  cerr << "\t(" << l << ',' << m << ")->(" << l-2 << ',' << m << ") val = " ;
	  val[ci] = CGpem(l,m)*CGmem((l-1),(m+1))*iY20N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l-2),(m));
	  rp++;ci++;
	}
      }
      if(lmYindex((l),(m))>= 0) {
   	  val00 = toast::complex(iY00N/iY20N,0);
	  if(lmYindex((l-1),(m+1)) >= 0)
	    val00 += CGpem(l,m)*CGmep((l-1),(m+1));
	  if(lmYindex((l+1),(m+1)) >= 0)
	    val00 += CGpep(l,m)*CGmem((l+1),(m+1));
	  //	  cerr << "\t\t" << iY00N << "\t" << CGpem(l,m)*CGmep((l-1),(m+1))*iY20N << "\t" << CGpep(l,m)*CGmem((l+1),(m+1))*iY20N<< "\t" << CGmep(l,m)*CGpem((l+1),(m-1))*iY20N<< "\t" << CGmem(l,m)*CGpep((l-1),(m-1))*iY20N << endl;
	  cerr << "\t(" << l << ',' << m << ")->(" << l << ',' << m <<") val = " ;
	  val[ci] = val00*iY20N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l),(m));
	  rp++;ci++;
	}

      if(l <nso-2) {        // raise by 2

      	  cerr << "\t(" << l << ',' << m << ")->(" << l+2 << ',' << m <<") val = " ;
	  val[ci] = CGpep(l,m)*CGmep((l+1),(m+1))*iY20N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l+2),(m));
	  rp++;ci++;


      }
      //     for(int vv = 0; vv < ci; vv++)
      //	cerr << " " << val[vv];
    }
    rowptr[nsp] = rp;
    Y20YYC.Initialise(rowptr,colidx,val);
    cerr << "\nY20YY struture \n" << Y20YYC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;

  /*---------------------------- Y2p1YY integral -------------------------*/
    cout << "\nBuilding Y2p1YY integral \n";
    CCompRowMatrix Y2p1YYC(nsp,nsp);
   /* set sparsity structure for Y2m1YYC */
    nnz = 3*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    const double iY2p1N = ((Y1p1N*Y10N)/Y22[1][2]);

    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      cerr << "row " << k << ":(" << l << ',' << m << ")" << endl;
     
      if(l>1) { // lowering by 2
	if(lmYindex((l-2),(m+1)) >= 0 ) {
	  cerr << "\t(" << l << ',' << m << ")->(" << l-2 << ',' << m+1 << ") val = " ;
	  val[ci] = CGpem(l,m)*CGzm((l-1),(m+1))*iY2p1N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l-2),(m+1));
	  rp++;ci++;
	}
      }
      if(lmYindex((l),(m+1))>= 0) {
   	  val00 = toast::complex(0,0);
	  if(lmYindex((l-1),(m+1)) >= 0)
	    val00 += CGpem(l,m)*CGzp((l-1),(m+1));
	  if(lmYindex((l+1),(m+1)) >= 0)
	    val00 += CGpep(l,m)*CGzm((l+1),(m+1));
	  //	  cerr << "\t\t" << iY00N << "\t" << CGpem(l,m)*CGmep((l-1),(m+1))*iY20N << "\t" << CGpep(l,m)*CGmem((l+1),(m+1))*iY20N<< "\t" << CGmep(l,m)*CGpem((l+1),(m-1))*iY20N<< "\t" << CGmem(l,m)*CGpep((l-1),(m-1))*iY20N << endl;
	  cerr << "\t(" << l << ',' << m << ")->(" << l << ',' << m+1 <<") val = " ;
	  val[ci] = val00*iY2p1N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l),(m+1));
	  rp++;ci++;
	}

      if(l <nso-2) {        // raise by 2

      	  cerr << "\t(" << l << ',' << m << ")->(" << l+2 << ',' << m+1 <<") val = " ;
	  val[ci] = CGpep(l,m)*CGzp((l+1),(m+1))*iY2p1N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l+2),(m+1));
	  rp++;ci++;

      }
      //     for(int vv = 0; vv < ci; vv++)
      //	cerr << " " << val[vv];
    }
    rowptr[nsp] = rp;
    Y2p1YYC.Initialise(rowptr,colidx,val);
    cerr << "\nY2p1YY struture \n" << Y2p1YYC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;  
    /*---------------------------- Y2p2YY integral -------------------------*/
    cout << "\nBuilding Y2p2YY integral \n";
    CCompRowMatrix Y2p2YYC(nsp,nsp);
   /* set sparsity structure for Y22YYC */
    nnz = 3*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    const double iY2p2N = (SQR(Y1p1N)/Y22[2][2]);

    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      cerr << "row " << k << ":(" << l << ',' << m << ")" << endl;
     
      if(l>1) { // lowering by 2
	if(lmYindex((l-2),(m+2)) >= 0 ) {
	  cerr << "\t(" << l << ',' << m << ")->(" << l-2 << ',' << m+2 << ") val = " ;
	  val[ci] = CGpem(l,m)*CGpem((l-1),(m+1))*iY2p2N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l-2),(m+2));
	  rp++;ci++;
	}
      }
      if(lmYindex((l),(m+2))>= 0) {
   	  val00 = toast::complex(0,0);
	  if(lmYindex((l-1),(m+1)) >= 0)
	    val00 += CGpem(l,m)*CGpep((l-1),(m+1));
	  if(lmYindex((l+1),(m+1)) >= 0)
	    val00 += CGpep(l,m)*CGpem((l+1),(m+1));
	  cerr << "\t(" << l << ',' << m << ")->(" << l << ',' << m+2 <<") val = " ;
	  val[ci] = val00*iY2p2N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l),(m+2));
	  rp++;ci++;
	}

      if(l <nso-2) {        // raise by 2

	  cerr << "\t(" << l << ',' << m << ")->(" << l+2 << ',' << m+2 <<") val = " ;
	  val[ci] = CGpep(l,m)*CGpep((l+1),(m+1))*iY2p2N;
	  cerr << val[ci] << endl;
	  colidx[ci] = lmYindex((l+2),(m+2));
	  rp++;ci++;


      }
      //     for(int vv = 0; vv < ci; vv++)
      //	cerr << " " << val[vv];
    }
    rowptr[nsp] = rp;
    Y2p2YYC.Initialise(rowptr,colidx,val);
    cerr << "\nY2p2YY struture \n" << Y2p2YYC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;
 
 
    /*---------------------------- zzYY integral --------------------------*/
    /*--------------- this one is exact, for comparison -------------------*/
    cout << "\nBuilding RRzz\n";
    CCompRowMatrix YYzzC(nsp,nsp);
   /* set sparsity structure for YYzC */
    nnz = 3*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      //      cerr << "\n(l,m) : (" << l << ',' << m << ") ind(l-1,m) : " << lmYindex((l-1),m) << " ind(l+1,m) : " << lmYindex((l+1),m);

      if(l >1&& abs(m) < l-1) { // lowering by 2
	val[ci] = CGzzm(l,m);
	colidx[ci] = lmYindex((l-2),m);
	rp++;ci++;
      }
      val00 = CGzz1(l,m); // raise and lower
//      cerr << "\nzz1 : (l,m) : (" << l << ',' << m << ") val00 = " << val00;
      if(l >0&& abs(m) < l) {// lower and raise
	val00 += CGzz0(l,m);
//	cerr << "\nzz0 : (l,m) : (" << l << ',' << m << ") val00 = " << val00;
      }

      val[ci] = val00;
      colidx[ci] = lmYindex(l,m);
      rp++;ci++;
      if(l <nso-2) {        // raise by 2
	val[ci] = CGzzp(l,m);
	colidx[ci] = lmYindex((l+2),m);
	rp++;ci++;
      }
    }
    rowptr[nsp] = rp;
    YYzzC.Initialise(rowptr,colidx,val);
    //    cerr << "\nYYzz struture \n" << YYzzC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;
    
 
    tmp = RY*YYzzC*RYT;
    //  cerr << tmp;
    
 
    //    CCompRowMatrix tmp = RY*YYzzC*RYT; // not this!
    
    //    Aintcc = RealMat(tmp);
    cout << "\nBuilt RRzz OK\n";
    cerr <<  RealMat(tmp) << endl;

/*------------ now form all the real matrices by combination and rotation */
    const double Y2N = sqrt(4*M_PI/5);
    RCompRowMatrix sxx(nsp,nsp);
    tmp = RY*((Y2m2YYC+Y2p2YYC)*toast::complex(sqrt(1.0/6.0),0) - Y20YYC*toast::complex(1.0/3.0))*RYT*toast::complex(Y2N,0) +idnspC*toast::complex(1.0/3.0);
    sxx = RealMat(tmp);
    Aintscsc = Chop(sxx);
    cout << "\nBuilt RRxx by Y2 method\n";
    cerr << Aintscsc << endl;

    RCompRowMatrix syy(nsp,nsp);
    tmp = RY*((Y2m2YYC+Y2p2YYC)*toast::complex(-sqrt(1.0/6.0),0) - Y20YYC*toast::complex(1.0/3.0))*RYT*toast::complex(Y2N,0) +idnspC*toast::complex(1.0/3.0);
    syy = RealMat(tmp);
    Aintssss = Chop(syy);
    cout << "\nBuilt RRyy by Y2 method\n";
    cerr << Aintssss << endl;

    RCompRowMatrix szz(nsp,nsp);
    tmp = RY*Y20YYC*RYT*toast::complex(2*Y2N/3,0)+idnspC*toast::complex(1.0/3.0);
    szz = RealMat(tmp);
    Aintcc = Chop(szz);
    cout << "\nBuilt RRzz by Y2 method\n";
    cerr << Aintcc << endl;

    RCompRowMatrix sxy(nsp,nsp);

    tmp = RY*((Y2m2YYC-Y2p2YYC)*toast::complex(0,sqrt(1.0/6.0)))*RYT*toast::complex(Y2N,0);
    sxy = RealMat(tmp);
    Aintscss = Chop(sxy);
    cout << "\nBuilt RRxy by Y2 method\n";
    cerr << Aintscss << endl;

   RCompRowMatrix sxz(nsp,nsp);
   //   cout << "testing Y21 matrices\n(Y2m1YYC-Y2p2YYC) : " << Sparse2Dense(Y2m1YYC-Y2p1YYC) <<endl;
    tmp = RY*(Y2m1YYC-Y2p1YYC)*RYT*toast::complex(Y2N*sqrt(1.0/6.0),0);
    sxz = RealMat(tmp);
    Aintscc = Chop(sxz);
    cout << "\nBuilt RRxz by Y2 method\n";
    cerr << Aintscc << endl;

   RCompRowMatrix syz(nsp,nsp);
    tmp = RY*((Y2m1YYC+Y2p1YYC)*toast::complex(0,sqrt(1.0/6.0)) )*RYT*toast::complex(Y2N,0);
    syz = RealMat(tmp);
    Aintssc = Chop(syz);
    cout << "\nBuilt RRyz by Y2 method\n";
    cerr << Aintssc << endl;
}
