/***************************************************************************
 * genmat_angint_3D_PN.cc         Simon Arridge                20.07.10    *
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


 
void genmat_angint_3D_PN(RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Anvec, CCompRowMatrix& RY, const int sporder)
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

  Aint.New(nsp,nsp);   // integrals of products of spherical harmonics 

  Aintsc.New(nsp,nsp); // integrals of products of spherical harmonics with x

  Aintss.New(nsp,nsp); // integrals of products of spherical harmonics with y

  Aintc.New(nsp,nsp);  // integrals of products of spherical harmonics with z

  Anvec.New(nsp,nsp);  // and this one...

  /*--------- first the trivial identity matrix for Aint -----------*/
  cout << "Building Aint nsp = "<< nsp << endl;
    nnz = nsp;
    rowptr = new int [nsp+1];
    colidx = new int [nnz];
    rval = new double [nnz];
    for (i = 0; i < nsp; i++){
 //      cerr << "i " << i;
      rowptr[i] = rp++;
      rval[ci] = 1.0;
      colidx[ci++] = i;
 //      cerr << " " << rowptr[i] << " " << colidx[i] << " " << val[i] <<endl;
    }
    rowptr[nsp] = rp;
    Aint.Initialise(rowptr,colidx,rval);
    //    cerr << "\nAint struture \n" << Aint << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] rval;
  /*------------- the strange one.. Anvec. It seems to be same as Aint  --*/
    Anvec = Aint;
  /*--------- Rotation matrix from complex to real spherical harmonics  --*/

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
    //    cerr << "rotation \n" << Sparse2Dense(RY) << endl;
    CCompRowMatrix RYT = transp(RY);
    for (i = 0; i < nsp; i++)
      for(j = 0; j < nsp; j++)
	if(RY.Exists(i,j))
	  RYT(i,j) = conj(RYT(i,j));
    //    cerr << "conjg transp rotation \n" << Sparse2Dense(RYT) << endl;
    

    /*---------------------------- z YY integral --------------------------*/
    cout << "\nBuilding RRz\n";
    CCompRowMatrix YYzC(nsp,nsp);
   /* set sparsity structure for YYzC */
    nnz = 2*nsp-2;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz]; 
    val = new toast::complex[nnz];
    rp = 0; ci = 0;
    for(int k = 0; k < nsp; k++) {
      rowptr[k] = rp;
      int l = lind(k);
      int m = mind(k);
      //      cerr << "\n(l,m) : (" << l << ',' << m << ") ind(l-1,m) : " << lmYindex((l-1),m) << " ind(l+1,m) : " << lmYindex((l+1),m);

      if(l >0&& abs(m) < l) {
	//	cerr << "\nci " << ci << " " << (colidx[ci] = lmYindex((l-1),m) );
	val[ci] = CGzm(l,m);
	colidx[ci] = lmYindex((l-1),m);
	rp++;ci++;
      }
      if(l <nso-1) {
	//        cerr << "\nci " << ci << " " << (colidx[ci] = lmYindex((l+1),m) );
	val[ci] = CGzp(l,m);
	colidx[ci] = lmYindex((l+1),m);
	rp++;ci++;
      }
    }
    rowptr[nsp] = rp;
    YYzC.Initialise(rowptr,colidx,val);
    //    cerr << "\nYYz struture \n" << YYzC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;
    
    CCompRowMatrix tmp1(nsp,nsp);
    tmp1 = RY*YYzC*RYT;
    //  cerr << tmp1;
    
    /*
    tmp1 = RY;
    tmp1 = tmp1*YYzC;
    tmp1 = tmp1*RYT;
    */
    //    CCompRowMatrix tmp1 = RY*YYzC*RYT; // not this!
    
    Aintc = RealMat(tmp1);
    cout << "\nBuilt RRz OK\n";
    tmp1.New(nsp,nsp);
    //    cerr << Aintc << endl;

    /*---------------------------- x YY integral --------------------------*/
    cout << "Building RRx\n";
    CCompRowMatrix YYxC(nsp,nsp);
  /* set sparsity structure for YYzC */
    nnz = 4*nsp;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz];
    val = new toast::complex[nnz];
    rp = 0; ci = 0;

    for(int k = 0; k < nsp; k++) {
      int l = lind(k);
      int m = mind(k);
      rowptr[k] = rp;
      if(l >0) {
	if(m > -l+1){
	  //       	  cerr << "\n(1)k : " << k << " (l,m) : (" << l-1 << ',' << m-1 << ")\t->  "<< lmYindex((l-1),(m-1));
	  val[ci] = CGmem(l,m)/2;
	  colidx[ci++] = lmYindex((l-1),(m-1));rp++;
	}
	if(m < l-1) {
	  //       	  cerr << "\n(2)k : " << k << " (l,m) : (" << l-1 << ',' << m+1 << ")\t->  "<< lmYindex((l-1),(m+1));
	  val[ci] = CGpem(l,m)/2;
	  colidx[ci++] = lmYindex((l-1),(m+1));rp++;
	}
      }
      if(l <nso-1) {
	if(m > -l-1) {
	  //       	cerr << "\n(3)k : " << k << " (l,m) : (" << l+1 << ',' << m-1 << ")\t->  "<< lmYindex((l+1),(m-1));
	  val[ci] = CGmep(l,m)/2;
	  colidx[ci++] = lmYindex((l+1),(m-1));rp++;
	}
      	if(m < l+1) {
	  //       	cerr << "\n(4)k : " << k << " (l,m) : (" << l+1 << ',' << m+1 << ")\t->  "<< lmYindex((l+1),(m+1));
	  val[ci] = CGpep(l,m)/2;
	  colidx[ci++] = lmYindex((l+1),(m+1));rp++;
	}
      }
    }
    rowptr[nsp] = rp;
    YYxC.Initialise(rowptr,colidx,val);
    //    cerr << "\nYYx struture \n" << YYxC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;
    CCompRowMatrix tmp2(nsp,nsp);
    tmp2 = RY*YYxC*RYT;
    //    cerr << "\nRotated OK\n" << tmp2;
    Aintsc = RealMat(tmp2);
    cerr << "\nBuilt RRx OK\n";
    //    tmp2.New(nsp,nsp);
    //    cerr << Aintsc << endl;

    /*---------------------------- y YY integral --------------------------*/
    cout << "\nBuilding RRy\n";

    CCompRowMatrix YYyC(nsp,nsp);
    nnz = 4*nsp;// this is a slight over estimate
    rowptr = new int [nsp+1];
    colidx = new int [nnz];
    val = new toast::complex[nnz];
    rp = 0; ci = 0;

   for(int k = 0; k < nsp; k++) {
      int l = lind(k);
      int m = mind(k);
      rowptr[k] = rp;

      if(l >0) {
	if(m > -l+1){
	  val[ci] = toast::complex(0,CGmem(l,m)/2);
	  colidx[ci++] = lmYindex((l-1),(m-1));rp++;
	}
	if(m < l-1){
	  val[ci] = toast::complex(0,-CGpem(l,m)/2);
	  colidx[ci++] = lmYindex((l-1),(m+1));rp++;
	}
      }
      if(l <nso-1) {
	if(m > -l-1){
	  val[ci] = toast::complex(0,CGmep(l,m)/2);
	  colidx[ci++] = lmYindex((l+1),(m-1));rp++;
	}
      	if(m < l+1) {
	  val[ci] = toast::complex(0,-CGpep(l,m)/2);
	  colidx[ci++] = lmYindex((l+1),(m+1));rp++;
	}
      }
    }
    rowptr[nsp] = rp;
    YYyC.Initialise(rowptr,colidx,val);
    //    cerr << "\nYYy struture \n" << YYyC << endl;
    delete [] colidx;
    delete [] rowptr;
    delete [] val;

    CCompRowMatrix tmp3 = RY*YYyC*RYT;
    Aintss  = RealMat(tmp3);
    cout << "\nBuilt RRy OK\n";
    //    cerr << Aintss << endl;
}
