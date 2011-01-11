//==========================================================================
// genmat_3DBman_newF.cc              S.Arridge                     27.11.06
//
//==========================================================================

#include <stdio.h>
#include <stdlib.h>
#include <stream.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include "toast.h"
#include "rte3D.h"
    using namespace toast;


// generates some sparse matrices of "system" size.
// sigmat is an input array of (mua+mus) parameters per element 
// sigma is array of angular scattering kernels, per element

void genmat_3DBman_newF(RCompRowMatrix& A3, RCompRowMatrix& A4, const Mesh& mesh, const RDenseMatrix* sigma, const RVector& sigmatot, const int SN, const RCompRowMatrix& Aint, const RCompRowMatrix& Aintsc, const RCompRowMatrix& Aintss, const RCompRowMatrix& Aintc, const RCompRowMatrix* SP, const RCompRowMatrix* SPdx, const RCompRowMatrix* SPdy, const RCompRowMatrix* SPdz, const RCompRowMatrix& Anvec, const RCompRowMatrix& Anvec_sc, const RCompRowMatrix& Anvec_ss, const RCompRowMatrix& Anvec_c, const ScatKernType sktyp)
{
   int sysdim = mesh.nlen();       // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;   // full size of angles X space nodes

   A3.New(fullsysdim,fullsysdim);

   cout << "calculating A3 integrals\n";

   //   RCompRowMatrix A3_rte(fullsysdim,fullsysdim) ;
   //   RCompRowMatrix A3_sdm(fullsysdim,fullsysdim) ;
   RCompRowMatrix spatA3_rte(sysdim,sysdim) ;
   RCompRowMatrix spatA3_sdmx(sysdim,sysdim) ;
   RCompRowMatrix spatA3_sdmy(sysdim,sysdim) ;
   RCompRowMatrix spatA3_sdmz(sysdim,sysdim) ;
   
   int  nel = mesh.elen();
   for (int el = 0; el < nel; el++) {
     if(!(el*100%nel) )
       cout << el*100/nel << "%\n";
     spatA3_rte = spatA3_rte + SP[el]*sigmatot[el];
     spatA3_sdmx = spatA3_sdmx + SPdx[el]*sigmatot[el];
     spatA3_sdmy = spatA3_sdmy + SPdy[el]*sigmatot[el];
     spatA3_sdmz = spatA3_sdmz + SPdz[el]*sigmatot[el];
   }

   A3 = kron(spatA3_rte,Aint) + kron(spatA3_sdmx,Aintsc) + kron(spatA3_sdmy,Aintss) + kron(spatA3_sdmz,Aintc);
   // clear space
   //   A3_rte.New(0,0);
   //   A3_sdm.New(0,0);
   spatA3_rte.New(0,0);
   spatA3_sdmx.New(0,0);
   spatA3_sdmy.New(0,0);
   spatA3_sdmz.New(0,0);
  
   cout << "\t- A3 finished\n";

   // A4 matrix. Identify 4 different cases for speed...

   A4.New(fullsysdim,fullsysdim);

   //   ScatKernType sktyp =  MUSINHOMOG_GENERAL;
   cout << "Scatter kernel type is " << sktyp << endl;

   switch (sktyp) {
     case MUSHOMOG_G0 :
       { // compound is needed to get scopes..
       RCompRowMatrix apu1(SN,SN);
       RCompRowMatrix apu1sc(SN,SN);
       RCompRowMatrix apu1ss(SN,SN);
       RCompRowMatrix apu1c(SN,SN);
       for(int bb = 0; bb < SN; bb++){ // not efficient !
	 RCompRowMatrix acol = Anvec.Subcols(bb,bb+1);
	 RCompRowMatrix acolsc = Anvec_sc.Subcols(bb,bb+1);
	 RCompRowMatrix acolss = Anvec_ss.Subcols(bb,bb+1);
	 RCompRowMatrix acolc  = Anvec_c.Subcols(bb,bb+1);
	 for(int cc = 0; cc < SN; cc++){
	   RCompRowMatrix bcol = Anvec.Subcols(cc,cc+1);
	   bcol.Transpone();// row vector of column

	   apu1 = apu1 + kron(acol,bcol);
	   apu1sc = apu1sc + kron(acolsc,bcol);
	   apu1ss = apu1ss + kron(acolss,bcol);
	   apu1c  = apu1c  + kron(acolc,bcol);
	 }
       }
       //       cout << "apu1 " << apu1 << endl;
       //       cout << "apu1sc " << apu1sc << endl;
       //       cout << "apu1ss " << apu1ss << endl;

       RCompRowMatrix SPS(sysdim,sysdim);
       RCompRowMatrix SPSdx(sysdim,sysdim);
       RCompRowMatrix SPSdy(sysdim,sysdim); 
       RCompRowMatrix SPSdz(sysdim,sysdim); 

       for (int el = 0; el < nel; el++) {
	 if(!(el*100%nel) )
	   cout << el*100/nel << "%\n";
	 //	 cout << SP[el] << endl;
       	 SPS = SPS + SP[el];
	 //	 cout << SPdx[el] << endl ;
       	 SPSdx = SPSdx + SPdx[el];
	 //       	 cout << SPdy[el] << endl;
       	 SPSdy = SPSdy + SPdy[el];
	 //       	 cout << SPdz[el] << endl;
       	 SPSdz = SPSdz + SPdz[el];
       }
       //      cout << "apu1 \n " << apu1 << endl;
       //       cout << "apu1sc \n " << apu1sc << endl;
       //       cout << "apu1ss \n " << apu1ss << endl;
       //       cout << "About to call kron \n";
       cout << "completed angular integrals : about to assemble A4\n";
       A4 = kron(SPS,apu1) + kron(SPSdx,apu1sc) + kron(SPSdy,apu1ss) + kron(SPSdz,apu1c);
       //       cout << "About to scale A4\n";
       A4 = A4*(-sigma[0](0,0)); // 
       }
       break;
     case MUSINHOMOG_G0 :
       break;
     case MUSHOMOG_GCONST :
      { // compound is needed to get scopes..
 
       RCompRowMatrix SPS(sysdim,sysdim);
       RCompRowMatrix SPSdx(sysdim,sysdim);
       RCompRowMatrix SPSdy(sysdim,sysdim); 
       RCompRowMatrix SPSdz(sysdim,sysdim); 

       for (int el = 0; el < nel; el++) {
	 if(!(el*100%nel) )
	   cout << el*100/nel << "%\n";
       	 SPS = SPS + SP[el];
       	 SPSdx = SPSdx + SPdx[el];
       	 SPSdy = SPSdy + SPdy[el];
       	 SPSdy = SPSdz + SPdz[el];
       }
       RVector apujjvec(SN*SN), apuscjjvec(SN*SN), apussjjvec(SN*SN), apucjjvec(SN*SN);
       for(int jj = 0; jj < SN; jj++) {
	   RCompRowMatrix acol  = Anvec.Subcols(jj,jj+1);
	   RCompRowMatrix asccol = Anvec_sc.Subcols(jj,jj+1);
	   RCompRowMatrix asscol = Anvec_ss.Subcols(jj,jj+1);
	   RCompRowMatrix accol  = Anvec_c.Subcols(jj,jj+1);

	   RCompRowMatrix apuii = kron(acol,Anvec);
	   RCompRowMatrix apuscii = kron(asccol,Anvec);
	   RCompRowMatrix apussii = kron(asscol,Anvec);
	   RCompRowMatrix apucii  = kron(asscol,Anvec);


	   RVector sigcol(SN); // assign a column of sigma[el];
	   //	   cout << "sigma[el] " << sigma[el] << endl;
     	   for (int ii = 0; ii < SN; ii++) {
	     //	     cout << sigma[el](ii,jj);
	     sigcol[ii] = sigma[0](ii,jj);
	   }
	   apujjvec += apuii*sigcol;
	   apuscjjvec += apuscii*sigcol;
	   apussjjvec += apussii*sigcol;
	   apucjjvec  += apucii*sigcol;

	 }
	 // do a "reshape "
	 RCompRowMatrix apujj_reshape(SN,SN);	 
	 int* arowptr;
	 if( !(arowptr = new int [SN+1]))
	   cerr << "Memory Allocation error arowptr = new int\n";
	 int* acolidx;
	 if( !(acolidx = new int [SN*SN]))
	   cerr << "Memory Allocation error acolidx = new int\n";
	 int rp = 0;
	 for(int ii = 0, cp = 0; ii < SN; ii++) {
	   arowptr[ii] = rp; rp += SN;
	   for(int kk = 0; kk < SN; kk++)
	     acolidx[cp++] = kk;
	 }
	 arowptr[SN] = rp;
	 apujj_reshape.Initialise(arowptr,acolidx);
	 RCompRowMatrix apuscjj_reshape = apujj_reshape;
	 RCompRowMatrix apussjj_reshape = apujj_reshape;
	 RCompRowMatrix apucjj_reshape  = apujj_reshape;
	 for (int ii = 0; ii < SN; ii++){
	   for(int kk = 0; kk < SN; kk++) {
	     apujj_reshape(ii,kk) = apujjvec[ii*SN + kk];
	     apuscjj_reshape(ii,kk) = apuscjjvec[ii*SN + kk];
	     apussjj_reshape(ii,kk) = apussjjvec[ii*SN + kk];
	     apucjj_reshape(ii,kk)  = apucjjvec[ii*SN + kk];
	   }
	 }
	 cout << "completed angular integrals : about to assemble A4\n";
      	 A4 = A4 - kron(SPS,apujj_reshape);
	 A4 = A4 - kron(SPSdx,apuscjj_reshape);
	 A4 = A4 - kron(SPSdy,apussjj_reshape);
	 A4 = A4 - kron(SPSdz,apucjj_reshape);

	 delete [] arowptr;
	 delete [] acolidx;
       }
       break;
     case MUSINHOMOG_GENERAL :
       { // compound is needed to get scopes..

       for (int el = 0; el < nel; el++) {
	 if(!(el*100%nel) )
	   cout << el*100/nel << "%\n";

	 //	 RCompRowMatrix apujj, apuscjj, apusjj;
	 RVector apujjvec(SN*SN), apuscjjvec(SN*SN), apussjjvec(SN*SN), apucjjvec(SN*SN);
	 for(int jj = 0; jj < SN; jj++) {
	   //#ifdef USECOLS
	   RCompRowMatrix acol  = Anvec.Subcols(jj,jj+1);
	   RCompRowMatrix asccol = Anvec_sc.Subcols(jj,jj+1);
	   RCompRowMatrix asscol = Anvec_ss.Subcols(jj,jj+1);
	   RCompRowMatrix accol  = Anvec_c.Subcols(jj,jj+1);
	   //#else
	   /*
	   RCompRowMatrix acol  = Anvec.Subrows(jj,jj+1);
	   acol.Transpone();
	   RCompRowMatrix asccol = Anvec_c.Subrows(jj,jj+1);
	   asccol.Transpone();
	   RCompRowMatrix asscol = Anvec_s.Subrows(jj,jj+1);
	   asscol.Transpone();
	   */
	   //#endif
	   //	   RCompRowMatrix sigcol = sigma[el].Subcols(jj,jj+1);

	   //      	   cout << "acol " << acol << endl;
	   //       	   cout << "asccol " << asccol << endl;
	   //       	   cout << "asscol " << asscol << endl;

	   RCompRowMatrix apuii = kron(acol,Anvec);
	   RCompRowMatrix apuscii = kron(asccol,Anvec);
	   RCompRowMatrix apussii = kron(asscol,Anvec);
	   RCompRowMatrix apucii  = kron(accol,Anvec);
	   //       	   cout << "apuii " << apuii << endl;
	   //       	   cout << "apuscii " << apuscii << endl;
	   //       	   cout << "apussii " << apussii << endl;
	

	   //	   apujj = apujj + apuii*sigcol;
	   //	   apuscjj = apuscjj + apuscii*sigcol;
	   //	   apussjj = apussjj + apuscii*sigcol;

	   RVector sigcol(SN); // assign a column of sigma[el];
	   //	   cout << "sigma[el] " << sigma[el] << endl;
     	   for (int ii = 0; ii < SN; ii++) {
	     //	     cout << sigma[el](ii,jj);
	     sigcol[ii] = sigma[el](ii,jj);
	   }
	   apujjvec += apuii*sigcol;
	   apuscjjvec += apuscii*sigcol;
	   apussjjvec += apussii*sigcol;
	   apucjjvec  += apucii*sigcol;
	   //	   cout << "apujjvec " << apujjvec << endl;
	   //	   cout << "apujjvec " << apujjvec << endl;
	   //	   cout << "apujjvec " << apujjvec << endl;
	 }
	 // do a "reshape "
	 RCompRowMatrix apujj_reshape(SN,SN);	 
	 int* arowptr;
	 if( !(arowptr = new int [SN+1]))
	   cerr << "Memory Allocation error arowptr = new int\n";
	 int* acolidx;
	 if( !(acolidx = new int [SN*SN]))
	   cerr << "Memory Allocation error acolidx = new int\n";
	 int rp = 0;
	 for(int ii = 0, cp = 0; ii < SN; ii++) {
	   arowptr[ii] = rp; rp += SN;
	   for(int kk = 0; kk < SN; kk++)
	     acolidx[cp++] = kk;
	 }
	 arowptr[SN] = rp;
	 apujj_reshape.Initialise(arowptr,acolidx);
	 RCompRowMatrix apuscjj_reshape = apujj_reshape;
	 RCompRowMatrix apussjj_reshape = apujj_reshape;
	 RCompRowMatrix apucjj_reshape  = apujj_reshape;
	 for (int ii = 0; ii < SN; ii++){
	   for(int kk = 0; kk < SN; kk++) {
	     apujj_reshape(ii,kk) = apujjvec[ii*SN + kk];
	     apuscjj_reshape(ii,kk) = apuscjjvec[ii*SN + kk];
	     apussjj_reshape(ii,kk) = apussjjvec[ii*SN + kk];
	     apucjj_reshape(ii,kk)  = apucjjvec[ii*SN + kk];
	   }
	 }
	 //	 cout << "apujj_reshape " << apujj_reshape << endl;
	 //	 cout << "apuscjj_reshape " << apuscjj_reshape << endl;
	 //	 cout << "apussjj_reshape " << apussjj_reshape << endl;
	 //	 RCompRowMatrix tmp4_rte = kron(SP[el],apujj_reshape);
	 //	 RCompRowMatrix tmp4_sdm_x = kron(SPdx[el],apuscjj_reshape);
	 //	 RCompRowMatrix tmp4_sdm_y = kron(SPdy[el],apussjj_reshape);
	 //     	 A4 = A4 - tmp4_rte - tmp4_sdm_x - tmp4_sdm_y;
	 cout << "completed angular integrals : about to assemble A4\n";
      	 A4 -= kron(SP[el],apujj_reshape);
	 A4 -= kron(SPdx[el],apuscjj_reshape);
	 A4 -= kron(SPdy[el],apussjj_reshape);
	 A4 -= kron(SPdz[el],apucjj_reshape);
	 //	 cout << "tmp4_rte "   << tmp4_rte;
	 //	 cout << "tmp4_sdm_x " << tmp4_sdm_x;
	 //	 cout << "tmp4_sdm_y " << tmp4_sdm_y;
	 //	 cout << "A4 " << A4;
	 delete [] arowptr;
	 delete [] acolidx;
       } // end loop on elements
       } // end of this switch case
       break;
   } // end of switch
   cout << "\t- A4 finished\n";

}









