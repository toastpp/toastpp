//==========================================================================
// precompute_3DBman_FEM.cc         S.Arridge                     21.02.06
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
// delta is an input array of parameters pertinent to the streamline diffusion
// functions.

void precompute_3DBman_FEM(const Mesh& mesh,  const RVector& delta, const Mesh& S2mesh, const double w, const double c, CCompRowMatrix& A0, RCompRowMatrix& A1, RCompRowMatrix& A2,RCompRowMatrix& b1, RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix* &SP, RCompRowMatrix* &SPdx, RCompRowMatrix* &SPdy,  RCompRowMatrix* &SPdz, RCompRowMatrix& Anvec, RCompRowMatrix& Anvec_sc, RCompRowMatrix& Anvec_ss, RCompRowMatrix& Anvec_c)
{

   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   //   const int& SE =  S2mesh.elen();       // number of spherical elements.
   const int sysdim = mesh.nlen();       // dimensions are size of nodes.
   const int fullsysdim = sysdim*SN;     // full size of angles X space nodes

   cout << "calculating spatial integrals\n";

   RCompRowMatrix Sint, Sgrad, Sx, Sy, Sz;

   genmat_spatint_nobf_3D(mesh, Sint, Sgrad, Sx, Sy, Sz, SP);
   
   cout << "\tfinished conventional spatial integrals over domain\n";
   RCompRowMatrix Sdx, Sdy, Sdz, Sdxx, Sdxy, Sdyx, Sdyy,  Sdxz, Sdzx,  Sdyz, Sdzy, Sdzz;

   genmat_spatint_sdm_nobf_3D(mesh, delta, Sdx, Sdy, Sdz, Sdxx, Sdxy, Sdyx, 
		     Sdyy,  Sdxz, Sdzx, Sdyz, Sdzy, Sdzz, SPdx, SPdy, SPdz);
   cout << "\tfinished streamline diffusion integrals over domain\n";

   cout << "calculating angular integrals\n";
 
   RCompRowMatrix Aintscsc,  Aintscss, Aintscc,  Aintssss,  Aintssc,  Aintcc;

   genmat_angint_3D(Aint, Aintsc, Aintss, Aintc, Anvec, S2mesh);
   genmat_angint_sdm_3D(Aintscsc,  Aintscss,   Aintscc,  Aintssss,  Aintssc,  Aintcc,  Anvec_sc, Anvec_ss,  Anvec_c, S2mesh);

   cout << "calculating boundary and source matrices A2 and b1\n";

   genmat_boundint_3D(A2, mesh, S2mesh );
   A2.Shrink();

   genmat_sourceint_3D(b1, mesh, S2mesh); 
   b1.Shrink();


   cout << "completing A0 and A1 matrices\n";

   RCompRowMatrix A0_rte = kron(Sint, Aint);
   RCompRowMatrix A0_sdm = kron(Sdx, Aintsc);
   A0_sdm += kron(Sdy, Aintss);
   A0_sdm +=  A0_sdm+ kron(Sdz, Aintc);
   // clear space
   Sint.New(0,0);
   Sgrad.New(0,0);   // Sgrad never seems to be used ... WHY not?
   Sdx.New(0,0);
   Sdy.New(0,0);
   Sdz.New(0,0);
   complex I1(0,1); // imaginary number
   A0 = cplx(A0_rte);
   A0 += cplx(A0_sdm);
   A0 *= (I1*(w/c));

   // clear space
   A0_rte.New(0,0);
   A0_sdm.New(0,0);
   
   cout << "\t- A0 finished\n";
   
   RCompRowMatrix mA1_rte = kron(Sx, Aintsc);
   mA1_rte +=   kron(Sy, Aintss);
   mA1_rte +=   kron(Sz, Aintc);
   RCompRowMatrix A1_sdm = kron(Sdxx, Aintscsc);
   A1_sdm += kron(Sdxy, Aintscss);
   A1_sdm += kron(Sdyx, Aintscss);
   A1_sdm += kron(Sdyy, Aintssss);
   A1_sdm += kron(Sdxz, Aintscc) ;
   A1_sdm += kron(Sdzx, Aintscc);
   A1_sdm += kron(Sdyz, Aintssc);
   A1_sdm += kron(Sdzy, Aintssc);
   A1_sdm += kron(Sdzz, Aintcc) ;

   A1 = A1_sdm - mA1_rte;
   // clear space
   mA1_rte.New(0,0);
   A1_sdm.New(0,0);
   Sx.New(0,0);

   Sy.New(0,0);
   Sz.New(0,0);
   Sdxx.New(0,0);
   Sdxy.New(0,0);
   Sdyx.New(0,0);
   Sdyy.New(0,0);
   Sdxz.New(0,0);
   Sdzx.New(0,0);
   Sdyz.New(0,0);
   Sdzy.New(0,0);
   Sdzz.New(0,0);
   
   cout << "\t- A1 finished\n";
  
}









