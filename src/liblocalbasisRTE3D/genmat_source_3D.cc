//==========================================================================
// genmat_source_3D.cc                  S.Arridge                   27.11.06
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


void genmat_source_3D(CCompRowMatrix* & Source, const Mesh& mesh,  const Mesh& S2mesh, const int* Nsource, const int ns, const CCompRowMatrix& b1)
  /*      
    Function generates the source values vector for FEM of the radiative 
    transfer equation
  */
{
   const int& SN =  S2mesh.nlen();   // dimensions are size of nodes.
   int sysdim = mesh.nlen();         // dimensions are size of nodes.
   int fullsysdim = sysdim*SN;       // full size of angles X space nodes

   if( !(Source = new  CCompRowMatrix [ns]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix\n";

   CCompRowMatrix Svec;
   for (int i = 0; i < ns; i++) {
     genmat_sourcevalvector_cos_3D(Svec, mesh, S2mesh, Nsource[i]);
     Source[i].New(fullsysdim,1);
     b1.AB(Svec,Source[i]);        // Toast weirdness
   }
}












