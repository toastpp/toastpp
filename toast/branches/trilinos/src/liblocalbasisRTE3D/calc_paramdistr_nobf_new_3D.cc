//==========================================================================
// calc_paramdistr_nobf_new_3D.cc         S.Arridge                 27.11.06
//
//==========================================================================

#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "toast.h"
#include "rte3D.h"

// generates some parameter distributions

void calc_paramdistr_nobf_new_3D(RVector& sigma, RVector& sigmatot, RVector& intst, const Mesh& mesh, const RVector& muscat, const RVector& muabs, const RVector& g,const Mesh& S2mesh)
{

   const int& SN =  S2mesh.nlen();       // dimensions are size of nodes.
   const int& SE =  S2mesh.elen();       // number of spherical elements.
   const int sysdim = mesh.nlen();       // dimensions are size of nodes.
   const int fullsysdim = sysdim*SN;     // full size of angles X space nodes

   // data for normalisation
   intst.New(SN);
   intst = 4*M_PI/SN ;  // just a hack value for now.
   // the normalised phase function multiplied by sigmascat

     sigmatot = muscat + muabs;

   RVector shati(3), shatj(3); // direction cosines
   for(int el = 0; el < mesh.elen(); el++) {
     //sigma[el].New(SN,SN);
     double ghg = g[el];           // anisotropy for this element
     for(int i = 0; i < SN; i++) {
       shati = S2mesh.nlist[i] ;       // do we really need to copy it ?
       double ileni = 1.0/length(shati); // normalise  just to be sure
       shati *= ileni;
       for(int j = 0; j < SN; j++) {
         shatj = S2mesh.nlist[j] ;       // do we really need to copy it ?
         double ilenj = 1.0/length(shatj); // normalise just to be sure
         shatj *= ilenj;
	 // evaluate Henyey-Greenstein function 3D version
	 double denom = 1.0;//(1 + ghg*ghg - 2*ghg*(shati&shatj));
	 sigma[el] = muscat[el]*1/(4*M_PI*denom*sqrt(denom));
       }
     }
     if(!(el*100%mesh.elen()) )
       std::cout << el*100/mesh.elen() << "%\n";
   } // end loop on elements
}









