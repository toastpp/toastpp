/***************************************************************************
 * test_edd.cc                       Simon Arridge          16.10.06       *
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
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
//#include <toast.h>
#include <sys/times.h>
#include <sys/param.h>
#include "stoastlib.h"
//#include "fwdsolver.h"
//#include "solution.h"
//#include "source.h"
//#include "pparse.h"
#include <ilutoast.h>

using namespace toast;
void Project (const QMMesh &mesh, int q, const CCompRowMatrix &mvec,
	      const CVector &phi, CVector &proj);
void ProjectAll (const QMMesh &mesh, 
    const CCompRowMatrix &mvec, const CVector *dphi, CVector &proj);
void OpenNIM (const char *nimname, const char *meshname, int size);
void WriteNIM (const char *nimname, const RVector &img, int size, int no);
bool ReadNim (char *nimname, RVector &img);
void WritePGM (const RVector &img, const IVector &gdim, char *fname);
void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname);
void WriteData (const RVector &data, char *fname);
void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname);
void SelectMesh (ParamParser &pp, char *meshname, QMMesh &mesh);
void SelectIntfMesh (ParamParser &pp, char *meshname, Mesh &mesh);
void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);
void SelectData (ParamParser &pp, int nqm, Measurement &dtype, double &freq);
void SelectInitialElParams (ParamParser &pp, Mesh &mesh);
double rfresnel_cos(const double c1, const double n1, const double n2);
double calcbeta(double n1, double n2);

#define MAXREGION 100
#define SOURCE_GAUSSIAN 0
#define SOURCE_COSINE 1
#define GMRESSOLVE

// =========================================================================

SourceMode srctp = SRCMODE_NEUMANN;   // source type
double lin_tol;
char *colormap = "../../scales/fire2.pal";

QMMesh qmmesh;
NodeList &nlist=qmmesh.nlist;
ElementList &elist=qmmesh.elist;
ParameterList &plist = qmmesh.plist;
Mesh intfmesh;                  // mesh holding interface surface

// error handler for FE library routines *************************************

void LocalErrorhandler (char *msg)
{
    cerr << "\nforward_refind (PID " << getpid() << ")\n" << msg << endl << flush;
    cerr << "Aborted.\n";
    //    logfile << msg << endl << "Aborted." << endl;
    exit (1);
}

// main routine **************************************************************

int main (int argc, char *argv[])
{

    char cbuf[256], meshname[256], imeshname[256], qmname[256] ;
    double elk_ij, elb_ij, ela_ij;
    int el, nodel, i, j, k,is, js;
    double freq, omega;    // modulation frequency (MHz and cycles/ps)
    int qprof, mprof;      // source/measurement profile
                           // (0=point, 1=Gaussian, 2=cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]
    Measurement datatype;
    int nQ, nM, nQM;
    CCompRowMatrix qvec, mvec;
    int idx, ofs, resettp, cmd;
    struct tms tm;
    double t0, t1;
    Point bbmin, bbmax;
    //    logfile << "Reading mesh" << endl;
    /*
    cout << "Reading mesh " << argv[1] << endl;
    ifstream ifs;
    ifs.open (argv[1]);
    xASSERT (ifs.is_open(), Mesh file not found.);
    ifs >> qmmesh;
    xASSERT (ifs.good(), Problem reading mesh.);
    ifs.close ();
    */

    ParamParser pp;
    SelectMesh (pp, meshname, qmmesh);

    int dimension = nlist[0].Dim();
    for (i = 1; i < nlist.Len(); i++)
	xASSERT(nlist[i].Dim() == dimension, Inconsistent node dimensions.);
    xASSERT(dimension >= 2 && dimension <= 3, Mesh dimension must be 2 or 3.);
    //    elist[0]->Initialise(nlist);
    qmmesh.Setup();

    double c = 0.3;     
    SelectSourceProfile (pp, qprof, qwidth, srctp);
    SelectMeasurementProfile (pp, mprof, mwidth);
    int dim = qmmesh.Dimension();
    int n = qmmesh.nlen();
    int elen = qmmesh.elen();
    nQ = qmmesh.nQ;
    nM = qmmesh.nM;
    nQM = qmmesh.nQM;
    qmmesh.BoundingBox (bbmin, bbmax);

    cout << "* " << qmmesh.elen() << " elements, " << qmmesh.nlen()
	 << " nodes\n";
    cout << "nQ " << nQ << " nM "<< nM  << endl;

    SelectData (pp, nQM, datatype, freq);
    //    FWS.datatype = datatype;
    omega = freq * 2.0*M_PI*1e-6;

    // build the source vectors
    qvec.New (nQ, n);
    LOGOUT1_INIT_PROGRESSBAR ("Source vectors", 50, nQ);
    for (i = 0; i < nQ; i++) {
	CVector q(n);
	switch (qprof) {
	case 0:
	    SetReal (q, QVec_Point (qmmesh, qmmesh.Q[i], srctp));
	    break;
	case 1:
	    SetReal (q, QVec_Gaussian (qmmesh, qmmesh.Q[i], qwidth, srctp));
	    break;
	case 2:
	    SetReal (q, QVec_Cosine (qmmesh, qmmesh.Q[i], qwidth, srctp));
	    break;
	case 3:
	  //	    q = CompleteTrigSourceVector (qmmesh, i);
	    break;
	}
	qvec.SetRow (i, q);

	cerr << "max=" << vmax(q) << endl;

	LOGOUT1_PROGRESS(i);
    }

    // build the measurement vectors
    mvec.New (nM, n);
    RVector c2avals(n); 
    for (el = 0; el < qmmesh.elen(); el++) { // fix boundary node c2s vals
      for (i = 0; i < elist[el]->nNode(); i++) {
  	int is;
	if ((is = elist[el]->Node[i]) >= n) continue;
	if(nlist[is].isBnd())
	  c2avals[is] =qmmesh.plist[el].C2A();
      }
    }    
//    cout << "c2avals " << c2avals << endl;
    LOGOUT1_INIT_PROGRESSBAR ("Meas. vectors", 50, nM);
    for (i = 0; i < nM; i++) {
	CVector m(n);
	switch (mprof) {
	case 0:
	    SetReal (m, QVec_Point (qmmesh, qmmesh.M[i], SRCMODE_NEUMANN));
	    break;
	case 1:
	    SetReal (m, QVec_Gaussian (qmmesh, qmmesh.M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case 2:
	    SetReal (m, QVec_Cosine (qmmesh, qmmesh.M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	case 3:
	  //	    m = CompleteTrigSourceVector (qmmesh, i);
	    break;
	}
	for (j = 0; j < n; j++) m[j] *= c2avals[j];//qmmesh.plist[j].C2A();
	mvec.SetRow (i, m);
	LOGOUT1_PROGRESS(i);
    }
    //    cout << "Measurement vectors \n" << mvec << endl;
     // reset initial parameter estimates
    //    Solution sol(OT_NPARAM, n);
    SelectInitialElParams (pp, qmmesh);
    int    er,nr;
    IVector nelr(MAXREGION);
    RVector erind(MAXREGION);     // record the refractive indices
    for(i = 0; i < qmmesh.elen(); i++)  {
	er = qmmesh.elist[i]->Region();
	if (er >= 0 && er < MAXREGION) nelr[er]++;
	erind[er] = plist[i].N();
    }
    for (nr = i = 0; i < MAXREGION; i++)
	if (nelr[i]) {
	  nr++;
        }
    cout << nr << " distinct regions with refractive indices " ;
    for(i = 0; i < nr ; i++)
      cout << erind[i] << " ";
    cout << endl;

    /* get the parameters */
    RVector mua( plist.Mua());
    RVector diff(plist.Kappa());
    RVector rind(plist.N());
    RVector c2a(plist.C2A());
    RVector cspeed(plist.C()); // this is the speed including refind : 0.3/N


    //    cout << "mua    " << mua << endl;
    //    cout << "diff   " << diff << endl;
    //    cout << "refind " << rind << endl;
    //    cout << "c2a" << c2a;
    //    cout << "speed " << cspeed << endl;
 


    IVector interfaceindex(qmmesh.nlen());

    // list the interface nodes
    int ninterface = 0;

    cout << "\nInterface nodes : ";
    for(i = 0, k=0; i< qmmesh.nlen();i++) {
      if(nlist[i].isInternalInterface()) {
	interfaceindex[i] =  k++;
	//	cout << i << " ";
	ninterface++;
      }
    }
    cout << "\nTotal : " << ninterface << endl;
    ninterface = ninterface/2 ; // careful!
    cout << "Unique Total : " << ninterface << endl;
    IVector interfacelist(ninterface);
    for(i = 0, k=0; i< qmmesh.nlen()-ninterface;i++) {
      if(nlist[i].isInternalInterface()) {
	interfacelist[k++] = i;
      }
    }
    for(i = qmmesh.nlen()-ninterface, k=0; i<qmmesh.nlen(); i++) {
      if(nlist[i].isInternalInterface()) {
	interfaceindex[i] = k++;
      }
    }
    
    cout << "Interface List " << interfacelist << endl;
    cout << "Interface Index " << interfaceindex << endl;


    // read in the surface mesh

    SelectIntfMesh (pp, imeshname, intfmesh);
    intfmesh.Setup();
    cout << "* " << intfmesh.elen() << " elements, " << intfmesh.nlen()
	 << " nodes\n";
    for(i =0; i < intfmesh.nlen(); i++) {
      cout << i << " " << intfmesh.nlist[i] << "\t" << nlist[interfacelist[i]] << endl ;
    }

    //************ system matrices *********************************
    cout << "Setting up matrix structure\n";

    int sysdim = qmmesh.nlen();       // dimensions are size of nodes.
    cout << "sysdim " << sysdim << endl;
// reset storage allocation from mesh neighbour list
    int *rowptr, *colidx, nzero;
    qmmesh.SparseRowStructure (rowptr, colidx, nzero);
    //    CCompRowMatrix AA(sysdim, sysdim);
    //    AA.Initialise (rowptr, colidx);
    cout << "Got sparse stucture and initialised. Non zeros : " << nzero <<endl;

    //** now do the coupling matrices mesh-interface and interface-interface

    int isysdim = intfmesh.nlen();       // dimensions are size of nodes.
    cout << "isysdim " << isysdim << endl;
// reset storage allocation from mesh neighbour list
    int *BBrowptr, *BBcolidx, BBnzero;
    intfmesh.SparseRowStructure (BBrowptr, BBcolidx, BBnzero);
;
    cout <<"BB rows : ";
    for(i =0 ;i <= isysdim;i++)
      cout << BBrowptr[i] << " ";
    cout <<"\nBB cols : ";
    for(i =0 ;i < BBnzero;i++)
      cout << (BBcolidx[i]) << " ";
    cout << endl;

    CCompRowMatrix BB(isysdim, isysdim);
    BB.Initialise (BBrowptr, BBcolidx);
    //    cout << iBB << endl;
    cout << "Got sparse stucture and initialised. Non zeros : " << BBnzero <<endl;



    int *ABrowptr, *ABcolidx, ABnzero;
    if( !(ABrowptr = new int [sysdim +1])) // holds sparse strutucre 
     cerr << "Memory Allocation error ABrowptr = new int\n";
    if( !(ABcolidx = new int [nzero]))     // maximum size 
     cerr << "Memory Allocation error ABcolindx = new int\n";

    // BB matrix has same column structure. The row structure just misses out
    // the empty rows of AB
    int rp = 0, cp = 0;
    for(i = 0, k=0; i < sysdim-ninterface; i++){
       ABrowptr[i] = cp;
       if(!nlist[i].isInternalInterface()) continue;
       int ki = interfaceindex[i];       // which row of BB ?
       for(j = BBrowptr[ki]; j < BBrowptr[ki+1]; j++) {
	   ABcolidx[cp++] = BBcolidx[j]; //
        }
    } // end loop on nodes
    // now the extra nodes
    for(i = sysdim-ninterface; i < sysdim; i++){
       ABrowptr[i] = cp;
       if(!nlist[i].isInternalInterface()) continue;
       int ki = interfaceindex[i];       // which row of BB ?
       for(j = BBrowptr[ki]; j < BBrowptr[ki+1]; j++) {
	   ABcolidx[cp++] = BBcolidx[j]; //
       }
    } // end loop on nodes


    ABnzero = ABrowptr[sysdim]  = cp; // final pointer
    cout <<"AB rows : ";
    for(i =0 ;i <= sysdim;i++)
      cout << ABrowptr[i] << " ";
    cout <<"\nAB cols : ";
    for(i =0 ;i < ABnzero;i++)
      cout << (ABcolidx[i]) << " ";
    cout << endl;

    CCompRowMatrix AB(sysdim, isysdim);
    AB.Initialise (ABrowptr, ABcolidx);
    //    cout << iBB << endl;
    cout << "Got sparse stucture and initialised. Non zeros : " << ABnzero <<endl;

    // *** add the matrices together (the hard way !)
    int fullnzero =  nzero+2*ABnzero + BBnzero;
    cout << "Allocating dummy arrays of size " << fullnzero << endl;
    int *dummyrow, *dummycol; 
    if( !(dummyrow = new int [fullnzero])) // holds row indices struture 
     cerr << "Memory Allocation error dummyrow = new int\n";
    if( !(dummycol = new int [fullnzero])) // holds col indices struture 
      cerr << "Memory Allocation error dummycol = new int\n";
    toast::complex *dummydata;
    if( !(dummydata = new toast::complex [fullnzero])) // holds data  
     cerr << "Memory Allocation error dummydata = new int\n";    
    // insert the structure into dummy rows 
    cout << "Building global matrix structure \n";
    for(i = 0, k=0; i < sysdim; i++){
      for(j = rowptr[i]; j < rowptr[i+1]; j++){
	dummyrow[k] =i; dummycol[k] = colidx[j];
	k++;
      }  
      for(j = ABrowptr[i]; j < ABrowptr[i+1]; j++){
	dummyrow[k] =i; 
	dummycol[k] = sysdim + ABcolidx[j];
	k++;
	dummycol[k] =i; 
	dummyrow[k] = sysdim + ABcolidx[j]; // it's the transpose!
	k++;
      }      
    }
    for(i = 0; i < ninterface; i++){ 
      for(j = BBrowptr[i]; j < BBrowptr[i+1]; j++){
	dummyrow[k] = sysdim + i; 
	dummycol[k] = sysdim + BBcolidx[j];
	k++;
	//	dummycol[k] = sysdim +i; 
	//	dummyrow[k] = sysdim + ABcolidx[j]; // it's the transpose!
	//	k++;
      }      
    }

    cout << "Total non zeros " << k << " estimate was " << fullnzero << endl;
    CCoordMatrix KKcd(sysdim+ ninterface,sysdim+ ninterface,fullnzero,dummyrow,dummycol,dummydata);
    cout << "Created CoordMatrix structure\n";

    for(i = 0; i < sysdim; i++){
      for(j = ABrowptr[i]; j < ABrowptr[i+1]; j++){
        k = ABcolidx[j];
	//	if(i < 836) 
	KKcd(i,sysdim + k) = 1;//AB(i,k);
	//	  	cout << i << " " << sysdim + k << " OK ";
	//      if(i < 836)
	KKcd(sysdim + k,i) = 1;// AB(i,k); // Transpose
	//	        cout << sysdim + k << " " << i << " " << AB(i,k) << " ";
      }      
      //            cout << endl;
    }
    //***************************** something goes wrong before the end ****
    cout << "finished copying AB\n";
    for(i = 0; i < ninterface; i++){
      for(j = BBrowptr[i]; j < BBrowptr[i+1]; j++){
        k = BBcolidx[j];
	KKcd(sysdim+i,sysdim + k) = 1;//BB(i,k);
	//	cout << BB(i,k) << " ";
      }      
      //      cout << endl;
    }
    cout << "finished copying BB\n";

    cout << "Got CoordMatrix structure\n";
    //    cout << KKcd << endl;
    CCompRowMatrix KK(sysdim+ ninterface,sysdim+ ninterface);
    KK = KKcd;
    cout << "Converted to Comp Row matrix structure\n";
    //    KK.merge(AA); 
    //    cout << "Merged OK\n";
    cout << "KK has " << KK.rowptr[sysdim+ ninterface] << " non zeros\n";
    /*
    cout << KK << endl;
    for(i = 0; i < sysdim+ninterface; i++) {
      for(j=0; j < sysdim+ninterface; j++){
       	cout << KK.Get(i,j) << " ";
      }
      cout << endl;
    }
    */   
    // release memory - this goes wrong!
    //   AA.New(0,0);
    //    AB.New(0,0);
    //    BB.New(0,0);
    //    KKcd.New(0,0);
    delete [] dummyrow;
    delete [] dummycol;
    delete [] dummydata;

    //**************** compute the matrix elements **********************
    cout << "calculating finite element integrals\n";

    double gamma = (1.2*1.2)/(1.4*1.4); // n2 < n1
    double homg = 0.5*(1-gamma); // half one minus gamma
    double beta = -2;    // arbitrary value
    if(nr > 0){
      gamma = erind[1]/erind[0];
      gamma = gamma*gamma;
      homg = 0.5*(1 - gamma);
      beta = -calcbeta(erind[0],erind[1]);
      cout << "Setting gamma " << gamma  << " beta " << beta << endl;
    }
    for (el = 0; el < qmmesh.elen(); el++) {
        int er = elist[el]->Region();
	//        cout << el << " region " << er << endl;
	nodel = elist[el]->nNode();

// now determine the element integrals
	for (i = 0; i < nodel; i++) {
	  if ((is = elist[el]->Node[i]) >= sysdim) continue;
	  for (j = 0; j < nodel; j++) {
      	    if ((js = elist[el]->Node[j]) >= sysdim) continue;
	    elb_ij = elist[el]->IntFF (i, j);
	    elk_ij = elist[el]->IntDD (i, j);
	    ela_ij = (elist[el]->HasBoundarySide() ? 
		      elist[el]->BndIntFF (i, j) : 0.0);
	    //	    cout << elb_ij << " " << elk_ij << " "  << ela_ij << endl; 
	    KK(is,js) += toast::complex(c*mua[el],-omega*rind[el])*elb_ij + 
	      c*diff[el]*elk_ij  + c2a[el]*ela_ij;

	  }
	}
    } // end element loop
    cout << "Done element integrals\n";
    // now add the interface integrals
    for (el = 0; el < intfmesh.elen(); el++) {
      //        cout << el  << endl;
	nodel = intfmesh.elist[el]->nNode();

// now determine the element integrals
	int ki, ki2, kj, kj2;
	for (i = 0; i < nodel; i++) {
	  if ((is = intfmesh.elist[el]->Node[i]) >= isysdim) continue;
	  for (j = 0; j < nodel; j++) {
      	    if ((js = intfmesh.elist[el]->Node[j]) >= isysdim) continue;
	    double ifint = intfmesh.elist[el]->IntFF (i, j);
 
	    ki = interfacelist[is];
	    kj = interfacelist[js];
	    ki2 =  sysdim-ninterface+is;
	    kj2 =  sysdim-ninterface+js;

	      KK(ki,sysdim+js) +=  -gamma*ifint;
	      KK(sysdim+js,ki) +=  -gamma*ifint;
	      KK(ki2,sysdim+js) +=  ifint;
	      KK(sysdim+js,ki2) +=  ifint;
      	      KK(is+sysdim,js+sysdim) += beta*ifint;
	  }
	  // following should be the same
//  cout << is << " " << interfaceindex[ki] << " " << interfaceindex[ki2] << endl;
	}
    } // end element loop

    cout << "finished KK matrix\n";

    /*
    cout << KK << endl;
    for(i = 0; i < sysdim+ninterface; i++) {
      for(j=0; j < sysdim+ninterface; j++){
       	cout << KK.Get(i,j) << " ";
      }
      cout << endl;
    }
     */

    //********** calculate radiance ***********************************

    cout << "calculating the radiance\n";

    CVector RHS(sysdim+ninterface);
    //      cout << "Built RHS\n";
    CVector x(sysdim+ninterface);
    //      cout << "Built x\n";
    CVector J(ninterface);
    //      cout << "Built J\n";
    CVector * dphi = new CVector [nQ];

    double tol = 1e-10;
#ifdef GMRESSOLVE
    //    CPrecon_Diag * KKCP = new  CPrecon_Diag;
    CPrecon_DILU * KKCP = new  CPrecon_DILU;

    KKCP->Reset(&KK);
#else
    ILUSetPermtype (ILUPERM_RCM);
    double droptol = 1e-3;
    int ierr;
#endif
    //      cout << "Reset Preconditioner\n";

    ofstream osPhi("Phi.sol");
    ofstream osJ("J.sol");

    for (j = 0; j < nQ ; j++) {
      cout << "start source " << j << endl;
      dphi[j].New(sysdim);
      cout << "Built dphi\n";
      RHS = toast::complex(0,0);
      for(i = 0; i < sysdim; i++) 
	RHS[i] = qvec.Get(j,i); //) << " ";
      //      BiCGSTAB (KK, RHS, x, tol,KKCP, 100);
      cout << "\ncalling GMRES\n";
#ifdef GMRESSOLVE
      GMRES (KK, RHS, x, tol,KKCP, 500);
#else
      ierr = ILUSolve (KK, RHS, x, tol, droptol, 500);
#endif
      cout << "finished source " << j << endl;
      for(i = 0; i < sysdim; i++) 
	dphi[j][i] = x[i];
      for(i = 0; i < ninterface; i++)
	J[i] = x[sysdim+i];
      cout << "max=" << vmax(dphi[j]) << endl;

      osPhi << "Phi " << j << "\n" << dphi[j] << endl;
      osJ << "J " << j << "\n" << J << endl;
    }
    osPhi.close();
    osJ.close();

    // output fields as NIM files
    cout << "Output nodal fields (1/0)? " << flush;
    cin >> cmd;
    if (cmd) {
        switch (datatype) {
	case MEAS_FRE_FIM:
	    OpenNIM ("phi_re.nim", meshname, sysdim);
	    OpenNIM ("phi_im.nim", meshname, sysdim);
	    for (i = 0; i < nQ; i++) {
	        WriteNIM ("phi_re.nim", Re(dphi[i]), sysdim, i);
		WriteNIM ("phi_im.nim", Im(dphi[i]), sysdim, i);
	    }
	    cout << "  Real field written to phi_re.nim" << endl;
	    cout << "  Imag field written to phi_im.nim" << endl;
	    break;
	case MEAS_FMOD_FARG:
	    OpenNIM ("phi_lnmod.nim", meshname, sysdim);
	    OpenNIM ("phi_arg.nim", meshname, sysdim);
	    for (i = 0; i < nQ; i++) {
	        WriteNIM ("phi_lnmod.nim", LogMod(dphi[i]), sysdim, i);
		WriteNIM ("phi_arg.nim", Arg(dphi[i]), sysdim, i);
	    }
	    cout << "  Log Mod field written to phi_lnmod.nim" << endl;
	    cout << "  Arg field written to phi_arg.nim" << endl;
	    break;
	}
    }
    // output data files
    
    cout << "Output data files (1/0)? " << flush;
    cin >> cmd;
    if (cmd) {
	CVector proj(nQM);
	ProjectAll (qmmesh, mvec, dphi, proj);
	//	cout << proj;
	switch (datatype) {
	case MEAS_FRE_FIM:
	    WriteData (Re(proj), "fre.fem");
	    WriteData (Im(proj), "fim.fem");
	    break;
	case MEAS_FMOD_FARG:
	    WriteData (LogMod(proj), "fmod.fem");
	    WriteData (Arg(proj), "farg.fem");
	    WriteDataBlock (qmmesh, LogMod(proj), "fmod.dat");
	    WriteDataBlock (qmmesh, Arg(proj), "farg.dat");
	    break;
	}     
    }
    delete []rowptr;
    delete []colidx;  
    delete []ABrowptr;
    delete []ABcolidx;   
    delete []BBrowptr;
}
//---------------------------------------------------------------------------
// generate a projection from a field - complex case

void Project (const QMMesh &mesh, int q, const CCompRowMatrix &mvec,
    const CVector &phi, CVector &proj)
{
    int i, m;
    for (i = 0; i < mesh.nQMref[q]; i++) {
	m = mesh.QMref[q][i];
	proj[i] = dot (phi, mvec.Row(m));
    }
}

void ProjectAll (const QMMesh &mesh, 
    const CCompRowMatrix &mvec, const CVector *dphi, CVector &proj)
{
    int i, len, ofs = 0;

    for (i = 0; i < mesh.nQ; i++) {
	len = mesh.nQMref[i];
	CVector proj_i (proj, ofs, len);
	Project (mesh, i, mvec, dphi[i], proj_i);
	ofs += len;
    }
}
void OpenNIM (const char *nimname, const char *meshname, int size)
{
    ofstream ofs (nimname);
    ofs << "NIM" << endl;
    ofs << "Mesh = " << meshname << endl;
    ofs << "SolutionType = N/A" << endl;
    ofs << "ImageSize = " << size << endl;
    ofs << "EndHeader" << endl;
}

void WriteNIM (const char *nimname, const RVector &img, int size, int no)
{
    ofstream ofs (nimname, ios::app);
    ofs << "Image " << no << endl;
    for (int i = 0; i < size; i++)
        ofs << img[i] << ' ';
    ofs << endl;
}

bool ReadNim (char *nimname, RVector &img)
{
    char cbuf[256];
    int i, imgsize = 0;

    ifstream ifs (nimname);
    if (!ifs.getline (cbuf, 256)) return false;
    if (strcmp (cbuf, "NIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return false;
    img.New(imgsize);
    do {
        ifs.getline (cbuf, 256);
    } while (strncasecmp (cbuf, "Image", 5));
    for (i = 0; i < imgsize; i++)
        ifs >> img[i];
    return true;
}

void WritePGM (const RVector &img, const IVector &gdim, char *fname)
{
    int i, ii, dim = gdim[0]*gdim[1];
    double imgmin = 1e100, imgmax = -1e100;
    unsigned char *pixmap = new unsigned char[dim];

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    for (i = 0; i < dim; i++) {
        if (img[i+ii] < imgmin) imgmin = img[i+ii];
	if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    double scale = 256.0/(imgmax-imgmin);
    for (i = 0; i < dim; i++) {
        int v = (int)((img[i+ii]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P5" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++) ofs << pixmap[i];
    ofs << endl;
    delete []pixmap;
}

void WritePPM (const RVector &img, const IVector &gdim,
    double *scalemin, double *scalemax, char *fname)
{
    typedef struct {
        unsigned char r,g,b;
    } RGB;

    int i, ii, dim = gdim[0]*gdim[1];
    unsigned char *pixmap = new unsigned char[dim];
    double imgmin, imgmax, scale;

    static RGB colmap[256];
    static bool have_colmap = false;

    if (!have_colmap) {
        int r, g, b;
        ifstream ifs (colormap);
	for (i = 0; i < 256; i++) {
	    ifs >> r >> g >> b;
	    colmap[i].r = (unsigned char)r;
	    colmap[i].g = (unsigned char)g;
	    colmap[i].b = (unsigned char)b;
	}
	have_colmap = true;
    }

    if (gdim.Dim() > 2)
        ii = (gdim[2]/2)*gdim[0]*gdim[1]; // simply plot middle layer
    else
        ii = 0;

    // rescale image
    if (scalemin) {
        imgmin = *scalemin;
    } else {
	for (i = 0, imgmin = 1e100; i < dim; i++)
	    if (img[i+ii] < imgmin) imgmin = img[i+ii];
    }
    if (scalemax) {
        imgmax = *scalemax;
    } else {
      for (i = 0, imgmax = -1e100; i < dim; i++)
	  if (img[i+ii] > imgmax) imgmax = img[i+ii];
    }
    scale = 256.0/(imgmax-imgmin);

    for (i = 0; i < dim; i++) {
        int v = (int)((img[i+ii]-imgmin)*scale);
	if      (v < 0  ) v = 0;
	else if (v > 255) v = 255;
	pixmap[i] = (unsigned char)v;
    }
    ofstream ofs(fname);
    ofs << "P6" << endl;
    ofs << "# CREATOR: TOAST (test_grid)" << endl;
    ofs << gdim[0] << ' ' << gdim[1] << endl;
    ofs << 255 << endl;
    for (i = 0; i < dim; i++)
        ofs << colmap[pixmap[i]].r
	    << colmap[pixmap[i]].g
	    << colmap[pixmap[i]].b;
    ofs << endl;

    delete []pixmap;
}

void WriteData (const RVector &data, char *fname)
{
    ofstream ofs (fname);
    ofs << setprecision(14);
    ofs << data << endl;
}

void WriteDataBlock (const QMMesh &mesh, const RVector &data, char *fname)
{
    int q, m, i;
    ofstream ofs (fname);
    for (m = i = 0; m < mesh.nM; m++) {
	for (q = 0; q < mesh.nQ; q++) {
	    if (mesh.Connected (q,m)) ofs << data[i++];
	    else                      ofs << '-';
	    ofs << (q == mesh.nQ-1 ? '\n' : '\t');
	}
    }   
}
void SelectMesh (ParamParser &pp, char *meshname, QMMesh &mesh)
{
    char qmname[256];

    if (!pp.GetString ("MESHFILE", meshname)) {
        cout << "\nMesh file name:\n>> ";
	cin >> meshname;
    }
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();

    if (!pp.GetString ("QMFILE", qmname)) {
        cout << "\nQM file name:\n>> ";
	cin >> qmname;
    }
    ifstream qmf (qmname);
    mesh.LoadQM (qmf);

    // write back
    pp.PutString ("MESHFILE", meshname);
    pp.PutString ("QMFILE", qmname);
}
void SelectIntfMesh (ParamParser &pp, char *meshname, Mesh &mesh)
{
    char qmname[256];

    if (!pp.GetString ("MESHFILE", meshname)) {
        cout << "\nMesh file name:\n>> ";
	cin >> meshname;
    }
    ifstream ifs (meshname);
    ifs >> mesh;
    mesh.Setup();


    // write back
    pp.PutString ("MESHFILE", meshname);
}
// ============================================================================

void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp)
{
    char cbuf[256];
    int cmd;

    bool typeok = false;
    if (pp.GetString ("SOURCETYPE", cbuf)) {
	if (!strcasecmp (cbuf, "NEUMANN")) {
	    srctp = SRCMODE_NEUMANN;
	    typeok = true;
	} else if (!strcasecmp (cbuf, "ISOTROPIC")) {
	    srctp = SRCMODE_ISOTROPIC;
	    typeok = true;
	}
    }
    while (!typeok) {
	cout << "\nSource type:\n";
	cout << "(1) Neumann boundary source\n";
	cout << "(2) Isotropic point source\n";
	cout << "[1|2] >> ";
	cin  >> cmd;
	switch (cmd) {
	    case 1: srctp = SRCMODE_NEUMANN;   typeok = true; break;
	    case 2: srctp = SRCMODE_ISOTROPIC; typeok = true; break;
	}
    }
    pp.PutString ("SOURCETYPE",
        srctp == SRCMODE_NEUMANN ? "NEUMANN" : "ISOTROPIC");

    qtype = -1;
    if (pp.GetString ("SOURCEPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    qtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    qtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    qtype = 2;
	}
    }
    while (qtype < 0) {
	cout << "\nSource profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> qtype;
	qtype -= 1;
    }
    if (qtype > 0 && !pp.GetReal ("SOURCEWIDTH", qwidth)) {
	switch (qtype) {
	case 1:
	    cout << "\nSource 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nSource support radius [mm]:\n>> ";
	    break;
	}
	cin >> qwidth;
    }
    switch (qtype) {
    case 0:
	pp.PutString ("SOURCEPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("SOURCEPROFILE", "GAUSSIAN");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    case 2:
	pp.PutString ("SOURCEPROFILE", "COSINE");
	pp.PutReal ("SOURCEWIDTH", qwidth);
	break;
    }
}

// ============================================================================

void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth)
{
    char cbuf[256];
    mtype = -1;
    if (pp.GetString ("MEASUREMENTPROFILE", cbuf)) {
	if (!strcasecmp (cbuf, "POINT")) {
	    mtype = 0;
	} else if (!strcasecmp (cbuf, "GAUSSIAN")) {
	    mtype = 1;
	} else if (!strcasecmp (cbuf, "COSINE")) {
	    mtype = 2;
	}
    }
    while (mtype < 0) {
	cout << "\nMeasurement profile:\n";
	cout << "(1) Point\n";
	cout << "(2) Gaussian\n";
	cout << "(3) Cosine\n";
	cout << "[1|2|3] >> ";
	cin  >> mtype;
	mtype -= 1;
    }
    if (mtype > 0 && !pp.GetReal ("MEASUREMENTWIDTH", mwidth)) {
	switch (mtype) {
	case 1:
	    cout << "\nMeasurement 1/e radius [mm]:\n>> ";
	    break;
	case 2:
	    cout << "\nMeasurement support radius [mm]:\n>> ";
	    break;
	}
	cin >> mwidth;
    }
    switch (mtype) {
    case 0:
	pp.PutString ("MEASUREMENTPROFILE", "POINT");
	break;
    case 1:
	pp.PutString ("MEASUREMENTPROFILE", "GAUSSIAN");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    case 2:
	pp.PutString ("MEASUREMENTPROFILE", "COSINE");
	pp.PutReal ("MEASUREMENTWIDTH", mwidth);
	break;
    }
}

// ============================================================================

void SelectData (ParamParser &pp, int nqm, Measurement &dtype, double &freq)
{
    char cbuf[256];
    bool def = false;
    RVector pdata;

    if (pp.GetString ("DATATYPE", cbuf)) {
	if (!strcasecmp (cbuf, "REALIMAG")) {
	    dtype = MEAS_FRE_FIM, def = true;
	} else if (!strcasecmp (cbuf, "MODARG")) {
	    dtype = MEAS_FMOD_FARG, def = true;
	} else if (!strcasecmp (cbuf, "INTENSITY")) {
	    dtype = MEAS_INTENSITY, def = true;
	}
    }
    while (!def) {
	int cmd;
	cout << "\nSelect data type:\n";
	cout << "(1) Complex intensity (real+imag)\n";
	cout << "(2) Complex intensity (mod+arg)\n";
	cout << "(3) CW intensity\n";
	cout << "[1|2|3] >> ";
	cin >> cmd;
	switch (cmd) {
	case 1: dtype = MEAS_FRE_FIM,   def = true; break;
	case 2: dtype = MEAS_FMOD_FARG, def = true; break;
	case 3: dtype = MEAS_INTENSITY, def = true; break;
	}
    }
    switch (dtype) {
    case MEAS_FRE_FIM:
	pp.PutString ("DATATYPE", "REALIMAG");
	break;
    case MEAS_FMOD_FARG:
	pp.PutString ("DATATYPE", "MODARG");
	break;
    case MEAS_INTENSITY:
	pp.PutString ("DATATYPE", "INTENSITY");
	break;
    }

    if (dtype != MEAS_INTENSITY) {
	if (!pp.GetReal ("FREQ", freq) || freq < 0.0) do {
	    cout << "\nSource modulation frequency [MHz]:\n>> ";
	    cin >> freq;
	} while (freq < 0.0);
	pp.PutReal ("FREQ", freq);
    } else freq = 0.0;
}
// ============================================================================

int ScanElRegions (const Mesh &mesh, int *nregnode)
{
    int i, reg, nreg;
    for (i = 0; i < MAXREGION; i++) nregnode[i] = 0;
    for (i = 0; i < mesh.elen(); i++) {
	reg = mesh.elist[i]->Region();
	if (reg >= 0 && reg < MAXREGION) nregnode[reg]++;
    }
    for (nreg = i = 0; i < MAXREGION; i++)
	if (nregnode[i]) nreg++;
    return nreg;
}
// ============================================================================

void SelectInitialElParams (ParamParser &pp,  Mesh &mesh)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    RVector param[3];
    int i, j, k, n, p, nreg, nregnode[MAXREGION];
    const char *resetstr[3] = {"RESET_MUA", "RESET_MUS", "RESET_N"};
    const ParameterType prmtp[3] = {PRM_MUA, PRM_MUS, PRM_N};
    for (p = 0; p < 3; p++) {

	param[p].New(mesh.elen());
	if (pp.GetString (resetstr[p], cbuf)) {
	    pp.PutString (resetstr[p], cbuf);
	    if (!strcasecmp (cbuf, "MESH")) {
		param[p] = mesh.plist.Param(prmtp[p]);
	    } else if (!strncasecmp (cbuf, "HOMOG", 5)) {
		sscanf (cbuf+5, "%lf", &prm);
		param[p] = prm;
	    } else if (!strncasecmp (cbuf, "REGION_HOMOG", 12)) {
		valstr = strtok (cbuf+12, " \t");
		for (n = 0; n < MAXREGION && valstr; n++) {
		    sscanf (valstr, "%lf", reg_prm+n);
		    valstr = strtok (NULL, " \t");
		}
		nreg = ScanElRegions (mesh, nregnode);
		for (i = k = 0; k < n && i < MAXREGION; i++) {
		    if (nregnode[i]) {
			for (j = 0; j < mesh.elen(); j++)
			    if (mesh.elist[j]->Region() == i)
				param[p][j] = reg_prm[k];
			k++;
		    }
		}	     
	    } /*	    else if (!strncasecmp (cbuf, "NIM", 3)) {
		ReadNim (cbuf+4, param[p]);
		} */
	} else {
	    cout << "\nSelect initial distribution for " << resetstr[p]
		 << endl;
	    cout << "(1) Use values stored in mesh\n";
	    cout << "(2) Global homogeneous\n";
	    cout << "(3) Homogeneous in regions\n";
	    //	    cout << "(4) Nodal image file (NIM)\n";
	    cout << "[1|2|3] >> ";
	    cin >> resettp;
	    switch (resettp) {
	    case 1:
		param[p] = mesh.plist.Param(prmtp[p]);
		strcpy (cbuf, "MESH");
		break;
	    case 2:
		cout << "\nGlobal value:\n>> ";
		cin >> prm;
		param[p] = prm;
		sprintf (cbuf, "HOMOG %f", prm);
		break;
	    case 3:
		nreg = ScanElRegions (mesh, nregnode);
		strcpy (cbuf, "REGION_HOMOG");
		cout << "\nFound " << nreg << " regions\n";
		for (i = 0; i < MAXREGION; i++) {
		    if (nregnode[i]) {
			cout << "Value for region " << i << " (" << nregnode[i]
			     << " elements):\n>> ";
			cin >> prm;
			sprintf (cbuf+strlen(cbuf), " %f", prm);
			for (j = 0; j < mesh.elen(); j++)
			    if (mesh.elist[j]->Region() == i)
				param[p][j] = prm;
		    }
		}
		break;
		/*
	    case 4:
		cout << "\nNIM file name:\n>> ";
		strcpy (cbuf, "NIM ");
		cin >> (cbuf+4);
		ReadNim (cbuf+4, param[p]);
		break;
		*/
	    }
	    pp.PutString (resetstr[p], cbuf);
	}
    }
    /*
    mesh.plist.SetParam (OT_CMUA,   param[0]*c0/param[2]);
    mesh.plist.SetParam (OT_CKAPPA, c0/(3.0*param[2]*(param[0]+param[1])));
    mesh.plist.SetParam (OT_N, param[2]);
    for (i = 0; i < param[OT_C2A].Dim(); i++)
	param[OT_C2A][i] = c0/(2*param[2][i]*A_Keijzer(param[OT_C2A][i]));
    mesh.plist.SetParam (OT_C2A, param[OT_C2A]);
    */
    mesh.plist.SetMua (param[0]);
    mesh.plist.SetMus (param[1]);
    mesh.plist.SetN (param[2]);
}

double calcbeta(double n1, double n2) 
{
  double RU21=0, RJ12=0,RJ21=0;
  double cinc = 1.0/1000.0;
  for (double c = 0.5*cinc; c < 1 ; c += cinc)
  {
    RU21 += (1 - rfresnel_cos(c,n2,n1))*c*cinc;
    RJ12 += 3*(1 -rfresnel_cos(c,n1,n2))*c*c*cinc;
    RJ21 += 3*(1 -rfresnel_cos(c,n2,n1))*c*c*cinc;
  }
  cout << "RU21 " << RU21 << " RJ12 " << RJ12 << " RJ21 " << RJ21 << endl; 
  return (2-RJ12 - RJ21)/RU21;
}

double rfresnel_cos(const double c1, const double n1, const double n2)
  {              
  // {              }
  //
  // c1 : cos of incident angle t1 in medium 1 (zero is normal indcidence), 
  // refractive index n1 in medium 1,
  // refractive index n2 in medium 2,
  // output is reflection coefficient
  //
  double tc1 = -1, tc2 = -1; // "critical angles" in the two media.
  double t1 = acos(c1);
  double s1 = sin(t1);
  double s2,t2,c2,R;
  if(n1 < n2)  
    tc2 = asin(n1/n2);
  else
    tc1 = asin(n2/n1);

  if( n1* s1/n2 > 1) // total internal reflection
  {              
    R = 1; 
  }
  else if( n1 == n2) // total internal reflection
  {              
    R = 0; 
  }
  else
  {
    s2 = s1*n1/n2;
    if (s2 <= 1)
    {
        t2 = asin(s2);
        c2 = cos(t2);
    }
    else
    {
        s2 = 1;
        c2 = 0;
    }
    double rparallel =  (n1*c2-n2*c1)/(n1*c2+n2*c1) ;
    //    tparallel =  (2*n1*c1)/(n1*c2+n2*c1) ;
    double rperp =  (n1*c1-n2*c2)/(n1*c1+n2*c2) ;
    //    tperp =  (2*n1*c1)/(n1*c1+n2*c2) ;
    R = 0.5*(rparallel*rparallel + rperp*rperp);
  }
  return R;
}


