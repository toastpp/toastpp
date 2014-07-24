#include "stoastlib.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include "slu_zdefs.h"
#include "supermatrix.h"
//#include "zsp_defs.h"
#include "fwdsolver.h"
#include "source.h"
//#include <pparse.h>
//#include <solution.h>
#include <mathlib.h>
#include <felib.h>


using namespace std;
using namespace toast;

#define MAXREGION 100

#define SOURCE_GAUSSIAN 0
#define SOURCE_COSINE 1


double lin_tol;
SourceMode srctp = SRCMODE_NEUMANN;   // source type
// local prototypes
CVector CompleteTrigSourceVector (const Mesh &mesh, int order);
void OpenNIM (const char *nimname, const char *meshname, int size);
void WriteNIM(const char *nimname, const RVector &img, int size, int no);
bool ReadNim (char *nimname, RVector &img);
void SelectSourceProfile (ParamParser &pp, int &qtype, double &qwidth,
    SourceMode &srctp);
void SelectMeasurementProfile (ParamParser &pp, int &mtype, double &mwidth);

void SelectInitialParams (ParamParser &pp, const Mesh &mesh, Solution &msol);
/*------------------------------------------------------------------*/
int main (int argc, char *argv[])
{
  /* defaults */
    QMMesh mesh;
    ifstream ifs("circle_small.msh");
    ifs >> mesh;
    ifs.close();
    ifstream qmf ("circle25_1x1.qm");
    mesh.LoadQM (qmf);
    qmf.close();

    mesh.Setup();


    mesh.InitSubdivisionSupport();
    mesh.PopulateNeighbourLists();
   
  /* set some parameters */

    const double c0 = 0.3;
    double tol = 1e-8;
    char meshname[256];
    double freq, omega;    // modulation frequency (MHz and cycles/ps)
    int qprof, mprof;      // source/measurement profile
                           // (0=point, 1=Gaussian, 2=cosine)
    double qwidth, mwidth; // source/measurement support radius [mm]

    ParamParser pp;
    if (argc > 1 && pp.Open (argv[1]))
        cout << "Reading parameters from " << argv[1] << endl;
    if (argc > 2) {
	pp.LogOpen (argv[2]);
	cout << "Writing log to " << argv[2] << endl;
    } else {
	pp.LogOpen ("test_adaptive.out");
	cout << "Writing log to test_adaptive.out" << endl;
    }

    int nQ, nM, nQM;
    CCompRowMatrix qvec, mvec;
    CVector *dphi;

    int i, j, cmd;
    double t0, dt;
    Point bbmin, bbmax;

    CFwdSolver FWS (pp);
    FWS.WriteParams (pp);
    lin_tol = FWS.GetLinSolverTol();

    //  SelectMesh (pp, meshname, mesh);
    SelectSourceProfile (pp, qprof, qwidth, srctp);
    SelectMeasurementProfile (pp, mprof, mwidth);
    int nnd = mesh.nlen();
    int dim = mesh.Dimension();
    nQ = mesh.nQ;
    nM = mesh.nM;
    nQM = mesh.nQM;
    mesh.BoundingBox (bbmin, bbmax);
    cout << nQ << " sources " << " nM " << " detectors " << endl ;
    cout << "Bounding box " << bbmin << "-" << bbmax << endl;

   omega = freq * 2.0*Pi*1e-6;

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New(nnd);

    // reset initial parameter estimates
    Solution sol(OT_NPARAM, nnd);
    SelectInitialParams (pp, mesh, sol);

    // build the source vectors
    qvec.New (nQ, nnd);
    LOGOUT1_INIT_PROGRESSBAR ("Source vectors", 50, nQ);
    for (i = 0; i < nQ; i++) {
	CVector q(nnd);
	switch (qprof) {
	case 0:
	    SetReal (q, QVec_Point (mesh, mesh.Q[i], srctp));
	    break;
	case 1:
	    SetReal (q, QVec_Gaussian (mesh, mesh.Q[i], qwidth, srctp));
	    break;
	case 2:
	    SetReal (q, QVec_Cosine (mesh, mesh.Q[i], qwidth, srctp));
	    break;
	case 3:
	    q = CompleteTrigSourceVector (mesh, i);
	    break;
	}
	qvec.SetRow (i, q);

	LOGOUT1_PROGRESS(i);
    }

    // build the measurement vectors
    mvec.New (nM, nnd);
    LOGOUT1_INIT_PROGRESSBAR ("Meas. vectors", 50, nM);
    for (i = 0; i < nM; i++) {
	CVector m(nnd);
	switch (mprof) {
	case 0:
	    SetReal (m, QVec_Point (mesh, mesh.M[i], SRCMODE_NEUMANN));
	    break;
	case 1:
	    SetReal (m, QVec_Gaussian (mesh, mesh.M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case 2:
	    SetReal (m, QVec_Cosine (mesh, mesh.M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	case 3:
	    m = CompleteTrigSourceVector (mesh, i);
	    break;
	}
	for (j = 0; j < nnd; j++) m[j] *= mesh.plist[j].C2A();
	mvec.SetRow (i, m);
	LOGOUT1_PROGRESS(i);
    }

    // allocate system matrix
    cout << endl << "Allocating system matrix" << endl;
    FWS.Allocate (mesh);

    cout << endl << endl << "----------" << endl;

    cout << "Assembling and pre-processing system matrix" << endl;
    t0 = clock();
    FWS.Reset (sol, omega);
    dt = ((double)clock()-t0)/(double)CLOCKS_PER_SEC;
    pp.PutReal("T(assembly): ", dt);
    cout << "T(assembly): " << dt << endl;

    // solve for fields
    cout << "Computing fields" << endl;
    IterativeSolverResult res;
    t0 = clock();
    FWS.CalcFields (qvec, dphi, &res);
    dt = ((double)clock()-t0)/(double)CLOCKS_PER_SEC;
    pp.PutReal("T(solve): ", dt);
    if (FWS.LinSolver() != LSOLVER_DIRECT) {
        cout << "Max iterations: " << res.it_count << endl;
	cout << "Max rel. error: " << res.rel_err << endl;
    }
    cout << "T(solve): " << dt << endl << endl;

    // output fields as NIM files
    cout << "Output nodal fields (1/0)? " << flush;
    cin >> cmd;
    if (cmd) {
        switch (FWS.GetDataScaling()) {
	case DATA_LIN:
	    OpenNIM ("phi_re.nim", meshname, nnd);
	    OpenNIM ("phi_im.nim", meshname, nnd);
	    for (i = 0; i < nQ; i++) {
	        WriteNIM ("phi_re.nim", Re(dphi[i]), nnd, i);
		WriteNIM ("phi_im.nim", Im(dphi[i]), nnd, i);
	    }
	    cout << "  Real field written to phi_re.nim" << endl;
	    cout << "  Imag field written to phi_im.nim" << endl;
	    break;
	case DATA_LOG:
	    OpenNIM ("phi_lnmod.nim", meshname, nnd);
	    OpenNIM ("phi_arg.nim", meshname, nnd);
	    for (i = 0; i < nQ; i++) {
	        WriteNIM ("phi_lnmod.nim", LogMod(dphi[i]), nnd, i);
		WriteNIM ("phi_arg.nim", Arg(dphi[i]), nnd, i);
	    }
	    cout << "  Log Mod field written to phi_lnmod.nim" << endl;
	    cout << "  Arg field written to phi_arg.nim" << endl;
	    break;
	}
    }

 /* construct an element wise error indicator */


    int  it, nel, subel;
    Point pt(2);
    pt[0] = pt[1] = 10.0;
    
    nel = mesh.elen();
    RVector eerr(nel);
    RVector verr(nel);
    RVector jerr(nel);
    RVector ephi(nel); // hold the average intensity in an element
    Point *elcen = new Point[nel];
    Element *ep;
      for (i=0; i < nel; i++ ){
      elcen[i] = Point(2);
      for (j = 0; j < mesh.elist[i]->nNode (); j++) {
	      int en = mesh.elist[i]->Node[j];
	      ephi[i] += mod(dphi[0][en]); // intensity at centre of this node
	      elcen[i] += mesh.nlist[en];
      }
      ephi[i] /= mesh.elist[i]->nNode ();
      elcen[i] /= mesh.elist[i]->nNode ();
    }
    
    for (i=0; i < nel; i++ ){
      cout <<"\nElement " << i << " centre " << elcen[i] << " size " << mesh.elist[i]->Size();
      //      mesh.plist.Param(PRM_MUA);
       for (j = 0; j < mesh.elist[i]->nNode (); j++) {
	    int en = mesh.elist[i]->Node[j];
     //	cout << " Node " << en << " mua " <<  mesh.plist[en].Mua() << " ";
	    verr[i] +=  mesh.plist[en].Mua()*mod(dphi[0][en])* mesh.elist[i]->IntF(j) ;
      }
      cout << " volume integral " << verr[i] << " average intensity " << ephi[i];
      for(j = 0; j <= mesh.elist[i]->nSide() ; j++){	
      	ep= mesh.elist[i]->SideNeighbour(j);
	if(ep != NULL) {// valid neighbour
	  //	  cout << "element number " << (ep - mesh.elist[0]);
	  cout << "\n" << *ep;
	  Point  encen(2);
	  double enphi= 0;
	  int je; // inefficient, but how to get el number from pointer?
	 
	    cout << "\n\t neighbour with " << ep->nNode() << " nodes";
	  
	  for(je = 0; je < ep->nNode(); je++) {
	      int en = ep->Node[je];
	      enphi += mod(dphi[0][en]); // intensity at centre of this node
	      encen += mesh.nlist[en];	    
	  }
	  
	  enphi /= ep->nNode ();
	  encen /= ep->nNode ();
	  double edist = elcen[i].Dist(encen);
	  double J = (enphi - ephi[i])/edist;
	  jerr[i] += J;
	  cout << " current " << J;
	}
	else
	  cout << "side " << j << " has no neighbour " << endl;
      } // end loop on element sides
      cout << "\t finished with element " << i << endl;
    } // finished loop on elements
    cout << "\nFinished volume error term\n";

    bool *refine = new bool[nel];
    refine[6] = refine[12] = 1;
    cout << "going to refine mesh\n";
    RefineTriangle3Mesh (&mesh, refine);
    delete []refine;   
    ofstream ofs1("refine_v1.msh");
    ofs1 << mesh;
   
    for (it = 0; it < 4; it++) {
	nel = mesh.elen();
	subel = mesh.ElFind (pt);
	bool *refine = new bool[nel];
	for (i = 0; i < nel; i++)
	    refine[i] = (i==subel /*|| i==2*/);

	RefineTriangle3Mesh (&mesh, refine);
	delete []refine;
        char num[20];// = {0,0,0,0,0};
       	sprintf(num,"%d",it);
        cout << "output level " << num  << endl;
	char name[20] ;
	strcpy(name,"refine");
	name[7] = name[8] = 0;
	strcat(name,num);
	strcat(name,".msh");
	//       cout << "output file name " << name << endl;
        ofstream ofs(name);
        ofs << mesh;
    }
    
    ofstream ofs("subdiv.msh");
    ofs << mesh;

    return 0;
}
/*============================================================*/
void Integrate_Lin_Cosine (double d, double a, double x0, double x1,
    double &int_cos_u0, double &int_cos_u1)
{
    double arg1 = 2.0*a / (Pi*Pi*(x0-x1));

    int_cos_u0 = arg1 * (-2*a * cos(Pi*(d-x0)/(2*a)) +
			 2*a * cos(Pi*(d-x1)/(2*a)) +
			 Pi * (x0-x1) * sin(Pi*(d-x0)/(2*a)));
    int_cos_u1 = arg1 * (2*a * cos(Pi*(d-x0)/(2*a)) -
			 2*a * cos(Pi*(d-x1)/(2*a)) -
			 Pi * (x0-x1) * sin(Pi*(d-x1)/(2*a)));
}

CVector CompleteTrigSourceVector (const Mesh &mesh, int order)
{
    // currently only works with 2D circular mesh centered at origin
    int el, sd, nnode, *node;
    int n = mesh.nlen();
    double phi0, phi1, a, f0, f1;
    Element *pel;
    CVector qvec (n);

    for (el = 0; el < mesh.elen(); el++) {
	pel = mesh.elist[el];
	xASSERT (pel->Type() == ELID_TRI3, Element type not supported);
	nnode = pel->nNode();
	node  = pel->Node;
	for (sd = 0; sd < pel->nSide(); sd++) {
	    if (!pel->IsBoundarySide (sd)) continue;
	    Node &nd0 = mesh.nlist[node[pel->SideNode (sd, 0)]];
	    Node &nd1 = mesh.nlist[node[pel->SideNode (sd, 1)]];
	    phi0 = atan2 (nd0[1], nd0[0]);
	    phi1 = atan2 (nd1[1], nd1[0]);

	    if (fabs (phi0-phi1) > Pi) {
		if (phi1 > phi0) phi0 += 2.0*Pi;
		else             phi1 += 2.0*Pi;
	    }
	    if (order) {
		a    = 2.0*Pi/4.0/order;
		Integrate_Lin_Cosine (0, a, phi0, phi1, f0, f1);
	    } else {
		f0 = f1 = 0.0;
	    }
	    f0 += fabs (phi1-phi0);
	    f1 += fabs (phi1-phi0);
	    qvec[node[pel->SideNode(sd,0)]] += f0;
	    qvec[node[pel->SideNode(sd,1)]] += f1;
	}
    }
    return qvec;
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
      ofs << setprecision(10) << img[i] << ' ';
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

int ScanRegions (const Mesh &mesh, int *nregnode)
{
    int i, reg, nreg;
    for (i = 0; i < MAXREGION; i++) nregnode[i] = 0;
    for (i = 0; i < mesh.nlen(); i++) {
	reg = mesh.nlist[i].Region();
	if (reg >= 0 && reg < MAXREGION) nregnode[reg]++;
    }
    for (nreg = i = 0; i < MAXREGION; i++)
	if (nregnode[i]) nreg++;
    return nreg;
}

// ============================================================================

void SelectInitialParams (ParamParser &pp, const Mesh &mesh, Solution &msol)
{
    char cbuf[256], *valstr;
    int resettp = 0;
    double prm, reg_prm[MAXREGION];
    RVector param[3];
    int i, j, k, n, p, nreg, nregnode[MAXREGION];
    const char *resetstr[3] = {"RESET_MUA", "RESET_MUS", "RESET_N"};
    const ParameterType prmtp[3] = {PRM_MUA, PRM_MUS, PRM_N};
    for (p = 0; p < 3; p++) {

	param[p].New(mesh.nlen());
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
		nreg = ScanRegions (mesh, nregnode);
		for (i = k = 0; k < n && i < MAXREGION; i++) {
		    if (nregnode[i]) {
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = reg_prm[k];
			k++;
		    }
		}	     
	    } else if (!strncasecmp (cbuf, "NIM", 3)) {
		ReadNim (cbuf+4, param[p]);
	    }
	} else {
	    cout << "\nSelect initial distribution for " << resetstr[p]
		 << endl;
	    cout << "(1) Use values stored in mesh\n";
	    cout << "(2) Global homogeneous\n";
	    cout << "(3) Homogeneous in regions\n";
	    cout << "(4) Nodal image file (NIM)\n";
	    cout << "[1|2|3|4] >> ";
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
		nreg = ScanRegions (mesh, nregnode);
		strcpy (cbuf, "REGION_HOMOG");
		cout << "\nFound " << nreg << " regions\n";
		for (i = 0; i < MAXREGION; i++) {
		    if (nregnode[i]) {
			cout << "Value for region " << i << " (" << nregnode[i]
			     << " nodes):\n>> ";
			cin >> prm;
			sprintf (cbuf+strlen(cbuf), " %f", prm);
			for (j = 0; j < mesh.nlen(); j++)
			    if (mesh.nlist[j].Region() == i)
				param[p][j] = prm;
		    }
		}
		break;
	    case 4:
		cout << "\nNIM file name:\n>> ";
		strcpy (cbuf, "NIM ");
		cin >> (cbuf+4);
		ReadNim (cbuf+4, param[p]);
		break;
	    }
	    pp.PutString (resetstr[p], cbuf);
	}
    }
    msol.SetParam (OT_CMUA,   param[0]*c0/param[2]);
    msol.SetParam (OT_CKAPPA, c0/(3.0*param[2]*(param[0]+param[1])));
    msol.SetParam (OT_N, param[2]);
    for (i = 0; i < param[OT_C2A].Dim(); i++)
	param[OT_C2A][i] = c0/(2*param[2][i]*A_Keijzer(param[OT_C2A][i]));
    msol.SetParam (OT_C2A, param[OT_C2A]);
}
