/***************************************************************************
 * PN_varOrd.cc             Surya Mohan           11/01/2011         *
 *                                                                         *
 ***************************************************************************/
#include <pthread.h>
#include "PN_spatint.h"
#include "PN_angint.h"
#include "PN_boundint.h"
#include "quad_unitSphere.h"
#include "toast_io.h"
#include "phaseFunc.h"
#include "source.h"
#include "detector.h"

using namespace toast;
using namespace std;
/*Data context class containing necessary objects required by GMRES implicit solver*/
class MyDataContext{
public:
CCompRowMatrix SAint_co, SAintsc_co, SAintss_co, SAintc_co;
RCompRowMatrix SAintscss_co, SAintscc_co, SAintssc_co;
RCompRowMatrix Sdxx, Sdyy, Sdzz;
RCompRowMatrix SPS, SPSdx, SPSdy, SPSdz;
RCompRowMatrix Aint, Aintsc, Aintss, Aintc;
RCompRowMatrix Aintscsc,  Aintscss, Aintscc, Aintssss,  Aintssc, Aintcc;
RCompRowMatrix apu1, apu1sc, apu1ss, apu1c;
RCompRowMatrix A2, b1;
int mfc_count;
int spatN, maxSphOrder, maxAngN;
IVector node_angN, offset, lowfi_sphOrder, hifi_sphOrder;
const toast::complex *saint_coval, *saintsc_coval, *saintss_coval, *saintc_coval;
const double *saintscss_coval, *saintscc_coval, *saintssc_coval; 

const double *sdxxval, *sdyyval, *sdzzval; 
const double *spsval, *spsdxval, *spsdyval, *spsdzval; 

const double *aintval, *aintscval, *aintssval, *aintcval, *apu1val, *apu1scval, *apu1ssval, *apu1cval;
const double *aintscscval, *aintscssval, *aintssssval, *aintsccval, *aintsscval, *aintccval;
const double *a2val;
RVector ref;
double ref_out;
int BCType;
RDenseMatrix *Ylm;
double w, c;
int *xrowptr, *xcolidx;
CCompRowMatrix augA;
void initialize(QMMesh &spatMesh, const IVector& hifi_sphOrder, RVector &delta, RVector &mua, RVector &mus, RVector &ref, const double ref_out, const char BCType, const double g, toast::complex (*phaseFunc)(const double g, const double costheta), double w, double c, const RDenseMatrix& pts, const RVector &wts)
{ 
	this->w = w;
	this->c = c;
	this->hifi_sphOrder = hifi_sphOrder;
	this->ref = ref;	
	this->ref_out = ref_out;
	this->BCType = BCType; 
	maxSphOrder = vmax(hifi_sphOrder);

	/*Precomputing spherical harmonics over the quadarture points 
	which are required by the boundary integrals*/
	Ylm = new RDenseMatrix[maxSphOrder +1];
	for(int l=0; l<=maxSphOrder; l++)
		Ylm[l].New(2*l+1, pts.nRows());
	sphericalHarmonics(maxSphOrder, pts.nRows(), pts, Ylm);
	
	spatN = spatMesh.nlen();
	lowfi_sphOrder.New(spatN);
	for(int i=0; i < spatN; i++)
		lowfi_sphOrder[i] = 1;
	
	/*Computing the angular degrees of freedom at each spatial node*/
	node_angN.New(spatN);
	for(int i=0; i < spatN; i++)
	{
		node_angN[i] = 0;
		for(int l=0; l <= hifi_sphOrder[i]; l++)
			node_angN[i] += (2*l+1);
	}
	maxAngN = vmax(node_angN);

	/*Computes the start location for a given spatial node 
	in the system matrix or solution vector*/
	offset.New(spatN);
	offset[0] = 0;
	for(int i=1; i < spatN; i++)
		offset[i] = offset[i-1] + node_angN[i-1];

 
	/*Allocating memory for Ax where A is an angular matrix and x is the solution vector*/
	CCompRowMatrix Sint, Sdx, Sdy, Sdz;
	RCompRowMatrix Sx, Sy, Sz, Sdxy, Sdyx, Sdzx, Sdxz, Sdyz, Sdzy;
	RCompRowMatrix spatA3_rte, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz;	
	cout<<"Generating spatial integrals ..."<<endl;
	gen_spatint(spatMesh, mua, mus, ref, delta, w, c, Sint, Sdx, Sdy, Sdz, Sx, Sy, Sz, Sdxx, Sdxy, Sdyx, Sdyy, Sdxz, Sdzx, Sdyz, Sdzy, Sdzz, spatA3_rte, spatA3_sdmx, spatA3_sdmy, spatA3_sdmz, SPS, SPSdx, SPSdy, SPSdz);
	SAint_co = Sint + cplx(spatA3_rte); SAintsc_co = Sdx + cplx(spatA3_sdmx - Sx); SAintss_co = Sdy + cplx(spatA3_sdmy - Sy);
	SAintc_co = Sdz + cplx(spatA3_sdmz - Sz); SAintscss_co = Sdxy + Sdyx; SAintscc_co = Sdxz + Sdzx; SAintssc_co = Sdyz + Sdzy;
	
	/*Dereferencing the val pointers of spatial matrices*/
	saint_coval = SAint_co.ValPtr(); saintsc_coval = SAintsc_co.ValPtr(); saintss_coval = SAintss_co.ValPtr();
	saintc_coval = SAintc_co.ValPtr(); saintscss_coval = SAintscss_co.ValPtr(); saintscc_coval = SAintscc_co.ValPtr(); 
	saintssc_coval = SAintssc_co.ValPtr(); sdxxval = Sdxx.ValPtr(); sdyyval = Sdyy.ValPtr(); sdzzval = Sdzz.ValPtr(); 
	spsval = SPS.ValPtr(); spsdxval = SPSdx.ValPtr(); spsdyval = SPSdy.ValPtr(); spsdzval = SPSdz.ValPtr(); 
	
	cout<<"Generating angular integrals ..."<<endl;
      	genmat_angint(maxSphOrder, maxAngN, Aint, Aintsc, Aintss, Aintc, Aintscsc, Aintscss, Aintscc,  Aintssss, Aintssc, Aintcc);
	
	cout<<"Generating phase integrals ..."<<endl;	
	genmat_apu(phaseFunc, g, maxAngN, maxSphOrder, apu1, apu1sc, apu1ss, apu1c);

	//Writing Sint Aint and apu1 for Jacobian Computation

	cout<<"Generating boundary integrals ..."<<endl;//slow process
	genmat_boundint(spatMesh, ref, ref_out, BCType, hifi_sphOrder, node_angN, offset, pts, wts, Ylm, A2, b1);
		

	/*Preparing angular integrals for computing Kronecker products implicitly*/	
        apu1.Transpone(); apu1sc.Transpone(); apu1ss.Transpone(); apu1c.Transpone(); 	
	Aint.Transpone(); Aintsc.Transpone(); Aintss.Transpone(); Aintc.Transpone();   
	Aintscsc.Transpone(); Aintscss.Transpone(); Aintscc.Transpone(); Aintssc.Transpone(); 
	Aintcc.Transpone(); Aintssss.Transpone();

        /*Dereferencing the val pointers of angular matrices*/ 
	aintval = Aint.ValPtr(); aintscval = Aintsc.ValPtr(); aintssval = Aintss.ValPtr(); aintcval = Aintc.ValPtr();
	apu1val = apu1.ValPtr(); apu1scval = apu1sc.ValPtr(); apu1ssval = apu1ss.ValPtr(); apu1cval = apu1c.ValPtr();
	aintscscval = Aintscsc.ValPtr(); aintscssval = Aintscss.ValPtr(); aintssssval = Aintssss.ValPtr(); aintsccval = Aintscc.ValPtr();
	aintsscval = Aintssc.ValPtr(); aintccval = Aintcc.ValPtr();

	a2val = A2.ValPtr();
	initializePrecon(spatMesh, g, phaseFunc, pts, wts);


        /*The solution vector is stored in a matrix form for which memory is being allocated*/
	xrowptr = new int[spatN + 1];
	
	xrowptr[0]=0;
	for(int i=1;i < spatN+1; i++)
		xrowptr[i] = xrowptr[i-1] + node_angN[i-1];
	
	xcolidx = new int[xrowptr[spatN]];
	int k=0;
	for(int i = 0; i < spatN; i++)
	{
		for(int j=0; j < node_angN[i]; j++)
		{
			xcolidx[k] = j;
			k++;
		}
	 }

}
void initializePrecon(QMMesh &spatMesh, const double g, toast::complex (*phaseFunc)(const double g, const double costheta), const RDenseMatrix& pts, const RVector &wts)

{
	RCompRowMatrix Aint_lowfi, Aintsc_lowfi, Aintss_lowfi, Aintc_lowfi;
	RCompRowMatrix Aintscsc_lowfi,  Aintscss_lowfi, Aintscc_lowfi, Aintssss_lowfi,  Aintssc_lowfi, Aintcc_lowfi;
	RCompRowMatrix apu1_lowfi, apu1sc_lowfi, apu1ss_lowfi, apu1c_lowfi;
	RCompRowMatrix A2_lowfi, b1_lowfi;
	CCompRowMatrix A;
	int spatN, max_lowfi_sphOrder, lowfi_maxAngN;
	IVector lowfi_node_angN, lowfi_offset;
 
	max_lowfi_sphOrder = vmax(lowfi_sphOrder);

	spatN = spatMesh.nlen();

	/*Computing the angular degrees of freedom at each spatial node*/
	lowfi_node_angN.New(spatN);
	for(int i=0; i < spatN; i++)
	{
		lowfi_node_angN[i] = 0;
		for(int l=0; l <= lowfi_sphOrder[i]; l++)
			lowfi_node_angN[i] += (2*l+1);
	}
	lowfi_maxAngN = vmax(lowfi_node_angN);

	IVector idx2node(sum(lowfi_node_angN));
	int idx = 0;
	for(int i=0; i < spatN; i++)
	{
		for(int j = 0; j < lowfi_node_angN[i]; j++)
			idx2node[idx++] = i;	
	}

	/*Computes the start location for a given spatial node 
	in the system matrix or solution vector*/
	lowfi_offset.New(spatN);
	lowfi_offset[0] = 0;
	for(int i=1; i < spatN; i++)
		lowfi_offset[i] = lowfi_offset[i-1] + lowfi_node_angN[i-1];

      	genmat_angint(max_lowfi_sphOrder, lowfi_maxAngN, Aint_lowfi, Aintsc_lowfi, Aintss_lowfi, Aintc_lowfi, Aintscsc_lowfi, Aintscss_lowfi, Aintscc_lowfi,  Aintssss_lowfi, Aintssc_lowfi, Aintcc_lowfi);
	
	genmat_apu(phaseFunc, g, lowfi_maxAngN, max_lowfi_sphOrder, apu1_lowfi, apu1sc_lowfi, apu1ss_lowfi, apu1c_lowfi);

	genmat_boundint(spatMesh, ref, ref_out, BCType, lowfi_sphOrder, lowfi_node_angN, lowfi_offset, pts, wts, Ylm, A2_lowfi, b1_lowfi);
	b1_lowfi.New(0, 0);
	
	A = kron(SAint_co, cplx(Aint_lowfi)) + kron(SAintsc_co, cplx(Aintsc_lowfi)) + kron(SAintss_co, cplx(Aintss_lowfi));
	A = A + kron(SAintc_co, cplx(Aintc_lowfi)) + cplx(kron(Sdxx, Aintscsc_lowfi) + kron(SAintscss_co, Aintscss_lowfi));
	A = A + cplx(kron(Sdyy, Aintssss_lowfi) + kron(SAintscc_co, Aintscc_lowfi) + kron(SAintssc_co, Aintssc_lowfi));
	A = A + cplx(kron(Sdzz, Aintcc_lowfi));
	A = A - cplx(kron(SPS, apu1_lowfi) + kron(SPSdx, apu1sc_lowfi) + kron(SPSdy, apu1ss_lowfi) + kron(SPSdz, apu1c_lowfi));  
	A = A + cplx(A2_lowfi);	

	const int *rowptr, *colidx;
	A.GetSparseStructure(&rowptr, &colidx);

	int lowfi_sysdim = sum(lowfi_node_angN);
	int hifi_sysdim = sum(node_angN);
	int *augrowptr, *augcolidx;
	toast::complex *val;
	
	augrowptr = new int[hifi_sysdim + 1];
	augrowptr[0] = 0;
	for(int i=0; i < spatN; i++)
	{
		for(int j=0; j < lowfi_node_angN[i]; j++)
			augrowptr[offset[i]+j+1] = augrowptr[offset[i]+j]+(rowptr[lowfi_offset[i]+j+1]-rowptr[lowfi_offset[i]+j]);
		
		for(int j=lowfi_node_angN[i]; j < node_angN[i]; j++)
			augrowptr[offset[i]+j+1] = augrowptr[offset[i]+j] + 1;
	}

	augcolidx = new int[augrowptr[hifi_sysdim]];
	val = new toast::complex[augrowptr[hifi_sysdim]];
	toast::complex *aval = A.ValPtr();
	toast::complex a0;
	int idxaug=0;
	idx=0;
	for(int i=0; i < spatN; i++)
	{
		for(int j=0; j < lowfi_node_angN[i]; j++)
		{
			for(int k=rowptr[lowfi_offset[i]+j]; k<rowptr[lowfi_offset[i]+j+1]; k++)
			{
				int node_num = idx2node[colidx[idx]];
				augcolidx[idxaug] = offset[node_num] + colidx[idx] - lowfi_offset[node_num];//colidx[idx];
				val[idxaug] = aval[idx];
				idxaug++;
				idx++;
			}
			
		}
		for(int j=lowfi_node_angN[i]; j < node_angN[i]; j++)
		{	
			augcolidx[idxaug] = offset[i] + j;
			a0 = SAint_co.Get(i, i)*Aint.Get(j, j) + SAintsc_co.Get(i, i)*Aintsc.Get(j, j) + SAintss_co.Get(i, i)*Aintss.Get(j, j);
			a0 += SAintc_co.Get(i, i)*Aintc.Get(j, j) + Sdxx.Get(i, i)*Aintscsc.Get(j, j) + SAintscss_co.Get(i, i)*Aintscss.Get(j, j);
			a0 += Sdyy.Get(i, i)*Aintssss.Get(j, j) + SAintscc_co.Get(i, i)*Aintscc.Get(j, j) + SAintssc_co.Get(i, i)*Aintssc.Get(j, j);
			a0 += Sdzz.Get(i, i)*Aintcc.Get(j, j) + A2.Get(offset[i] + j, offset[i] + j); 
			a0 -= SPS.Get(i, i)*apu1.Get(j, j);
			a0 -= SPSdx.Get(i, i)*apu1sc.Get(j, j);
			a0 -= SPSdy.Get(i, i)*apu1ss.Get(j, j);
			a0 -= SPSdz.Get(i, i)*apu1c.Get(j, j);   
			val[idxaug] = a0;         
			idxaug++;
		}
	}
	augA.New(hifi_sysdim, hifi_sysdim);
	augA.Initialise(augrowptr, augcolidx, val);
	
}


~MyDataContext()
{
  delete []Ylm;
  delete []xrowptr;
  delete []xcolidx;

}

};// end MyDataContext class

// =========================================================================
// global parameters
int NUM_THREADS = 1;
pthread_mutex_t sid_mutex = PTHREAD_MUTEX_INITIALIZER;
CCompRowMatrix *Xmat;
CVector *result;
CDenseMatrix *Aintx, *Aintscx, *Aintssx, *Aintcx;
CDenseMatrix *apu1x, *apu1scx, *apu1ssx, *apu1cx;
CDenseMatrix *Aintscscx, *Aintscssx, *Aintssssx, *Aintsccx, *Aintsscx, *Aintccx;
MyDataContext ctxt;
int nQ;
CPrecon_ILU * AACP;
CVector *RHS, *Phi;
int sources_processed = -1;
double tol = 1e-09; 
QMMesh qmmesh;
NodeList &nlist=qmmesh.nlist;
ElementList &elist=qmmesh.elist;
inline CVector matrixFreeCaller(const CVector& x, void * context);
void *solveForSource(void *source_idx);
CVector getDiag(void * context);

int main (int argc, char *argv[])
{
    char cbuf[200];
    int el;
    int ind = 1;
    char mesh_fname[400], QM_fname[400], prm_fname[400], sphorder_fname[400], dirs_fname[400];
    bool specify_QM = false, specify_mesh = false, specify_prm = false, specify_nodal_sphorder = false, specify_prefix = false;
    double freq = 100, g=0;
    int srctp = 1;
    RVector *dirVec;
    char file_extn[200];
    int nM=1;
    RCompRowMatrix qvec, mvec;
    double w, c = 0.3, ref_out = 1.0;
    char BCType = 'v';
    RVector muscat;
    RVector muabs;
    RVector ref;
    IVector sphOrder;
	
    while(ind<argc)
    {
	string str = argv[ind];
	string::iterator it;
	for(it = str.begin(); it < str.end(); it++)
		*it = tolower(*it);
	if(!str.compare(string("-nt")))
	{
		NUM_THREADS = atoi(argv[ind+1]);
		cout<<"Number of threads used for this computation: "<<NUM_THREADS<<endl;
	}
	if(!str.compare(string("-mesh")))
	{
		strcpy(mesh_fname, argv[ind+1]);
		ReadMesh(mesh_fname, qmmesh);
		muabs.New(qmmesh.elen());
		muscat.New(qmmesh.elen());
		ref.New(qmmesh.elen());
		sphOrder.New(qmmesh.nlen());
		specify_mesh = true;
	}
	if(!str.compare(string("-qm")))
	{
		specify_QM = true;
		strcpy(QM_fname, argv[ind+1]);
		ReadQM(QM_fname, qmmesh, nQ, nM,  qvec, mvec);
	
	}
	if(!str.compare(string("-g")))
	{
		g = strtod(argv[ind+1], NULL);
		xASSERT(g>=-1 && g<=1, "g should strictly belong to [-1 1]");
		cout<< "g (WARNING !! This g is considered constant throughout the domain): "<< g<<endl;
    		
	}
	if(!str.compare(string("-freq")))
	{
		freq = strtod(argv[ind+1], NULL);
		cout << "Frequency: "<<freq<< " MHz"<<endl;
		w = freq * 2.0*M_PI*1e-6;

	}
	if(!str.compare(string("-sap")))
	{
		srctp = atoi(argv[ind+1]);
		xASSERT(srctp ==0 || srctp == 1 || srctp == 2, "Unknown source angular profile");
		cout << "Angular profile of the source (0. Directed 1. Cosine 2. Uncollided line source ): "<<srctp<<endl;
    	}
	if(!str.compare(string("-prefix")))
	{
		strcpy(file_extn, argv[ind+1]);
		cout<<"Output files' prefix: "<<file_extn<<endl;
		specify_prefix = true;
	}
	if(!str.compare(string("-bctype")))
	{
		BCType = argv[ind+1][0];
		BCType = tolower(BCType);
		xASSERT(BCType == 'v' || BCType == 'r', "Boundary condition type should be either v (vacuum) or r (reflection)");
		cout<<(BCType == 'v' ? "Vacuum" : "Reflection")<<" boundary conditions being applied"<<endl;
	}
	if(!str.compare(string("-ri_external")))
	{
		ref_out = strtod(argv[ind+1], NULL);
		cout<<"Refractive index for the external medium: "<<ref_out<<endl;
	}
	if(!str.compare(string("-prm")))
	{
		strcpy(prm_fname, argv[ind+1]);
		cout<<"Reading parameters from "<<prm_fname<<endl;
		ReadParams(prm_fname, muabs, muscat, ref);
		cout<<"Min and max values for absorption, scattering and refractive index: ("<<vmin(muabs)<<", "<<vmax(muabs)<<")  "<<"("<<vmin(muscat)<<", "<<vmax(muscat)<<")  "<<"("<<vmin(ref)<<", "<<vmax(ref)<<")  "<<endl; 
		specify_prm = true; 
	}
	if(!str.compare(string("-nodal_sphorder")))
	{
		strcpy(sphorder_fname, argv[ind+1]);
		cout<<"Reading nodal spherical harmonic order from "<<sphorder_fname<<endl;
		ReadSphOrder(sphorder_fname, sphOrder);
		cout<<"Min and max values for spherical harmonic orders : ("<<vmin(sphOrder)<<", "<<vmax(sphOrder)<<")"<<endl;
		specify_nodal_sphorder = true; 
	}

	ind++;
    }
    xASSERT(specify_mesh, "A valid mesh file should be specified.");
    xASSERT(specify_prm, "A valid text file containing the optical parameters (mua, mus, n) for each spatial element should be specified.");
    xASSERT(specify_nodal_sphorder, "A valid text file containing spherical harmonic order for each node should be specified.");
    xASSERT(specify_prefix, "Output files prefix needs to be specified.");
  
    IVector Nsource;
    IVector Ndetector;
    bool specify_nq = false, specify_nm = false;
    if(!specify_QM) 
    { // point sources and detectors 
	ind = 1;
	while(ind < argc)
	{
		string str = argv[ind];
		string::iterator it;
		for(it = str.begin(); it < str.end(); it++)
			*it = tolower(*it);
	 	if(!str.compare(string("-nq")))
		{
			nQ = atoi(argv[ind+1]);
			Nsource.New(nQ);
			cout<<"Number of sources: "<<nQ<<endl;
			for(int i=0; i<nQ; i++)
			{
				cout<<"Entered a node number for source number "<<i+1<<endl;
				cin >> Nsource[i];
				xASSERT(Nsource[i]>=0 && Nsource[i] < qmmesh.nlen(), "This node does not exist.");
				cout<<nlist[Nsource[i]]<<endl;
			}
			specify_nq = true;
			break;
		}
		ind++;
	 }
	 xASSERT(specify_nq, "-nQ should be specified as QM file has not been specified.");
	 ind=1;
	 while(ind < argc)
	 {
		string str = argv[ind];
		string::iterator it;
		for(it = str.begin(); it < str.end(); it++)
			*it = tolower(*it);
	 	if(!str.compare(string("-nm")))
		{
			nM = atoi(argv[ind+1]);
			Ndetector.New(nM);
			cout<<"Number of detectors: "<<nM<<endl;
			for(int i=0; i<nM; i++)
			{
				cout<<"Entered a node number for detector number "<<i+1<<endl;
				cin >> Ndetector[i];
				xASSERT(Ndetector[i]>=0 && Ndetector[i] < qmmesh.nlen(), "This node does not exist.");
				cout<<nlist[Ndetector[i]]<<endl;
			}
			specify_nm = true;
			break;
		}
		ind++;
	  }
	  xASSERT(specify_nm, "-nM should be specified as QM file has not been specified.");
	}

	dirVec = new RVector[nQ];
	for(int i=0; i < nQ; i++) dirVec[i].New(qmmesh.Dimension());

 	ind=1;
	while(ind < argc)
	{
		string str = argv[ind];
		string::iterator it;
		for(it = str.begin(); it < str.end(); it++)
			*it = tolower(*it);

		if(!str.compare(string("-dir")))
		{
			strcpy(dirs_fname, argv[ind+1]);
			ReadDirections(dirs_fname, nQ, dirVec);
			for(int i=0; i < nQ; i++)
			{
				dirVec[i] = dirVec[i]*1/(double)length(dirVec[i]); // normalize the direction vector just in case
				cout<< "The direction vector for the source for directed or uncollided case (default: Element normal pointing inwards): "<<dirVec[i]<<endl;
			}
    			 
		}
		ind++;
	}

	if(nQ < NUM_THREADS)
		NUM_THREADS = nQ;	

    //***** parameters 
    int numpts;
    numpts = 110;
    RDenseMatrix pts(numpts, 3); 
    RVector wts(numpts);
    for(int i=0; i < numpts; i++){
	wts[i] = wts17[i]; 
	for(int j=0; j < 3; j++)
		pts(i, j) = pts17[i][j];
    }	

    //smoothing parameter for the streamline diffusion modification
    RVector delta(qmmesh.elen());
    double min = 1e20, max = -1e20;
    for(el = 0; el <  qmmesh.elen(); el++)
    { 
	 double sval =  elist[el]->Size()/((muscat[el]));
	 delta[el] = sval;
    }
    cout<<"Min and max values of delta: "<<vmin(delta)<<" "<<vmax(delta)<<endl;
     ctxt.initialize(qmmesh, sphOrder, delta, muabs, muscat, ref, ref_out, BCType, g, &phaseFunc, w, c, pts, wts);

     
    CVector proj(nQ*nM); // projection data
    RCompRowMatrix *Source, *Detector;
    RCompRowMatrix *b2;
    if( !(Source = new  RCompRowMatrix [nQ]))
          cerr << "Memory Allocation error Source = new  RCompRowMatrix[nQ]\n";
    if( !(b2 = new  RCompRowMatrix [nQ]))
          cerr << "Memory Allocation error b2 = new  RCompRowMatrix[nQ]\n";
    if( !(Detector = new  RCompRowMatrix [nM]))
          cerr << "Memory Allocation error Detector = new  RCompRowMatrix[nM]\n";

    cout<<"Computing the source and detector vectors ..."<<endl;
   if(!specify_QM)	
   {
      genmat_source(ctxt.hifi_sphOrder, ctxt.node_angN, ctxt.offset, Source, qmmesh, Nsource, nQ, dirVec, srctp, pts, wts, ctxt.b1, ctxt.Ylm);
      genmat_detector(ctxt.hifi_sphOrder, ctxt.node_angN, ctxt.offset, Detector, qmmesh, Ndetector, nM, pts, wts, ctxt.Ylm);
      genmat_intsource(ctxt.hifi_sphOrder, ctxt.node_angN, ctxt.offset, b2, qmmesh, Nsource, nQ, dirVec, srctp, pts, wts, ctxt.Ylm, delta);
      if(srctp == 2) genmat_intsourceuncollided(ctxt.hifi_sphOrder, ctxt.node_angN, ctxt.offset, b2, qmmesh, Nsource, nQ, dirVec, pts, wts, muabs, muscat, delta, g);
    }
    else 
    {
      genmat_toastsource(ctxt.hifi_sphOrder,  ctxt.node_angN, ctxt.offset, Source, qmmesh, qvec, nQ, dirVec, srctp, pts, wts, ctxt.b1, ctxt.Ylm);
      genmat_toastdetector(ctxt.hifi_sphOrder,  ctxt.node_angN, ctxt.offset, Detector, qmmesh, mvec, nM, pts, wts, ctxt.Ylm);
    }

    cout << "calculating the radiance\n";
    int sysdim = qmmesh.nlen();
    int fullsysdim = sum(ctxt.node_angN);
    AACP = new  CPrecon_ILU;
    //CVector idiag = getDiag(&ctxt);
    //AACP->ResetFromDiagonal(idiag);
    char partition[] = "amd";
    AACP->Reset(ctxt.augA, 1, partition, 0.1, 5, 5);
   
    RHS = new CVector[nQ];
    Phi = new CVector [nQ];
    for (int j = 0; j < nQ ; j++) {   
      Phi[j].New(fullsysdim);
      RHS[j].New(fullsysdim);
       for(int i = 0; i < fullsysdim; i++)
	 RHS[j][i] = Source[j].Get(i,0) + b2[j].Get(i, 0);
    }

    CVector *detect = new CVector[nM];
    for(int j=0; j<nM; j++)
    {
       detect[j].New(fullsysdim);
	for(int i =0; i < fullsysdim; i++)
		detect[j][i] = Detector[j].Get(i, 0);
    }	
    int num_thread_loops = (int)ceil((double)nQ/NUM_THREADS);
    pthread_t thread[NUM_THREADS];
   
//GMRES(&matrixFreeCaller, &ctxt, RHS, Phi[0], tol, AACP, 100);
    int *thread_id = new int[NUM_THREADS];
    result = new CVector[NUM_THREADS];
    Xmat = new CCompRowMatrix[NUM_THREADS];
    Aintx = new CDenseMatrix[NUM_THREADS];
    Aintscx = new CDenseMatrix[NUM_THREADS];    
    Aintssx = new CDenseMatrix[NUM_THREADS];
    Aintcx = new CDenseMatrix[NUM_THREADS];    
    Aintscscx = new CDenseMatrix[NUM_THREADS];
    Aintscssx = new CDenseMatrix[NUM_THREADS];    
    Aintssssx = new CDenseMatrix[NUM_THREADS];
    Aintsccx = new CDenseMatrix[NUM_THREADS]; 
    Aintsscx = new CDenseMatrix[NUM_THREADS];
    Aintccx = new CDenseMatrix[NUM_THREADS];    
    apu1x = new CDenseMatrix[NUM_THREADS];
    apu1scx = new CDenseMatrix[NUM_THREADS]; 
    apu1ssx = new CDenseMatrix[NUM_THREADS];
    apu1cx = new CDenseMatrix[NUM_THREADS];  
    for(int i=0; i < NUM_THREADS; i++)
    {
	result[i].New(sum(ctxt.node_angN));
	Xmat[i].New(ctxt.spatN, ctxt.maxAngN);
	Xmat[i].Initialise(ctxt.xrowptr, ctxt.xcolidx);
	Aintx[i].New(ctxt.spatN, ctxt.maxAngN); Aintscx[i].New(ctxt.spatN, ctxt.maxAngN); Aintssx[i].New(ctxt.spatN, ctxt.maxAngN); 
	Aintcx[i].New(ctxt.spatN, ctxt.maxAngN); apu1x[i].New(ctxt.spatN, ctxt.maxAngN); apu1scx[i].New(ctxt.spatN, ctxt.maxAngN); 
	apu1ssx[i].New(ctxt.spatN, ctxt.maxAngN); apu1cx[i].New(ctxt.spatN, ctxt.maxAngN);Aintscscx[i].New(ctxt.spatN, ctxt.maxAngN); 
	Aintscssx[i].New(ctxt.spatN, ctxt.maxAngN); Aintssssx[i].New(ctxt.spatN, ctxt.maxAngN); Aintsccx[i].New(ctxt.spatN, ctxt.maxAngN);
	Aintsscx[i].New(ctxt.spatN, ctxt.maxAngN); Aintccx[i].New(ctxt.spatN, ctxt.maxAngN);
	thread_id[i] = i;
    }  
    //clock_t start = clock();
    time_t start, end;
    time(&start);
    for(int i=0; i < num_thread_loops; i++)
    {
	for(int j=0; j < NUM_THREADS; j++)
	{
		pthread_create(&thread[j], NULL, &solveForSource, (void *)&thread_id[j]);
			
	}
	for(int j=0; j<NUM_THREADS; j++)
	   pthread_join(thread[j], NULL);

    }
   time(&end);
   // clock_t end = clock();
    cout<<"The solver took "<<difftime(end, start)<<" seconds"<<endl;
    cout<<"Writing output files ..."<<endl;

    char fphi_re[300], fphi_im[300];
    strcpy(fphi_re, file_extn); strcpy(fphi_im, file_extn); 
    strcat(fphi_re, "_Phi_re.sol"); strcat(fphi_im, "_Phi_im.sol"); 
    
    char fRHS_re[300], fRHS_im[300];
    strcpy(fRHS_re, file_extn); strcpy(fRHS_im, file_extn);
    strcat(fRHS_re, "_RHS_re.sol"); strcat(fRHS_im, "_RHS_im.sol");
    
    FILE *osRHS_re, *osRHS_im, *osPhi_re, *osPhi_im;
    osRHS_re = fopen(fRHS_re, "w"); osRHS_im = fopen(fRHS_im, "w");
    osPhi_re = fopen(fphi_re, "w"); osPhi_im = fopen(fphi_im, "w");
    
    char fDet_re[300], fDet_im[300];
    strcpy(fDet_re, file_extn); strcpy(fDet_im, file_extn);
    strcat(fDet_re, "_Det_re.sol"); strcat(fDet_im, "_Det_im.sol");
    FILE *osDet_re, *osDet_im;
    osDet_re = fopen(fDet_re, "w"); osDet_im = fopen(fDet_im, "w");

    CVector *Phisum = new CVector [nQ];
    for (int j = 0; j < nQ ; j++) {   
      Phisum[j].New(sysdim);
      for(int i =0; i<fullsysdim; i++)
      {
	fprintf(osRHS_re, "%12e ", RHS[j][i].re);
	fprintf(osRHS_im, "%12e ", RHS[j][i].im);
	fprintf(osPhi_re, "%12e ", Phi[j][i].re);
	fprintf(osPhi_im, "%12e ", Phi[j][i].im);
	}
      fprintf(osRHS_re, "\n");fprintf(osRHS_im, "\n");
      fprintf(osPhi_re, "\n");fprintf(osPhi_im, "\n");
      for (int k = 0; k < sysdim; k++)
      {
	 Phisum[j][k] += Phi[j][ctxt.offset[k]]*sqrt(4*M_PI);
      }
      
    }
    fclose(osPhi_re);fclose(osPhi_im);
    fclose(osRHS_re);fclose(osRHS_im);

    for (int j = 0; j < nM ; j++) {   
	for(int i =0; i<fullsysdim; i++)
      	{
    		fprintf(osDet_re, "%12e ", detect[j][i].re);
		fprintf(osDet_im, "%12e ", detect[j][i].im);
	}
	fprintf(osDet_re, "\n");fprintf(osDet_im, "\n");
    }
    fclose(osDet_re);fclose(osDet_im);
    
    char flnmod[300], farg[300], ftime[300];
    strcpy(flnmod, file_extn); strcpy(farg, file_extn);strcpy(ftime, file_extn);
    strcat(flnmod, "_re.nim"); strcat(farg, "_im.nim"); strcat(ftime, "_time.txt");
    OpenNIM (flnmod, argv[1], sysdim);
    OpenNIM (farg, argv[1], sysdim);
    for (int i = 0; i < nQ; i++) {
	   WriteNIM (flnmod, Re(Phisum[i]), sysdim, i);
	   WriteNIM (farg, Im(Phisum[i]), sysdim, i);
    }
    cout << "  Log Mod field written to "<< flnmod << endl;
    cout << "  Arg field written to "<<farg << endl;

    FILE *fid;
    fid = fopen(ftime, "w");
    fprintf(fid, "Time taken by solver: %f\n", difftime(end, start));
    fclose(fid);
    
    delete []Phi;
    delete []RHS;
    delete []Phisum;
    delete []Source;
    delete []Detector;
    delete []b2;
    delete []detect;
    delete []result;
    delete []Xmat;
    delete []Aintx;
    delete []Aintscx;
    delete []Aintssx;
    delete []Aintcx;
    delete []Aintscscx;
    delete []Aintscssx;
    delete []Aintssssx;
    delete []Aintsccx;
    delete []Aintsscx;
    delete []Aintccx;
    delete []apu1x;
    delete []apu1scx;
    delete []apu1ssx;
    delete []apu1cx;
    delete []dirVec;
    delete []thread_id;
}

void *solveForSource(void *thread_id)
{
	//int *sid = (int *)source_idx;
        int source_num;
	int thread_num;
	pthread_mutex_lock(&sid_mutex);
	int *tid = (int *)thread_id;
	if(sources_processed < nQ-1)
	{
		sources_processed++;
		source_num = sources_processed;
	        thread_num = *tid;	
		pthread_mutex_unlock(&sid_mutex);
		cout<<"Processing source number: "<<source_num<< "  "<<thread_num<<"  "<<*tid<<endl;

	}
	else 
	{
	    pthread_mutex_unlock(&sid_mutex);
	    return 0;
	}
	 GMRES(&matrixFreeCaller, &thread_num, RHS[source_num], Phi[source_num], tol, AACP, (int)ceil(100.0/NUM_THREADS));
	
}

inline CVector matrixFreeCaller(const CVector& x, void *tid)
{
	int *thread_num = (int *)tid;  
	int dof = sum(ctxt.node_angN);
	int spatN = ctxt.spatN;
	int maxAngN = ctxt.maxAngN;
	
	int tidx = *thread_num;
	/*Implict Kronecker product implementation
	*	(S \circplus A)x = Sx_{r}A^{T}
	* where 'x_{r}' is a matrix resulting form reshaping of 'x'. 
	*/

	/*Reshaping 'x' to 'x_{r}'*/
	toast::complex *xmatval = Xmat[tidx].ValPtr();
	
	memcpy (xmatval, x.data_buffer(), dof*sizeof(toast::complex));
	int i, j, k, m, ra, ra1, ra2, rb, rb1, rb2;
    	int nr = Xmat[tidx].nRows();
    	int nc = ctxt.Aint.nCols();
   
	toast::complex *aintxval = Aintx[tidx].data_buffer(); toast::complex *aintscxval = Aintscx[tidx].data_buffer(); toast::complex *aintssxval = Aintssx[tidx].data_buffer();
	toast::complex *aintcxval = Aintcx[tidx].data_buffer(); toast::complex *apu1xval = apu1x[tidx].data_buffer(); toast::complex *apu1scxval = apu1scx[tidx].data_buffer();  
	toast::complex *apu1ssxval = apu1ssx[tidx].data_buffer(); toast::complex *apu1cxval = apu1cx[tidx].data_buffer(); 
	toast::complex *aintscscxval = Aintscscx[tidx].data_buffer(); toast::complex *aintscssxval = Aintscssx[tidx].data_buffer();  
	toast::complex *aintssssxval = Aintssssx[tidx].data_buffer(); toast::complex *aintsccxval = Aintsccx[tidx].data_buffer();  
	toast::complex *aintsscxval = Aintsscx[tidx].data_buffer();  toast::complex *aintccxval = Aintccx[tidx].data_buffer();  

	 
	/*Intialize Ax's to zero where A is the angular matrix*/	
	memset(aintxval, 0, dof*sizeof(toast::complex)); memset(aintscxval, 0, dof*sizeof(toast::complex)); 
	memset(aintssxval, 0, dof*sizeof(toast::complex)); memset(aintcxval, 0, dof*sizeof(toast::complex)); 
	memset(aintscscxval, 0, dof*sizeof(toast::complex)); memset(aintscssxval, 0, dof*sizeof(toast::complex)); 
	memset(aintsccxval, 0, dof*sizeof(toast::complex)); memset(aintsscxval, 0, dof*sizeof(toast::complex));
	memset(aintssssxval, 0, dof*sizeof(toast::complex)); memset(aintccxval, 0, dof*sizeof(toast::complex)); 
	memset(apu1xval, 0, dof*sizeof(toast::complex)); memset(apu1scxval, 0, dof*sizeof(toast::complex)); 
	memset(apu1ssxval, 0, dof*sizeof(toast::complex)); memset(apu1cxval, 0, dof*sizeof(toast::complex));
	/*Dereference the val pointers of Ax's*/	
	
	/*Computing x_{r}A^{T}*/
	toast::complex xval;
    	for (i = 0; i < nr; i++) {
    		ra1 = Xmat[tidx].rowptr[i];
		ra2 = Xmat[tidx].rowptr[i+1];
		for (ra = ra1; ra < ra2; ra++) {
			j = Xmat[tidx].colidx[ra];
			xval = xmatval[ra];

	    		rb1 = ctxt.Aint.rowptr[j];
	    		rb2 = ctxt.Aint.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aint.colidx[rb];
				aintxval[i*maxAngN + k] += xval*ctxt.aintval[rb];		
	    		}

			rb1 = ctxt.Aintsc.rowptr[j];
	    		rb2 = ctxt.Aintsc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintsc.colidx[rb];
				aintscxval[i*maxAngN + k] += xval*ctxt.aintscval[rb];		
	    		}

			rb1 = ctxt.Aintss.rowptr[j];
	    		rb2 = ctxt.Aintss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintss.colidx[rb];
				aintssxval[i*maxAngN + k] += xval*ctxt.aintssval[rb];		
	    		}
			
			rb1 = ctxt.Aintc.rowptr[j];
	    		rb2 = ctxt.Aintc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintc.colidx[rb];
				aintcxval[i*maxAngN + k] += xval*ctxt.aintcval[rb];		
	    		}

			rb1 = ctxt.apu1.rowptr[j];
	    		rb2 = ctxt.apu1.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.apu1.colidx[rb];
				apu1xval[i*maxAngN + k] += xval*ctxt.apu1val[rb];		
	    		}

			rb1 = ctxt.apu1sc.rowptr[j];
	    		rb2 = ctxt.apu1sc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.apu1sc.colidx[rb];
				apu1scxval[i*maxAngN + k] += xval*ctxt.apu1scval[rb];		
	    		}

			rb1 = ctxt.apu1ss.rowptr[j];
	    		rb2 = ctxt.apu1ss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.apu1ss.colidx[rb];
				apu1ssxval[i*maxAngN + k] += xval*ctxt.apu1ssval[rb];		
	    		}
			
			rb1 = ctxt.apu1c.rowptr[j];
	    		rb2 = ctxt.apu1c.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.apu1c.colidx[rb];
				apu1cxval[i*maxAngN + k] += xval*ctxt.apu1cval[rb];		
	    		}

			rb1 = ctxt.Aintscsc.rowptr[j];
	    		rb2 =ctxt.Aintscsc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintscsc.colidx[rb];
				aintscscxval[i*maxAngN + k] += xval*ctxt.aintscscval[rb];		
	    		}

			rb1 = ctxt.Aintscss.rowptr[j];
	    		rb2 = ctxt.Aintscss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintscss.colidx[rb];
				aintscssxval[i*maxAngN + k] += xval*ctxt.aintscssval[rb];		
	    		}

			rb1 = ctxt.Aintssss.rowptr[j];
	    		rb2 = ctxt.Aintssss.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintssss.colidx[rb];
				aintssssxval[i*maxAngN + k] += xval*ctxt.aintssssval[rb];		
	    		}
			
			rb1 = ctxt.Aintscc.rowptr[j];
	    		rb2 = ctxt.Aintscc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintscc.colidx[rb];
				aintsccxval[i*maxAngN + k] += xval*ctxt.aintsccval[rb];		
	    		}
			
			rb1 = ctxt.Aintssc.rowptr[j];
	    		rb2 = ctxt.Aintssc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintssc.colidx[rb];
				aintsscxval[i*maxAngN + k] += xval*ctxt.aintsscval[rb];		
	    		}
			
			rb1 = ctxt.Aintcc.rowptr[j];
	    		rb2 = ctxt.Aintcc.rowptr[j+1];
	    		for (rb = rb1; rb < rb2; rb++) {
	        		k = ctxt.Aintcc.colidx[rb];
				aintccxval[i*maxAngN + k] += xval*ctxt.aintccval[rb];		
	    		}


		}
    }
    /*Computing S(x_{r}A^{T})*/
    int scol;
    for(int is = 0; is < spatN; is++)
    {
	for(int ia = 0; ia < ctxt.node_angN[is]; ia++)
	{
		toast::complex temp(0, 0);
		for(int js = ctxt.SAint_co.rowptr[is]; js < ctxt.SAint_co.rowptr[is+1]; js++)
		{
			scol = ctxt.SAint_co.colidx[js];
			temp += ctxt.saint_coval[js]*aintxval[scol*maxAngN + ia];
			temp += ctxt.saintsc_coval[js]*aintscxval[scol*maxAngN + ia];
			temp += ctxt.saintss_coval[js]*aintssxval[scol*maxAngN + ia];
			temp += ctxt.saintc_coval[js]*aintcxval[scol*maxAngN +  ia];
			temp += aintscscxval[scol*maxAngN + ia]*ctxt.sdxxval[js];
			temp += aintscssxval[scol*maxAngN + ia]*ctxt.saintscss_coval[js];
			temp += aintssssxval[scol*maxAngN + ia]*ctxt.sdyyval[js];
			temp += aintsccxval[scol*maxAngN + ia]*ctxt.saintscc_coval[js];
			temp += aintsscxval[scol*maxAngN + ia]*ctxt.saintssc_coval[js];
			temp += aintccxval[scol*maxAngN + ia]*ctxt.sdzzval[js];
			temp -= apu1xval[scol*maxAngN + ia]*ctxt.spsval[js];
			temp -= apu1scxval[scol*maxAngN + ia]*ctxt.spsdxval[js];
			temp -= apu1ssxval[scol*maxAngN +  ia]*ctxt.spsdyval[js];
			temp -= apu1cxval[scol*maxAngN +  ia]*ctxt.spsdzval[js];
		}
		result[tidx][ctxt.offset[is] + ia]  = temp;
	} 
	}
	/*Computing A_{2}x explicitly where A_{2} is the matrix resulting from boundary*/
  	int r, i2;
    	toast::complex br;

    	for (r = i = 0; r < ctxt.A2.nRows();) {
		i2 = ctxt.A2.rowptr[r+1];
		for (br = toast::complex(0, 0); i < i2; i++)
	    		br += x[ctxt.A2.colidx[i]]*ctxt.a2val[i];
		result[tidx][r++] += br;
    	}

    	return result[tidx];
}
/** Computes the diagonal of the system matrix which is required for preconditioning
**/
CVector getDiag(void * context)
{
    MyDataContext *ctxt = (MyDataContext*) context;
    int spatN = ctxt->spatN; 
    int nDim = sum(ctxt->node_angN);
    CVector result(nDim);
    int arow, brow;
    toast::complex a;
    double  coeff = ctxt->w/ctxt->c;

    for(int i = 0; i < spatN; i++)
    {
	for(int j = 0; j < ctxt->node_angN[i]; j++)
	{
		a = ctxt->SAint_co.Get(i, i)*ctxt->Aint.Get(j, j);
		a += ctxt->SAintsc_co.Get(i, i)*ctxt->Aintsc.Get(j, j);
		a += ctxt->SAintss_co.Get(i, i)*ctxt->Aintss.Get(j, j);
		a += ctxt->SAintc_co.Get(i, i)*ctxt->Aintc.Get(j, j);
		a +=ctxt->SAintscss_co.Get(i, i)*ctxt->Aintscss.Get(j, j);
		a += ctxt->SAintscc_co.Get(i, i)*ctxt->Aintscc.Get(j, j);
		a += ctxt->SAintssc_co.Get(i, i)*ctxt->Aintssc.Get(j, j);
		a += ctxt->Sdxx.Get(i, i)*ctxt->Aintscsc.Get(j, j);
		a +=ctxt->Sdyy.Get(i, i)*ctxt->Aintssss.Get(j, j);	
		a += ctxt->Sdzz.Get(i, i)*ctxt->Aintcc.Get(j, j);
		a -= ctxt->SPS.Get(i, i)*ctxt->apu1.Get(j, j);
		a -= ctxt->SPSdx.Get(i, i)*ctxt->apu1sc.Get(j, j);
		a -= ctxt->SPSdy.Get(i, i)*ctxt->apu1ss.Get(j, j);
		a -= ctxt->SPSdz.Get(i, i)*ctxt->apu1c.Get(j, j);
		a += ctxt->A2.Get(ctxt->offset[i]+j, ctxt->offset[i]+j);

	
		result[ctxt->offset[i] + j] = a;
	}
     }
     return result;
}


//=========================================================================


// ============================================================================

