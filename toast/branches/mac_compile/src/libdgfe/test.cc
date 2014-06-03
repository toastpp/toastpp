// =========================================================================
// toastMakeNonconformingMesh
// Creates a TOAST mesh from vertex coordinates and element connectivity
// data.
//
// RH parameters:
//     1: vertex array (double, n x 2 or n x 3, where n is number of nodes)
//     2: element index list (integer, m x v, where m is number of
//        elements, and v is max. number of nodes in an element). The list is
//        1-based.
//     3: element type list (integer, m x 1). Supported types:
//           3: 4-noded tetrahedra
//     4: array of nodal parameter values (double, n x (1-3)), optional
// LH parameters:
//     1: mesh handle (pointer)
// =========================================================================
#include <alloc.h>
#include <iostream.h>
#include <fstream>
#include <iomanip.h>
#include "mathlib.h"
#include "nonconformingMesh.h"
#include "dgfwdsolver.h"
#include "dgSource.h"
#include "util.h"

using namespace std;
void CalcSysmat (NonconformingMesh *mesh, RVector &mua, RVector &mus, RVector &ref,
		 double freq, bool elbasis, TCompRowMatrix<toast::complex> res)
{
    int n = (elbasis ? mesh->elen() : mesh->nlen());

    // Set optical coefficients
    Solution sol (OT_NPARAM, n);
    sol.SetParam (OT_CMUA, mua*c0/ref);
    sol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (int i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    sol.SetParam (OT_C2A, c2a);
    // Create forward solver to initialise system matrix
    CDGFwdSolver DGFWS(LSOLVER_DIRECT, 1e-10);
    //DGFWS.SetDataScaling (DATA_LOG);
    double omega = freq * 2.0*Pi*1e-6;
   cout<<"Memory allocation begin ..."<<endl;

    DGFWS.Allocate (*mesh);
    cout<<"Memory allocation done ..."<<endl;
    DGFWS.AssembleSystemMatrix (sol, omega, elbasis);

    // Return system matrix to MATLAB
    res = *DGFWS.F;
}

int main()
{
    int i, j, k;
    int nvtx, nel, dim, nnd0;
    double *vtx;
    float temp;
    int *idx, *etp;
    /*size_t nvtx = mxGetM(prhs[0]);
    size_t nel  = mxGetM(prhs[1]);
    size_t dim  = mxGetN(prhs[0]);
    size_t nnd0 = mxGetN(prhs[1]);
    double *vtx = mxGetPr (prhs[0]);
    double *idx = mxGetPr (prhs[1]);
    double *etp = mxGetPr (prhs[2]);*/
    FILE *f_coords, *f_connec, *f_eltype, *f_qdg;
    f_coords = fopen("temp_coords.txt", "r");
    if(f_coords)
    {
	fscanf(f_coords, "%d", &nvtx);
	fscanf(f_coords, "%d", &dim);
	int k=0;
	vtx = (double *) malloc(nvtx*dim*sizeof(double));
	for(int i=0; i<nvtx; i++)
		for(int j=0; j<dim; j++)
		{fscanf(f_coords, "%f", &temp); vtx[k] = (double)temp; k++;}
	fclose(f_coords);
    }
    else
    {
 	cerr << "mesh_coords.txt not found ..."<<endl;
	exit(1);
     }
   
   f_connec = fopen("temp_connec.txt", "r");
    if(f_connec)
    {
	fscanf(f_connec, "%d", &nel);
	fscanf(f_connec, "%d", &nnd0);
	printf("%d %d\n", nel, nnd0);
	int k=0;
	idx = (int *) malloc(nel*nnd0*sizeof(int));
	for(int i=0; i<nel; i++)
		for(int j=0; j<nnd0; j++)
			{fscanf(f_connec, "%d", &idx[k]);k++;} 	
	fclose(f_connec);
	}
    else
    {
 	cerr << "mesh_connec.txt not found ..."<<endl;
	exit(1);
     }

    /*f_eltype = fopen("tet2_eltype.txt", "r");
    if(f_eltype)
    {*/
	etp = (int *)malloc(nel*sizeof(int));
	for(int i=0; i<nel; i++)
		etp[i] = 3;
	/*	fscanf(f_eltype, "%d", &etp[i]); 	
	fclose(f_eltype);
	}
    else
    {
 	cerr << "mesh_eltype.txt not found ..."<<endl;
	exit(1);
     }*/

    CVector qDG((dim+1)*nel);
	
    /*f_qdg = fopen("qDG.txt", "r");
    if(f_qdg)
    {
	for(int i=0; i<(dim+1)*nel; i++)
		{fscanf(f_coords, "%f", &temp); qDG[i] = (double)temp;}
	fclose(f_qdg);
    }
    else
    {
 	cerr << "qDG.txt not found ..."<<endl;
	exit(1);
     }
  */

 printf("Read all files ...\n");

    NonconformingMesh *mesh = new NonconformingMesh();

    // create node list
    mesh->nlist.New ((int)nvtx);
    for (i = 0; i < (int)nvtx; i++) {
	mesh->nlist[i].New((int)dim);
	mesh->nlist[i].SetBndTp (BND_NONE); // don't know
    }
    for (i = k = 0; i < (int)nvtx; i++) {
	for (j = 0; j < (int)dim; j++) {
	    mesh->nlist[i][j] = vtx[k++];
	}
    }

    // create element list
    Element *el, **list = new Element*[nel];
    for (i = 0; i < (int)nel; i++) {
	int eltp = (int)(etp[i]+0.5);
	switch (eltp) {
	case ELID_TRI3OLD:
	    list[i] = new Triangle3old;
	    break;
	case ELID_TET4:
	    list[i] = new Tetrahedron4;
	    break;
	case ELID_WDG6:
	    list[i] = new Wedge6;
	    break;
	case ELID_VOX8:
	    list[i] = new Voxel8;
	    break;
	case ELID_TRI6:
	    list[i] = new Triangle6;
	    break;
	case ELID_TET10:
	    list[i] = new Tetrahedron10;
	    break;
	case ELID_TRI6_IP:
	    list[i] = new Triangle6_ip;
	    break;
	case ELID_TRI10:
	    list[i] = new Triangle10;
	    break;
	case ELID_TRI10_IP:
	    list[i] = new Triangle10_ip;
	    break;
	case ELID_TET10_IP:
	    list[i] = new Tetrahedron10_ip;
	    break;
	case ELID_PIX4:
	    list[i] = new Pixel4;
	    break;
	case ELID_TRI3:
	    list[i] = new Triangle3;
	    break;
	case ELID_TRI3D3:
	    list[i] = new Triangle3D3;
	    break;
	case ELID_TRI3D6:
	    list[i] = new Triangle3D6;
	    break;
	default:
	    printf ("Element type not supported!\n");
	    list[i] = 0;
	    break;
	}
    }
    mesh->elist.Clear();
    mesh->elist.AppendList ((int)nel, list);
    delete []list;
    printf("Updated eltype ...\n");

    for (i = k = 0; i < nel; i++) {
	for (j = 0; j < nnd0; j++) {
	    if (el = mesh->elist[i]) {
		if (j < el->nNode())
		    el->Node[j] = (int)(idx[k]-0.5);
	    }
	    k++;
	}
    }

    // set up mesh
    mesh->MarkBoundary();
    mesh->Setup();
    mesh->SetupEdgeTables();
    printf("Set up done ..."); 


    printf ("Mesh: %d nodes, %d elements, dimension %d number of edges %d\n", nvtx, nel, dim, mesh->iedgelen());

   
   // CVector phiDG((dim+1)*nel);
   // DGFWS.SetLinSolver("GMRES", 1e-10);
    //DGFWS.CalcField(qDG, phiDG);
    /*std::ostream f_phidg;
    f_phidg.open("phiDG.txt", ios::out);
    f_phidg << phiDG;
    f_phidg.close();	 */

    /*for(int i=0; i<(dim+1)*nel; i++)
    cout << phiDG[i]<<endl;
    cout<<vmax(phiDG)<<endl;*/
    /*Node *node = new Node(3, BND_NONE);
    (*node)[0] = 4.796393300000000e-01;  (*node)[1]= 3.800521400000000e-01; (*node)[2]= 1.954439200000000e-13;
    cout<<"Exists test  "<<mesh->nlist.Exists(*node, 1e-06)<<endl;
    cout<<mesh->nlist[49][0]<<"  "<< mesh->nlist[49][1]<<"  "<<mesh->nlist[49][2]<<endl;
    cout<<(*node)[0]<<"  "<<(*node)[1]<<"  "<<(*node)[2]<<endl;
*/
    /*ofstream os("cube.msh");
    os<<*mesh;*/

    /*int len = mesh->nlen();
    int *perm;
	perm = new int[len];
    for (int i = 0; i < len; i++) perm[i] = i;

     int res = Optimise_MMD (*mesh, perm, 0, len);*/
    std::vector< std::vector<int> > ::iterator it1, it2;
    std::vector<short> :: iterator it3;
    std::vector<int>::iterator it4;
    std::vector<int> vec1, vec2;


    //mesh->RefineTetElem(0);
   /* for(it1 = mesh->iedge_elist.begin(), it2 = mesh->iedge_nlist.begin(), it3 = mesh->iedge_state.begin(); it1 != mesh->iedge_elist.end(); it1++, it2++, it3++)
    {
	vec1 = *it1; vec2 = *it2;
	cout<<vec1[0]<<"  "<<vec1[1]<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	getchar();
	}
    getchar();

    for(it4 = mesh->bedge_elist.begin(), it2 = mesh->bedge_nlist.begin(), it3 = mesh->bedge_state.begin(); it4 != mesh->bedge_elist.end(); it4++, it2++, it3++)
    {
	vec2 = *it2;
	cout<<*it4<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	getchar();
	}
    getchar();
*/
   //mesh->RefineTetElem(1);
    /*for(it1 = mesh->iedge_elist.begin(), it2 = mesh->iedge_nlist.begin(), it3 = mesh->iedge_state.begin(); it1 != mesh->iedge_elist.end(); it1++, it2++, it3++)
    {
	vec1 = *it1; vec2 = *it2;
	cout<<vec1[0]<<"  "<<vec1[1]<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	}
    getchar();

    for(it4 = mesh->bedge_elist.begin(), it2 = mesh->bedge_nlist.begin(), it3 = mesh->bedge_state.begin(); it4 != mesh->bedge_elist.end(); it4++, it2++, it3++)
    {
	vec2 = *it2;
	cout<<*it4<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	}
    getchar();

*/

/*   for(int i=0; i< nel; i++){
	mesh->RefineTetElem(i);
   for(it1 = mesh->iedge_elist.begin(), it2 = mesh->iedge_nlist.begin(), it3 = mesh->iedge_state.begin(); it1 != mesh->iedge_elist.end(); it1++, it2++, it3++)
    {
	vec1 = *it1; vec2 = *it2;
	cout<<vec1[0]<<"  "<<vec1[1]<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	}
    getchar();}

    
    nel = mesh->elist.Len();
    cout<<"Elements "<<endl;
    for(int i=0; i<nel; i++)
	cout<<i<<" "<<mesh->elist[i]->Node[0]<<"  "<<mesh->elist[i]->Node[1]<<"  "<<mesh->elist[i]->Node[2]<<"  "<<mesh->elist[i]->Node[3]<<endl;
getchar();
     for(int i=0; i< nel; i++){
	mesh->RefineTetElem(i);
 for(it1 = mesh->iedge_elist.begin(), it2 = mesh->iedge_nlist.begin(), it3 = mesh->iedge_state.begin(); it1 != mesh->iedge_elist.end(); it1++, it2++, it3++)
    {
	vec1 = *it1; vec2 = *it2;
	cout<<vec1[0]<<"  "<<vec1[1]<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	}
    getchar();}*/

 /*   for(it4 = mesh->bedge_elist.begin(), it2 = mesh->bedge_nlist.begin(), it3 = mesh->bedge_state.begin(); it4 != mesh->bedge_elist.end(); it4++, it2++, it3++)
    {
	vec2 = *it2;
	cout<<*it4<<"  "<<vec2[0]<<"  "<<vec2[1]<<"  "<<vec2[2]<<"  "<< *it3<<endl;
	}
    getchar();*/

 /*nel = mesh->elist.Len();
 cout<<"Elements "<<endl;
    for(int i=0; i<nel; i++)
	cout<<i<<"  "<<mesh->elist[i]->Node[0]<<"  "<<mesh->elist[i]->Node[1]<<"  "<<mesh->elist[i]->Node[2]<<"  "<<mesh->elist[i]->Node[3]<<endl;
getchar();
*/


    cout<<"Updated number of elements  "<<nel<<endl;
    RVector mua(nel), mus(nel), ref(nel);
    double freq=100;
    TCompRowMatrix<toast::complex> S;
    for(int i=0; i<nel; i++){ mua[i] = 1.0; mus[i]=0.005; ref[i]=1.4;}
    //CalcSysmat(mesh, mua, mus, ref, freq, true, S);
   
    DGQVec_Gaussian(*mesh, mesh->nlist[1], 2.0, SRCMODE_ISOTROPIC);  
    // Set optical coefficients
    Solution sol (OT_NPARAM, nel);
    sol.SetParam (OT_CMUA, mua*c0/ref);
    sol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    sol.SetParam(OT_N, ref);
    RVector c2a(nel);
    for (int i = 0; i < nel; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    sol.SetParam (OT_C2A, c2a);
    // Create forward solver to initialise system matrix
    CDGFwdSolver DGFWS(LSOLVER_ITERATIVE, 1e-10);
    //DGFWS.SetDataScaling (DATA_LOG);
    double omega = freq * 2.0*Pi*1e-6;
    cout<<"Memory allocation begin ..."<<endl;

    DGFWS.Allocate (*mesh);
    cout<<"Memory allocation done ..."<<endl;
    DGFWS.AssembleSystemMatrix (sol, omega, true);

    /*for(int i=0; i<10; i++)
	for(int j=0; j<10; j++)
		cout<<DGFWS.F->Get(i, j)<<" ";
		cout<<endl;  */
    /*for(int i =0 ; i<10; i++)
    { 
	std::set<elems>::iterator it1;
	std::set<nodes>::iterator it2;
	it1 = mesh->iedgelist.elems.begin();
	it2 = mesh->iedgelist.nodes.begin();
	std::vector<int> vec1 = *it1;
	std::vector<int> vec2 = *it2;
        int e = vec1.back();
	vec1.pop_back();
	int el = vec1.back();
	vec1.pop_back();

	int n1 = vec2.back(); vec2.pop_back();
	int n2 = vec2.back(); vec2.pop_back();

	mexPrintf("%d %d %d %d\n",e, el, n1, n2);
	it1++;
	it2++;
    }*/
}


