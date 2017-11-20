#include "matlabtoast.h"
#include "mexutil.h"
#include "gmsh_interface.h"

GmshInterface *gmsh_interface = 0;

GModel *GetGModel(const mxArray *idx)
{
    if (!mxIsUint64(idx))
        mexErrMsgTxt("GetGModel: Invalid handle format (expected uint64)");
    uint64_T modelidx = *(uint64_T*)mxGetData(idx);
    GModel *model = (GModel*)Handle2Ptr(modelidx);
    if (!model)
        mexErrMsgTxt("GetGModel: Index out of range");

    return model;
}

void MatlabToast::CreateGModel(int nlhs, mxArray *plhs[], int nrhs,
			       const mxArray *prhs[])
{
    if (!gmsh_interface)
	gmsh_interface = new GmshInterface;
    GModel *model = gmsh_interface->Create_Model();

    if (model)
	mexLock();
    
    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(model);
}

void MatlabToast::LoadGModel(int nlhs, mxArray *plhs[], int nrhs,
			       const mxArray *prhs[])
{
    if (!gmsh_interface)
	gmsh_interface = new GmshInterface;
    GModel *model = GetGModel(prhs[0]);

    char modelname[256];
    if (mxIsChar(prhs[1]))
	mxGetString(prhs[1], modelname, 256);
    else
	mexErrMsgTxt("LoadGModel: Argument 2: file name expected.");
    
    gmsh_interface->LoadModel(model, modelname);
}

void MatlabToast::MeshGModel(int nlhs, mxArray *plhs[], int nrhs,
			       const mxArray *prhs[])
{
    if (!gmsh_interface)
	gmsh_interface = new GmshInterface;
    GModel *model = GetGModel(prhs[0]);
    int dim;
    double minlength = -1.0, maxlength = -1.0;
    
    if (mxIsNumeric(prhs[1]))
	dim = (int)mxGetScalar(prhs[1]);
    else
	mexErrMsgTxt("MeshGModel: Argument 2: scalar expected (mesh dimension).");
    
    if (dim < 2 || dim > 3)
	mexErrMsgTxt("MeshGModel: Argument 2: mesh dimension 2 or 3 expected.");

    if (nrhs > 2) {
	if (mxIsNumeric(prhs[2]))
	    minlength = mxGetScalar(prhs[2]);
	else
	    mexErrMsgTxt("MeshGModel: Argument 3: min length (scalar) expected.");
    }
    if (nrhs > 3) {
	if (mxIsNumeric(prhs[3]))
	    maxlength = mxGetScalar(prhs[3]);
	else
	    mexErrMsgTxt("MeshGModel: Argument 3: max length (scalar) expected.");
    }
    
    gmsh_interface->MeshModel(model, dim, minlength, maxlength);
}

void MatlabToast::RemeshGModel(int nlhs, mxArray *plhs[], int nrhs,
			       const mxArray *prhs[])
{
    if (!gmsh_interface)
	gmsh_interface = new GmshInterface;
    GModel *model = GetGModel(prhs[0]);

    double minlength;
    double maxlength;
    if (mxIsNumeric(prhs[1]))
	minlength = mxGetScalar(prhs[1]);
    else
	mexErrMsgTxt("RemeshGModel: Argument 2: scalar expected (minlength).");
    if (mxIsNumeric(prhs[2]))
	maxlength = mxGetScalar(prhs[2]);
    else
	mexErrMsgTxt("RemeshGModel: Argument 3: scalar expected (maxlength).");
    
    gmsh_interface->RemeshModel(model, minlength, maxlength);
}

void MatlabToast::gmshWriteMesh(int nlhs, mxArray *plhs[], int nrhs,
					  const mxArray *prhs[])
{
    if (!gmsh_interface)
	gmsh_interface = new GmshInterface;
    GModel *model = GetGModel(prhs[0]);

    char fname[256];
    if (mxIsChar(prhs[1]))
	mxGetString(prhs[1], fname, 256);
    else
	mexErrMsgTxt("gmshWriteMesh: Argument 2: file name expected.");

    gmsh_interface->WriteMesh(model, fname);
}

void MatlabToast::GetToastMeshFromGModel(int nlhs, mxArray *plhs[], int nrhs,
					  const mxArray *prhs[])
{
    if (!gmsh_interface)
	gmsh_interface = new GmshInterface;
    GModel *model = GetGModel(prhs[0]);

    Mesh *mesh = gmsh_interface->Mesh_From_gMesh(model);
    mesh->Setup();

    plhs[0] = mxCreateNumericMatrix (1, 1, mxUINT64_CLASS, mxREAL);
    uint64_T *ptr = (uint64_T*)mxGetData (plhs[0]);
    *ptr = Ptr2Handle(mesh);
}
