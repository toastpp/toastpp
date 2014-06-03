// -*- C++ -*-
#include <Python.h>
#include <numpy/arrayobject.h>
#include "felib.h"
#include "stoastlib.h"
#include <fstream>

#include "../common/objmgr.h"

// ===========================================================================
// A python wrapper object that handled deallocation of C++ allocated memory
// for setting up a scipy CSR matrix

typedef struct {
    PyObject_HEAD
    int *rp;
    int *ci;
    void *val;
} CSR_mem_wrapper;

// ===========================================================================
// Globals

ObjectManager<Mesh> g_meshmgr;
ObjectManager<Raster> g_rastermgr;
ObjectManager<Regularisation> g_regmgr;

// ===========================================================================
// Helper functions

// Copy a row or column vector from python to toast
RVector CopyVector (PyObject *pyvec)
{
    npy_intp *dims = PyArray_DIMS(pyvec);
    int dim = dims[0]*dims[1];
    // should check that one of the dimensions is 1

    RVector v(dim);
    memcpy (v.data_buffer(), PyArray_DATA(pyvec), dim*sizeof(double));
    return v;
}

// Copy a toast vector to a 1D python array object
void CopyVector (PyObject **pyvec, const RVector &vec)
{
    npy_intp dim = vec.Dim();
    *pyvec = PyArray_SimpleNew (1, &dim, NPY_DOUBLE);
    double *pydata = (double*)PyArray_DATA (*pyvec);
    const double *data = vec.data_buffer();
    memcpy (pydata, data, dim*sizeof(double));
}

// ===========================================================================

static PyObject *mesh_read (PyObject *self, PyObject *args)
{
    const char *meshname;
    
    if (!PyArg_ParseTuple (args, "s", &meshname))
        return NULL;
    
    QMMesh *mesh = new QMMesh;
    std::ifstream ifs (meshname);
    ifs >> *mesh;
    
    mesh->Setup();
    int hmesh = g_meshmgr.Add (mesh);

    return Py_BuildValue ("i", hmesh);
}

// ===========================================================================

static PyObject *mesh_write (PyObject *self, PyObject *args)
{
    const char *meshname;
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "is", &hmesh, &meshname))
        return NULL;

    if ((mesh = g_meshmgr.Get(hmesh))) {
        std::ofstream ofs (meshname);
	ofs << *mesh;
    }

    Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *mesh_clear (PyObject *self, PyObject *args)
{
    int hmesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;

    g_meshmgr.Delete (hmesh);
    Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *mesh_data (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;

    int i, j;
    int nlen = mesh->nlen();
    int elen = mesh->elen();
    int dim  = mesh->Dimension();
    npy_intp node_dims[2] = {nlen, dim};

    PyObject *nodelist = PyArray_SimpleNew (2, node_dims, NPY_DOUBLE);
    double *v, *vtx_data = (double*)PyArray_DATA (nodelist);

    for (i = 0, v = vtx_data; i < nlen; i++)
        for (j = 0; j < dim; j++)
	    *v++ =  mesh->nlist[i][j];

    // max number of nodes per element
    int nnd = mesh->elist[0]->nNode();
    for (i = 1; i < elen; i++)
	nnd = std::max (nnd, mesh->elist[i]->nNode());
    npy_intp el_dims[2] = {elen, nnd};

    PyObject *idx = PyArray_SimpleNew (2, el_dims, PyArray_INT32);
    int *e, *el_data = (int*)PyArray_DATA (idx);

    // element index list
    // (0-based; value -1 indicates unused matrix entry)
    for (i = 0, e = el_data; i < elen; i++) {
        for (j = 0; j < mesh->elist[i]->nNode(); j++)
	    *e++ = mesh->elist[i]->Node[j];
	for (; j < nnd; j++)
	    *e++ = -1;
    }

    // element type list (see element.h)
    PyObject *eltp = PyArray_SimpleNew (1, el_dims, PyArray_INT32);
    int *et, *etp_data = (int*)PyArray_DATA (eltp);

    for (i = 0, et = etp_data; i < elen; i++) {
        *et++ = mesh->elist[i]->Type();
    }

    return Py_BuildValue ("OOO", nodelist, idx, eltp);
}

// ===========================================================================

static PyObject *toast_surf_data (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;
    
    int i, j, k;
    int nlen = mesh->nlen();
    int dim  = mesh->Dimension();
    int nbnd = mesh->nlist.NumberOf (BND_ANY);
    int *bndidx = new int[nlen];
    npy_intp node_dims[2] = {nbnd, dim};

    double *v, *vtx_data = new double[nbnd*dim];
    for (i = 0, v = vtx_data; i < nlen; i++)
        if (mesh->nlist[i].isBnd())
	    for (j = 0; j < dim; j++)
	        *v++ = mesh->nlist[i][j];

    // vertex coordinate list
    PyObject *vtx = PyArray_SimpleNewFromData (2, node_dims, PyArray_FLOAT64,
        vtx_data);

    for (j = k = 0; j < nlen; j++)
	bndidx[j] = (mesh->nlist[j].isBnd() ? k++ : -1);
    
    // boundary element index list
    // note: this currently assumes that all elements contain the
    // same number of vertices!

    int nnd = 0, nface, sd, nn, nd, bn, *bndellist, *bndsdlist;
    nface = mesh->BoundaryList (&bndellist, &bndsdlist);
    for (j = 0; j < nface; j++)
	nnd = std::max (nnd, mesh->elist[bndellist[j]]->nSideNode(bndsdlist[j]));
    npy_intp face_dims[2] = {nface,nnd};

    int *id, *idx_data = new int[nface*nnd];
    for (j = 0, id = idx_data; j < nface; j++) {
        Element *pel = mesh->elist[bndellist[j]];
	sd = bndsdlist[j];
	nn = pel->nSideNode (sd);
	for (i = 0; i < nnd; i++) {
	    if (i < nn) {
		nd = pel->Node[pel->SideNode (sd, i)];
		bn = bndidx[nd];
	    } else bn = -1;
	    *id++ = bn;
	}
    }
    PyObject *face = PyArray_SimpleNewFromData (2, face_dims, PyArray_INT32,
        idx_data);

    // generate nodal permutation index list
    int *p, *perm_data = new int[nbnd];
    for (i = 0, p = perm_data; i < nlen; i++) {
	if (bndidx[i] >= 0)
	    *p++ = i;
    }
    npy_intp perm_dim = nbnd;
    PyObject *perm = PyArray_SimpleNewFromData (1, &perm_dim, PyArray_INT32,
        perm_data);

    // cleanup
    delete []bndidx;    

    return Py_BuildValue ("OOO", vtx, face, perm);
}

// ===========================================================================

static PyObject *toast_element_data (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;
    
    // TO BE DONE
}

// ===========================================================================

static PyObject *toast_mesh_node_count (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;

    return Py_BuildValue ("i", mesh->nlen());
}

// ===========================================================================

static PyObject *toast_mesh_element_count (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;

    return Py_BuildValue ("i", mesh->elen());
}

// ===========================================================================

static PyObject *toast_mesh_dim (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;

    return Py_BuildValue ("i", mesh->Dimension());
}

// ===========================================================================

static PyObject *toast_mesh_bb (PyObject *self, PyObject *args)
{
    int hmesh;
    Mesh *mesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;
    if (!(mesh = g_meshmgr.Get (hmesh)))
        return NULL;

    int i, dim = mesh->Dimension();
    Point pmin(dim), pmax(dim);
    mesh->BoundingBox (pmin, pmax);

    double *bb_data = new double[dim*2];
    npy_intp bb_dims[2] = {dim,2};

    for (i = 0; i < dim; i++) {
        bb_data[i*2] = pmin[i];
	bb_data[i*2+1] = pmax[i];
    }
    PyObject *bb = PyArray_SimpleNewFromData (2, bb_dims, PyArray_FLOAT64,
        bb_data);

    return Py_BuildValue ("O", bb);
}

// ===========================================================================

static PyObject *toast_basis_pos (PyObject *self, PyObject *args)
{
    int hraster;
    Raster *raster;
    RDenseMatrix pos;

    if (!PyArg_ParseTuple (args, "i", &hraster))
        return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    raster->BasisVoxelPositions (pos);
    npy_intp pos_dims[2] = {pos.nRows(), pos.nCols()};

    PyObject *pypos = PyArray_SimpleNew (2, pos_dims, PyArray_FLOAT64);
    double *data = (double*)PyArray_DATA (pypos);
    memcpy (data, pos.ValPtr(), pos.nRows()*pos.nCols()*sizeof(double));

    return Py_BuildValue ("N", pypos);
}

// ===========================================================================

static PyObject *toast_sol_pos (PyObject *self, PyObject *args)
{
    int hraster;
    Raster *raster;
    RDenseMatrix pos;

    if (!PyArg_ParseTuple (args, "i", &hraster))
        return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    raster->SolutionVoxelPositions (pos);
    npy_intp pos_dims[2] = {pos.nRows(), pos.nCols()};

    PyObject *pypos = PyArray_SimpleNew (2, pos_dims, PyArray_FLOAT64);
    double *data = (double*)PyArray_DATA (pypos);
    memcpy (data, pos.ValPtr(), pos.nRows()*pos.nCols()*sizeof(double));

    return Py_BuildValue ("N", pypos);
}

// ===========================================================================

static PyObject *toast_make_mesh (PyObject *self, PyObject *args)
{
    PyObject *py_ndlist, *py_ellist, *py_eltp;
    if (!PyArg_ParseTuple (args, "OOO", &py_ndlist, &py_ellist, &py_eltp))
        return NULL;

    Mesh *mesh = new QMMesh;
    int i, j, k;

    // create node list
    npy_intp *dims = PyArray_DIMS(py_ndlist);
    double *vtx = (double*)PyArray_DATA(py_ndlist);
    int nvtx = dims[0];
    int dim  = dims[1];

    mesh->nlist.New (nvtx);
    for (i = 0; i < nvtx; i++) {
        mesh->nlist[i].New(dim);
	mesh->nlist[i].SetBndTp (BND_NONE);
    }
    for (i = k = 0; i < nvtx; i++)
        for (j = 0; j < dim; j++)
	    mesh->nlist[i][j] = vtx[k++];

    // create element list
    dims = PyArray_DIMS(py_ellist);
    int *idx = (int*)PyArray_DATA(py_ellist);
    int *etp = (int*)PyArray_DATA(py_eltp);
    int nel = dims[0];
    int nnd0 = dims[1];

    Element *el, **list = new Element*[nel];
    for (i = 0; i < nel; i++) {
        int eltp = etp[i];
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
	    std::cerr << "Element type not supported!" << std::endl;
	    list[i] = 0;
	    break;
	}
    }
    mesh->elist.SetList (nel, list);
    delete []list;

    for (i = k = 0; i < nel; i++) {
        for (j = 0; j < nnd0; j++) {
	  if ((el = mesh->elist[i])) {
	        if (j < el->nNode())
		    el->Node[j] = idx[k];
	    }
	    k++;
	}
    }

    // create dummy parameter list
    mesh->plist.New (nvtx);
    mesh->plist.SetMua (0.01);
    mesh->plist.SetMus (1);
    mesh->plist.SetN (1);

    // set up mesh
    mesh->Setup();
    int hmesh = g_meshmgr.Add (mesh);

    return Py_BuildValue ("i", hmesh);
}

// ===========================================================================

static PyObject *toast_make_raster (PyObject *self, PyObject *args)
{
  int hmesh, hraster;
    Mesh *mesh;
    RDenseMatrix *bb = 0;
    PyObject *py_size;

    if (!PyArg_ParseTuple (args, "iO", &hmesh, &py_size))
        return NULL;
    if (!(mesh = (Mesh*)g_meshmgr.Get (hmesh)))
        return NULL;

    npy_intp *dims = PyArray_DIMS(py_size);
    int dim = dims[0];

    if (dim != mesh->Dimension())
        return NULL;

    int *size = (int*)PyArray_DATA(py_size);
    IVector bdim(dim, size);
    
    std::cerr << "bdim=" << bdim << std::endl;
    std::cerr << size[0] << ", " << size[1] << std::endl;

    Raster *raster;
    raster = new Raster_Pixel (bdim, bdim, mesh, bb);
    hraster = g_rastermgr.Add (raster);

    return Py_BuildValue ("i", hraster);
}

// ===========================================================================

static PyObject *toast_clear_raster (PyObject *self, PyObject *args)
{
    int hraster;

    if (!PyArg_ParseTuple (args, "i", &hraster))
        return NULL;

    g_rastermgr.Delete (hraster);
    Py_RETURN_NONE;
}

// ===========================================================================

template<typename T>
void MapBasis (const Raster *raster, char srcid, char tgtid,
	       const TVector<T> &src, TVector<T> &tgt)
{
    switch (srcid) {
    case 'M':
        switch (tgtid) {
	case 'G':
	    raster->Map_MeshToGrid (src, tgt);
	    break;
	case 'B':
	    raster->Map_MeshToBasis (src, tgt);
	    break;
	case 'S':
	    raster->Map_MeshToSol (src, tgt);
	    break;
	default:
	    std::cerr << "toast.MapBasis: target id not recognised"
		      << std::endl;
	}
	break;
    case 'G':
        switch (tgtid) {
	case 'M':
	    raster->Map_GridToMesh (src, tgt);
	    break;
	case 'B':
	    raster->Map_GridToBasis (src, tgt);
	    break;
	case 'S':
	    raster->Map_GridToSol (src, tgt);
	    break;
	default:
	    std::cerr << "toast.MapBasis: target id not recognised"
		      << std::endl;
	}
	break;
    case 'B':
        switch (tgtid) {
	case 'G':
	    raster->Map_BasisToGrid (src, tgt);
	    break;
	case 'S':
	    raster->Map_BasisToSol (src, tgt);
	    break;
	case 'M':
	    raster->Map_BasisToMesh (src, tgt);
	    break;
	default:
	    std::cerr << "toast.MapBasis: target id not recognised"
		      << std::endl;
	}
	break;
    case 'S':
        switch (tgtid) {
	case 'B':
	    raster->Map_SolToBasis (src, tgt);
	    break;
	case 'G':
	    raster->Map_SolToGrid (src, tgt);
	    break;
	case 'M':
	    raster->Map_SolToMesh (src, tgt);
	    break;
	default:
	    std::cerr << "toast.MapBasis: target id not recognised"
		      << std::endl;
	}
	break;
    default:
        std::cerr << "toast.MapBasis: source id not recognised"
		  << std::endl;
    }

}

static PyObject *toast_map_basis (PyObject *self, PyObject *args)
{
    int hraster;
    Raster *raster;
    const char *mapstr;
    char srcid, tgtid;
    PyObject *py_srcvec;

    if (!PyArg_ParseTuple (args, "isO", &hraster, &mapstr, &py_srcvec))
        return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    if (strlen (mapstr) != 4) {
        std::cerr << "toast.MapBasis: mapping string not recognised"
		  << std::endl;
	return NULL;
    }
    if (!strncmp (mapstr+1, "->", 2)) {
        srcid = toupper(mapstr[0]);
	tgtid = toupper(mapstr[3]);
    } else if (!strncmp (mapstr+1, "<-", 2)) {
        srcid = toupper(mapstr[3]);
	tgtid = toupper(mapstr[0]);
    } else {
        std::cerr << "toast.MapBasis: mapping string not recognised"
		  << std::endl;
	return NULL;
    }
    
    npy_intp *dims = PyArray_DIMS(py_srcvec);
    int dtype = PyArray_TYPE(py_srcvec);

    int nsrc = dims[0];
    int ntgt;
    switch (tgtid) {
    case 'M': ntgt = raster->mesh().nlen(); break;
    case 'G': ntgt = raster->GLen(); break;
    case 'B': ntgt = raster->BLen(); break;
    case 'S': ntgt = raster->SLen(); break;
    default:
	std::cerr << "toast.MapBasis: target id not recognised"
		  << std::endl;
	return NULL;
    }
    npy_intp py_nsrc = nsrc, py_ntgt = ntgt;
    PyObject *py_tgtvec = PyArray_SimpleNew (1, &py_ntgt, dtype);

    switch (dtype) {
    case NPY_DOUBLE: {
	double *src_data = (double*)PyArray_DATA(py_srcvec);
	RVector tgt, src(nsrc, src_data);
	MapBasis (raster, srcid, tgtid, src, tgt);
	double *tgt_data = (double*)PyArray_DATA(py_tgtvec);
	memcpy (tgt_data, tgt.data_buffer(), ntgt*sizeof(double));
        } break;
    case NPY_CDOUBLE: {
	toast::complex *src_data = (toast::complex*)PyArray_DATA(py_srcvec);
	CVector tgt, src(nsrc, src_data);
	MapBasis (raster, srcid, tgtid, src, tgt);
	toast::complex *tgt_data = (toast::complex*)PyArray_DATA(py_tgtvec);
	memcpy (tgt_data, tgt.data_buffer(), ntgt*sizeof(toast::complex));
        } break;
#ifdef UNDEF
    case NPY_FLOAT: {
	float *src_data = (float*)PyArray_DATA(py_srcvec);
	FVector tgt, src(nsrc, src_data);
	MapBasis (raster, srcid, tgtid, src, tgt);
	float *tgt_data = (float*)PyArray_DATA(py_tgtvec);
	memcpy (tgt_data, tgt.data_buffer(), ntgt*sizeof(float));
        } break;
    case NPY_CFLOAT: {
	scomplex *src_data = (scomplex*)PyArray_DATA(py_srcvec);
	SCVector tgt, src(nsrc, src_data);
	MapBasis (raster, srcid, tgtid, src, tgt);
	scomplex *tgt_data = (scomplex*)PyArray_DATA(py_tgtvec);
	memcpy (tgt_data, tgt.data_buffer(), ntgt*sizeof(scomplex));
        } break;
#endif
    default:
	std::cerr << "toast.MapBasis: vector type not recognised"
		  << std::endl;
	return NULL;

    }

    return Py_BuildValue ("O", py_tgtvec);
}

// ===========================================================================

static PyObject *toast_readqm (PyObject *self, PyObject *args)
{
    const char *qmname;
    int hmesh;
    QMMesh *mesh;

    if (!PyArg_ParseTuple (args, "is", &hmesh, &qmname))
        return NULL;
    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;

    std::ifstream ifs(qmname);
    mesh->LoadQM (ifs);

    std::cout << "QM: " << mesh->nQ << " sources, " << mesh->nM
	      << " detectors, " << mesh->nQM << " measurements" << std::endl;

    Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_read_nim (PyObject *self, PyObject *args)
{
    const char *nimname;
    char cbuf[256];
    int i, j = 0, idx, imgsize = 0;

    if (!PyArg_ParseTuple (args, "si", &nimname, &idx))
        return NULL;
    
    std::ifstream ifs(nimname);
    if (!ifs.getline (cbuf, 256)) return NULL;
    if (strcmp (cbuf, "NIM") && strcmp (cbuf, "RIM")) return false;
    do {
        ifs.getline (cbuf, 256);
	if (!strncasecmp (cbuf, "ImageSize", 9))
	    sscanf (cbuf+11, "%d", &imgsize);
    } while (strcasecmp (cbuf, "EndHeader"));
    if (!imgsize) return NULL;

    double *img = new double[imgsize];
    for (;;) {
	do {
	    ifs.getline (cbuf, 256);
	} while (ifs.good() && strncasecmp (cbuf, "Image", 5));
	if (!ifs.good()) break;
	for (i = 0; i < imgsize; i++)
	    ifs >> img[i];
	if (j++ == idx) break;
    }
    
    npy_intp dims = imgsize;
    PyObject *nim = PyArray_SimpleNewFromData (1, &dims, PyArray_FLOAT64, img);

    return Py_BuildValue ("O", nim);
}

// ===========================================================================

static PyObject *toast_write_nim (PyObject *self, PyObject *args)
{
    const char *nimname, *meshname;
    PyObject *py_nim;

    if (!PyArg_ParseTuple (args, "ssO", &nimname, &meshname, &py_nim))
        return NULL;

    std::ofstream ofs(nimname);
    ofs << "NIM" << std::endl;
    ofs << "Mesh = " << meshname << std::endl;
    ofs << "SolutionType = N/A" << std::endl;
    
    npy_intp nd = PyArray_NDIM(py_nim);
    npy_intp n = 1;
    for (int i = 0; i < nd; i++)
        n *= PyArray_DIM(py_nim,i);
    ofs << "ImageSize = " << n << std::endl;

    ofs << "EndHeader" << std::endl;

    ofs << "Image 0" << std::endl;

    ofs.precision(12);
    ofs.setf (std::ios::scientific);
    for (int i = 0; i < n; i++)
	ofs << (*(double*)PyArray_GETPTR1(py_nim, i)) << ' ';
    ofs << std::endl;

    Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_sysmat_cw (PyObject *self, PyObject *args)
{
    bool elbasis = false;  // for now

    int i, hmesh;
    QMMesh *mesh;
    PyObject *py_mua, *py_mus, *py_ref;
    
    if (!PyArg_ParseTuple (args, "iOOO", &hmesh, &py_mua, &py_mus, &py_ref))
        return NULL;
    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;

    double *mua = (double*)PyArray_DATA(py_mua);
    double *mus = (double*)PyArray_DATA(py_mus);
    double *ref = (double*)PyArray_DATA(py_ref);
    int nlen = mesh->nlen();

    RVector prm(nlen);

    // Set optical coefficients
    Solution sol (OT_NPARAM, nlen);
    for (i = 0; i < nlen; i++) prm[i] = mua[i]*c0/ref[i];
    sol.SetParam (OT_CMUA, prm);
    for (i = 0; i < nlen; i++) prm[i] = c0/(3.0*ref[i]*(mua[i]+mus[i]));
    sol.SetParam (OT_CKAPPA, prm);
    for (int i = 0; i < nlen; i++) prm[i] = c0/(2.0*ref[i]*A_Keijzer(ref[i]));
    sol.SetParam (OT_C2A, prm);

    // Create forward solver to initialise system matrix
    RFwdSolver FWS (mesh, LSOLVER_ITERATIVE, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    
    FWS.Allocate ();
    FWS.AssembleSystemMatrix (sol, 0, elbasis);

    const idxtype *rowptr, *colidx;
    npy_intp nnz = FWS.F->GetSparseStructure (&rowptr, &colidx);
    const double *Fval = FWS.F->ValPtr();
    npy_intp nrp = nlen+1;

    // Allocate the numpy arrays for the CSR matrix
    PyObject *py_rp = PyArray_SimpleNew (1, &nrp, PyArray_INT32);
    PyObject *py_ci = PyArray_SimpleNew (1, &nnz, PyArray_INT32);
    PyObject *py_vl = PyArray_SimpleNew (1, &nnz, PyArray_FLOAT64);

    // Copy the data over
    int *rp = (int*)PyArray_DATA(py_rp);
    for (i = 0; i < nrp; i++) rp[i] = rowptr[i];
    int *ci = (int*)PyArray_DATA(py_ci);
    for (i = 0; i < nnz; i++) ci[i] = colidx[i];
    double *val = (double*)PyArray_DATA(py_vl);
    for (i = 0; i < nnz; i++) val[i] = Fval[i];

    return Py_BuildValue ("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

static PyObject *toast_sysmat (PyObject *self, PyObject *args)
{
    bool elbasis = false;  // for now

    int i, hmesh;
    double freq;
    QMMesh *mesh;
    PyObject *py_mua, *py_mus, *py_ref;
    
    if (!PyArg_ParseTuple (args, "iOOOd", &hmesh, &py_mua, &py_mus, &py_ref,
        &freq))
        return NULL;
    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;

    double omega = freq * 2.0*Pi*1e-6;
    double *mua = (double*)PyArray_DATA(py_mua);
    double *mus = (double*)PyArray_DATA(py_mus);
    double *ref = (double*)PyArray_DATA(py_ref);
    int nlen = mesh->nlen();

    RVector prm(nlen);

    // Set optical coefficients
    Solution sol (OT_NPARAM, nlen);
    for (i = 0; i < nlen; i++) prm[i] = mua[i]*c0/ref[i];
    sol.SetParam (OT_CMUA, prm);
    for (i = 0; i < nlen; i++) prm[i] = c0/(3.0*ref[i]*(mua[i]+mus[i]));
    sol.SetParam (OT_CKAPPA, prm);
    for (int i = 0; i < nlen; i++) prm[i] = c0/(2.0*ref[i]*A_Keijzer(ref[i]));
    sol.SetParam (OT_C2A, prm);

    // Create forward solver to initialise system matrix
    CFwdSolver FWS (mesh, LSOLVER_ITERATIVE, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    
    FWS.Allocate ();
    FWS.AssembleSystemMatrix (sol, omega, elbasis);

    const idxtype *rowptr, *colidx;
    npy_intp nnz = FWS.F->GetSparseStructure (&rowptr, &colidx);
    const toast::complex *Fval = FWS.F->ValPtr();
    npy_intp nrp = nlen+1;

    // Allocate the numpy arrays for the CSR matrix
    PyObject *py_rp = PyArray_SimpleNew (1, &nrp, NPY_INT);
    PyObject *py_ci = PyArray_SimpleNew (1, &nnz, NPY_INT);
    PyObject *py_vl = PyArray_SimpleNew (1, &nnz, NPY_CDOUBLE);

    // Copy the data over
    int *rp = (int*)PyArray_DATA(py_rp);
    for (i = 0; i < nrp; i++) rp[i] = rowptr[i];
    int *ci = (int*)PyArray_DATA(py_ci);
    for (i = 0; i < nnz; i++) ci[i] = colidx[i];
    toast::complex *val = (toast::complex*)PyArray_DATA(py_vl);
    for (i = 0; i < nnz; i++) val[i] = Fval[i];

    return Py_BuildValue ("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

void CalcQvec (const QMMesh *mesh, SourceMode qtype,
    SRC_PROFILE qprof, double qwidth, CCompRowMatrix *qvec) 
{
    int i, n, nQ;

    n = mesh->nlen();
    nQ = mesh->nQ;

    // build the source vectors
    qvec->New (nQ, n);

    for (i = 0; i < nQ; i++) {
	CVector q(n);
	switch (qprof) {
	case PROF_POINT:
	    SetReal (q, QVec_Point (*mesh, mesh->Q[i], qtype));
	    break;
	case PROF_GAUSSIAN:
	    SetReal (q, QVec_Gaussian (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COSINE:
	    SetReal (q, QVec_Cosine (*mesh, mesh->Q[i], qwidth, qtype));
	    break;
	case PROF_COMPLETETRIG:
	    std::cerr << "Not implemented" << std::endl;
	    //q = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	qvec->SetRow (i, q);
    }
}

static PyObject *toast_qvec (PyObject *self, PyObject *args, PyObject *keywds)
{
    int i, hmesh;
    const char *typestr = "Neumann";
    const char *profstr = "Gaussian";
    double qwidth = 1.0;
    QMMesh *mesh;

    static char *kwlist[] = {"mesh", "type", "shape", "width", NULL};

    if (!PyArg_ParseTupleAndKeywords (args, keywds, "i|ssd", kwlist, &hmesh,
        &typestr, &profstr, &qwidth))
        return NULL;
    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;

    SourceMode qtype = SRCMODE_NEUMANN;
    if      (!strcasecmp (typestr, "Neumann"))   qtype = SRCMODE_NEUMANN;
    else if (!strcasecmp (typestr, "Isotropic")) qtype = SRCMODE_ISOTROPIC;
    else    std::cerr << "toast.Qvec: Invalid source type" << std::endl;

    SRC_PROFILE qprof = PROF_POINT;
    if      (!strcasecmp (profstr, "Point"))     qprof = PROF_POINT;
    else if (!strcasecmp (profstr, "Gaussian"))  qprof = PROF_GAUSSIAN;
    else if (!strcasecmp (profstr, "Cosine"))    qprof = PROF_COSINE;
    else if (!strcasecmp (profstr, "TrigBasis")) qprof = PROF_COMPLETETRIG;
    else    std::cerr << "toast.Qvec: Invalid source profile" << std::endl;

    if (qprof != PROF_POINT) {
        if (qwidth <= 0.0)
	    std::cerr << "toast.Qvec: Invalid source width" << std::endl;
    }

    CCompRowMatrix qvec;
    CalcQvec (mesh, qtype, qprof, qwidth, &qvec);

    const idxtype *rowptr, *colidx;
    npy_intp nnz = qvec.GetSparseStructure (&rowptr, &colidx);
    const toast::complex *qval = qvec.ValPtr();
    npy_intp nrp = qvec.nRows()+1;

    // Allocate the numpy arrays for the CSR matrix
    PyObject *py_rp = PyArray_SimpleNew (1, &nrp, NPY_INT);
    PyObject *py_ci = PyArray_SimpleNew (1, &nnz, NPY_INT);
    PyObject *py_vl = PyArray_SimpleNew (1, &nnz, NPY_CDOUBLE);

    // Copy the data over
    int *rp = (int*)PyArray_DATA(py_rp);
    for (i = 0; i < nrp; i++) rp[i] = rowptr[i];
    int *ci = (int*)PyArray_DATA(py_ci);
    for (i = 0; i < nnz; i++) ci[i] = colidx[i];
    toast::complex *val = (toast::complex*)PyArray_DATA(py_vl);
    for (i = 0; i < nnz; i++) val[i] = qval[i];

    return Py_BuildValue ("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

void CalcMvec (const QMMesh *mesh, SRC_PROFILE mprof, double mwidth,
    CCompRowMatrix *mvec) 
{
    int n, nM;
    int i, j;

    n = mesh->nlen();
    nM = mesh->nM;

    // build the measurement vectors
    mvec->New (nM, n);
    for (i = 0; i < nM; i++) {
	CVector m(n);
	switch (mprof) {
	case PROF_GAUSSIAN:
	    SetReal (m, QVec_Gaussian (*mesh, mesh->M[i], mwidth,
				       SRCMODE_NEUMANN));
	    break;
	case PROF_COSINE:
	    SetReal (m, QVec_Cosine (*mesh, mesh->M[i], mwidth,
				     SRCMODE_NEUMANN));
	    break;
	case PROF_COMPLETETRIG:
	    std::cerr << "Not implemented" << std::endl;
	    //m = CompleteTrigSourceVector (*mesh, i);
	    break;
	}
	for (j = 0; j < n; j++) m[j] *= mesh->plist[j].C2A();
	mvec->SetRow (i, m);
    }
}

static PyObject *toast_mvec (PyObject *self, PyObject *args, PyObject *keywds)
{
    int i, hmesh;
    const char *profstr = "Gaussian";
    double mwidth = 1.0;
    QMMesh *mesh;

    static char *kwlist[] = {"mesh", "shape", "width", NULL};

    if (!PyArg_ParseTupleAndKeywords (args, keywds, "i|sd", kwlist, &hmesh,
        &profstr, &mwidth))
        return NULL;
    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;

    SRC_PROFILE mprof = PROF_POINT;
    if      (!strcasecmp (profstr, "Point"))     mprof = PROF_POINT;
    else if (!strcasecmp (profstr, "Gaussian"))  mprof = PROF_GAUSSIAN;
    else if (!strcasecmp (profstr, "Cosine"))    mprof = PROF_COSINE;
    else if (!strcasecmp (profstr, "TrigBasis")) mprof = PROF_COMPLETETRIG;
    else    std::cerr << "toast.Qvec: Invalid source profile" << std::endl;

    if (mprof != PROF_POINT) {
        if (mwidth <= 0.0)
	    std::cerr << "toast.Mvec: Invalid detector width" << std::endl;
    }

    CCompRowMatrix mvec;
    CalcMvec (mesh, mprof, mwidth, &mvec);

    const idxtype *rowptr, *colidx;
    npy_intp nnz = mvec.GetSparseStructure (&rowptr, &colidx);
    const toast::complex *mval = mvec.ValPtr();
    npy_intp nrp = mvec.nRows()+1;

    // Allocate the numpy arrays for the CSR matrix
    PyObject *py_rp = PyArray_SimpleNew (1, &nrp, NPY_INT);
    PyObject *py_ci = PyArray_SimpleNew (1, &nnz, NPY_INT);
    PyObject *py_vl = PyArray_SimpleNew (1, &nnz, NPY_CDOUBLE);

    // Copy the data over
    int *rp = (int*)PyArray_DATA(py_rp);
    for (i = 0; i < nrp; i++) rp[i] = rowptr[i];
    int *ci = (int*)PyArray_DATA(py_ci);
    for (i = 0; i < nnz; i++) ci[i] = colidx[i];
    toast::complex *val = (toast::complex*)PyArray_DATA(py_vl);
    for (i = 0; i < nnz; i++) val[i] = mval[i];

    return Py_BuildValue ("NNN", py_rp, py_ci, py_vl);
}

// ===========================================================================

void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol,
    PyObject **dfield, PyObject **afield)
{
    const double c0 = 0.3;
    int i, j, n, dim, nQ, nM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    slen = (raster ? raster->SLen() : n);

    CVector *dphi, *aphi;
    CFwdSolver FWS (mesh, solver, tol);

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    FWS.SetDataScaling (DATA_LOG);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);
    CVector sphi(slen);
    
    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);

    FWS.CalcFields (qvec, dphi);

    npy_intp dmx_dims[2] = {slen,nQ};
    PyObject *dmx = PyArray_SimpleNew (2, dmx_dims, NPY_CDOUBLE);
    toast::complex *dmx_ptr = (toast::complex*)PyArray_DATA (dmx);
    for (i = 0; i < nQ; i++) {
        toast::complex *dp = dmx_ptr + i;
	if (raster) raster->Map_MeshToSol (dphi[i], sphi);
	else        sphi = dphi[i];
	for (j = 0; j < slen; j++) {
	    *dp = sphi[j];
	    dp += nQ;
	}
    }
    delete []dphi;
    *dfield = dmx;

    // build adjoint field vectors
    if (afield) {
	aphi = new CVector[nM];
	for (i = 0; i < nM; i++) aphi[i].New (n);
	FWS.CalcFields (mvec, aphi);

	npy_intp amx_dims[2] = {slen,nM};
	PyObject *amx = PyArray_SimpleNew (2, amx_dims, NPY_CDOUBLE);
	toast::complex *amx_ptr = (toast::complex*)PyArray_DATA (amx);
	for (i = 0; i < nM; i++) {
	    toast::complex *ap = amx_ptr + i;
	    if (raster) raster->Map_MeshToSol (aphi[i], sphi);
	    else        sphi = aphi[i];
	    for (j = 0; j < slen; j++) {
		*ap = sphi[j];
		ap += nM;
	    }
	}
	delete []aphi;
	*afield = amx;
    }
}                                                    

static PyObject *toast_fields (PyObject *self, PyObject *args)
{
    int hmesh, hraster;
    double freq;
    QMMesh *mesh;
    Raster *raster;
    const char *modestr = "d";
    bool do_direct = false;
    bool do_adjoint = false;

    PyObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;
    PyObject *py_mvec_vl, *py_mvec_rp, *py_mvec_ci;
    PyObject *py_mua, *py_mus, *py_ref;

    if (!PyArg_ParseTuple (args, "iiOOOOOOOOOds",
	&hmesh, &hraster,
        &py_qvec_vl, &py_qvec_rp, &py_qvec_ci,
	&py_mvec_vl, &py_mvec_rp, &py_mvec_ci,
	&py_mua, &py_mus, &py_ref,
	&freq, &modestr))
        return NULL;

    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;
    if (hraster < 0)
	raster = NULL;
    else if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    for (const char *c = modestr; *c; c++) {
	if (toupper(*c) == 'D') do_direct = true;
	if (toupper(*c) == 'A') do_adjoint = true;
    }

    int n = mesh->nlen();
    int nQ = mesh->nQ;
    int nM = mesh->nM;

    int *qrowptr = (int*)PyArray_DATA(py_qvec_rp);
    int *qcolidx = (int*)PyArray_DATA(py_qvec_ci);
    toast::complex *qval = (toast::complex*)PyArray_DATA(py_qvec_vl);
    CCompRowMatrix qvec(nQ, n, qrowptr, qcolidx, qval);

    int *mrowptr = (int*)PyArray_DATA(py_mvec_rp);
    int *mcolidx = (int*)PyArray_DATA(py_mvec_ci);
    toast::complex *mval = (toast::complex*)PyArray_DATA(py_mvec_vl);
    CCompRowMatrix mvec(nM, n, mrowptr, mcolidx, mval);

    double *mua_ptr = (double*)PyArray_DATA(py_mua);
    RVector mua (n, mua_ptr);

    double *mus_ptr = (double*)PyArray_DATA(py_mus);
    RVector mus (n, mus_ptr);

    double *ref_ptr = (double*)PyArray_DATA(py_ref);
    RVector ref (n, ref_ptr);

    PyObject *dfield, *afield;

    CalcFields (mesh, raster, qvec, mvec, mua, mus, ref, freq,
		"DIRECT", 1e-12, &dfield, do_adjoint ? &afield : NULL);

    PyObject *ret;

    if (do_direct) {
	if (do_adjoint)
	    ret = Py_BuildValue ("OO", dfield, afield);
	else
	    ret = Py_BuildValue ("O", dfield);
    } else {
	if (do_adjoint)
	    ret = Py_BuildValue ("O", afield);
	else
	    ret = NULL;
    }
    Py_DECREF(dfield);
    if (do_adjoint) Py_DECREF(afield);

    if (ret)
	return ret;
    else
	Py_RETURN_NONE;
}

// ===========================================================================

// Calculate Jacobian from given direct and adjoint fields and boundary
// projection data

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CVector *dphi, const CVector *aphi,
    const CVector *proj, DataScale dscale, PyObject **res)
{
    int nQM, slen, ndat, nprm;
    nQM  = mesh->nQM;
    slen = (raster ? raster->SLen() : mesh->nlen());
    ndat = nQM * 2;
    nprm = slen * 2;

    RDenseMatrix J(ndat,nprm);

    ofstream ofs("dbg_p.dat");
    ofs << dphi[0] << endl;
    ofs.close();

    GenerateJacobian (raster, mesh, dphi, aphi, proj, dscale, J);

    npy_intp J_dims[2] = {ndat, nprm};
    PyObject *pyJ = PyArray_SimpleNew (2, J_dims, NPY_DOUBLE);
    double *pyJ_data = (double*)PyArray_DATA (pyJ);
    const double *J_data = J.ValPtr();
    memcpy (pyJ_data, J_data, ndat*nprm*sizeof(double));
    *res = pyJ;
}

// Calculate Jacobian from given optical parameters.
// This version calculates the direct and adjoint fields, and boundary
// projection data on the fly from the provided optical parameters

void CalcJacobian (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, PyObject **res)
{
    const double c0 = 0.3;
    int i, n, dim, nQ, nM, nQM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    nQM  = mesh->nQM;
    slen = (raster ? raster->SLen() : n);

    CVector *dphi, *aphi;
    CFwdSolver FWS (mesh, solver, tol);

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    FWS.SetDataScaling (DATA_LOG);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    aphi = new CVector[nM];
    for (i = 0; i < nM; i++) aphi[i].New (n);

    // Calculate direct and adjoint fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    FWS.CalcFields (mvec, aphi);

    // Calculate projections if required
    CVector *proj = 0;
    DataScale dscale = FWS.GetDataScaling();
    if (dscale == DATA_LOG) {
    	proj = new CVector(nQM);
	*proj = FWS.ProjectAll (mvec, dphi, DATA_LIN);
    	//ProjectAll (*mesh, FWS, mvec, dphi, *proj);
    }

    ofstream ofs("dbg_p.dat");
    ofs << *proj << endl;
    ofs.close();

    // Calculate Jacobian
    CalcJacobian (mesh, raster, dphi, aphi, proj, dscale, res);

    delete []dphi;
    delete []aphi;
    if (proj) delete proj;
}                                                                              

static PyObject *toast_jacobian (PyObject *self, PyObject *args)
{
    int hmesh, hraster, i, j;
    QMMesh *mesh;
    Raster *raster;
    PyObject *py_dphi, *py_aphi, *py_proj;

    if (!PyArg_ParseTuple (args, "iiOOO", &hmesh, &hraster,
			   &py_dphi, &py_aphi, &py_proj))
	return NULL;

    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    int n   = mesh->nlen();
    int nq  = mesh->nQ;
    int nm  = mesh->nM;
    int nqm = mesh->nQM;

    // set up direct fields
    toast::complex *dphi_ptr = (toast::complex*)PyArray_DATA(py_dphi);
    CVector *dphi = new CVector[nq];
    for (i = 0; i < nq; i++) {
	dphi[i].Relink (dphi_ptr + i*n, n);
    }

    // set up adjoint fields
    toast::complex *aphi_ptr = (toast::complex*)PyArray_DATA(py_aphi);
    CVector *aphi = new CVector[nm];
    for (i = 0; i < nm; i++) {
	aphi[i].Relink (aphi_ptr + i*n, n);
    }

    // copy projections
    toast::complex *proj_ptr = (toast::complex*)PyArray_DATA(py_proj);
    CVector proj(nqm, proj_ptr, SHALLOW_COPY);

    PyObject *J;
    CalcJacobian (mesh, raster, dphi, aphi, &proj, DATA_LOG, &J);

    // cleanup
    delete []dphi;
    delete []aphi;

    PyObject *res = Py_BuildValue ("O", J);
    Py_DECREF(J);
    return res;
}

// ===========================================================================
// Calculate the gradient of the forward operator

void AddDataGradient (QMMesh *mesh, Raster *raster, const CFwdSolver &FWS,
    const RVector &proj, const RVector &data, const RVector &sd, CVector *dphi,
    const CCompRowMatrix &mvec, RVector &grad)
{
    int i, j, q, m, n, idx, ofs_mod, ofs_arg;
    double term;
    int glen = raster->GLen();
    int slen = raster->SLen();
    int dim  = raster->Dim();
    const IVector &gdim = raster->GDim();
    const RVector &gsize = raster->GSize();
    const int *elref = raster->Elref();
    CVector wqa (mesh->nlen());
    RVector wqb (mesh->nlen());
    CVector dgrad (slen);
    ofs_mod = 0;
    ofs_arg = mesh->nQM;
    RVector grad_cmua (grad, 0, slen);
    RVector grad_ckappa (grad, slen, slen);

    for (q = 0; q < mesh->nQ; q++) {

	// expand field and gradient
        CVector cdfield (glen);
        CVector *cdfield_grad = new CVector[dim];
	raster->Map_MeshToGrid (dphi[q], cdfield);
	ImageGradient (gdim, gsize, cdfield, cdfield_grad, elref);

        n = mesh->nQMref[q];

	RVector y_mod (data, ofs_mod, n);
	RVector s_mod (sd, ofs_mod, n);
	RVector ypm_mod (proj, ofs_mod, n);
	RVector b_mod(n);
	b_mod = (y_mod-ypm_mod)/s_mod;

	RVector y_arg (data, ofs_arg, n);
	RVector s_arg (sd, ofs_arg, n);
	RVector ypm_arg (proj, ofs_arg, n);
	RVector b_arg(n);
	b_arg = (y_arg-ypm_arg)/s_arg;

	RVector ype(n);
	ype = 1.0;  // will change if data type is normalised

	CVector cproj(n);
	//Project_cplx (*mesh, q, dphi[q], cproj);
	cproj = ProjectSingle (mesh, q, mvec, dphi[q]);
	wqa = toast::complex(0,0);
	wqb = 0.0;

	for (m = idx = 0; m < mesh->nM; m++) {
	    if (!mesh->Connected (q, m)) continue;
	    const CVector qs = mvec.Row(m);
	    double rp = cproj[idx].re;
	    double ip = cproj[idx].im;
	    double dn = 1.0/(rp*rp + ip*ip);

	    // amplitude term
	    term = -2.0 * b_mod[idx] / (ype[idx]*s_mod[idx]);
	    wqa += qs * toast::complex (term*rp*dn, -term*ip*dn);

	    // phase term
	    term = -2.0 * b_arg[idx] / (ype[idx]*s_arg[idx]);
	    wqa += qs * toast::complex (-term*ip*dn, -term*rp*dn);

	    //wqb += Re(qs) * (term * ypm[idx]);
	    idx++;
	}

	// adjoint field and gradient
	CVector wphia (mesh->nlen());
	FWS.CalcField (wqa, wphia);

	CVector cafield(glen);
	CVector *cafield_grad = new CVector[dim];
	raster->Map_MeshToGrid (wphia, cafield);
	ImageGradient (gdim, gsize, cafield, cafield_grad, elref);

	// absorption contribution
	raster->Map_GridToSol (cdfield * cafield, dgrad);
	grad_cmua -= Re(dgrad);

	// diffusion contribution
	// multiply complex field gradients
	CVector gk(glen);
	for (i = 0; i < glen; i++)
	    for (j = 0; j < dim; j++)
	        gk[i] += cdfield_grad[j][i] * cafield_grad[j][i];
	raster->Map_GridToSol (gk, dgrad);
	grad_ckappa -= Re(dgrad);

	delete []cdfield_grad;
	delete []cafield_grad;

	ofs_mod += n; // step to next source
	ofs_arg += n;
    }
}

void GetGradient (QMMesh *mesh, Raster *raster, CFwdSolver &FWS,
    const RVector &mua, const RVector &mus, const RVector &ref, double freq,
    const RVector &data, const RVector &sd,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    RVector &grad)
{
    const double c0 = 0.3;
    int i, n = mesh->nlen();
    int nQ = mesh->nQ;
    CVector *dphi;
    RVector proj(data.Dim());

    // Solution in mesh basis
    Solution msol(OT_NPARAM, n);
    msol.SetActive (OT_CMUA, true);
    msol.SetActive (OT_CKAPPA, true);

    // Set optical coefficients
    msol.SetParam (OT_CMUA, mua*c0/ref);
    msol.SetParam (OT_CKAPPA, c0/(3.0*ref*(mua+mus)));
    RVector c2a(n);
    for (i = 0; i < n; i++)
	c2a[i] = c0/(2.0*ref[i]*A_Keijzer (ref[i]));
    msol.SetParam (OT_C2A, c2a);

    double omega = freq * 2.0*Pi*1e-6; // convert from MHz to rad

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);

    // Calculate fields
    FWS.Allocate ();
    FWS.Reset (msol, omega);
    FWS.CalcFields (qvec, dphi);
    proj = FWS.ProjectAll_real (mvec, dphi);

    AddDataGradient (mesh, raster, FWS, proj, data, sd, dphi, mvec, grad);

    delete []dphi;
}

static PyObject *toast_gradient (PyObject *self, PyObject *args)
{
    int hmesh, hraster;
    double freq;
    QMMesh *mesh;
    Raster *raster;
    PyObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;
    PyObject *py_mvec_vl, *py_mvec_rp, *py_mvec_ci;
    PyObject *py_mua, *py_mus, *py_ref;
    PyObject *py_data, *py_sd;

    const char *solver = "direct"; // for now
    const double tol = 1e-12; // for now

    if (!PyArg_ParseTuple (args, "iiOOOOOOOOOdOO", &hmesh, &hraster,
			   &py_qvec_vl, &py_qvec_rp, &py_qvec_ci,
			   &py_mvec_vl, &py_mvec_rp, &py_mvec_ci,
			   &py_mua, &py_mus, &py_ref,
			   &freq, &py_data, &py_sd))
	return NULL;

    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
	return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
	return NULL;

    int n = mesh->nlen();
    int nQ = mesh->nQ;
    int nM = mesh->nM;
    int nQM = mesh->nQM;

    int *qrowptr = (int*)PyArray_DATA(py_qvec_rp);
    int *qcolidx = (int*)PyArray_DATA(py_qvec_ci);
    toast::complex *qval = (toast::complex*)PyArray_DATA(py_qvec_vl);
    CCompRowMatrix qvec(nQ, n, qrowptr, qcolidx, qval, SHALLOW_COPY);

    int *mrowptr = (int*)PyArray_DATA(py_mvec_rp);
    int *mcolidx = (int*)PyArray_DATA(py_mvec_ci);
    toast::complex *mval = (toast::complex*)PyArray_DATA(py_mvec_vl);
    CCompRowMatrix mvec(nM, n, mrowptr, mcolidx, mval, SHALLOW_COPY);

    double *mua_ptr = (double*)PyArray_DATA(py_mua);
    RVector mua (n, mua_ptr, SHALLOW_COPY);

    double *mus_ptr = (double*)PyArray_DATA(py_mus);
    RVector mus (n, mus_ptr, SHALLOW_COPY);

    double *ref_ptr = (double*)PyArray_DATA(py_ref);
    RVector ref (n, ref_ptr, SHALLOW_COPY);

    double *data_ptr = (double*)PyArray_DATA(py_data);
    RVector data (nQM*2, data_ptr, SHALLOW_COPY);

    double *sd_ptr = (double*)PyArray_DATA(py_sd);
    RVector sd (nQM*2, sd_ptr, SHALLOW_COPY);

    CFwdSolver FWS (mesh, solver, tol);
    FWS.SetDataScaling (DATA_LOG);

    npy_intp grad_dim = raster->SLen()*2;
    PyObject *py_grad = PyArray_SimpleNew (1, &grad_dim, NPY_DOUBLE);
    double *py_grad_data = (double*)PyArray_DATA (py_grad);
    memset (py_grad_data, 0, grad_dim*sizeof(double));
    RVector grad (grad_dim, py_grad_data, SHALLOW_COPY);
    GetGradient (mesh, raster, FWS, mua, mus, ref, freq, data, sd,
		 qvec, mvec, grad);

    PyObject *res = Py_BuildValue ("O", py_grad);
    Py_DECREF (py_grad);

    return res;
}

static PyObject *toast_krylov (PyObject *self, PyObject *args)
{
    PyObject *py_x, *py_J;

    if (!PyArg_ParseTuple (args, "OO", &py_x, &py_J))
	return NULL;

    npy_intp *J_dims = PyArray_DIMS(py_J);
    std::cerr << "Jm=" << J_dims[0] << ", Jn=" << J_dims[1] << std::endl;

    Py_RETURN_NONE;
}

// ===========================================================================
// ===========================================================================
// Regularisation methods

static PyObject *toast_regul (PyObject *self, PyObject *args, PyObject *keywds)
{
    Regularisation *reg = 0;
    const char *regtype;
    int hraster, hreg;
    double tau, beta = 1;
    void *kapref = 0;
    bool istensor = false;      // reference diffusivity in tensor format?
    Raster *raster;
    PyObject *py_x;

    static char *kwlist[] = {"regtype", "raster", "x", "tau", "beta", NULL};

    if (!PyArg_ParseTupleAndKeywords (args, keywds, "siOd|d", kwlist, &regtype,
        &hraster, &py_x, &tau, &beta))
	return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    npy_intp *x_dims = PyArray_DIMS(py_x);
    int xlen = x_dims[0]*x_dims[1];

    RVector x0(xlen);
    memcpy (x0.data_buffer(), PyArray_DATA(py_x), xlen*sizeof(double));

    if (!strcasecmp (regtype, "TK0")) {
	RVector xs (x0.Dim());
	xs = 1;
	reg = new Tikhonov0 (tau, &x0, &xs);
    } else if (!strcasecmp (regtype, "TK1")) {
	reg = new TK1 (tau, &x0, raster, kapref, istensor);
    } else if (!strcasecmp (regtype, "TV")) {
	reg = new TV (tau, beta, &x0, raster, 0, false);
    }
    
    hreg = g_regmgr.Add (reg);
    return Py_BuildValue ("i", hreg);
}

// ===========================================================================

static PyObject *toast_regul_value (PyObject *self, PyObject *args)
{
    int hreg;
    Regularisation *reg;
    PyObject *py_x;

    if (!PyArg_ParseTuple (args, "iO", &hreg, &py_x))
	return NULL;
    if (!(reg = (Regularisation*)g_regmgr.Get (hreg)))
	return NULL;

    npy_intp *x_dims = PyArray_DIMS(py_x);
    int xlen = x_dims[0]*x_dims[1];

    RVector x(xlen);
    memcpy (x.data_buffer(), PyArray_DATA(py_x), xlen*sizeof(double));

    double rval = reg->GetValue (x);
    return Py_BuildValue ("d", rval);
}

// ===========================================================================

static PyObject *toast_regul_grad (PyObject *self, PyObject *args)
{
    int hreg;
    Regularisation *reg;
    PyObject *py_x;

    if (!PyArg_ParseTuple (args, "iO", &hreg, &py_x))
	return NULL;
    if (!(reg = (Regularisation*)g_regmgr.Get (hreg)))
	return NULL;

    npy_intp *x_dims = PyArray_DIMS(py_x);
    int xlen = x_dims[0]*x_dims[1];

    RVector x(xlen);
    memcpy (x.data_buffer(), PyArray_DATA(py_x), xlen*sizeof(double));

    RVector grad(xlen);
    grad = reg->GetGradient (x);

    npy_intp grad_dim = grad.Dim();
    PyObject *py_grad = PyArray_SimpleNew (1, &grad_dim, NPY_DOUBLE);
    double *py_grad_data = (double*)PyArray_DATA (py_grad);
    const double *grad_data = grad.data_buffer();
    memcpy (py_grad_data, grad_data, grad.Dim()*sizeof(double));

    PyObject *res = Py_BuildValue ("O", py_grad);
    Py_DECREF (py_grad);
    return res;
}

// ===========================================================================

static PyObject *toast_regul_hdiag (PyObject *self, PyObject *args)
{
    int hreg;
    Regularisation *reg;
    PyObject *py_x, *py_hdiag;

    if (!PyArg_ParseTuple (args, "iO", &hreg, &py_x))
	return NULL;
    reg = (Regularisation*)g_regmgr.Get (hreg);
    // Note: we allow reg=0 as a "null" regularisation

    RVector x = CopyVector (py_x);
    RVector diag(x.Dim());
    if (reg) diag = reg->GetHessianDiag (x);
    CopyVector (&py_hdiag, diag);

    PyObject *res = Py_BuildValue ("O", py_hdiag);
    Py_DECREF (py_hdiag);
    return res;
}

// ===========================================================================

static PyObject *toast_test (PyObject *self, PyObject *args)
{
    npy_intp dmx_dims[2] = {100,10000};
    PyObject *dmx = PyArray_SimpleNew (2, dmx_dims, NPY_CDOUBLE);

    PyObject *ret = Py_BuildValue ("O", dmx);
    Py_DECREF(dmx);
    return ret;
}

// ===========================================================================

static PyMethodDef ToastMethods[] = {
    {"ReadMesh", mesh_read, METH_VARARGS, "Read a Toast mesh from file"},
    {"WriteMesh", mesh_write, METH_VARARGS, "Write a Toast mesh to file"},
    {"ClearMesh", mesh_clear, METH_VARARGS, "Delete a mesh from memory"},
    {"MeshData", mesh_data, METH_VARARGS, "Extract node and element data from a mesh"},
    {"SurfData", toast_surf_data, METH_VARARGS, "Extract surface node and face data from a mesh"},
    {"ElementData", toast_element_data, METH_VARARGS, "Extract node data for a mesh element"},
    {"MeshNodeCount", toast_mesh_node_count, METH_VARARGS, "Return the number of mesh nodes"},
    {"MeshElementCount", toast_mesh_element_count, METH_VARARGS, "Return the number of mesh elements"},
    {"MeshDim", toast_mesh_dim, METH_VARARGS, "Return the mesh dimension (2 or 3)"},
    {"MeshBB", toast_mesh_bb, METH_VARARGS, "Return the corner coordinates of the mesh bounding box"},
    {"RasterBasisPoints", toast_basis_pos, METH_VARARGS, "Return the positions of the basis points in a matrix"},
    {"RasterSolutionPoints", toast_sol_pos, METH_VARARGS, "Return the positions of the solution points in a matrix"},
    {"MakeMesh", toast_make_mesh, METH_VARARGS, "Create a Toast mesh from node and element data"},
    {"MakeRaster", toast_make_raster, METH_VARARGS, "Create a mesh-to-raster mapper"},
    {"ClearRaster", toast_clear_raster, METH_VARARGS, "Delete a raster object from memory"},
    {"MapBasis", toast_map_basis, METH_VARARGS, "Map a field from one basis to another"},
    {"ReadQM", toast_readqm, METH_VARARGS, "Load a QM description file into a mesh"},
    {"ReadNim", toast_read_nim, METH_VARARGS, "Read a nodal image file"},
    {"WriteNim", toast_write_nim, METH_VARARGS, "Write a nodal image to file"},
    {"Sysmat", toast_sysmat, METH_VARARGS, "Create a complex sparse system matrix from optical parameters"},
    {"Sysmat_CW", toast_sysmat_cw, METH_VARARGS, "Create a real sparse system matrix from optical parameters"},
    {"Qvec", (PyCFunction)toast_qvec, METH_VARARGS | METH_KEYWORDS, "Construct an array of source vectors from QM information"},
    {"Mvec", (PyCFunction)toast_mvec, METH_VARARGS | METH_KEYWORDS, "Construct an array of measurement vectors from QM information"},
    {"Fields", toast_fields, METH_VARARGS, "Calculate direct and adjoint fields"},
    {"Jacobian", toast_jacobian, METH_VARARGS, "Calculate the Jacobian matrix of the DOT forward operator"},
    {"Gradient", toast_gradient, METH_VARARGS, "Calculate the gradient of the DOT forward operator from the optical parameters"},
    {"Krylov", toast_krylov, METH_VARARGS, "Solve Hx = y with implicit Hessian H = J^T J"},
    {"Regul", (PyCFunction)toast_regul, METH_VARARGS | METH_KEYWORDS, "Return a regularisation object"},
    {"RegValue", toast_regul_value, METH_VARARGS, "Returns the regularisation value R(x) for a parameter vector x"},
    {"RegGradient", toast_regul_grad, METH_VARARGS, "Returns the regularisation gradient for a parameter vector x"},
    {"RegHDiag", toast_regul_hdiag, METH_VARARGS, "Returns the diagonal of the Hessian of the regularisation operator"},

    {"Test", toast_test, METH_VARARGS, "A dummy test function"},
    {NULL, NULL, 0, NULL}
};

// ===========================================================================

PyMODINIT_FUNC
inittoastmod(void)
{
    (void)Py_InitModule ("toastmod", ToastMethods);
    import_array();
}
