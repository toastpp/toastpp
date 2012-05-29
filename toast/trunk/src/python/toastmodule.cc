// -*- C++ -*-
#include <Python.h>
#include <numpy/arrayobject.h>
#include "felib.h"
#include "stoastlib.h"
#include <fstream>

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

class MeshManager {
public:
    MeshManager ();
    int Add (Mesh *mesh);
    Mesh *Get (int idx) const;
    bool Delete (int idx);
    void Clear ();

private:
    Mesh **list;
    int nlist;
    int nmesh;
} g_meshmgr;

MeshManager::MeshManager ()
{
    nlist = nmesh = 0;
}

int MeshManager::Add (Mesh *mesh)
{
    int i;

    // find first empty slot
    for (i = 0; i < nlist; i++) {
        if (!list[i]) break;
    }

    // no available slot - extend list
    if (i == nlist) {
        Mesh **tmp = new Mesh*[nlist+1];
	if (nlist) {
	    memcpy (tmp, list, nlist*sizeof(Mesh*));
	    delete []list;
	}
	list = tmp;
	nlist++;
    }

    list[i] = mesh;
    nmesh++;
    return i;
}

Mesh *MeshManager::Get (int idx) const
{
    Mesh *mesh = (idx >= 0 && idx < nlist ? list[idx] : 0);
    if (!mesh) std::cerr << "Not a valid mesh handle" << std::endl;
    return mesh;
}

bool MeshManager::Delete (int idx)
{
    Mesh *mesh = (idx >= 0 && idx < nlist ? list[idx] : 0);
    if (mesh) {
        delete mesh;
	list[idx] = 0;
	nmesh--;
	return true;
    } else {
        std::cerr << "Not a valid mesh handle" << std::endl;
	return false;
    }
    
}

void MeshManager::Clear ()
{
    for (int i = 0; i < nlist; i++) {
        if (list[i]) {
  	    delete list[i];
	    nmesh--;
	}
    }
    if (nlist) {
        delete []list;
	nlist = 0;
    }
}

// ===========================================================================

class RasterManager {
public:
    RasterManager ();
    int Add (Raster *raster);
    Raster *Get (int idx) const;
    bool Delete (int idx);
    void Clear ();

private:
    Raster **list;
    int nlist;
    int nraster;
} g_rastermgr;

RasterManager::RasterManager ()
{
    nlist = nraster = 0;
}

int RasterManager::Add (Raster *raster)
{
    int i;

    // find first empty slot
    for (i = 0; i < nlist; i++) {
        if (!list[i]) break;
    }

    // no available slot - extend list
    if (i == nlist) {
        Raster **tmp = new Raster*[nlist+1];
	if (nlist) {
	    memcpy (tmp, list, nlist*sizeof(Raster*));
	    delete []list;
	}
	list = tmp;
	nlist++;
    }

    list[i] = raster;
    nraster++;
    return i;
}

Raster *RasterManager::Get (int idx) const
{
    Raster *raster = (idx >= 0 && idx < nlist ? list[idx] : 0);
    if (!raster) std::cerr << "Not a valid raster handle" << std::endl;
    return raster;
}

bool RasterManager::Delete (int idx)
{
    Raster *raster = (idx >= 0 && idx < nlist ? list[idx] : 0);
    if (raster) {
        delete raster;
	list[idx] = 0;
	nraster--;
	return true;
    } else {
        std::cerr << "Not a valid raster handle" << std::endl;
	return false;
    }
    
}

void RasterManager::Clear ()
{
    for (int i = 0; i < nlist; i++) {
        if (list[i]) {
  	    delete list[i];
	    nraster--;
	}
    }
    if (nlist) {
        delete []list;
	nlist = 0;
    }
}

// ===========================================================================

static PyObject *toast_read_mesh (PyObject *self, PyObject *args)
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

static PyObject *toast_write_mesh (PyObject *self, PyObject *args)
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

static PyObject *toast_clear_mesh (PyObject *self, PyObject *args)
{
    int hmesh;

    if (!PyArg_ParseTuple (args, "i", &hmesh))
        return NULL;

    g_meshmgr.Delete (hmesh);
    Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_mesh_data (PyObject *self, PyObject *args)
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

    double *v, *vtx_data = new double[nlen*dim];
    for (i = 0, v = vtx_data; i < nlen; i++)
        for (j = 0; j < dim; j++)
	    *v++ =  mesh->nlist[i][j];

    PyObject *nodelist = PyArray_SimpleNewFromData (2, node_dims,
        PyArray_FLOAT64, vtx_data);

    // max number of nodes per element
    int nnd = mesh->elist[0]->nNode();
    for (i = 1; i < elen; i++)
	nnd = std::max (nnd, mesh->elist[i]->nNode());
    npy_intp el_dims[2] = {elen, nnd};

    // element index list
    // (0-based; value -1 indicates unused matrix entry)

    int *e, *el_data = new int[elen*nnd];
    for (i = 0, e = el_data; i < elen; i++) {
        for (j = 0; j < mesh->elist[i]->nNode(); j++)
	    *e++ = mesh->elist[i]->Node[j];
	for (; j < nnd; j++)
	    *e++ = -1;
    }
    PyObject *idx = PyArray_SimpleNewFromData (2, el_dims,
        PyArray_INT32, el_data);

    // element type list (see element.h)
    int *et, *etp_data = new int[elen];
    for (i = 0, et = etp_data; i < elen; i++) {
        *et++ = mesh->elist[i]->Type();
    }
    PyObject *eltp = PyArray_SimpleNewFromData (1, el_dims,
        PyArray_INT32, etp_data);

    // note: not allowed to delete the vertex array here - shallow copy
    // so where are we supposed to delete vtx_data, el_data and etp_data ???
    // for now, it's a memory leak

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
    memcpy (data, pos.valptr(), pos.nRows()*pos.nCols()*sizeof(double));

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
    memcpy (data, pos.valptr(), pos.nRows()*pos.nCols()*sizeof(double));

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
    }
    
    npy_intp *dims = PyArray_DIMS(py_srcvec);
    double *src_data = (double*)PyArray_DATA(py_srcvec);
    int nsrc = dims[0];

    RVector src (nsrc, src_data);
    RVector tgt;

    srcid='M';
    tgtid='B';

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
	    return NULL;
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
	    return NULL;
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
	    return NULL;
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
	    return NULL;
	}
	break;
    default:
        std::cerr << "toast.MapBasis: source id not recognised"
		  << std::endl;
	return NULL;
    }

    npy_intp ntgt = tgt.Dim();
    double *tgt_data = new double[ntgt];
    memcpy (tgt_data, tgt.data_buffer(), ntgt*sizeof(double));
    PyObject *py_tgtvec = PyArray_SimpleNewFromData (1, &ntgt, PyArray_FLOAT64,
        tgt_data);

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
    RFwdSolver FWS (LSOLVER_ITERATIVE, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    
    FWS.Allocate (*mesh);
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
    CFwdSolver FWS (LSOLVER_ITERATIVE, 1e-10);
    FWS.SetDataScaling (DATA_LOG);
    
    FWS.Allocate (*mesh);
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
    int i, j, k, idx, ofs, resettp, cmd;

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

#ifdef UNDEF
void CalcFields (QMMesh *mesh, Raster *raster,
    const CCompRowMatrix &qvec, const CCompRowMatrix &mvec,
    const RVector &mua, const RVector &mus, const RVector &ref,
    double freq, char *solver, double tol, mxArray **dfield, mxArray **afield)
{
    const double c0 = 0.3;
    int i, j, idx, n, dim, nQ, nM, nQM, slen;

    n    = mesh->nlen();
    dim  = mesh->Dimension();
    nQ   = mesh->nQ;
    nM   = mesh->nM;
    slen = (raster ? raster->SLen() : n);

    CVector *dphi, *aphi;
    CFwdSolver FWS (solver, tol);

//  if (FWS.LinSolver() == LSOLVER_DIRECT)
//	mexPrintf ("Using direct solver\n");
//  else
//	mexPrintf ("Using iterative solver, tol=%f\n", FWS.ItSolverTol());

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
    FWS.Allocate (*mesh);
    FWS.Reset (msol, omega);
    CVector sphi(slen);

    // build the field vectors
    dphi = new CVector[nQ];
    for (i = 0; i < nQ; i++) dphi[i].New (n);
    FWS.CalcFields (qvec, dphi);
    mxArray *dmx = mxCreateDoubleMatrix (slen, nQ, mxCOMPLEX);
    double *pr = mxGetPr (dmx);
    double *pi = mxGetPi (dmx);

    for (i = idx = 0; i < nQ; i++) {
	if (raster) raster->Map_MeshToSol (dphi[i], sphi);
	else        sphi = dphi[i];
	for (j = 0; j < slen; j++) {
	    pr[idx] = sphi[j].re;
	    pi[idx] = sphi[j].im;
	    idx++;
	}
    }
    delete []dphi;
    *dfield = dmx;

    // build adjoint field vectors
    if (afield) {
	aphi = new CVector[nM];
	for (i = 0; i < nM; i++) aphi[i].New (n);
	FWS.CalcFields (mvec, aphi);
	mxArray *amx = mxCreateDoubleMatrix (slen, nM, mxCOMPLEX);
	pr = mxGetPr (amx);
	pi = mxGetPi (amx);

	for (i = idx = 0; i < nM; i++) {
	    if (raster) raster->Map_MeshToSol (aphi[i], sphi);
	    else        sphi = aphi[i];
	    for (j = 0; j < slen; j++) {
		pr[idx] = sphi[j].re;
		pi[idx] = sphi[j].im;
		idx++;
	    }
	}
	delete []aphi;
	*afield = amx;
    }
}                                                    
#endif

static PyObject *toast_fields (PyObject *self, PyObject *args)
{
#ifdef UNDEF
    int hmesh, hraster;
    QMMesh *mesh;
    Raster *raster;
    PyObject *py_qvec_vl, *py_qvec_rp, *py_qvec_ci;

    if (!PyArg_ParseTuple (args, "iiOOO", &hmesh, &hraster, &py_qvec_vl,
        &py_qvec_rp, &py_qvec_ci))
        return NULL;

    if (!(mesh = (QMMesh*)g_meshmgr.Get (hmesh)))
        return NULL;
    if (!(raster = (Raster*)g_rastermgr.Get (hraster)))
        return NULL;

    int *rowptr = (int*)PyArray_DATA(py_qvec_rp);
    int *colidx = (int*)PyArray_DATA(py_qvec_ci);
    double *val = (double*)PyArray_DATA(py_qvec_vl);

    //RCompRowMatrix qvec(qrows, qcols, rowptr, colidx, val);
#endif
    Py_RETURN_NONE;
}

// ===========================================================================

static PyObject *toast_test (PyObject *self, PyObject *args)
{
    PyObject *py_test;

    if (!PyArg_ParseTuple (args, "O", &py_test))
        return NULL;

    if (PyTuple_Check(py_test))
      std::cout << "check passed" << std::endl;
    else
      std::cout << "check failed" << std::endl;

    Py_RETURN_NONE;
}

// ===========================================================================

static PyMethodDef ToastMethods[] = {
    {"ReadMesh", toast_read_mesh, METH_VARARGS, "Read a Toast mesh from file"},
    {"WriteMesh", toast_write_mesh, METH_VARARGS, "Write a Toast mesh to file"},
    {"ClearMesh", toast_clear_mesh, METH_VARARGS, "Delete a mesh from memory"},
    {"MeshData", toast_mesh_data, METH_VARARGS, "Extract node and element data from a mesh"},
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
    {"Test", toast_test, METH_VARARGS, "A dummy test function"},
    {NULL, NULL, 0, NULL}
};

// ===========================================================================

PyMODINIT_FUNC
initctoast(void)
{
    (void)Py_InitModule ("ctoast", ToastMethods);
    import_array();
}
