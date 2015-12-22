#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

// =========================================================================
// class Raster_CPixel_Tree
// =========================================================================

Raster_CPixel_Tree::Raster_CPixel_Tree (const IVector &_bdim,
    const IVector &_gdim, Mesh *mesh, RDenseMatrix *bb, double _map_tol)
  : Raster2 (_bdim, _gdim, mesh, bb, _map_tol)
{
    // currently the following conditions apply:
    // Dimensions _bdim must be the same in all directions, and must be
    // powers of 2 (_gdim is ignored)

    std::cerr << "raster_cpixel_tree constructor" << std::endl;

    int res = _bdim[0];
    int treedepth = 0;
    while (res > 1) {
	res >>= 1;
	treedepth++;
    }
    tree = new BasisTree (this, treedepth);
    blen = 0;
    slen = 0;
    
    std::cerr << "before finalisetree" << std::endl;
    FinaliseTree();
}

// =========================================================================

void Raster_CPixel_Tree::FinaliseTree ()
{
    std::cerr << "FinaliseTree" << std::endl;

    if (blen) {
	delete []node;
    }
    blen = tree->NumEndNodes();
    slen = blen;
    glen = blen;
    node = new TreeNode*[blen];
    tree->EnumerateEndNodes (node);
    for (int i = 0; i < blen; i++)
	node[i]->linidx = i;

    //Init();
    if (Bvv) delete Bvv;
    Bvv = CreateBvv();
    if (Buv) delete Buv;
    Buv = CreateBuv();
    if (map_tol) {
	std::cerr << "reset precon" << std::endl;
	if (Bvv_precon) delete Bvv_precon;
	Bvv_precon = new RPrecon_IC; Bvv_precon->Reset (Bvv);
    } else {
	std::cerr << "reset cholesky" << std::endl;
	if (Bvv_Cholesky_L) delete Bvv_Cholesky_L;
	if (Bvv_Cholesky_d) delete Bvv_Cholesky_d;
	int *rowptr, *colidx;
	Bvv->SymbolicCholeskyFactorize (rowptr, colidx);
	Bvv_Cholesky_L = new RCompRowMatrix(blen, blen, rowptr, colidx);
	delete []rowptr;
	delete []colidx;
	Bvv_Cholesky_d = new RVector(blen);
    }
    // may need to think about reconstructing D
    
}

// =========================================================================

double Raster_CPixel_Tree::Value_nomask (const Point &p, int i, bool is_solidx)
    const
{
    int bi = (is_solidx ? sol2basis[i] : i);
    int iz = (dim < 3 ? 0 : bi / (bdim[0]*bdim[1]));
    bi -= iz*bdim[0]*bdim[1];
    int iy = bi/bdim[0];
    int ix = bi - iy*bdim[0];

    int px = (int)((p[0]-bbmin[0])/bbsize[0]*bdim[0]);
    int py = (int)((p[1]-bbmin[1])/bbsize[1]*bdim[1]);
    if (px != ix || py != iy) return 0.0;
    if (dim == 3) {
	int pz = (int)((p[2]-bbmin[2])/bbsize[2]*bdim[2]);
	if (pz != iz) return 0.0;
    }
    return 1.0;
}

// =========================================================================

RDenseMatrix Raster_CPixel_Tree::SupportArea (int idx)
{
    dASSERT(idx >= 0 && idx < blen, "Index out of range");

    RDenseMatrix sa(2,2);
    TreeNode *n = node[idx];
    sa(0,0) = n->xmin;
    sa(1,0) = n->xmax;
    sa(0,1) = n->ymin;
    sa(1,1) = n->ymax;
    return sa;
}

// =========================================================================

void Raster_CPixel_Tree::Refine (int idx)
{
    dASSERT(idx >= 0 && idx < blen, "Index out of range");

    node[idx]->Refine();
    FinaliseTree();    
}

// =========================================================================

void Raster_CPixel_Tree::Refine (int *idx, int nidx)
{
    for (int i = 0; i < nidx; i++) {
	dASSERT(idx[i] >= 0 && idx[i] < blen, "Index out of range");
	node[idx[i]]->Refine();
    }
    FinaliseTree();    
}

// =========================================================================

RCompRowMatrix *Raster_CPixel_Tree::CreateBvv () const
{
    // construct a sparse diagonal matrix
    idxtype *rowptr = new idxtype[blen+1];
    idxtype *colidx = new idxtype[blen];
    double *val = new double[blen];
    rowptr[0] = 0;
    for (int i = 0; i < blen; i++) {
	rowptr[i+1] = rowptr[i]+1;
	colidx[i] = i;
	val[i] = node[i]->Size();
    }
    RCompRowMatrix *bvv = new RCompRowMatrix(blen,blen,rowptr,colidx,val);
    delete []rowptr;
    delete []colidx;
    delete []val;

    return bvv;
}

// =========================================================================

RCompRowMatrix *Raster_CPixel_Tree::CreateBuv () const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    case ELID_TRI6:
    case ELID_TRI10:
	return CreateBuv_tri();
	return 0;
    case ELID_TET4:
	//return CreateBuv_tet4();
	return 0;
    default:
	xERROR("Raster_CPixel: Unsupported element type");
	return 0;
    }
}

// ==========================================================================
// Adds contribution from single element to system matrix

void Raster_CPixel_Tree::AddToElMatrix (int el,
    RGenericSparseMatrix &M, const RVector *pxcoeff, int mode) const
{
    switch (meshptr->elist[0]->Type()) {
    case ELID_TRI3:
    //case ELID_TRI6:
    //case ELID_TRI10:
	//AddToElMatrix_tri (el, M, pxcoeff, mode);
	break;
    case ELID_TET4:
	//AddToElMatrix_tet (el, M, pxcoeff, mode);
	break;
    default:
	xERROR("Raster_Pixel2: Unsupported element type");
    }
}



// =========================================================================
// class BasisTree
// =========================================================================

TreeNode::TreeNode (const BasisTree *Tree, TreeNode *Parent, int subidx)
    : tree(Tree), parent(Parent), idx(subidx)
{
    linidx = -1;

    for (int i = 0; i < 4; i++)
	child[i] = NULL;

    if (parent) {
	lvl = parent->lvl+1;
	xmin = (idx & 1 ? (parent->xmin + parent->xmax) * 0.5 : parent->xmin);
	xmax = (idx & 1 ? parent->xmax : (parent->xmin + parent->xmax) * 0.5);
	ymin = (idx < 2 ? parent->ymin : (parent->ymin + parent->ymax) * 0.5);
	ymax = (idx < 2 ? (parent->ymin + parent->ymax) * 0.5 : parent->ymax);
    } else {
	lvl = 0;
	RDenseMatrix bb = tree->raster->BoundingBox();
	xmin = bb(0,0);
	xmax = bb(1,0);
	ymin = bb(0,1);
	ymax = bb(1,1);
    }
}

void TreeNode::EnsureLevel (int minlvl)
{
    if (lvl < minlvl) {
	Refine(); // make sure subnodes exist
	for (int i = 0; i < 4; i++) {
	    child[i]->EnsureLevel (minlvl);
	    // recursively refine to desired level
	}
    }
}

void TreeNode::DeleteBranch (int subidx)
{
    if (child[subidx]) {
	for (int i = 0; i < 4; i++)
	    child[subidx]->DeleteBranch (i);
	delete child[subidx];
	child[subidx] = NULL;
    }
}

void TreeNode::Refine ()
{
    for (int i = 0; i < 4; i++) {
	if (!child[i])
	    child[i] = new TreeNode (tree, this, i);
    }
}

const TreeNode *TreeNode::FindNode (const Point &p)
{
    if (p[0] < xmin || p[0] > xmax || p[1] < ymin || p[1] > ymax)
	return NULL;
    
    int idx = (p[0] < (xmin+xmax)*0.5 ? 0:1);
    if (p[1] >= (ymin+ymax)*0.5) idx += 2;

    return (child[idx] ? child[idx]->FindNode(p) : this);
}

double TreeNode::Size () const
{
    return (xmax-xmin)*(ymax-ymin);
}

int TreeNode::NumEndNodes () const
{
    int nnode = 0;
    for (int i = 0; i < 4; i++)
	if (child[i]) nnode += child[i]->NumEndNodes();

    if (!nnode) nnode = 1;
    return nnode;
}

void TreeNode::EnumerateEndNodes (TreeNode **node, int *idx)
{
    bool is_end = true;
    for (int i = 0; i < 4; i++)
	if (child[i]) {
	    child[i]->EnumerateEndNodes (node, idx);
	    is_end = false;
	}
    if (is_end) {
	node[*idx] = this;
	(*idx)++;
    }
}

void TreeNode::EnumerateOverlappingEndNodes (TreeNode **node, int *idx,
    int maxidx, double _xmin, double _ymin, double _xmax, double _ymax)
{
    if (_xmin >= xmax || _xmax <= xmin || _ymin >= ymax || _ymax <= ymin)
	return;

    bool is_end = true;
    for (int i = 0; i < 4; i++)
	if (child[i]) {
	    child[i]->EnumerateOverlappingEndNodes (node, idx, maxidx,
	        _xmin, _ymin, _xmax, _ymax);
	    is_end = false;
	}
    if (is_end) {
	if (*idx < maxidx)
	    node[*idx] = this;
	(*idx)++;
    }
}

// =========================================================================

BasisTree::BasisTree (const Raster_CPixel_Tree *Raster, int startres)
    : raster(Raster)
{
    root = new TreeNode (this);
    root->EnsureLevel (startres);
}

const TreeNode *BasisTree::FindNode (const Point &p)
{
    return root->FindNode (p);
}

int BasisTree::NumEndNodes () const
{
    return root->NumEndNodes();
}

void BasisTree::EnumerateEndNodes (TreeNode **node)
{
    int idx = 0;
    root->EnumerateEndNodes (node, &idx);
}

int BasisTree::EnumerateOverlappingEndNodes (TreeNode **node, int maxidx,
    double xmin, double ymin, double xmax, double ymax)
{
    int idx = 0;
    root->EnumerateOverlappingEndNodes (node, &idx, maxidx,
        xmin, ymin, xmax, ymax);
    return idx;
}
