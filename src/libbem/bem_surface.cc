#define BEMLIB_IMPLEMENTATION

#include "bemlib.h"
#include "bem_surface.h"
#include "bem_tri6.h"

using namespace std;

BEM_Surface::BEM_Surface (RDenseMatrix &N, IDenseMatrix &E)
{
	nNd = N.nRows();
	int i, j, dim = N.nCols();

	// Create list of nodes from node matrix N
	Nd = new Point3D[nNd];
	for (i = 0; i < nNd; i++) {
		for (j = 0; j < dim; j++)
			Nd[i][j] = N(i,j);
	}

	// Find element type (currently this only allows
	// 3-noded and 6-noded triangles, depending on number
	// of elements in the list)
	enum {TRI3, TRI6} ElType;
	switch (E.nCols()) {
		case 6: ElType = TRI6; break;
		default: xERROR("Element type not supported"); exit(0);
	}

	// Create list of elements from element index matrix E
	nEl = E.nRows();
	El = new BEM_Element*[nEl];
	for (i = 0; i < nEl; i++) {
		switch (ElType) {
			case TRI6:
				El[i] = new BEM_Triangle6 (this, E.Row(i)); 
				break;
			default:
				xERROR ("This should not have happened");
				break;
		}
	}
}

BEM_Surface::BEM_Surface (std::istream &is)
{
	Mesh mesh;
	is >> mesh;

	nNd = mesh.nlen();
	int i, j, dim = mesh.Dimension();

	Nd = new Point3D[nNd];
	for (i = 0; i < nNd; i++) {
		Nd[i].New(dim);
		for (j = 0; j < dim; j++)
			Nd[i][j] = mesh.nlist[i][j];
	}

	nEl = mesh.elen();
	El = new BEM_Element*[nEl];
	for (i = 0; i < nEl; i++) {
		switch (mesh.elist[i]->Type()) {
			case ELID_TRI3D6:
				El[i] = new BEM_Triangle6 (this, IVector(6, mesh.elist[i]->Node));
				break;
			default:
				xERROR("Element type not supported");
				break;
		}
	}
}

BEMLIB std::ostream& operator<< (std::ostream& o, BEM_Surface &surf)
{
	int j,k;
	o << "NodeList " << surf.nNd << endl;
	for (j = 0; j < surf.nNd; j++) {
		for (k = 0; k < 3; k++) {
			o << surf.Nd[j][k] << ' ';
		}
		o << endl;
	}
	o << "\nElementList " << surf.nEl << endl;
	for (j = 0; j < surf.nEl; j++) {
		for (k = 0; k < surf.El[j]->nNode(); k++)
			o << surf.El[j]->NodeIndex(k) << ' ';
		o << endl;
	}
	return o;
}

CVector BEM_Surface::Integrate_Nonsingular (BEM_Kernel *kernel, const Point3D &point, int el, bool invert)
{
	dASSERT(el >= 0 && el < nEl, "Element index out of range");
	return El[el]->Integrate_Nonsingular (kernel, point, invert);
}

CVector BEM_Surface::Integrate (BEM_Kernel *kernel, int nd, int el, bool invert)
{
	dASSERT(el >= 0 && el < nEl, "Element index out of range");
	dASSERT(nd >= 0 && nd < nNd, "Node index out of range");

	bool isSingular = false;
	BEM_Element *pel = El[el];
	int i;
	
	for (i = 0; i < pel->nNode(); i++)
		if (pel->NodeIndex(i) == nd) {
			isSingular = true;
			break;
		}

	if (isSingular)
		return pel->Integrate_Singular (kernel, nd, Nd[nd], invert);
	else
		return pel->Integrate_Nonsingular (kernel, Nd[nd], invert);
}