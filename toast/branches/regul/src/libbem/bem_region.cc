#define BEMLIB_IMPLEMENTATION

#include "bemlib.h"
#include "bem_region.h"
#include "bem_surface.h"

using namespace std;

BEM_Region::BEM_Region (BEM_Surface *outer, BEM_Region *parent)
{
	nChildren = 0;
	outerSurf = outer;
	Parent = parent;
	kernel = new BEM_Kernel_Helmholtz;
}

BEM_Region::~BEM_Region ()
{
	if (nChildren) delete []Child;
	delete kernel;
}

std::complex<double> BEM_Region::WaveNumberMinus () const
{
	double kappa = 1.0/(3.0*(mua+mus)); // diffusion coefficient [mm]
	double c = 0.3e12/ref;              // speed of light [mm/s]
	double omega = 2.0*pi*freq;         // angular velocity [cycles/s]
	return sqrt(std::complex<double> (mua/kappa, -omega/(c*kappa)));
}

std::complex<double> BEM_Region::WaveNumber () const // the wavenumber with plus
{
	double kappa = 1.0/(3.0*(mua+mus)); // diffusion coefficient [mm]
	double c = 0.3e12/ref;              // speed of light [mm/s]
	double omega = 2.0*pi*freq;         // angular velocity [cycles/s]
	return sqrt(std::complex<double> (mua/kappa, +omega/(c*kappa)));
}


void BEM_Region::ResetOuterSurface (BEM_Surface *outer)
{
	outerSurf = outer;
	if (Parent)
		Parent->NotifyChildSurface (this);
}

int BEM_Region::AddChild (BEM_Region *child)
{
	int i;

	BEM_Region **tmp = new BEM_Region*[nChildren+1];
	if (nChildren) { // copy existing list
		for (i = 0; i < nChildren; i++)
			tmp[i] = Child[i];
		delete []Child;
	}
	Child = tmp;
	Child[nChildren] = child;
	return nChildren++;
}

void BEM_Region::DeleteChild (int idx)
{
	dASSERT(idx >= 0 && idx < nChildren, "Index out of range");

	BEM_Region **tmp;
	int i, j;

	if (nChildren > 1) {
		tmp = new BEM_Region*[nChildren-1];
		for (i = j = 0; i < nChildren; i++)
			if (i != idx) tmp[j++] = Child[i];
		delete []Child;
	} else
		tmp = NULL;
	Child = tmp;
	nChildren--;
}

void BEM_Region::NotifyChildSurface (BEM_Region *child)
{
	// Add initialisation code here if required
}


BEM_Surface *BEM_Region::GetInnerSurface (int idx)
{
	dASSERT(idx >= 0 && idx < nChildren, "Index out of range");
	return Child[idx]->outerSurf;
}

void BEM_Region::ConstructRegionMatrix (CDenseMatrix &A, CDenseMatrix &B)
{
	int i, j, k, nd, el;
	int colofs, rowofs = 0;

	
	// Collect outer and inner surfaces in a single list
	// for eaase of addressing
	BEM_Surface **surf = new BEM_Surface*[nChildren+1];
	surf[0] = outerSurf;
	for (i = 0; i < nChildren; i++)
		surf[i+1] = GetInnerSurface(i);
	int nsurf = nChildren + 1;

	int totnode = 0;
	for (i = 0; i < nsurf; i++)
		totnode += surf[i]->nNodes();
	A.New (totnode, totnode);
	B.New (totnode, totnode);

	// reset the wavenumber on the kernel
	kernel->SetWavenumber (WaveNumber());

	bool invert;
	bool integrate_local;
	CVector integral;
	for (i = 0; i < nsurf; i++) {
		Point3D *nlist = surf[i]->NodeList();
		for (nd = 0; nd < surf[i]->nNodes(); nd++) {
			//cerr << nd << endl;
			colofs = 0;
			for (j = 0; j < nsurf; j++) {
				invert = (j > 0);
				integrate_local = (i == j);
				for (el = 0; el < surf[j]->nElements(); el++) {
					//cerr<<el<<endl;
					if (integrate_local){
						integral = surf[j]->Integrate (kernel, nd, el, invert);
						
					}
					else{
						integral = surf[j]->Integrate_Nonsingular (kernel, nlist[nd], el, invert);
					}
					int nnd = surf[j]->Element(el)->nNode();
					for (k = 0; k < nnd; k++) {
						int ndidx = surf[j]->Element(el)->NodeIndex(k);
						A(rowofs+nd,colofs+ndidx) += integral[k];
						B(rowofs+nd,colofs+ndidx) += integral[k+nnd];
						
					}

				} // end loop el
				//colofs += surf[j]->nNodes();
			} // end loop element nsurf
		} // end loop nd
		rowofs += surf[i]->nNodes();
	} // end loop node nsurf

	delete []surf;
}
