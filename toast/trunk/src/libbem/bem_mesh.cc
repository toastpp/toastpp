#define BEMLIB_IMPLEMENTATION

#include "bemlib.h"
#include "bem_mesh.h"
#include "bem_region.h"

BEM_Mesh::BEM_Mesh ()
{
	Region = NULL;
	nRegion = 0;
	Surf = NULL;
	nSurf = 0;
}

int BEM_Mesh::AddSurface (BEM_Surface *surf)
{
	int i;

	BEM_Surface **tmp = new BEM_Surface*[nSurf+1];
	if (nSurf) {
		for (i = 0; i < nSurf; i++)
			tmp[i] = Surf[i];
		delete []Surf;
	}
	Surf = tmp;
	Surf[nSurf++] = surf;
	return nSurf;
}

void BEM_Mesh::AddRegion (BEM_Region *region, double mua, double mus, double ref, BEM_Region *parent)
{
	int i;

	BEM_Region **tmp = new BEM_Region*[nRegion+1];
	if (nRegion) {
		for (i = 0; i < nRegion; i++)
			tmp[i] = Region[i];
		delete []Region;
	}
	Region = tmp;
	Region[nRegion++] = region;

	if (parent)
		parent->AddChild (region);

}
