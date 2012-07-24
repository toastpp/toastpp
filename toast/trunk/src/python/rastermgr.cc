#include "felib.h"
#include "stoastlib.h"
#include "rastermgr.h"

// ===========================================================================

RasterManager::RasterManager ()
{
    nlist = nraster = 0;
}

// ===========================================================================

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

// ===========================================================================

Raster *RasterManager::Get (int idx) const
{
    Raster *raster = (idx >= 0 && idx < nlist ? list[idx] : 0);
    if (!raster) std::cerr << "Not a valid raster handle" << std::endl;
    return raster;
}

// ===========================================================================

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

// ===========================================================================

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
