#include "felib.h"
#include "meshmgr.h"

// ===========================================================================

MeshManager::MeshManager ()
{
    nlist = nmesh = 0;
}

// ===========================================================================

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

// ===========================================================================

Mesh *MeshManager::Get (int idx) const
{
    Mesh *mesh = (idx >= 0 && idx < nlist ? list[idx] : 0);
    if (!mesh) std::cerr << "Not a valid mesh handle" << std::endl;
    return mesh;
}

// ===========================================================================

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

// ===========================================================================

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
