#include "objmgr.h"

// ===========================================================================

template<typename ObjTp>
ObjectManager<ObjTp>::ObjectManager ()
{
    nlist = nobj = 0;
}

// ===========================================================================

template<typename ObjTp>
ObjectManager<ObjTp>::~ObjectManager ()
{
    Clear ();
}

// ===========================================================================

template<typename ObjTp>
int ObjectManager<ObjTp>::Add (ObjTp *obj)
{
    int i;

    // find first empty slot
    for (i = 0; i < nlist; i++) {
        if (!list[i]) break;
    }

    // no available slot - extend list
    if (i == nlist) {
        ObjTp **tmp = new ObjTp*[nlist+1];
	if (nlist) {
	    memcpy (tmp, list, nlist*sizeof(ObjTp*));
	    delete []list;
	}
	list = tmp;
	nlist++;
    }

    list[i] = obj;
    nobj++;
    return i+1;
    
}

// ===========================================================================

template<typename ObjTp>
ObjTp *ObjectManager<ObjTp>::Get (int idx) const
{
    ObjTp *obj = (idx > 0 && idx <= nlist ? list[idx-1] : 0);
    if (!obj) std::cerr << "Not a valid object handle" << std::endl;
    return obj;
}

// ===========================================================================

template<typename ObjTp>
bool ObjectManager<ObjTp>::Delete (int idx)
{
    ObjTp *obj = (idx > 0 && idx <= nlist ? list[idx-1] : 0);
    if (obj) {
        delete obj;
	list[idx-1] = 0;
	nobj--;
	return true;
    } else {
        std::cerr << "Not a valid object handle: " << idx << std::endl;
	return false;
    }
}

// ===========================================================================

template<typename ObjTp>
void ObjectManager<ObjTp>::Clear ()
{
    for (int i = 0; i < nlist; i++) {
        if (list[i]) {
  	    delete list[i];
	    nobj--;
	}
    }
    if (nlist) {
        delete []list;
	nlist = 0;
    }
}
