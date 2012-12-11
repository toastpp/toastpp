#ifndef __PYTOAST_OBJMGR_H
#define __PYTOAST_OBJMGR_H

template<typename ObjType>
class ObjectManager {
public:
    ObjectManager ();
    ~ObjectManager ();
    int Add (ObjType *obj);
    ObjType *Get (int idx) const;
    bool Delete (int idx);
    void Clear ();

private:
    ObjType **list;
    int nlist;
    int nobj;
};

#include "objmgr_imp.hpp"

#endif // !__PYTOAST_OBJMGR_H
