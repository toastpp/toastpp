#ifndef __PYTOAST_MESHMGR_H
#define __PYTOAST_MESHMGR_H

/**
 * \brief Manages a list of persistent active meshes
 */
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
};

#endif // !__PYTOAST_MESHMGR_H
