#ifndef __MESHMGR_H
#define __MESHMGR_H

/**
 * \brief Manages a list of persistent active meshes.
 *
 * Meshes are accessed by a handle, consisting of the 0-based index of their
 * position in the mesh.
 * This class is used by Matlab and python script interfaces to keep track of
 * meshes created on the C++ side.
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

#endif // !__MESHMGR_H
