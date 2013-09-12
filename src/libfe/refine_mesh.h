#ifndef __REFINE_MESH_H
#define __REFINE_MESH_H

FELIB void RefineTriangle3Mesh (Mesh *mesh, bool *elrefine);
FELIB void RefineTetrahedron4Mesh (Mesh *mesh, bool *elrefine);

FELIB void Mesh_SplitTriangle3 (Mesh *mesh, int el);

FELIB void JiggleMesh (Mesh *mesh, double scale = 0.1, int iterations = 1);

#endif // !__REFINE_MESH_H
