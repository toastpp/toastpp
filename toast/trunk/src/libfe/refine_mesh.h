#ifndef __REFINE_MESH_H
#define __REFINE_MESH_H

void RefineTriangle3Mesh (Mesh *mesh, bool *elrefine);
void RefineTetrahedron4Mesh (Mesh *mesh, bool *elrefine);

void JiggleMesh (Mesh *mesh, double scale = 0.1, int iterations = 1);

#endif // !__REFINE_MESH_H
