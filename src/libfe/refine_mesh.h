#ifndef __REFINE_MESH_H
#define __REFINE_MESH_H

void RefineTriangle3Mesh (Mesh *mesh, bool *elrefine);
void RefineTetrahedron4Mesh (Mesh *mesh, bool *elrefine);

#endif // !__REFINE_MESH_H
