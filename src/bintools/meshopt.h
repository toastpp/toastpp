// -*-C++-*-
// ============================================================================
// meshopt
// Mesh optimisation utility routines
// ============================================================================

#ifndef __MESHOPT_H
#define __MESHOPT_H

void Reorder (Mesh &mesh, int *perm);
int SortBndToEnd (Mesh &mesh, int *perm);
int Optimise_MinBandwidth (Mesh &mesh, int *perm, int ofs, int len);
int Optimise_MMD (Mesh &mesh, int *perm, int ofs, int len);
int Optimise_Tinney2 (Mesh &mesh, int *perm, int ofs, int len);

#endif // !__MESHOPT_H
