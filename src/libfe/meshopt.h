// -*-C++-*-
// ============================================================================
// meshopt
// Mesh optimisation utility routines
// ============================================================================

#ifndef __MESHOPT_H
#define __MESHOPT_H

FELIB int SortBndToEnd (Mesh &mesh, int *perm);
FELIB int Optimise_MinBandwidth (Mesh &mesh, int *perm, int ofs, int len);
FELIB int Optimise_MMD (Mesh &mesh, int *perm, int ofs, int len);
FELIB int Optimise_Tinney2 (Mesh &mesh, int *perm, int ofs, int len);

#endif // !__MESHOPT_H
