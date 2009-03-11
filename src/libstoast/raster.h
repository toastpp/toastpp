// -*-C++-*-
// =========================================================================
// This class represents the intermediate high-resolution pixel grid
// between the mesh basis and the "user" basis
//
// The following bases are used:
//   * FEM mesh basis ("mesh")
//   * optionally: High-res intermediate pixel grid ("grid")
//   * Low-res user pixel grid ("basis")
//   * Low-res user pixel basis excluding pixels without mesh support ("sol")
//
// It provides methods to map between the forward and inverse bases
// =========================================================================

#ifndef __RASTER_H
#define __RASTER_H

#include "mathlib.h"
#include "felib.h"
#include "solution.h"

// =========================================================================
#define MGM_DIRECTMAP     1
#define MGM_TRANSPOSE     2
#define MGM_PSEUDOINVERSE 3

#define MAP_GRIDTOMESH    MGM_DIRECTMAP
// Mapping strategy for grid->mesh basis mapping
// MGM_DIRECTMAP: Use linear interpolation of grid points at mesh nodes
// MGM_TRANSPOSE: Use transpose of mesh->grid mapping matrix
// MGM_PSEUDOINVERSE: Use pseudoinverse of mesh->grid matrix, solve with CG
// =========================================================================
#ifndef SQR
#define SQR(_X) ((_X)*(_X))
#endif
static int POW2[] = {1,2,4,8,16,32,64,128,256,512,1024,2096,4192};

class STOASTLIB Raster {
public:
    Raster (IVector &_gdim, Mesh *mesh);
    // gdim: dimensions of solution pixel basis
    // mesh: FEM basis

    Raster (IVector &_gdim, IVector &_bdim, Mesh *mesh);
    // gdim: dimensions of intermediate hires grid
    // bdim: dimensions of solution pixel basis
    // mesh: FEM basis

    ~Raster();

    const int Dim()  const { return dim; }
    const int GLen() const { return glen; }
    const int BLen() const { return blen; }
    const int SLen() const { return slen; }
    const IVector &GDim() const { return gdim; }
    const IVector &BDim() const { return bdim; }
    const RVector &GSize() const { return gsize; }
    const Mesh &mesh() const { return *meshptr; }
    const int *Elref() const { return elref; }
    const int *BElref() const { return belref; }
    const int *BasisToSol() const { return basis2sol.data_buffer(); }
    const RGenericSparseMatrix &Mesh2GridMatrix()  const { return *B; }
    const RGenericSparseMatrix &Grid2MeshMatrix()  const { return *BI; }
    const RGenericSparseMatrix &Grid2BasisMatrix() const { return *B2; }
    const RGenericSparseMatrix &Basis2GridMatrix() const { return *B2I; }
    const RGenericSparseMatrix &Mesh2BasisMatrix() const;
    const RGenericSparseMatrix &Basis2MeshMatrix() const;

    void GetPixelCoords (int i, IVector &pix) const;
    // convert solution index i into an (x,y[,z]) pixel coordinate

    const int Basis2Sol (int i) const { return basis2sol[i]; }
    const int Sol2Basis (int j) const { return sol2basis[j]; }
    // convert between BDim and SDim basis indices

    const RVector &BasisSupportArea () const { return bsupport; }
    // Returns vector of length blen which basis support weights
    // Range: 0 (no mesh support) to 1 (full mesh support)

    void NeighbourGraph (int *&rowptr, int *&colidx, int &nzero) const;
    void NeighbourGraph (ICompRowMatrix &NG) const;
    IVector NeighbourShift (const ICompRowMatrix &NG, int i, int j) const;
    RVector RNeighbourShift (const ICompRowMatrix &NG, int i, int j) const;

    void MapGraph2Sol (const int *browptr, const int *bcolidx,
        const int bnzero, int *&srowptr, int *&scolidx, int &snzero) const;

    void BuildRegularRestriction (const IVector &gdim, RCompRowMatrix &R);

    // Mapping between bases: 1. Raw real vectors
    inline void Map_MeshToGrid (const RVector &mvec, RVector &gvec) const
    { B->Ax (mvec, gvec); }
    inline void Map_GridToMesh (const RVector &gvec, RVector &mvec) const
    {
#if MAP_GRIDTOMESH == MGM_PSEUDOINVERSE
	Map_GridToMeshPI (gvec, mvec);
#else
	BI->Ax (gvec, mvec);
#endif
    }
    inline void Map_GridToBasis (const RVector &gvec, RVector &bvec) const
    { if (bCoarseBasis) B2->Ax (gvec, bvec); else bvec = gvec; }
    inline void Map_BasisToGrid (const RVector &bvec, RVector &gvec) const
    { if (bCoarseBasis) B2I->Ax (bvec, gvec); else gvec = bvec; }
    inline void Map_MeshToBasis (const RVector &mvec, RVector &bvec) const
    { if (bCoarseBasis) C->Ax (mvec, bvec); else Map_MeshToGrid (mvec, bvec); }
    inline void Map_BasisToMesh (const RVector &bvec, RVector &mvec) const
    { if (bCoarseBasis) CI->Ax (bvec, mvec); else Map_GridToMesh (bvec, mvec);}
    inline void Map_BasisToSol (const RVector &bvec, RVector &svec) const
    { D->Ax (bvec, svec); }
    inline void Map_SolToBasis (const RVector &svec, RVector &bvec) const
    { D->ATx (svec, bvec); }
    void Map_GridToSol  (const RVector &gvec, RVector &svec) const;

    void Map_MeshToSol  (const RVector &mvec, RVector &svec) const;
    void Map_SolToMesh  (const RVector &svec, RVector &mvec) const;
    inline void Map_SolToGrid  (const RVector &svec, RVector &gvec) const
    { RVector bvec(blen); D->ATx (svec, bvec); Map_BasisToGrid (bvec, gvec); }

    void Map_GridToMeshPI (const RVector &gsol, RVector &msol) const;
    // map from grid to mesh using pseudoinverse of B

    // Mapping between bases: 2. Raw complex vectors
    inline void Map_MeshToGrid  (const CVector &mvec, CVector &gvec) const
    { ((RCompRowMatrix*)B)->Ax_cplx (mvec, gvec); }
    inline void Map_GridToMesh  (const CVector &gvec, CVector &mvec) const
    {   ((RCompRowMatrix*)BI)->Ax_cplx (gvec, mvec); }
    inline void Map_GridToBasis (const CVector &gvec, CVector &bvec) const
    {   if (bCoarseBasis) ((RCompRowMatrix*)B2)->Ax_cplx (gvec, bvec);
	else bvec = gvec;
    }
    inline void Map_BasisToGrid (const CVector &bvec, CVector &gvec) const
    {   if (bCoarseBasis) ((RCompRowMatrix*)B2I)->Ax_cplx (bvec, gvec);
	else gvec = bvec;
    }
    inline void Map_MeshToBasis (const CVector &mvec, CVector &bvec) const
    {   if (bCoarseBasis) ((RCompRowMatrix*)C)->Ax_cplx (mvec, bvec);
	else Map_MeshToGrid (mvec, bvec);
    }
    inline void Map_BasisToMesh (const CVector &bvec, CVector &mvec) const
    {   if (bCoarseBasis) ((RCompRowMatrix*)CI)->Ax_cplx (bvec, mvec);
	else Map_GridToMesh (bvec, mvec);
    }
    inline void Map_BasisToSol (const CVector &bvec, CVector &svec) const
    { D->Ax_cplx (bvec, svec); }
    inline void Map_SolToBasis (const CVector &svec, CVector &bvec) const
    { D->ATx_cplx (svec, bvec); }
    void Map_MeshToSol (const CVector &mvec, CVector &svec) const;
    void Map_SolToMesh (const CVector &svec, CVector &mvec) const;
    void Map_GridToSol (const CVector &gvec, CVector &svec) const;
    inline void Map_SolToGrid  (const CVector &svec, CVector &gvec) const
    { CVector bvec(blen); D->ATx_cplx (svec, bvec); Map_BasisToGrid (bvec, gvec); }


    // Mapping between bases: 3. Solutions
    void Map_MeshToGrid (const Solution &msol, Solution &gsol,
        bool mapall = false) const;
    void Map_GridToMesh (const Solution &gsol, Solution &msol,
	bool mapall = false) const;
    void Map_GridToBasis (const Solution &gsol, Solution &bsol,
	bool mapall = false) const;
    void Map_BasisToGrid (const Solution &bsol, Solution &gsol,
	bool mapall = false) const;
    void Map_MeshToBasis (const Solution &msol, Solution &bsol,
	bool mapall = false) const;
    void Map_BasisToMesh (const Solution &bsol, Solution &msol,
	bool mapall = false) const;

    void Map_BasisToSol (const Solution &bsol, Solution &ssol,
	bool mapall = false) const;
    void Map_SolToBasis (const Solution &ssol, Solution &bsol,
	bool mapall = false) const;
    void Map_MeshToSol (const Solution &msol, Solution &ssol,
	bool mapall = false) const;
    void Map_SolToMesh (const Solution &ssol, Solution &msol,
	bool mapall = false ) const;

    // Mapping from active parameter vector to solution
    void Map_ActiveSolToMesh (const RVector &asol, Solution &msol) const;

    RVector SmoothImage(const RVector &x, double sd) const;
    RVector * ImageJet(const RVector &x, double sd, bool *iflags) const;
    RVector ImageGradient(const RVector &x, double sd) const;

private:
    bool bCoarseBasis;  // lo-res solution basis supported?
    IVector gdim;       // dimensions of intermediate grid
    IVector paddim;     // dimensions padded to a power of 2
    IVector bdim;       // dimensions of user basis
    Mesh *meshptr;      // mesh pointer
    Point bbmin, bbmax; // mesh bounding box
    RVector gsize;      // extents of bounding box
    int dim;            // dimension
    int glen, blen;     // vector length in intermediate and user grid
    int padlen;         // vector length of padded grid
    int slen;           // vector length of solution basis
    int *elref;         // fine grid (gdim) ->element reference list
    int *belref;        // user basis (bdim) -> element reference list
    RGenericSparseMatrix *B, *BI;   // map matrices mesh->grid and grid->mesh
    RGenericSparseMatrix *B2, *B2I; // map matrix grid->basis and basis->grid
    RCompRowMatrix *C;              // B2*B : map mesh->basis
    RCompRowMatrix *CI;             // BI*B2I: map basis->mesh
    RCompRowMatrix *D;              // map basis->solution (limited support)
    RCompRowMatrix BB; // experimental
    RPreconditioner *pBTB;
    bool BBover;

    RVector bsupport;
    // Returns vector of length blen which basis support weights
    // Range: 0 (no mesh support) to 1 (full mesh support)
    // currently only binary states (0 or 1) are supported

    IVector basis2sol;
    // index vector of length blen containing the index of the solution
    // vector for each pixel in the user grid. For pixel (i,j,k) the solution
    // index is given by basis2sol[i+bdim[0]*j+bdim[0]*bdim[1]*k]
    // Value -1 means pixel not present in solution vector (no mesh support)

    IVector sol2basis;
    // index vector of length slen for mapping from solution to rectangular
    // basis.
};

#endif // !__RASTER_H
