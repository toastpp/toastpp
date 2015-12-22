// -*-C++-*-

#ifndef __RASTER_CPX_TREE_H
#define __RASTER_CPX_TREE_H

/**
 * \file Defines class Raster_CPixel_Tree (piecewise constant pixel/voxel
 *   image basis with quad/octree adaptivity)
 */

#include "raster2.h"

class STOASTLIB BasisTree;
class STOASTLIB Raster_CPixel_Tree;

// =========================================================================

class STOASTLIB TreeNode {
public:
    TreeNode (const BasisTree *Tree, TreeNode *Parent=NULL, int subidx=0);
    
    // make sure the node expands to at least level minlvl, adding nodes
    // as required
    void EnsureLevel (int minlvl);

    // delete the child branch subidx
    void DeleteBranch (int subidx);

    // Add subnodes, if not already present
    void Refine ();

    const TreeNode *FindNode (const Point &p);

    double Size() const;

    int NumEndNodes () const;

    void EnumerateEndNodes (TreeNode **node, int *idx);

    void EnumerateOverlappingEndNodes (TreeNode **node, int *idx,
	int maxidx, double _xmin, double _ymin, double _xmax, double _ymax);

    const BasisTree *tree;
    TreeNode *parent;
    TreeNode *child[4]; // do quadtree for now
    int lvl;
    int idx; // child subindex (0-3)
    int linidx; // end node basis index (-1 for non-endnodes)
    double xmin, ymin, xmax, ymax; // pixel extents
};

// =========================================================================

class STOASTLIB BasisTree {
public:
    BasisTree (const Raster_CPixel_Tree *Raster, int startres=0);

    // returns the furthest node containing point p
    const TreeNode *FindNode (const Point &p);

    // number of nodes at the branch ends
    int NumEndNodes() const;

    void EnumerateEndNodes (TreeNode **node);

    // return a list of end nodes overlapping the specified rectangle
    // idx: list to be filled
    // maxidx: length of idx
    // xmin, ymin, xmax, ymax: area to check for overlap
    // return value is number of overlapping nodes found
    // if return value > maxidx, the list is truncated. Caller should
    // re-allocate to required size and call again
    int EnumerateOverlappingEndNodes (TreeNode **node, int maxidx,
        double xmin, double ymin, double xmax, double ymax);

    const Raster_CPixel_Tree *raster;
    TreeNode *root;
};

// =========================================================================
/**
 * \brief Basis for 2D and 3D image representation with a regular piecewise
 *   constant pixel/voxel grid with quad/octree adaptivity.
 */

class STOASTLIB Raster_CPixel_Tree: public Raster2 {
    friend class Raster2;
    // to make constructor accessible to Raster2 factory

public:
    /**
     * \brief Value of basis function b_i at point p
     * This does not check for mesh support
     * \sa Value
     */
    double Value_nomask (const Point &p, int i, bool is_solidx=true) const;

    /**
     * \brief Support area for an individual basis function
     * \param idx basis index (0...blen-1)
     * \return Matrix containing pixel bounding box
     */
    RDenseMatrix SupportArea (int idx);

    void Refine (int idx);
    void Refine (int *idx, int nidx);

    /**
     * \brief Single-element system matrix assembly
     *
     * Assemble single-element contribution for element "el" into global
     * system matrix M, where coefficients (where applicable) are given in
     * pixel basis.
     * \param el element index (>= 0)
     * \param M global system matrix
     * \param pxcoeff pointer to coefficient vector in pixel basis (only
     *   required for integration modes involving a parameter distribution)
     * \param mode integration type index
     */
    void AddToElMatrix (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;


protected:
    /**
     * \brief Pixel basis constructor
     * \param _bdim basis dimensions (the length of this vector is 2 or 3,
     *   corresponding to the mesh dimension).
     * \param _gdim Not used by this basis type. Must be identical to _bdim.
     * \param mesh pointer to mesh instance
     * \param bb bounding box for the basis representation (size 2x3 or 3x3)
     * \param _map_tol Tolerance limit for least-squares solver in
     *   Map_MeshToBasis and Map_BasisToMesh
     * \note If the bounding box is not provided, it is calculated from the
     *   mesh geometry.
     * \note The mesh pointer must remain valid for the lifetime of the basis
     *   instance.
     */
    Raster_CPixel_Tree (const IVector &_bdim, const IVector &_gdim, Mesh *mesh,
        RDenseMatrix *bb=0, double _map_tol=1e-10);

    RCompRowMatrix *CreateBvv () const;
    RCompRowMatrix *CreateBuv () const;
    RCompRowMatrix *CreateBuv_tri () const;
    RCompRowMatrix *CreateBuv_tet4 () const;

    void AddToElMatrix_tri (int el, RGenericSparseMatrix &M,
        const RVector *pxcoeff, int mode) const;

private:
    int SutherlandHodgman (int el, TreeNode *node, Point *clip_poly,
        int npoly) const;

    void FinaliseTree();

    BasisTree *tree;
    TreeNode **node;   ///< list of nodes for each basis index
};

#endif // !__RASTER_CPX_H
