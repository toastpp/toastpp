// -*-C++-*-
// ==========================================================================
// Module libfe
// File mesh.h
// Declaration of class Mesh
// ==========================================================================
#ifndef __NONCONFORMINGMESH_H
#define __NONCONFORMINGMESH_H
#include <vector>
#include<set>
#include "felib.h"

#define MIN_INTERIOR_EDGES 100000
#define MIN_BOUNDARY_EDGES 50000
#define EDGE_ACTIVE 1
#define EDGE_INACTIVE 0
#define EDGE_GARBAGE -1
// ==========================================================================
// class NonconformingMesh
// ==========================================================================
/**
 * \brief Finite-element nonconforming mesh management
 *
 * Defines node coordinates (via NodeList), element connectivity (via
 * ElementList) interior edge lists (via iedge_elist, iedge_nlist), boundary edge list (via bedge_elist, bedge_nlist) and provides general mesh-related methods.
 * This version also defines node-based parameter coefficient sets for
 * optical tomography (via ParameterList), but this should be removed.
 */

class FELIB NonconformingMesh: public Mesh{
public:

     /** contains the list of elements involved in each of the interior edges
      */
     std::vector< std::vector<int> > iedge_elist; 

     /** contains the list of nodes involved in each of the interior edges
     */
     std::vector< std::vector<int> > iedge_nlist; 

    /** contains the list of nodes involved in each of the boundary edges
     */
    std::vector< std::vector<int> > bedge_nlist; 

    /** contains the list of elements corresponding to each boundary edge
     */
    std::vector<int> bedge_elist; 

    /** List of 'state' of each interior edge. The possible 'state' are:
     * EDGE_ACTIVE => It is a physical edge and should be accounted while computing edge integrals in DG formulation
     * EDGE_INACTIVE => The edge is no longer valid but exists to aid in mesh refinement. NOTE: DO NOT DELETE THESE ENTRIES
     * EDGE_GARBAGE => The edge is no longer valid and is marked for 'GarbageCollection()'
     */
    std::vector<short>  iedge_state;

    /** List of 'state' of each boundary edge. The possible 'state' are:
     * EDGE_ACTIVE => It is a physical edge and should be accounted while computing boundary integrals
     * EDGE_GARBAGE => The edge is no longer valid and is marked for 'GarbageCollection()'
     */
    std::vector<short>  bedge_state;	

    /**
     * \brief Constructs an empty mesh.
     * \note The new mesh has no nodes, elements or parameters.
     * \note The nlist, elist and plist data members must be initialised
     *   to define an actual mesh geometry.
     */
    NonconformingMesh ();

    // destructor
    virtual ~NonconformingMesh ();

    /** Calls the Mesh::Setup() first and then
     *	populates the interior edge and boundary edge data structures and sets them to 'EDGE_ACTIVE' state
     *  NOTE: Call this function before performing mesh refinement or other operations on the mesh
     */	
    void SetupEdgeTables();
    
    // Computes 'rowptr', 'colidx' data structures required to initilaise the system matrix in CRS format
    void SparseRowStructure (int *&rowptr, int *&colidx, int &nzero) const; 
   
    // Outputs the number of interior edges
    int iedgelen() const { return iedge_elist.size(); }

    // Outputs the number of boundary edges
    int bedgelen() const { return bedge_elist.size(); }

    /** Refines a given tetrahedral element
     *  Input:
     *		el -> Global element number to be refined (0 based)
     */			
    int RefineTetElem(const int el);

private:
    // flag to keep track of whether Setup() has been called
    bool is_set_up;

    // Adds a given interior edge
    void PushTriIFace(const int e1, const int e2, const int n1, const int n2, const int n3);

    // Adds a given boundary edge
    void PushTriBFace(const int e1, const int n1, const int n2, const int n3);

    // Finds the edge/face on element number 'el' with nodes  'n1', 'n2', 'n3' and removes it. 
    int FindTrashTriFace(const int el, const int n1, const int n2, const int n3);

    //Goes through each interior and boundary edge and removes them if their 'state' is set to 'EDGE_GARBAGE'
    void GarbageCollection();			
};


/** Computes the terms of the form \langle \phi_{i}^{e} \phi_{j}^{el} \rangle_{\gamma_{l}},  \langle \phi_{i}^{e} \nabla\phi_{j}^{el} \rangle_{\gamma_{l}}
 *  where \phi dentoes the shape functions, \nabla \phi denotes the gradient of shape function 
 *  'e' and 'el' denote two elements sharing an edge/face \gamma_{l}
 *  Inputs:
 *      mesh -> nonconforming mesh object
 *      e -> element number 'e'
 *      el -> element number 'el'
 *      *edge_coord -> Vectors of edge (\gamma_{l}) coordinates
 *  Outputs:
 *       eelbndFD -> 2 x (3x3)(2D) or 3 x (4x4)(3D)
 *       eelbndFF -> 3 x3 (2D) or 4 x 4 (3D) 	
 *  NOTE: Currently only supported for triangular and tetrahedral elements
 */ 
void computeFD_FFeel(NonconformingMesh &mesh, int e, int el, RVector *edge_coord, RDenseMatrix *eelbndFD, RDenseMatrix &eelbndFF);
void computeFFeel(NonconformingMesh &mesh, int e, int el, RVector *edge_coord, RDenseMatrix &eelbndFF);
void computeDDeel(NonconformingMesh &mesh, int e, int el, RVector *edge_coord, RDenseMatrix &eelDD, RVector &normal);

/** Computes the normal to a given edge/face
 * Inputs:
 *      mesh -> nonconforming mesh object
 *      e    -> element number 'e'
 *      *nds -> node numbers of the desired edge
 * Outputs:
 *      normal
 *  NOTE: Currently only supported for triangular and tetrahedral elements
 */	
bool computeEdgeGlobalNormal(NonconformingMesh &mesh, int e, int eside, RVector &normal);

// Add a component to element matrix 'M', given 'mesh' and 'el'
// Element matrix type is defined by 'mode' (see mesh.h)
// nodal or element coefficients are given by 'coeff'
void DGAddToElMatrix (NonconformingMesh &mesh, int el, CGenericSparseMatrix &M, const RVector *coeff, int mode);

// Add a component to system matrix 'M', given 'mesh'
// Element matrix type is defined by 'mode' (see mesh.h)
// nodal coefficients are given by 'coeff'
void DGAddToSysMatrix (NonconformingMesh &mesh, RGenericSparseMatrix &M, const RVector *coeff, int mode);

void DGAddToSysMatrix (NonconformingMesh &mesh, CGenericSparseMatrix &M, const RVector *coeff, int mode);
void DGAddToSysMatrix (NonconformingMesh &mesh, CGenericSparseMatrix &M, const double coeff, int mode);

// Adds interior edge contributions to the system matrix 
void DGAddToSysMatrixInteriorEdgeCont(NonconformingMesh &mesh, CGenericSparseMatrix &M, const RVector *coeff1, const RVector *coeff2);

// Adds interior edge contributions to the system matrix 
void DGAddToSysMatrixBoundaryEdgeCont(NonconformingMesh &mesh, CGenericSparseMatrix &M, const RVector *coeff);

// find an std::vector in std:vector<std::vector>
std::vector< std::vector<int> >::iterator findVecVec(std::vector< std::vector<int> > *source, std::vector<int> target);

double computeBeta(double n1, double n2);
#endif //!_NONCONFORMINGMESH_H
