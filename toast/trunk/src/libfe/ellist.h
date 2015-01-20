// -*-C++-*-
// ==========================================================================
// Module libfe
// File ellist.h
// Declaration of class ElementList
// ==========================================================================

#ifndef __ELLIST_H
#define __ELLIST_H

#define BW_TOTAL 0	// bandwidth type ids
#define BW_INTERNAL 1
#define BW_AUTO 2

// ==========================================================================
// class ElementList

class FELIB ElementList {
public:
    ElementList () { Length=0; /*Segmented = false;*/ isKappa = false; };
    ~ElementList ();
    // constructor, destructor

    PElement& operator[] (int rec) const;
    // returns reference to an element in the list

    ElementList &operator= (const ElementList &elist);
    // make *this a copy of elist

    void New (int length);
    // create new (empty) list of with 'length' entries

    void Clear() { New(0); }
    // empty list

    int Len (void) const { return Length; };
    // Length of the list

    void Insert (PElement pel, int pos);
    // inserts the element pointed to by pel into the element list at position
    // pos (0 <= pos <= Len).
    // typical usage: elist.Insert (new Triangle3, 5);

    void Append (PElement pel);
    // appends the element pointer pel to the end of the current element list
    // typical usage: elist.Append (new Triangle3);

    void AppendList (int no, PElement *ppel);
    // Appends an array of element pointers to the end of the current element
    // list. The pointers in the array are expected to point to initialised
    // elements. The array itself (but not the elements) can be deleted after
    // the call

    void SetList (int no, PElement *ppel);
    // Replaces the current list with ppel, length no
    // The array must not be deleted by the calling function

    void Delete (int rec);
    // removes record 'rec' from the list. 0 <= rec < Len

    void RemoveNodeRef (int nd);
    // Updates the node references in all elements under the assuption that
    // node nd has been removed from the node list. If there are still
    // references to nd when RemoveNodeRef is called, these will afterwards be
    // references to the following node (which is then assigned the identifier
    // nd)

    Element *SideNeighbour (int el, int side);
    // Return neighbour element of 'el' at 'side', or NULL if side has
    // no neighbour, i.e. if the side is part of the mesh surface
    // Note: This function requires a preceeding call to
    // Mesh::PopulateNeighbourLists
    // Note: This replaces 'EdgeAdjacentElement'

    int EdgeAdjacentElement (int el, int side) const;
    // returns the number of the element adjacent to 'el' at side 'side'
    // returns -1 if no element is found, i.e. if 'side' is a boundary

    int VertexAdjacentElement (int el, int nd, int min_el=0);
    // Returns the number of an element sharing the node 'nd' with 'el'.
    // nd is the local node number: 0 <= nd < nNode
    // Searching starts from element number 'min_el' upwards
    // Returns -1 if no matching element >= min_el is found

    friend std::istream& operator>> (std::istream& i, ElementList& nlist);
    friend std::ostream& operator<< (std::ostream& o, ElementList& elist);
    // read/write element list from/to streams

    //bool Segmented;
    // TRUE if elements contain region information

    bool isKappa;
    // TRUE if Coeff[MUS] contains c*kappa for each element in the list

protected:

    PElement *List;
    // array of pointers to elements

    int Length;
    // length of list

};

// ==========================================================================
// inline functions

#ifndef FEM_DEBUG
inline PElement& ElementList::operator[] (int rec) const
{
    return List[rec];
}
#endif // !FEM_DEBUG

#endif // !__ELLIST_H
