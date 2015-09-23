// -*-C++-*-
// ==========================================================================
// Module libfe
// File ndlist.h
// Declaration of class NodeList
// ==========================================================================

#ifndef __NDLIST_H
#define __NDLIST_H

// ==========================================================================
// class NodeList

class FELIB NodeList {
public:
    // constructors, destructor
    NodeList ();
    NodeList (int length);
    NodeList (const NodeList &nl);
    ~NodeList ();

    void New (int length);
    // create new (empty) list of with 'length' entries

    void Append (int number);
    // append 'number' empty entries to the list

    void Clear ();
    // clear node list

    void SetList (int no, Node *nds);
    // replaces the node list with array nds of size no, discarding prev list
    // nds should not be deallocated by the caller

    void Remove (int nd);
    // removes node 'nd' from the list. All node numbers higher than 'nd' will
    // be decremented by 1.
    // WARNING: the element list will be invalid after the node is removed!
    // remove all references to the deleted node from the element list and run
    // ElementList::RemoveNodeRef(nd) to update the element list.

    int Len (void) const { return size; };
    // returns number of nodes in the list

    int Exists (const Node &node, double rad = 1e-8) const;
    // returns the number of the node in the list that coincides with 'node'
    // within radius 'rad', or -1 if no node is found
    // If more than one node lies within 'rad' of 'node' then the first in
    // the list is returned

    NodeList &operator= (const NodeList &nl);
    // make *this a copy of nl

    int TotBnd (void) const;
    // returns the number of boundary nodes in the list  OLD!

    int NumberOf (BYTE bndtype) const;
    // returns the number of nodes of type bndtype in the list, where bndtype is
    // any of the following: BND_NONE, BND_DIRICHLET, BND_NEUMANN, BND_CAUCHY,
    // BND_ANY

#ifdef FEM_DEBUG
    Node& operator[] (int rec) const;
#else
    Node& operator[] (int rec) const { return list[rec]; }
#endif
    // returns reference to a node in the list

    void Swap (int nd1, int nd2);
    // exchange elements el1 and el2 in the list

    // friends
    friend std::istream& operator>> (std::istream& is, NodeList& nlist);
    friend std::ostream& operator<< (std::ostream& os, NodeList& nlist);

private:
  int size;	// list length
  Node *list;	// array of nodes
  int dofnod;	// degrees of freedom per node; always 1 so far

};

// ==========================================================================
// friend prototypes

std::istream& operator>> (std::istream& is, NodeList& nlist);
std::ostream& operator<< (std::ostream& os, NodeList& nlist);

#endif // !__NDLIST_H
