// -*-C++-*-
// ==========================================================================
// Module libfe
// File node.h
// Declaration of class Node
//
// Inheritance:
// ------------
// RVector ----> Point ----> Node
// ==========================================================================

#ifndef __NODE_H
#define __NODE_H

#define BND_NONE      0
#define BND_DIRICHLET 3
#define BND_NEUMANN   5
#define BND_ROBIN     (BND_DIRICHLET & BND_NEUMANN)
#define BND_ANY       (BND_DIRICHLET | BND_NEUMANN)
#define BND_INTERNAL  2    // physical boundary for extrapolated b.c.
#define INTERNAL      BND_NONE
#define XLAYER_INTERNAL 4  // internal nodes in extrapolation band
#define BND_VOID      8

// ==========================================================================
// prototypes

class Node;

FELIB void Swap (Node &n1, Node &n2);
FELIB double Dist (const Node &n1, const Node &n2);
FELIB std::istream& operator>> (std::istream& is, Node& nd);
FELIB std::ostream& operator<< (std::ostream& os, Node& nd);


// ==========================================================================
// class Node

class FELIB Node : public Point {
public:
    // constructors
    Node ();
    Node (int dim, char _bndtp = BND_NONE);
    Node (const Node& nd);

    void Copy (const Node &nd);
    // as operator= but reallocates *this, so dimensions need not be the same

    // assignment operators
    Node &operator= (const Node &nd);

    Node &operator= (const Point &pt);

    // relational operators
    bool operator== (const Node& nd) const;
    bool operator!= (const Node& nd) const;

    friend FELIB void Swap (Node &n1, Node &n2);
    // Swap two nodes. Assumes same dimension

    friend FELIB double Dist (const Node &n1, const Node &n2);
    // distance between two nodes

    // boundary type-related functions
    char BndTp() const { return bndtp; }
    bool isInternalInterface() const { return bndtp == 2; } // label 'I'
    bool isBnd() const { return bndtp != 0 &&  bndtp != 2; }
    bool isAnyBnd() const { return bndtp != 0 ; }
    void SetBndTp (char _bndtp) { bndtp = _bndtp; }
    double Phi () const { return phi; }

    // region-related functions
    int Region() const { return region; }
    void SetRegion (int _region) { region = _region; }

    // I/O
    friend FELIB std::istream& operator>> (std::istream& i, Node& nd);
    friend FELIB std::ostream& operator<< (std::ostream& o, Node& nd);

protected:
    char bndtp;         // id for boundary type
    int region;         // region number (-1 = no region)
    double phi;         // Dirichlet condition (assumes bndtp=BND_DIRICHLET)
};


// ==========================================================================
// external variables, defined in node.cc

#ifndef __NODE_CC
extern char BndId[8];
// character identifier for boundary type, defined in node.cc
#endif // !__NODE_CC

#endif // !__NODE_H

