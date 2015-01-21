// ==========================================================================
// Module libfe
// File ellist.cc
// Definition of class ElementList
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "arch.h"

using namespace std;

ElementList::~ElementList ()
{
    Clear();
}

void ElementList::New (int length)
{
    int i;

    if (Length) {
        for (i = 0; i < Length; i++)
	    if (List[i]) delete List[i];
        delete []List;
    }
    if ((Length = length)) {
        List = new PElement[Length];
	for (i = 0; i < Length; i++) List[i] = 0;
    } else {
        List = 0;
    }
}

void ElementList::Insert (PElement pel, int pos)
{
    dASSERT(pos >= 0 && pos <= Length, "Index out of range.");
    PElement *TmpList = new PElement[Length+1];
    memcpy (TmpList, List, pos * sizeof (PElement));
    memcpy (TmpList+pos+1, List+pos, (Length-pos) * sizeof (PElement));
    TmpList[pos] = pel;
    if (Length) delete []List;
    List = TmpList;
    Length++;
}

void ElementList::Append (PElement pel)
{
    Insert (pel, Length);
}

void ElementList::Delete (int rec)
{
    dASSERT(rec >= 0 && rec < Length, "Index out of range.");
    delete List[rec];
    PElement *TmpList = new PElement[--Length];
    memcpy (TmpList, List, rec * sizeof (PElement));
    memcpy (TmpList+rec, List+rec+1, (Length-rec) * sizeof (PElement));
    delete []List;
    List = TmpList;
}

ElementList &ElementList::operator= (const ElementList &elist)
{
    New (elist.Len());
    for (int i = 0; i < Length; i++)
	List[i] = elist[i]->Copy();
    isKappa = elist.isKappa;
    return *this;
}

void ElementList::SetList (int no, PElement *ppel)
{
    Clear(); // clear previous list
    List = new PElement[no];
    memcpy (List, ppel, no*sizeof(PElement));
    //List = ppel;
    Length = no;
}

void ElementList::AppendList (int no, PElement *ppel)
{
    PElement *TmpList = new PElement[Length+no];
    memcpy (TmpList, List, Length * sizeof (PElement));
    memcpy (TmpList+Length, ppel, no * sizeof (PElement));
    if (Length) delete []List;
    List = TmpList;
    Length += no;
}

void ElementList::RemoveNodeRef (int nd)
{
    for (int el = 0; el < Length; el++)
	for (int n = 0; n < List[el]->nNode(); n++)
	    if (List[el]->Node[n] > nd) List[el]->Node[n]--;
}

Element *ElementList::SideNeighbour (int el, int side)
{
    dASSERT(el >= 0 && el < Length, "Element index out of range");
    dASSERT(side >= 0 && side < List[el]->nSide(), "Side index out of range");

    Element *pel = List[el];
    dASSERT(pel->sdnbhr, "Element neighour lists not available");

    return pel->sdnbhr[side];
}

int ElementList::EdgeAdjacentElement (int el, int side) const
{
    const int max_sidenode = 10;
    int i, i1, i2, nn, nd[max_sidenode];

    dASSERT(side >= 0 && side < List[el]->nSide(), "Index out of range.");
    nn = List[el]->nSideNode (side);
    dASSERT(nn <= max_sidenode, "Nodes per side limit exceeded");
    for (i = 0; i < nn; i++)
        nd[i] = List[el]->Node[List[el]->SideNode (side, i)];

    // search from el outward to improve efficiency
    for (i1 = el-1, i2 = el+1; i1 >= 0 || i2 < Length; i1--, i2++) {
        if (i1 >= 0 && (List[i1]->IsSide (nn, nd) >= 0)) return i1;
	if (i2 < Length && (List[i2]->IsSide (nn, nd) >= 0)) return i2;
    }
    return -1;
}

int ElementList::VertexAdjacentElement (int el, int nd, int min_el)
{
    int n, i;

    dASSERT(nd >= 0 && nd < List[el]->nNode(), "Index out of range.");
    n = List[el]->Node[nd];
    for (i = min_el; i < Length; i++)
	if (i != el && List[i]->IsNode (n)) return i;
    return -1;
}

#ifdef FEM_DEBUG // otherwise inline
PElement& ElementList::operator[] (int rec) const
{
    dASSERT(rec >= 0 && rec < Length, "Index out of range.");
    return List[rec];
}
#endif // FEM_DEBUG

istream& operator>> (istream& i, ElementList& elist)
{
    char cbuf[200], cid;
    int id, el;

    while (i.getline (cbuf, 200) && strncasecmp (cbuf, "ElementList", 11));
    if (!i) return i;	// not found
    elist.Clear();	// clear previous list
    std::istringstream iss(cbuf+11);
    iss >> elist.Length;
    elist.List = new PElement[elist.Length];
    for (el=0; el<elist.Length; el++) {
	do { i >> cid; } while (cid==' ');
	if (cid>='a' && cid<='z') id = 1+(cid-'a'); // new type descriptor
	else                      id = 1+(cid-'1'); // old type descriptor
	switch (id) {
  	    case ELID_TRI3:    elist.List[el] = new Triangle3; break;
	    case ELID_TRI3OLD: elist.List[el] = new Triangle3old; break;
	    case ELID_TRI6:    elist.List[el] = new Triangle6; break;
//	    case ELID_RCT4: elist.List[el] = new Rectangle4; break;
	    case ELID_TET4:    elist.List[el] = new Tetrahedron4; break;
	    case ELID_TET10:   elist.List[el] = new Tetrahedron10; break;
	    case ELID_WDG6:    elist.List[el] = new Wedge6; break;
	    case ELID_VOX8:    elist.List[el] = new Voxel8; break;
   	    case ELID_PIX4:    elist.List[el] = new Pixel4; break;
	    case ELID_TRI6_IP: elist.List[el] = new Triangle6_ip; break;
	    case ELID_TRI10:   elist.List[el] = new Triangle10; break;
	    case ELID_TRI10_IP: elist.List[el] = new Triangle10_ip; break;
	    case ELID_TET10_IP: elist.List[el] = new Tetrahedron10_ip; break;
	    case ELID_TRI3D3: elist.List[el] = new Triangle3D3; break;
	    case ELID_TRI3D6: elist.List[el] = new Triangle3D6; break;
	    case ELID_LINE2D2: elist.List[el] = new Line2D2; break;
	    default: xERROR("Unknown element type.");
	}
	i >> *elist.List[el];
    }
    return i;
}  

ostream& operator<< (ostream& o, ElementList& elist)
{
    o << "ElementList " << elist.Length << endl;
    for (register int el=0; el<elist.Length; el++) o << *elist.List[el];
    return o;
}

