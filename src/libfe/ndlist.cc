// ==========================================================================
// Module libfe
// File ndlist.cc
// Definition of class NodeList
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "arch.h"

using namespace std;

NodeList::NodeList ()
{
    size = 0;
    dofnod = 1;
    New (0);
}

NodeList::NodeList (int length)
{
    size = 0;
    dofnod = 1;
    New (length);
}

NodeList::NodeList (const NodeList &nl)
{
    size = 0;
    dofnod = nl.dofnod;
    New (nl.Len());
    for (int i = 0; i < size; i++)
	list[i].Copy (nl.list[i]);
}

NodeList::~NodeList ()
{
    New (0);
}

void NodeList::New (int length)
{
    if (size) delete []list;
    size = length;
    if (size) {
	list = new Node[size];
	dASSERT (list, "Memory allocation failed.");
    } else list = 0;
}

void NodeList::Clear()
{
    New (0);
}

void NodeList::Append (int number)
{
    Node *TmpList = new Node[size+number];
    dASSERT(TmpList, "Memory allocation failed.");
    for (int i = 0; i < size; i++) {
	TmpList[i].New (list[i].Dim());
	TmpList[i] = list[i];
    }
    if (size) delete []list;
    list = TmpList;
    size += number;
}

void NodeList::SetList (int no, Node *nds)
{
    Clear();
    list = nds;
    size = no;
}

void NodeList::Remove (int nd)
{
    dASSERT(nd >= 0 && nd < size, "Index out of range.");
    int i;
    Node *TmpList = new Node[size-1];
    dASSERT(TmpList, "Memory allocation failed.");
    for (i = 0; i < nd; i++) TmpList[i].Copy(list[i]);
    for (i = nd+1; i < size; i++) TmpList[i-1].Copy(list[i]);
    delete []list;
    list = TmpList;
    size--;
}

int NodeList::Exists (const Node &node, double rad) const
{
    for (int i = 0; i < size; i++)
        if (Dist (node, list[i]) < rad) return i;
    return -1;
}

NodeList &NodeList::operator= (const NodeList &nl)
{
    New (nl.Len());
    for (int i = 0; i < size; i++)
	list[i].Copy (nl.list[i]);
    dofnod = nl.dofnod;
    return *this;
}

int NodeList::TotBnd (void) const
{
    int n, bnd;

    for (bnd = n = 0; n < size; n++) if (list[n].isBnd()) bnd++;
    return bnd;
}

int NodeList::NumberOf (BYTE bndtype) const
{
    int n, bnd;

    if (bndtype != BND_ANY) {
	for (bnd = n = 0; n < size; n++) if (list[n].BndTp() == bndtype) bnd++;
    } else {
	for (bnd = n = 0; n < size; n++) if (list[n].isBnd()) bnd++;
    }
    return bnd;
}

void NodeList::Swap (int nd1, int nd2)
{
    dASSERT(nd1 >= 0 && nd1 < size && nd2 >= 0 && nd2 < size,
	"Index out of range.");
    ::Swap (list[nd1], list[nd2]);
}

#ifdef FEM_DEBUG
Node& NodeList::operator[] (int rec) const
{
    dASSERT(rec >= 0 && rec < size, "Index out of range.");
    return list[rec];
}
#endif

istream& operator>> (istream& is, NodeList& nlist)
{
    char cbuf[200];
    int length, i;

    do {
	is >> cbuf;
    } while (strcasecmp (cbuf, "NodeList"));
    is >> length >> nlist.dofnod;
    nlist.New (length);
    is.getline (cbuf, 200); // skip to eol
    for (i = 0; i < length; i++) is >> nlist.list[i];
    return is;
}

ostream& operator<< (ostream& os, NodeList& nlist)
{
    os << "NodeList " << nlist.size << " " << nlist.dofnod << endl;
    for (int i = 0; i < nlist.size; i++)
	os << nlist.list[i] << endl;
    return os;
}

