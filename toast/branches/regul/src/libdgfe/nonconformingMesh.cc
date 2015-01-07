// ==========================================================================
// Module libfe
// File nonconforminMesh.cc
// Definition of class nonconformingMesh
// ==========================================================================

#define FELIB_IMPLEMENTATION

#include <iostream>
#include <fstream>
#include <string.h>
#include <algorithm>
#include <memory>
#include "nonconformingMesh.h"
#include "mathlib.h"

using namespace std;

NonconformingMesh::NonconformingMesh()
{
    is_set_up = false;
}

NonconformingMesh::~NonconformingMesh()
{
}

//calls super.setup() and setsup iedgelist and bedgelist
void NonconformingMesh::SetupEdgeTables()
{
    int numside_node;
    // Sets up the mesh
    Mesh::Setup();        
    
    std::vector<int> vec1, vec2;
    vec1.reserve(2); vec2.reserve(elist[0]->nNode());
    int *nds = new int[10];
    //nds = (int *)malloc(10*sizeof(int));
   /* iedge_elist.reserve(3*elen());
    bedge_elist.reserve(elen());
    iedge_nlist.reserve(3*elen());
    bedge_nlist.reserve(elen());
    iedge_state.reserve(3*elen());
    bedge_state.reserve(elen());
    for(int i=0; i < 3*elen(); i++){
	iedge_elist[i].reserve(2);
	iedge_nlist[i].reserve(3);
     }

    for(int i=0; i < elen(); i++)
	bedge_nlist[i].reserve(3);	
*/
   	
    //setup the iedgelist and bedgelist
    for(int e=0; e < elen(); e++)
    {
//	if((e+1)%1000==0) cout<<"Elements: "<<e+1<<endl;
	for(int side=0; side < elist[e]->nSide(); side++) 
	{       numside_node = elist[e]->nSideNode(side);
		int el = elist.EdgeAdjacentElement(e, side);
		if(el != -1) 
		{
			if(e<el)
			{
				vec1.clear();vec2.clear();
				vec1.push_back(e);
				vec1.push_back(el);
	 	      			
				for(int node_num = 0 ; node_num < numside_node; node_num++)
					nds[node_num] = elist[e]->Node[elist[e]->SideNode(side, node_num)];
				std::sort(nds, nds+numside_node);
			
				for(int node_num =0 ; node_num < numside_node; node_num++)
					vec2.push_back(nds[node_num]);
				iedge_state.push_back(EDGE_ACTIVE);
				iedge_elist.insert(iedge_elist.end(), vec1);
				iedge_nlist.insert(iedge_nlist.end(), vec2);
			}
		}
		else
		{
			vec2.clear();
			bedge_elist.insert(bedge_elist.end(), e);
			
			for(int node_num = 0 ; node_num < numside_node; node_num++)
				nds[node_num] = elist[e]->Node[elist[e]->SideNode(side, node_num)];
			std::sort(nds, nds+numside_node);
			
			for(int node_num =0 ; node_num < numside_node; node_num++)
				vec2.push_back(nds[node_num]);
			bedge_state.push_back(EDGE_ACTIVE);	
			bedge_nlist.insert(bedge_nlist.end(), vec2);
		}
	}	
    }
    delete []nds;	
    is_set_up = true;
}

// Hack: assumes all elements in the mesh are of the same type 
void NonconformingMesh::SparseRowStructure (int *&rowptr, int *&colidx, int &nzero) const
{
    // M.S. 1.10.99: Check that this works for higher-order element types
    // (TRI6 and TET10)
   int dim = elist[0]->nNode(); // Hack: assumes all elements in the mesh are of the same type 
   std::vector< std::vector<int> > :: const_iterator it;
   
   std::vector<int>  ed;
   int e, el;
   std::vector<int> *temp;
   temp =  new std::vector<int>[dim*elen()];	
    for(int i=0; i<dim*elen(); i++)
	temp[i].reserve((dim+1)*dim);

    
   
    for(int elem=0; elem<elen(); elem++){
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			temp[dim*elem+i].push_back(dim*elem+j);
		}
	}
    }
    int count=0;;
     for (it = iedge_elist.begin(); it != iedge_elist.end(); it++, count++) {
	if(iedge_state[count] == EDGE_ACTIVE){
	ed = *it;
	e = ed[0]; el = ed[1];
	for(int i=0; i<dim; i++){
		for(int j=0; j<dim; j++){
			temp[dim*e+i].push_back(dim*el+j);
			temp[dim*el+i].push_back(dim*e+j);
		}
	}
     }}
     
     rowptr = new int[dim*elen()+1];
     rowptr[0] = 0;
     rowptr[1] = temp[0].size();
     for(int i=1; i<dim*elen(); i++){
	rowptr[i+1] = rowptr[i] + temp[i].size();
	sort(temp[i].begin(), temp[i].end());
	}

     nzero = rowptr[dim*elen()];	
     colidx = new int[nzero];
     int k=0;
     std::vector<int> :: iterator it2;
     for(int i=0; i<dim*elen(); i++){
	for(it2 = temp[i].begin(); it2 != temp[i].end(); it2++){
		colidx[k] = *it2;
		k++;
	}
     }
    delete []temp;	
}

int NonconformingMesh::RefineTetElem(const int el){
 int nel, nnode, niedges, nbedges;
 nel =  elen(); nnode = nlen(); niedges = iedgelen(); nbedges = bedgelen();
 Node *node = new Node(3, BND_NONE);
 int nds[6];

//cout<<"Refining:  "<<el<<endl;

 /*  OBTAIN ALL THE UNIQUE LINES IN THE TETRAHEDRA */
 int pedges[][2] ={elist[el]->Node[elist[el]->SideNode(0, 0)], elist[el]->Node[elist[el]->SideNode(0, 1)], // 0-1 
		    elist[el]->Node[elist[el]->SideNode(0, 0)], elist[el]->Node[elist[el]->SideNode(0, 2)], // 0-2
		    elist[el]->Node[elist[el]->SideNode(0, 1)], elist[el]->Node[elist[el]->SideNode(0, 2)], // 1-2	
		    elist[el]->Node[elist[el]->SideNode(2, 0)], elist[el]->Node[elist[el]->SideNode(2, 2)], // 0-3
		    elist[el]->Node[elist[el]->SideNode(2, 1)], elist[el]->Node[elist[el]->SideNode(2, 2)], // 2-3	
                    elist[el]->Node[elist[el]->SideNode(1, 2)], elist[el]->Node[elist[el]->SideNode(1, 1)]}; //1-3
/*cout<<pedges[0][0] << "  "<<pedges[0][1]<<endl;
cout<<pedges[1][0] << "  "<<pedges[1][1]<<endl;
cout<<pedges[2][0] << "  "<<pedges[2][1]<<endl;
cout<<pedges[3][0] << "  "<<pedges[3][1]<<endl;
cout<<pedges[4][0] << "  "<<pedges[4][1]<<endl;
cout<<pedges[5][0] << "  "<<pedges[5][1]<<endl;
getchar();
*/


 /* START ADDING NODES */
 int curr_node = nnode, indx;
 for(int edge = 0; edge < 6; edge++){ (*node)[0]=0.0; (*node)[1]=0.0; (*node)[2]=0.0;
 	for(int nedgenode=0; nedgenode < 2; nedgenode++)
		for(int dim =0 ; dim < 3 ; dim++){ indx = pedges[edge][nedgenode];
			(*node)[dim] += nlist[indx][dim]/2.0;}
	nds[edge] = nlist.Exists(*node, 1e-6); // check if the node already exists before adding it
	if( nds[edge] == -1){
		nds[edge] = curr_node;
		nlist.Append(1);
		nlist[curr_node].New(3);
		nlist[curr_node].SetBndTp(BND_NONE);
		nlist[curr_node] = (*node);
		curr_node++;
	}
		
}
/*END ADDING NODES*/

/*START ADDING ELEMENTS */
int curr_nel = nel;

elist[el]->Node[0] = pedges[0][0]; elist[el]->Node[1] = nds[0]; elist[el]->Node[2] = nds[1]; elist[el]->Node[3] = nds[3]; //0-4-5-7
elist[el]->Initialise(nlist);
if(elist[el]->Size() < 0){  elist[el]->Node[2] = nds[3]; elist[el]->Node[3] = nds[1]; 
elist[el]->Initialise(nlist);}

for(int i=0; i < 7; i++) elist.Append(new Tetrahedron4);
elist[curr_nel]->Node[0] = pedges[0][1]; elist[curr_nel]->Node[1] = nds[0]; elist[curr_nel]->Node[2] = nds[2]; elist[curr_nel]->Node[3] = nds[5];//1-4-6-9
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0){ elist[curr_nel]->Node[2] = nds[5]; elist[curr_nel]->Node[3] = nds[2]; 
elist[curr_nel]->Initialise(nlist);}
curr_nel++;

elist[curr_nel]->Node[0] = pedges[1][1]; elist[curr_nel]->Node[1] = nds[1]; elist[curr_nel]->Node[2] = nds[2]; elist[curr_nel]->Node[3] = nds[4];//2-5-6-8
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0) {elist[curr_nel]->Node[2] = nds[4]; elist[curr_nel]->Node[3] = nds[2];
elist[curr_nel]->Initialise(nlist);}
curr_nel++;

elist[curr_nel]->Node[0] = pedges[3][1]; elist[curr_nel]->Node[1] = nds[3]; elist[curr_nel]->Node[2] = nds[4]; elist[curr_nel]->Node[3] = nds[5];//3-7-8-9
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0){ elist[curr_nel]->Node[2] = nds[5]; elist[curr_nel]->Node[3] = nds[4];
elist[curr_nel]->Initialise(nlist);}
curr_nel++;

elist[curr_nel]->Node[0] = nds[2]; elist[curr_nel]->Node[1] = nds[3]; elist[curr_nel]->Node[2] = nds[4]; elist[curr_nel]->Node[3] = nds[5]; //6-7-8-9
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0){ elist[curr_nel]->Node[2] = nds[5]; elist[curr_nel]->Node[3] = nds[4];
elist[curr_nel]->Initialise(nlist);}
curr_nel++;

elist[curr_nel]->Node[0] = nds[0]; elist[curr_nel]->Node[1] = nds[2]; elist[curr_nel]->Node[2] = nds[3]; elist[curr_nel]->Node[3] = nds[5]; //4-6-7-9
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0) { elist[curr_nel]->Node[2] = nds[5]; elist[curr_nel]->Node[3] = nds[3];
elist[curr_nel]->Initialise(nlist);}
curr_nel++;

elist[curr_nel]->Node[0] = nds[0]; elist[curr_nel]->Node[1] = nds[1]; elist[curr_nel]->Node[2] = nds[2]; elist[curr_nel]->Node[3] = nds[3];//4-5-6-7
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0) {elist[curr_nel]->Node[2] = nds[3]; elist[curr_nel]->Node[3] = nds[2];
elist[curr_nel]->Initialise(nlist);}
curr_nel++;

elist[curr_nel]->Node[0] = nds[1]; elist[curr_nel]->Node[1] = nds[2]; elist[curr_nel]->Node[2] = nds[3]; elist[curr_nel]->Node[3] = nds[4]; //5-6-7-8
elist[curr_nel]->Initialise(nlist);
if(elist[curr_nel]->Size() < 0){ elist[curr_nel]->Node[2] = nds[4]; elist[curr_nel]->Node[3] = nds[3];
elist[curr_nel]->Initialise(nlist);}
curr_nel++;
/*END ADDING ELEMENTS */

std::vector< std::vector<int> > :: iterator it, itn1, ite1;
std::vector<int>::iterator itb1;
std::vector<short>::iterator itisA, itisB;
std::vector<int> ed, vec2, vec1;
int el2;

/* START ADDING FACES */

/* ADDING NEW INTERNAL FACES WHICH IS TRIVIAL */
PushTriIFace(el, nel+5, nds[0], nds[1], nds[3]); // edge 4-5-7
PushTriIFace(nel, nel+4, nds[0], nds[2], nds[5]); //edge 4-6-9
PushTriIFace(nel+1, nel+6, nds[1], nds[2], nds[4]); //edge 5-6-8
PushTriIFace(nel+2, nel+3, nds[3], nds[4], nds[5]); //edge 7-8-9
PushTriIFace(nel+3, nel+4, nds[2], nds[3], nds[5]); //edge 6-7-9
PushTriIFace(nel+4, nel+5, nds[0], nds[2], nds[3]); //edge 4-6-7
PushTriIFace(nel+5, nel+6, nds[1], nds[2], nds[3]);//edge 5-6-7
PushTriIFace(nel+3, nel+6, nds[2], nds[3], nds[4]);//edge 6-7-8
//cout<<"Added trivial edges ..."<<endl;

/*UPDATING EXISTING FACES WHICH IS TRICKY*/
// 0-1-2 to replace with 4 edges 0-4-5, 4-5-6, 1-4-6, 2-5-6
vec2.push_back(pedges[0][0]); vec2.push_back(pedges[0][1]); vec2.push_back(pedges[1][1]); sort(vec2.begin(), vec2.end()); //old edge 0-1-2
itn1 = findVecVec(&iedge_nlist, vec2);
if(itn1 != iedge_nlist.end()) // old edge 0-1-2 is an internal edge
{
	for(it = iedge_nlist.begin(), ite1 = iedge_elist.begin(), itisA = iedge_state.begin(); it != itn1; it++, ite1++, itisA++); // found e1-e2 in old 1-2-3
	ed = *ite1;
	el2 = ((ed[0] == el) ? ed[1] : ed[0]);
	
	if(*itisA == EDGE_ACTIVE) // 0-1-2 is 'active'
	{
		*itisA = EDGE_INACTIVE; // set edge 0-1-2 to 'inactive'
		PushTriIFace(el, el2, pedges[0][0], nds[0], nds[1]); // insert edge 0-4-5
		PushTriIFace(nel+5, el2, nds[0], nds[1], nds[2]); // insert edge 4-5-6 
		PushTriIFace(nel, el2, pedges[0][1], nds[0], nds[2]); //insert edge 1-4-6 
	      	PushTriIFace(nel+1, el2, pedges[1][1], nds[1], nds[2]); //insert edge 2-5-6 
	}
	else //0-1-2 is 'inactive'
	{
		*itisA = EDGE_GARBAGE; // mark old edge 0-1-2 for garbage collection
		/*insert edge 0-4-5 */
		el2 = FindTrashTriFace(el, pedges[0][0], nds[0], nds[1]);
		PushTriIFace(el, el2, pedges[0][0], nds[0], nds[1]);//push the new 0-4-5

 		/*insert edge 4-5-6 */
		el2 = FindTrashTriFace(el, nds[0], nds[1], nds[2]);
		PushTriIFace(nel+5, el2, nds[0], nds[1], nds[2]);//push the new 0-4-5

		/*insert edge 1-4-6 */
		el2 = FindTrashTriFace(el, pedges[0][1], nds[0], nds[2]);
		PushTriIFace(nel, el2, pedges[0][1], nds[0], nds[2]);//push the new 1-4-6

		/*insert edge 2-5-6 */
		el2 = FindTrashTriFace(el, pedges[1][1], nds[1], nds[2]);           
		PushTriIFace(nel+1, el2, pedges[1][1], nds[1], nds[2]);//push the new 2-5-6
	}
}
else // 0-1-2 is a boundary edge
{
	itn1 = findVecVec(&bedge_nlist, vec2);
	for(it = bedge_nlist.begin(), itisB = bedge_state.begin(); it != itn1; it++, itisB++);
	*itisB = EDGE_GARBAGE;
	PushTriBFace(el, pedges[0][0], nds[0], nds[1]);// push the new 0-4-5
	PushTriBFace(nel+5, nds[0], nds[1], nds[2]);//push the new 4-5-6
	PushTriBFace(nel, pedges[0][1], nds[0], nds[2]);//push the new 1-4-6
	PushTriBFace(nel+1, pedges[1][1], nds[1], nds[2]); // push the new 2-5-6
	
}


// 0-1-3 to replace with 4 edges 0-4-7, 1-4-9, 4-7-9, 3-7-9
vec2.clear();
vec2.push_back(pedges[0][0]); vec2.push_back(pedges[0][1]); vec2.push_back(pedges[3][1]); sort(vec2.begin(), vec2.end()); //old edge 0-1-3
itn1 = findVecVec(&iedge_nlist, vec2);

if(itn1 != iedge_nlist.end()) // old edge 0-1-3 is an internal edge
{
	for(it = iedge_nlist.begin(), ite1 = iedge_elist.begin(), itisA = iedge_state.begin(); it != itn1; it++, ite1++, itisA++); // found e1-e2 in old 1-2-3

	ed = *ite1;
	el2 = ((ed[0] == el) ? ed[1] : ed[0]);
	if(*itisA == EDGE_ACTIVE) // 0-1-3 is 'active'
	{
		*itisA = EDGE_INACTIVE; // set edge 0-1-2 to 'inactive' 
		PushTriIFace(el, el2, pedges[0][0], nds[0], nds[3]); // insert edge 0-4-7
		PushTriIFace(nel, el2, pedges[0][1], nds[0], nds[5]); // insert edge 1-4-9 
		PushTriIFace(nel+4, el2, nds[0], nds[3], nds[5]); //insert edge 4-7-9 
	      	PushTriIFace(nel+2, el2, pedges[3][1], nds[3], nds[5]); //insert edge 3-7-9 
	}
	else //0-1-3 is 'inactive'
	{
		*itisA = EDGE_GARBAGE; // mark old edge 0-1-3 for garbage collection
		/*insert edge 0-4-7 */
		el2 = FindTrashTriFace(el, pedges[0][0], nds[0], nds[3]);
		PushTriIFace(el, el2, pedges[0][0], nds[0], nds[3]);

 		/*insert edge 1-4-9 */
		el2 = FindTrashTriFace(el, pedges[0][1], nds[0], nds[5]);
		PushTriIFace(nel, el2, pedges[0][1], nds[0], nds[5]);

		/*insert edge 4-7-9 */
		el2 = FindTrashTriFace(el, nds[0], nds[3], nds[5]);
		PushTriIFace(nel+4, el2, nds[0], nds[3], nds[5]);

		/*insert edge 3-7-9 */
		el2 = FindTrashTriFace(el, pedges[3][1], nds[3], nds[5]);
		PushTriIFace(nel+2, el2, pedges[3][1], nds[3], nds[5]);
	}
}
else // 0-1-3 is a boundary edge
{
	itn1 = findVecVec(&bedge_nlist, vec2);
	for(it = bedge_nlist.begin(), itisB = bedge_state.begin(); it != itn1; it++, itisB++);
	*itisB = EDGE_GARBAGE;
	PushTriBFace(el, pedges[0][0], nds[0], nds[3]);// push the new 0-4-7
	PushTriBFace(nel, pedges[0][1], nds[0], nds[5]);//push the new 1-4-9
	PushTriBFace(nel+4, nds[0], nds[3], nds[5]);//push the new 4-7-9
	PushTriBFace(nel+2, pedges[3][1], nds[3], nds[5]); // push the new 2-5-6
	
}

// 0-2-3 to replace with 4 edges 0-5-7, 2-5-8, 5-7-8, 3-7-8
vec2.clear();
vec2.push_back(pedges[0][0]); vec2.push_back(pedges[1][1]); vec2.push_back(pedges[3][1]); sort(vec2.begin(), vec2.end()); //old edge 0-2-3
itn1 = findVecVec(&iedge_nlist, vec2);

if(itn1 != iedge_nlist.end()) // old edge 0-2-3 is an internal edge
{
		for(it = iedge_nlist.begin(), ite1 = iedge_elist.begin(), itisA = iedge_state.begin(); it != itn1; it++, ite1++, itisA++); // found e1-e2 in old 1-2-3
	ed = *ite1;
	el2 = ((ed[0] == el) ? ed[1] : ed[0]);
	if(*itisA == EDGE_ACTIVE) // 0-1-3 is 'active'
	{	
		*itisA = EDGE_INACTIVE; // set edge 0-2-3 to 'inactive' 
		PushTriIFace(el, el2, pedges[0][0], nds[1], nds[3]); // insert edge 0-5-7
		PushTriIFace(nel+1, el2, pedges[1][1], nds[1], nds[4]); // insert edge 2-5-8 
		PushTriIFace(nel+6, el2, nds[1], nds[3], nds[4]); //insert edge 5-7-8 
	      	PushTriIFace(nel+2, el2, pedges[3][1], nds[3], nds[4]); //insert edge 3-7-8 
	}
	else //0-2-3 is 'inactive'
	{
		*itisA = EDGE_GARBAGE; // mark old edge 0-2-3 for garbage collection
		/*insert edge 0-5-7 */
		el2 = FindTrashTriFace(el, pedges[0][0], nds[1], nds[3]);
		PushTriIFace(el, el2, pedges[0][0], nds[1], nds[3]);

 		/*insert edge 2-5-8 */
		el2 = FindTrashTriFace(el, pedges[1][1], nds[1], nds[4]);
		PushTriIFace(nel+1, el2, pedges[1][1], nds[1], nds[4]);

		/*insert edge 5-7-8 */
		el2 = FindTrashTriFace(el, nds[1], nds[3], nds[4]);
		PushTriIFace(nel+6, el2, nds[1], nds[3], nds[4]);

		/*insert edge 3-7-8 */
		el2 = FindTrashTriFace(el, pedges[3][1], nds[3], nds[4]);
		PushTriIFace(nel+2, el2, pedges[3][1], nds[3], nds[4]);
	}
}
else // 0-2-3 is a boundary edge
{
	itn1 = findVecVec(&bedge_nlist, vec2);
	for(it = bedge_nlist.begin(), itisB = bedge_state.begin(); it != itn1; it++, itisB++);
	*itisB = EDGE_GARBAGE;
	PushTriBFace(el, pedges[0][0], nds[1], nds[3]);// push the new 0-5-7
	PushTriBFace(nel+1, pedges[1][1], nds[1], nds[4]);//push the new 2-5-8
	PushTriBFace(nel+6, nds[1], nds[3], nds[4]);//push the new 5-7-8
	PushTriBFace(nel+2, pedges[3][1], nds[3], nds[4]); // push the new 3-7-8
	
}

// 1-2-3 to replace with 4 edges 1-6-9, 2-6-8, 3-8-9, 6-8-9
vec2.clear();
vec2.push_back(pedges[0][1]); vec2.push_back(pedges[1][1]); vec2.push_back(pedges[3][1]); sort(vec2.begin(), vec2.end()); //old edge 1-2-3
itn1 = findVecVec(&iedge_nlist, vec2);
if(itn1 != iedge_nlist.end()) // old edge 1-2-3 is an internal edge
{
	for(it = iedge_nlist.begin(), ite1 = iedge_elist.begin(), itisA = iedge_state.begin(); it != itn1; it++, ite1++, itisA++); // found e1-e2 in old 1-2-3

	
	ed = *ite1;
	el2 = ((ed[0] == el) ? ed[1] : ed[0]);
	if(*itisA == EDGE_ACTIVE) // 1-2-3 is 'active'
	{
		*(itisA) = EDGE_INACTIVE; // set edge 1-2-3 to 'inactive'
		PushTriIFace(nel, el2, pedges[0][1], nds[2], nds[5]); // insert edge 1-6-9
		PushTriIFace(nel+1, el2, pedges[1][1], nds[2], nds[4]); // insert edge 2-6-8
		PushTriIFace(nel+2, el2, pedges[3][1], nds[4], nds[5]); //insert edge 3-8-9
	      	PushTriIFace(nel+3, el2, nds[2], nds[4], nds[5]); //insert edge 6-8-9
	}
	else //1-2-3 is 'inactive'
	{
		*itisA = EDGE_GARBAGE; // mark old edge 1-2-3 for garbage collection
		/*insert edge 1-6-9 */
		el2 = FindTrashTriFace(el, pedges[0][1], nds[2], nds[5]);
		PushTriIFace(nel, el2, pedges[0][1], nds[2], nds[5]);

 		/*insert edge 2-6-8 */
		el2 = FindTrashTriFace(el, pedges[1][1], nds[2], nds[4]);
		PushTriIFace(nel+1, el2, pedges[1][1], nds[2], nds[4]);

		/*insert edge 3-8-9 */
		el2 = FindTrashTriFace(el, pedges[3][1], nds[4], nds[5]);
		PushTriIFace(nel+2, el2, pedges[3][1], nds[4], nds[5]);

		/*insert edge 6-8-9 */
		el2 = FindTrashTriFace(el, nds[2], nds[4], nds[5]);
		PushTriIFace(nel+3, el2, nds[2], nds[4], nds[5]);
	}
}
else // 1-2-3 is a boundary edge
{
	itn1 = findVecVec(&bedge_nlist, vec2);
	for(it = bedge_nlist.begin(), itisB = bedge_state.begin(); it != itn1; it++, itisB++);
	*itisB = EDGE_GARBAGE;
	PushTriBFace(nel, pedges[0][1], nds[2], nds[5]);// push the new 1-6-9
	PushTriBFace(nel+1, pedges[1][1], nds[2], nds[4]);//push the new 2-6-8
	PushTriBFace(nel+2, pedges[3][1], nds[4], nds[5]);//push the new 3-8-9
	PushTriBFace(nel+3, nds[2], nds[4], nds[5]); // push the new 6-8-9
	
}

GarbageCollection();
MarkBoundary();
Mesh::Setup();
/* END ADDING INTERNAL EDGES*/
    return 0;
}


void NonconformingMesh::PushTriIFace(const int e1, const int e2, const int n1, const int n2, const int n3)
{
	std::vector<int> vec1, vec2;
	vec1.push_back(e1); vec1.push_back(e2); sort(vec1.begin(), vec1.end());
	vec2.push_back(n1); vec2.push_back(n2); vec2.push_back(n3); sort(vec2.begin(), vec2.end());
	iedge_elist.insert(iedge_elist.end(), vec1); iedge_nlist.insert(iedge_nlist.end(), vec2); iedge_state.push_back(EDGE_ACTIVE);


}
void NonconformingMesh::PushTriBFace(const int e1, const int n1, const int n2, const int n3)
{
	std::vector<int> vec2;
	bedge_elist.insert(bedge_elist.end(), e1); 
	vec2.push_back(n1); vec2.push_back(n2); vec2.push_back(n3); sort(vec2.begin(), vec2.end());
	bedge_nlist.insert(bedge_nlist.end(), vec2); bedge_state.push_back(EDGE_ACTIVE);


}

int NonconformingMesh::FindTrashTriFace(const int el, const int n1, const int n2, const int n3)
{
	std::vector<int> vec2, ed;
	std::vector< std::vector<int> >::iterator it, itn1, ite1;
	std::vector<short>::iterator itisA;	

	vec2.push_back(n1); vec2.push_back(n2); vec2.push_back(n3); sort(vec2.begin(), vec2.end());
	itn1 = findVecVec(&iedge_nlist, vec2);
	for(it = iedge_nlist.begin(), ite1 = iedge_elist.begin(), itisA = iedge_state.begin(); it != itn1; it++, ite1++, itisA++);
	ed = *ite1;
	*itisA = EDGE_GARBAGE;
	return(((ed[0] == el) ? ed[1] : ed[0]));

}

void NonconformingMesh::GarbageCollection()
{
	std::vector<short>::iterator itstate;
	std::vector< std::vector<int> >::iterator itie, itn;
	std::vector<int>::iterator itbe;
	std::vector<int> vec;
	itstate = iedge_state.begin();
	itie = iedge_elist.begin();
	itn = iedge_nlist.begin();

	while(itstate != iedge_state.end())
	{
	  if(*itstate == EDGE_GARBAGE)
	  {
		iedge_elist.erase(itie);
		iedge_nlist.erase(itn);
		iedge_state.erase(itstate);
		itstate = iedge_state.begin();
		itie = iedge_elist.begin();
		itn = iedge_nlist.begin();
	  }	
	  else{
		itstate++; itie++; itn++;
	  }	
		
	} 
       
        itstate = bedge_state.begin();
	itbe = bedge_elist.begin(); 
	itn = bedge_nlist.begin();

	while(itstate != bedge_state.end())
	{
	  vec = *itn;
	  if(*itstate == EDGE_GARBAGE)
	  {
		bedge_elist.erase(itbe);
		bedge_nlist.erase(itn);
		bedge_state.erase(itstate);
	        itstate = bedge_state.begin();
		itbe = bedge_elist.begin(); 
		itn = bedge_nlist.begin();
	  }
	  else{
		itstate++; itbe++; itn++;
	  }	
	} 
	
}

// =========================================================================
// Nonmember functions
// =========================================================================

// Add a component to element matrix 'M', given 'mesh' and 'el'
// Element matrix type is defined by 'mode' (see mesh.h)
// nodal or element coefficients are given by 'coeff'
void DGAddToElMatrix (NonconformingMesh &mesh, int el, CGenericSparseMatrix &M,
    const RVector *coeff, int mode)
{
    int i, j, is, js, nnode;
    
    nnode = mesh.elist[el]->nNode();
    for (i = 0; i < nnode; i++) {
	is =  nnode*el + i;
	for (j = 0; j < nnode; j++) {
	    js = nnode*el + j;
	    double re = 0.0, im = 0.0;
	    switch (mode) {
	    case ASSEMBLE_FF:
		re = mesh.elist[el]->IntFF (i, j);
		break;
	    case ASSEMBLE_DD:
		re = mesh.elist[el]->IntDD (i, j);
		break;
	    case ASSEMBLE_PFF:
		re = mesh.elist[el]->IntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_PDD:
		re = mesh.elist[el]->IntPDD (i, j, *coeff);
		break;
	    case ASSEMBLE_PFF_EL:
		re = mesh.elist[el]->IntFF (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_PDD_EL:
		re = mesh.elist[el]->IntDD (i, j) * (*coeff)[el];
		break;
	    case ASSEMBLE_BNDPFF:
		re = mesh.elist[el]->BndIntPFF (i, j, *coeff);
		break;
	    case ASSEMBLE_BNDPFF_EL:
		re = mesh.elist[el]->BndIntFF (i, j)* (*coeff)[el];
		break;


	   }
	    M.Add (is, js, std::complex<double>(re,im));
	}
    }
}


// Add a component to system matrix 'M', given 'mesh'
// Element matrix type is defined by 'mode' (see mesh.h)
// nodal coefficients are given by 'coeff'
void DGAddToSysMatrix (NonconformingMesh &mesh, RGenericSparseMatrix &M,
    const RVector *coeff, int mode)
{
    int i, j, is, js, el, nnode;
    double entry;
     for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	for (i = 0; i < nnode; i++) {
	    is = nnode*el + i;
	    for (j = 0; j < nnode; j++) {
	        js = nnode*el + j;
	
	        switch (mode) {
		case ASSEMBLE_FF:
		    entry = mesh.elist[el]->IntFF (i, j);
		    break;
		case ASSEMBLE_DD:
		    entry = mesh.elist[el]->IntDD (i, j);
		    break;
		case ASSEMBLE_PFF:
		    entry = mesh.elist[el]->IntPFF (i, j, *coeff);
		    break;
		case ASSEMBLE_PDD:
		    entry = mesh.elist[el]->IntPDD (i, j, *coeff);
		    break;
		case ASSEMBLE_PFF_EL:
		    entry = mesh.elist[el]->IntFF (i, j) * (*coeff)[el];
		    break;
		case ASSEMBLE_PDD_EL:
		    entry = mesh.elist[el]->IntDD (i, j) * (*coeff)[el];
		    break;
		case ASSEMBLE_BNDPFF:
		    entry = mesh.elist[el]->BndIntPFF (i, j, *coeff);
		    break;
		case ASSEMBLE_BNDPFF_EL:
		    entry = mesh.elist[el]->BndIntFF (i, j)* (*coeff)[el];
		    break;

		}
		M.Add (is, js, entry);
	    }
	}
    }
}

void DGAddToSysMatrix (NonconformingMesh &mesh, CGenericSparseMatrix &M,
    const RVector *coeff, int mode)
{
    int el;
    for (el = 0; el < mesh.elen(); el++) {
	DGAddToElMatrix (mesh, el, M, coeff, mode);
    }
}

void DGAddToSysMatrix (NonconformingMesh &mesh, CGenericSparseMatrix &M,
    const double coeff, int mode)
{
    int i, j, is, js, el, nnode;
    for (el = 0; el < mesh.elen(); el++) {
        nnode = mesh.elist[el]->nNode();
	for (i = 0; i < nnode; i++) {
	    is = nnode*el +  i;
	    for (j = 0; j < nnode; j++) {
	        js = nnode*el + j;
		double re = 0.0, im = 0.0;
	        switch (mode) {
		case ASSEMBLE_CFF:
		    re = mesh.elist[el]->IntFF (i, j) * coeff;
		    break;
		case ASSEMBLE_CDD:
		    re = mesh.elist[el]->IntDD (i, j) * coeff;
		    break;
		case ASSEMBLE_iCFF:
		    im = mesh.elist[el]->IntFF (i, j) * coeff;
		    break;
		case ASSEMBLE_iCDD:
		    im = mesh.elist[el]->IntDD (i, j) * coeff;
		    break;
		}
		M.Add (is, js, std::complex<double>(re,im));
	   }
	}
    }  
}

/* Adds interior edge contributions to the system matrix */
void DGAddToSysMatrixInteriorEdgeCont(NonconformingMesh &mesh, CGenericSparseMatrix &M, const RVector *ck, const RVector *rf)
{
    std::vector< std::vector<int> > :: iterator it1, it2;
    std::vector<short>::iterator itstate;
    std::vector<int>  ed1, ed2;
    int e, el, dim, *nds;
    double edge_length;
    
    int etop = mesh.elist[0]->Type();
    if(etop == ELID_TRI3)//then assumes a triangular mesh	
	dim = 2;
    else if(etop == ELID_TET4)//assumes a tetrahedral mesh
	dim = 3;
    else
	xERROR("Unsupported mesh type\n");
    
    RVector edge_coord[3] = {RVector(dim), RVector(dim), RVector(dim)};
    RDenseMatrix eebndFF(dim+1, dim+1), eelbndFF(dim+1, dim+1), elebndFF(dim+1, dim+1), elelbndFF(dim+1, dim+1);
    RDenseMatrix eeFD(dim+1, dim+1), eelFD(dim+1, dim+1), eleFD(dim+1, dim+1), elelFD(dim+1, dim+1);
    RDenseMatrix eeDD(dim+1, dim+1), eelDD(dim+1, dim+1), eleDD(dim+1, dim+1), elelDD(dim+1, dim+1);
    RDenseMatrix eebndFD[3] = {RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1)};
    RDenseMatrix eelbndFD[3] = {RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1)};
    RDenseMatrix elebndFD[3] = {RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1)};
    RDenseMatrix elelbndFD[3] = {RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1), RDenseMatrix(dim+1, dim+1)};

    RVector normal(3);
    RVector ckappa, refind; 
    ckappa = *ck;
    refind = *rf;

    double gamma, beta=0;	
    nds = (int *)malloc(dim*sizeof(int));
    
    for (it1 = mesh.iedge_elist.begin(), it2 = mesh.iedge_nlist.begin(), itstate=mesh.iedge_state.begin(); it1 != mesh.iedge_elist.end(), it2 != mesh.iedge_nlist.end(), itstate!=mesh.iedge_state.end(); it1++, it2++, itstate++) {

	if(*itstate == EDGE_ACTIVE){

	ed1 = *it1; ed2 = *it2;
	e = ed1[0]; el = ed1[1];
	for(int i = 0; i<dim; i++){nds[i] = ed2[i];
		for(int j =0; j< dim; j++)
			edge_coord[i][j] = mesh.nlist[nds[i]][j];
                 }
	//cout << nds[0]<<"  "<<nds[1]<<"  "<<nds[2]<<endl;

	int eside = mesh.elist[e]->IsSide(dim, nds);
	//cout<<"Side: "<<eside<<endl;
	if(eside == -1){
		eside = mesh.elist[el]->IsSide(dim, nds);
		//cout<<"Modified eside: "<<eside<<endl;
		computeEdgeGlobalNormal(mesh, el, eside, normal); normal *= -1;
	}
	else
		computeEdgeGlobalNormal(mesh, e, eside, normal);

	//cout<<"Normal computed .. "<<normal[0]<<"  "<<normal[1]<<"  "<<normal[2]<<endl;
	/*for(int i=0; i<dim; i++)
		for(int j=0; j<dim; j++)
			edge_coord[i][j] = mesh.nlist[mesh.elist[e]->Node[mesh.elist[e]->SideNode(eside, i)]][j];*/
	/*int eside = mesh.elist[e]->IsSide(dim, nds);
	normal =  mesh.elist[e]->LNormal(eside);*/
	
        //cout<<"Edge coords assigned.."<<endl; 
	if(etop == ELID_TRI3)
		edge_length = l2norm(edge_coord[0] - (edge_coord[1]));
	else if (etop == ELID_TET4){ 
		double a, b, c, p;
		a = l2norm(edge_coord[0] - edge_coord[1]); b = l2norm(edge_coord[0] - edge_coord[2]); c = l2norm(edge_coord[1] - edge_coord[2]);
		p = 0.5*(a+b+c); 
		edge_length = sqrt(p*(p-a)*(p-b)*(p-c));
		}
	else 
		xERROR("Unsupported mesh type\n");


	computeFD_FFeel(mesh, e, e, edge_coord, eebndFD, eebndFF);
	computeFD_FFeel(mesh, e, el, edge_coord, eelbndFD, eelbndFF);
	computeFD_FFeel(mesh, el, e, edge_coord, elebndFD, elebndFF);
	computeFD_FFeel(mesh, el, el, edge_coord, elelbndFD, elelbndFF);
       	eeFD.Zero(); eelFD.Zero(); eleFD.Zero(); elelFD.Zero();
	eeDD.Zero(); eelDD.Zero(); 
        /*cout<<"relevant matrices computed "<<e<<"  "<<el<<endl;	
	cout<<eelbndFF.Get(0, 0)<< "  "<<eelbndFF.Get(0, 1)<< "  "<<eelbndFF.Get(0, 2)<< "  "<<eelbndFF.Get(0, 3)<<endl;
	cout<<eelbndFF.Get(1, 0)<< "  "<<eelbndFF.Get(1, 1)<< "  "<<eelbndFF.Get(1, 2)<< "  "<<eelbndFF.Get(1, 3)<<endl;
	cout<<eelbndFF.Get(2, 0)<< "  "<<eelbndFF.Get(2, 1)<< "  "<<eelbndFF.Get(2, 2)<< "  "<<eelbndFF.Get(2, 3)<<endl;
	cout<<eelbndFF.Get(3, 0)<< "  "<<eelbndFF.Get(3, 1)<< "  "<<eelbndFF.Get(3, 2)<< "  "<<eelbndFF.Get(3, 3)<<endl;
        getchar();

        cout<<elebndFF.Get(0, 0)<< "  "<<elebndFF.Get(0, 1)<< "  "<<elebndFF.Get(0, 2)<< "  "<<elebndFF.Get(0, 3)<<endl;
	cout<<elebndFF.Get(1, 0)<< "  "<<elebndFF.Get(1, 1)<< "  "<<elebndFF.Get(1, 2)<< "  "<<elebndFF.Get(1, 3)<<endl;
	cout<<elebndFF.Get(2, 0)<< "  "<<elebndFF.Get(2, 1)<< "  "<<elebndFF.Get(2, 2)<< "  "<<elebndFF.Get(2, 3)<<endl;
	cout<<elebndFF.Get(3, 0)<< "  "<<elebndFF.Get(3, 1)<< "  "<<elebndFF.Get(3, 2)<< "  "<<elebndFF.Get(3, 3)<<endl;
*/
	for(int i =0; i<dim; i++){
	 	for(int j=0; j < dim+1; j++){
			for(int k=0; k< dim+1; k++){ 
				eeFD.Set(j, k, eeFD.Get(j, k)+ normal[i]*eebndFD[i].Get(j, k)); 
				eelFD.Set(j, k, eelFD.Get(j, k)+ normal[i]*eelbndFD[i].Get(j, k));
				eleFD.Set(j, k, eleFD.Get(j, k)+ normal[i]*elebndFD[i].Get(j, k)); 
				elelFD.Set(j, k, elelFD.Get(j, k)+ normal[i]*elelbndFD[i].Get(j, k)); }}}

        double temp=0;
	gamma =  pow(refind[e], 2.0)/pow(refind[el], 2.0);
	//cout << gamma<<endl;
	computeDDeel(mesh, e, e, edge_coord, eeDD, normal);
	computeDDeel(mesh, e, el, edge_coord, eelDD, normal);
	computeDDeel(mesh, el, e, edge_coord, eleDD, normal);
	computeDDeel(mesh, el, el, edge_coord, elelDD, normal);

//	cout<<"normal: "<<normal[0]<<" "<<normal[1]<<"  "<<normal[2]<<endl;
	if(gamma == 1.0){
		for(int i=0; i< dim+1; i++){
			for(int j=0; j< dim+1; j++){
				temp =  0.5*ckappa[e]*(-eeFD.Get(i, j) - eeFD.Get(j, i)) + pow(10, 6)*eebndFF(i, j)/pow(edge_length, 3.0);
				M.Add((dim+1)*e+i, (dim+1)*e+j, temp);
			
				temp =  0.5*ckappa[el]*(elelFD.Get(i, j) + elelFD.Get(j, i)) + pow(10, 6)*elelbndFF(i, j)/pow(edge_length, 3.0);
				M.Add((dim+1)*el+i, (dim+1)*el+j, temp);
		        			
				temp =  0.5*(-ckappa[el]*eelFD.Get(i, j) + ckappa[e]*eleFD.Get(j, i)) -pow(10, 6)*eelbndFF(i, j)/pow(edge_length, 3.0);
		        	M.Add((dim+1)*e+i, (dim+1)*el+j, temp);

		      
				temp =  0.5*(ckappa[e]*eleFD.Get(i, j) - ckappa[el]*eelFD.Get(j, i)) -pow(10, 6)*elebndFF(i, j)/pow(edge_length, 3.0);
				M.Add((dim+1)*el+i, (dim+1)*e+j, temp);
			} 
		}
	}
	else{
		//beta = computeBeta(refind[e], refind[el]);
		
		  beta = 0;
	
		//cout << beta<<endl;
		for(int i=0; i< dim+1; i++){
			for(int j=0; j<dim+1; j++){

				temp =   0.5*ckappa[e]*(-eeFD.Get(i, j) - eeFD.Get(j, i))  - 0.5*beta*ckappa[e]*ckappa[e]*eeDD.Get(j, i) + pow(10, 6)*eebndFF(i, j)/pow(edge_length, 3.0) ;//+  pow(10, 6)*beta*ckappa[e]*eelFD.Get(i, j)/pow(edge_length, 3.0) ;
				M.Add((dim+1)*e+i, (dim+1)*e+j, temp);
			
				temp =   0.5*ckappa[el]*(elelFD.Get(i, j) + elelFD.Get(j, i))  + 0.5*(gamma-1)*ckappa[el]*elelFD.Get(j, i) + pow(10, 6)*gamma*elelbndFF(i, j)/pow(edge_length, 3.0);
				M.Add((dim+1)*el+i, (dim+1)*el+j, temp);
		        			
				temp =   0.5*(-ckappa[el]*eelFD.Get(i, j) + ckappa[e]*eleFD.Get(j, i)) + 0.5*(gamma-1)*ckappa[e]*eleFD.Get(j, i) -pow(10, 6)*gamma*eelbndFF(i, j)/pow(edge_length, 3.0);
		        	M.Add((dim+1)*e+i, (dim+1)*el+j, temp);

		      
				temp = 0.5*(ckappa[e]*eleFD.Get(i, j) - ckappa[el]*eelFD.Get(j, i))  - 0.5*beta*ckappa[e]*ckappa[el]*eleDD.Get(j, i)  - pow(10, 6)*elebndFF(i, j)/pow(edge_length, 3.0);// +  pow(10, 6)*beta*ckappa[e]*eleFD.Get(i, j)/pow(edge_length, 3.0) ;
				M.Add((dim+1)*el+i, (dim+1)*e+j, temp);

			} 
		}

	}
    }}
    free(nds);
}

void DGAddToSysMatrixBoundaryEdgeCont(NonconformingMesh &mesh, CGenericSparseMatrix &M, const RVector *coeff)
{
    std::vector<int>::iterator it1;
    std::vector< std::vector<int> > :: iterator it2;
    std::vector<short>::iterator itstate;
    std::vector<int>  ed2;
    int e, dim, *nds;
    double edge_length;
    
    int etop = mesh.elist[0]->Type();
    if(etop == ELID_TRI3)//then assumes a triangular mesh	
	dim = 2;
    else if(etop == ELID_TET4)//assumes a tetrahedral mesh
	dim = 3;
    else
	xERROR("Unsupported mesh type\n");
    
    RVector edge_coord[3] = {RVector(dim), RVector(dim), RVector(dim)};
    RDenseMatrix eebndFF(dim+1, dim+1);

    RVector c2a; 
    c2a = *coeff;
	
    nds = (int *)malloc(dim*sizeof(int));
    
    for (it1 = mesh.bedge_elist.begin(), it2 = mesh.bedge_nlist.begin(), itstate = mesh.bedge_state.begin(); it1 != mesh.bedge_elist.end(), it2 != mesh.bedge_nlist.end(), itstate != mesh.bedge_state.end(); it1++, it2++, itstate++) {
	if(*itstate == EDGE_ACTIVE){

		e = *it1; ed2 = *it2;
		//cout << e<<"  ";
		for(int i = 0; i < dim; i++){ nds[i] = ed2[i]; //cout<<nds[i]<<"  ";
			for(int j =0; j< dim; j++)
				edge_coord[i][j] = (double)mesh.nlist[nds[i]][j];
                 }
		//cout<<endl;
		computeFFeel(mesh, e, e, edge_coord, eebndFF);
	/*cout<<eebndFF.Get(0, 0)<< "  "<<eebndFF.Get(0, 1)<< "  "<<eebndFF.Get(0, 2)<< "  "<<eebndFF.Get(0, 3)<<endl;
	cout<<eebndFF.Get(1, 0)<< "  "<<eebndFF.Get(1, 1)<< "  "<<eebndFF.Get(1, 2)<< "  "<<eebndFF.Get(1, 3)<<endl;
	cout<<eebndFF.Get(2, 0)<< "  "<<eebndFF.Get(2, 1)<< "  "<<eebndFF.Get(2, 2)<< "  "<<eebndFF.Get(2, 3)<<endl;
	cout<<eebndFF.Get(3, 0)<< "  "<<eebndFF.Get(3, 1)<< "  "<<eebndFF.Get(3, 2)<< "  "<<eebndFF.Get(3, 3)<<endl;
*/
		for(int i=0; i< dim+1; i++)
			for(int j=0; j<dim+1; j++)
				M.Add((dim+1)*e+i, (dim+1)*e+j, eebndFF.Get(i, j)*c2a[e]);
	}
    }
    free(nds);
}

void computeFFeel(NonconformingMesh &mesh, int e, int el, RVector *edge_coord, RDenseMatrix &eelbndFF)
{
  int etp = mesh.elist[e]->Type(), eltp = mesh.elist[el]->Type();
  double x, y, z,temp1, temp2;
  if((etp == eltp) && (etp == ELID_TRI3)) // it is a triangle3
  {
	const double edge_length = l2norm(edge_coord[0] -edge_coord[1]);
	RVector fun1(3);
	RVector fun2(3);
	int nt = 10;
	const double t[] = {-0.973906528517172,0.973906528517172, -0.865063366688985, 0.865063366688985, -0.679409568299024, 0.679409568299024, -0.433395394129247, 0.433395394129247, -0.148874338981631, 0.148874338981631};
	const double W[] = {0.066671344308688, 0.066671344308688, 0.149451349150581, 0.149451349150581, 0.219086362515982, 0.219086362515982, 0.269266719309996, 0.269266719309996, 0.295524224714753, 0.295524224714753};
	eelbndFF.Zero();
 
	for(int i =0; i < nt; i++){
		x = 0.5*(edge_coord[1][0] - edge_coord[0][0])*t[i] + 0.5*(edge_coord[1][0] + edge_coord[0][0]);
		y = 0.5*(edge_coord[1][1] - edge_coord[0][1])*t[i] + 0.5*(edge_coord[1][1] + edge_coord[0][1]);
		
		Point2D pt(x, y);
 		fun1 = mesh.elist[e]->GlobalShapeF(mesh.nlist, pt);
		fun2 = mesh.elist[el]->GlobalShapeF(mesh.nlist, pt);		
			
		for(int j=0; j<3; j++){
		    for(int k=0; k<3; k++){
			temp2 =  eelbndFF.Get(j, k) + W[i]*fun1[j]*fun2[k];
			eelbndFF.Set(j, k, temp2);
		   }
		 }
	}

	eelbndFF = eelbndFF*0.5*edge_length;
   }
   else if((etp == eltp) && (etp == ELID_TET4))
   {
	double a, b, c, p;
	a = l2norm(edge_coord[0] -edge_coord[1]); b = l2norm(edge_coord[0] - edge_coord[2]); c = l2norm(edge_coord[1] - edge_coord[2]);
	p = 0.5*(a+b+c); 
	const double edge_length = sqrt(p*(p-a)*(p-b)*(p-c));
	RVector fun1(4);
	RVector fun2(4);

	int nt = 7;
	const double alpha1 = 0.0597158717, alpha2 = 0.7974269853, beta1 = 0.4701420641, beta2 = 0.1012865073;
	const double w1 = 0.225, w2 = 0.1323941527, w3 = 0.1259391805;
	const double t[][7] = {1/3.0, alpha1, beta1, beta1, alpha2, beta2, beta2,
			      1/3.0, beta1, alpha1, beta1, beta2, alpha2, beta2,
			      1/3.0, beta1, beta1, alpha1, beta2, beta2, alpha2};
	const double W[] = {w1, w2, w2, w2, w3, w3, w3};
	eelbndFF.Zero();
		
	for(int i =0; i < nt; i++){
		x = edge_coord[0][0]*t[0][i] + edge_coord[1][0]*t[1][i] + edge_coord[2][0]*t[2][i];
		y = edge_coord[0][1]*t[0][i] + edge_coord[1][1]*t[1][i] + edge_coord[2][1]*t[2][i];
		z = edge_coord[0][2]*t[0][i] + edge_coord[1][2]*t[1][i] + edge_coord[2][2]*t[2][i];

		Point3D pt(x, y, z);
 		fun1 = mesh.elist[e]->GlobalShapeF(mesh.nlist, pt);
		fun2 = mesh.elist[el]->GlobalShapeF(mesh.nlist, pt);		
		for(int j=0; j<4; j++){
			for(int k=0; k<4; k++){
				temp2 =  eelbndFF.Get(j, k) + W[i]*fun1[j]*fun2[k];
				eelbndFF.Set(j, k, temp2);}}
	 }
	 eelbndFF = eelbndFF*edge_length;
   }
	
}

void computeDDeel(NonconformingMesh &mesh, int e, int el, RVector *edge_coord, RDenseMatrix &eelbndDD, RVector &normal)
{
  int etp = mesh.elist[e]->Type(), eltp = mesh.elist[el]->Type();
  double x, y, z,temp1, temp2;
  if((etp == eltp) && (etp == ELID_TRI3)) // it is a triangle3
  {
	const double edge_length = l2norm(edge_coord[0] -edge_coord[1]);
	RVector fun1(3);
	RVector fun2(3);
	RDenseMatrix mat1(2, 3), mat2(2, 3);

	int nt = 10;
	const double t[] = {-0.973906528517172,0.973906528517172, -0.865063366688985, 0.865063366688985, -0.679409568299024, 0.679409568299024, -0.433395394129247, 0.433395394129247, -0.148874338981631, 0.148874338981631};
	const double W[] = {0.066671344308688, 0.066671344308688, 0.149451349150581, 0.149451349150581, 0.219086362515982, 0.219086362515982, 0.269266719309996, 0.269266719309996, 0.295524224714753, 0.295524224714753};
	eelbndDD.Zero();
 
	for(int i =0; i < nt; i++){
		x = 0.5*(edge_coord[1][0] - edge_coord[0][0])*t[i] + 0.5*(edge_coord[1][0] + edge_coord[0][0]);
		y = 0.5*(edge_coord[1][1] - edge_coord[0][1])*t[i] + 0.5*(edge_coord[1][1] + edge_coord[0][1]);
		
		Point2D pt(x, y);
		mat1 = mesh.elist[e]->GlobalShapeD(mesh.nlist, pt);
		mat2 = mesh.elist[el]->GlobalShapeD(mesh.nlist, pt);
		for(int j=0; j<3; j++)
		{
		  fun1[j] = mat1.Get(0, j)*normal[0] + mat1.Get(1, j)*normal[1];
	 	  fun2[j] = mat2.Get(0, j)*normal[0] + mat2.Get(1, j)*normal[1];
		}
		
			
		for(int j=0; j<3; j++){
		    for(int k=0; k<3; k++){
			temp2 =  eelbndDD.Get(j, k) + W[i]*fun1[j]*fun2[k];
			eelbndDD.Set(j, k, temp2);
		   }
		 }
	}

	eelbndDD = eelbndDD*0.5*edge_length;
   }
   else if((etp == eltp) && (etp == ELID_TET4))
   {
	double a, b, c, p;
	a = l2norm(edge_coord[0] -edge_coord[1]); b = l2norm(edge_coord[0] - edge_coord[2]); c = l2norm(edge_coord[1] - edge_coord[2]);
	p = 0.5*(a+b+c); 
	const double edge_length = sqrt(p*(p-a)*(p-b)*(p-c));
	RVector fun1(4);
	RVector fun2(4);
	RDenseMatrix mat1(3, 4), mat2(3, 4);


	int nt = 7;
	const double alpha1 = 0.0597158717, alpha2 = 0.7974269853, beta1 = 0.4701420641, beta2 = 0.1012865073;
	const double w1 = 0.225, w2 = 0.1323941527, w3 = 0.1259391805;
	const double t[][7] = {1/3.0, alpha1, beta1, beta1, alpha2, beta2, beta2,
			      1/3.0, beta1, alpha1, beta1, beta2, alpha2, beta2,
			      1/3.0, beta1, beta1, alpha1, beta2, beta2, alpha2};
	const double W[] = {w1, w2, w2, w2, w3, w3, w3};
	eelbndDD.Zero();
		
	for(int i =0; i < nt; i++){
		x = edge_coord[0][0]*t[0][i] + edge_coord[1][0]*t[1][i] + edge_coord[2][0]*t[2][i];
		y = edge_coord[0][1]*t[0][i] + edge_coord[1][1]*t[1][i] + edge_coord[2][1]*t[2][i];
		z = edge_coord[0][2]*t[0][i] + edge_coord[1][2]*t[1][i] + edge_coord[2][2]*t[2][i];

		Point3D pt(x, y, z);
		mat1 = mesh.elist[e]->GlobalShapeD(mesh.nlist, pt);
		mat2 = mesh.elist[el]->GlobalShapeD(mesh.nlist, pt);
		for(int j=0; j<4; j++)
		{
			fun1[j] = mat1.Get(0, j)*normal[0] +  mat1.Get(1, j)*normal[1] + mat1.Get(2, j)*normal[2];
			fun2[j] = mat2.Get(0, j)*normal[0] +  mat2.Get(1, j)*normal[1] + mat2.Get(2, j)*normal[2];

		}		

		for(int j=0; j<4; j++){
			for(int k=0; k<4; k++){
				temp2 =  eelbndDD.Get(j, k) + W[i]*fun1[j]*fun2[k];
				eelbndDD.Set(j, k, temp2);}}
	 }
	 eelbndDD = eelbndDD*edge_length;
   }
	
}


void computeFD_FFeel(NonconformingMesh &mesh, int e, int el, RVector *edge_coord, RDenseMatrix *eelbndFD, RDenseMatrix &eelbndFF)
{
  int etp = mesh.elist[e]->Type(), eltp = mesh.elist[el]->Type();
  double x, y, z,temp1, temp2;
  if((etp == eltp) && (etp == ELID_TRI3)) // it is a triangle3
  {
	const double edge_length = l2norm(edge_coord[0] -edge_coord[1]);
	RVector fun1(3);
	RDenseMatrix mat(2, 3);
	RVector fun2(3);
	int nt = 10;
	const double t[] = {-0.973906528517172,0.973906528517172, -0.865063366688985, 0.865063366688985, -0.679409568299024, 0.679409568299024, -0.433395394129247, 0.433395394129247, -0.148874338981631, 0.148874338981631};
	const double W[] = {0.066671344308688, 0.066671344308688, 0.149451349150581, 0.149451349150581, 0.219086362515982, 0.219086362515982, 0.269266719309996, 0.269266719309996, 0.295524224714753, 0.295524224714753};
	eelbndFD[0].Zero(); eelbndFD[1].Zero();
	eelbndFF.Zero();
 
	for(int i =0; i < nt; i++){
		x = 0.5*(edge_coord[1][0] - edge_coord[0][0])*t[i] + 0.5*(edge_coord[1][0] + edge_coord[0][0]);
		y = 0.5*(edge_coord[1][1] - edge_coord[0][1])*t[i] + 0.5*(edge_coord[1][1] + edge_coord[0][1]);
		
		Point2D pt(x, y);
 		fun1 = mesh.elist[e]->GlobalShapeF(mesh.nlist, pt);
		fun2 = mesh.elist[el]->GlobalShapeF(mesh.nlist, pt);		
		mat = mesh.elist[el]->GlobalShapeD(mesh.nlist, pt);
		
		for(int j=0; j<3; j++){
		    for(int k=0; k<3; k++){
			for(int dim=0; dim < 2; dim++){ 
			temp1 =  eelbndFD[dim].Get(j, k) + W[i]*fun1[j]*mat.Get(dim, k);
			eelbndFD[dim].Set(j, k, temp1);}
			
			temp2 =  eelbndFF.Get(j, k) + W[i]*fun1[j]*fun2[k];
			eelbndFF.Set(j, k, temp2);
		   }
		 }
	}

	eelbndFF = eelbndFF*0.5*edge_length;
	eelbndFD[0] = eelbndFD[0]*0.5*edge_length;
	eelbndFD[1] = eelbndFD[1]*0.5*edge_length;

        
   }
   else if((etp == eltp) && (etp == ELID_TET4))
   {
	double a, b, c, p;
	a = l2norm(edge_coord[0] -edge_coord[1]); b = l2norm(edge_coord[0] - edge_coord[2]); c = l2norm(edge_coord[1] - edge_coord[2]);
	p = 0.5*(a+b+c); 
	const double edge_length = sqrt(p*(p-a)*(p-b)*(p-c));
	RVector fun1(4);
	RDenseMatrix mat(3, 4);
	RVector fun2(4);

	int nt = 7;
	const double alpha1 = 0.0597158717, alpha2 = 0.7974269853, beta1 = 0.4701420641, beta2 = 0.1012865073;
	const double w1 = 0.225, w2 = 0.1323941527, w3 = 0.1259391805;
	const double t[][7] = {1/3.0, alpha1, beta1, beta1, alpha2, beta2, beta2,
			      1/3.0, beta1, alpha1, beta1, beta2, alpha2, beta2,
			      1/3.0, beta1, beta1, alpha1, beta2, beta2, alpha2};
	const double W[] = {w1, w2, w2, w2, w3, w3, w3};
	eelbndFD[0].Zero(); eelbndFD[1].Zero(); eelbndFD[2].Zero();
	eelbndFF.Zero();
		
	for(int i =0; i < nt; i++){
		x = edge_coord[0][0]*t[0][i] + edge_coord[1][0]*t[1][i] + edge_coord[2][0]*t[2][i];
		y = edge_coord[0][1]*t[0][i] + edge_coord[1][1]*t[1][i] + edge_coord[2][1]*t[2][i];
		z = edge_coord[0][2]*t[0][i] + edge_coord[1][2]*t[1][i] + edge_coord[2][2]*t[2][i];

		Point3D pt(x, y, z);
 		fun1 = mesh.elist[e]->GlobalShapeF(mesh.nlist, pt);
		fun2 = mesh.elist[el]->GlobalShapeF(mesh.nlist, pt);		
		mat = mesh.elist[el]->GlobalShapeD(mesh.nlist, pt);
		for(int j=0; j<4; j++){
			for(int k=0; k<4; k++){
				for(int dim =0 ; dim < 3; dim++){
					temp1 =  eelbndFD[dim].Get(j, k) + W[i]*fun1[j]*mat.Get(dim, k);
					eelbndFD[dim].Set(j, k, temp1);}
				
				temp2 =  eelbndFF.Get(j, k) + W[i]*fun1[j]*fun2[k];
				eelbndFF.Set(j, k, temp2);}}
	 }
	 eelbndFF = eelbndFF*edge_length;
	 eelbndFD[0] = eelbndFD[0]*edge_length;
	 eelbndFD[1] = eelbndFD[1]*edge_length;
	 eelbndFD[2] = eelbndFD[2]*edge_length;

        }
	
}

bool computeEdgeGlobalNormal(NonconformingMesh &mesh, int e, int eside, RVector &normal)
{
	int dim; 
	bool success_state;
	double temp;
        int etop = mesh.elist[e]->Type();
        if(etop == ELID_TRI3)//then assumes a triangular mesh	
	   dim = 2;
        else if(etop == ELID_TET4)//assumes a tetrahedral mesh
	  dim = 3;
        else
	xERROR("Unsupported mesh type\n");
	
//	cout<<nds[0]<<"  "<<nds[1]<<"  "<<nds[2]<<endl;
  	
//	cout<<eside<<endl;
        if(eside == -1){
	   return false;	
	}

	if(etop == ELID_TRI3){
		RDenseMatrix vertices =  mesh.elist[e]->Elgeom(mesh.nlist);
		double n1[3];
		double v1[3] = {vertices.Get(1, 0)-vertices(0, 0), vertices(1, 1)-vertices(0, 1), 0};
		double v2[3] = {vertices.Get(2, 0)-vertices(1, 0), vertices.Get(2, 1)-vertices.Get(1, 1), 0};
		double v3[3] = {vertices.Get(0, 0)-vertices(2, 0), vertices.Get(0, 1)-vertices.Get(2, 1), 0};
		n1[0] = v1[1]*v2[2] - v1[2]*v2[1]; n1[1] =v1[2]*v2[0] - v1[0]*v2[2]; n1[2] = v1[0]*v2[1] - v1[1]*v2[0];
		switch(eside){
			case 0:
				normal[0] = n1[1]*v1[2] - n1[2]*v1[1]; normal[1] =n1[2]*v1[0] - n1[0]*v1[2]; normal[2] = n1[0]*v1[1] - n1[1]*v1[0];
				break;
			case 1:
				normal[0] = n1[1]*v2[2] - n1[2]*v2[1]; normal[1] = n1[2]*v2[0] - n1[0]*v2[2]; normal[2] = n1[0]*v2[1] - n1[1]*v2[0];
				break;
			case 2:
				normal[0] = n1[1]*v3[2] - n1[2]*v3[1]; normal[1] = n1[2]*v3[0] - n1[0]*v3[2]; normal[2] = n1[0]*v3[1] - n1[1]*v3[0];
				break;
			default:
				return false;
		}
	}
	else{ //its a tetrahedra
		RDenseMatrix vertices(3, 3);

		int gnodes[3];
		gnodes[0] = mesh.elist[e]->Node[mesh.elist[e]->SideNode(eside, 0)];
		gnodes[1] = mesh.elist[e]->Node[mesh.elist[e]->SideNode(eside, 1)];
		gnodes[2] = mesh.elist[e]->Node[mesh.elist[e]->SideNode(eside, 2)];
		//cout<<gnodes[0]<<"  "<<gnodes[1]<<"  "<<gnodes[2]<<endl;
		vertices.Set(0, 0, mesh.nlist[gnodes[0]][0]); vertices.Set(0, 1, mesh.nlist[gnodes[0]][1]); vertices.Set(0, 2, mesh.nlist[gnodes[0]][2]); 
		vertices.Set(1, 0, mesh.nlist[gnodes[1]][0]); vertices.Set(1, 1, mesh.nlist[gnodes[1]][1]); vertices.Set(1, 2, mesh.nlist[gnodes[1]][2]); 
		vertices.Set(2, 0, mesh.nlist[gnodes[2]][0]); vertices.Set(2, 1, mesh.nlist[gnodes[2]][1]); vertices.Set(2, 2, mesh.nlist[gnodes[2]][2]); 
		
	/*	cout<<vertices.Get(0, 0) << "  "<<vertices.Get(0, 1)<< "  "<<vertices.Get(0, 2)<<endl;
		cout<<vertices.Get(1, 0) << "  "<<vertices.Get(1, 1)<< "  "<<vertices.Get(1, 2)<<endl;
		cout<<vertices.Get(2, 0) << "  "<<vertices.Get(2, 1)<< "  "<<vertices.Get(2, 2)<<endl;
*/
		double v1[3] = {vertices.Get(2, 0)-vertices.Get(0, 0), vertices.Get(2, 1)-vertices.Get(0, 1), vertices.Get(2, 2)-vertices.Get(0, 2)};
		double v2[3] = {vertices.Get(1, 0)-vertices.Get(0, 0), vertices.Get(1, 1)-vertices.Get(0, 1), vertices.Get(1, 2)-vertices.Get(0, 2)};
	
		normal[0] = v1[1]*v2[2] - v1[2]*v2[1]; normal[1] =v1[2]*v2[0] - v1[0]*v2[2]; normal[2] = v1[0]*v2[1] - v1[1]*v2[0];
	}

	temp =  sqrt((normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]));
	normal[0] /= temp; normal[1] /= temp; normal[2] /= temp;
	/*cout<<normal[0]<<"  "<<normal[1]<<"  "<<normal[2]<<endl;
	getchar();*/
	return true; 
}

std::vector< std::vector<int> >::iterator findVecVec(std::vector< std::vector<int> > *source, std::vector<int> target)
{
	std::vector< std::vector<int> >::iterator it;
	std::vector<int> vec;
	for(it = source->begin(); it < source->end(); it++)
	{
		vec = *it;
		if(*it == target){
			return(it);}
	}
	it=source->end();
	return(it);
}

double computeBeta(double n1, double n2)
{
	double newN1, newN2, beta; 
	newN1 = n1 > n2 ? n1 : n2;
	newN2 = n1 > n2 ? n2 : n1;
	
	double n, p, p1;

        n = newN1/newN2;
	p = sqrt((n-1)/(n+1)); p1 = p - 0.64;
	if(n == 1)
	        return(0); 
	else if(n > 1 && n < 1.82){
		beta = 25.6*pow(p, 3.0)*(1-(45.0/32.0)*p + (83.0/21.0)*pow(p, 2.0));
		return(beta);}
	else if(n>=1.82 && n <= 3.73){
	 	beta = 14.3*(1 + 9.35*p1 + 49.1*pow(p1, 2.0) + 327.0*pow(p1, 3.0) + 1800.0*pow(p1, 4.0));
		return(beta);}
	else{
		cout<<"beta could not be computed properly";
		return(0);
		}

}

