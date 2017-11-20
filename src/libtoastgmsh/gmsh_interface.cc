#define TOASTGMSHLIB_IMPLEMENTATION

#include "gmsh_interface.h"
#include "MElement.h"
#include "MTriangle.h"
#include "MTetrahedron.h"
#include "discreteRegion.h"

GmshInterface::GmshInterface()
{
    GmshInitialize(0, NULL);
}

GmshInterface::~GmshInterface()
{
  GmshFinalize();
}

GModel *GmshInterface::Create_Model()
{
    return new GModel();
}

GModel *GmshInterface::Load_Model(const char *fname)
{
    GModel *m = new GModel();
    if (!m->readGEO(fname))
	xERROR("Geometry file not found");
    return m;
}

void GmshInterface::LoadModel(GModel *model, const char *fname)
{
    if (!model->readGEO(fname))
	xERROR("Geometry file not found");
}

void GmshInterface::MeshModel(GModel *model, int dim,
    double minlength, double maxlength)
{
    if (maxlength < 0.0)
	maxlength = minlength;
    
    GmshSetOption("Mesh", "CharacteristicLengthExtendFromBoundary", 1.);
    GmshSetOption("Mesh", "OldRefinement", 1.);
    if (minlength > 0.0)
	GmshSetOption("Mesh", "CharacteristicLengthMin", minlength);
    if (maxlength > 0.0)
	GmshSetOption("Mesh", "CharacteristicLengthMax", maxlength);
    GmshSetOption("Mesh", "Optimize", 0.); // not yet: need boundary!

    model->mesh(dim);
}

void GmshInterface::RemeshModel(GModel *model, double minlength,
    double maxlength)
{
    if (maxlength < 0.0) maxlength = minlength;

    GmshSetOption("Mesh", "CharacteristicLengthExtendFromBoundary", 1.);
    GmshSetOption("Mesh", "OldRefinement", 1.);
    GmshSetOption("Mesh", "CharacteristicLengthMin", minlength);
    GmshSetOption("Mesh", "CharacteristicLengthMax", maxlength);
    GmshSetOption("Mesh", "Optimize", 0.); // not yet: need boundary!

    std::set<GFace*> faceset;
    for (GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it) {
	std::cerr << "Found face" << std::endl;
    	GFace *face = *it;
	faceset.insert(face);
    }
    model->classifyFaces(faceset);
    
    model->createTopologyFromMeshNew();
    model->makeDiscreteRegionsSimplyConnected();
    
    for (GModel::riter it = model->firstRegion(); it != model->lastRegion();
	 ++it) {
	std::cerr << "Found region" << std::endl;
	discreteRegion *r = dynamic_cast<discreteRegion*>(*it);
	if (r) {
	    std::cerr << "Found discrete region" << std::endl;
	    r->remesh();
	}
    }
}

void GmshInterface::WriteMesh(GModel *model, const char *fname)
{
    model->writeMSH(fname);
}

GModel *GmshInterface::Load_gMesh(const char *fname)
{
    GModel *m = new GModel();
    m->readMSH(fname);
    return m;
}

GModel *GmshInterface::gMesh_From_Model(GModel *model)
{
    model->mesh(3);
    return model;
}

QMMesh *GmshInterface::Mesh_From_gMesh(GModel *model)
{
    int i, j, k;

    QMMesh *mesh = new QMMesh;
    NodeList &nlist = mesh->nlist;
    ElementList &elist = mesh->elist;
    
    // retrieve info from gmsh
    int numVertices = model->indexMeshVertices(true, 0, false);

    std::cerr << "Vertices: " << numVertices << std::endl;

    std::vector<GEntity*> entities;
    model->getEntities(entities);

    // sanity check
    int nv = 0;
    for (i = 0; i < entities.size(); i++)
	nv += entities[i]->mesh_vertices.size();
    if (nv != numVertices)
	std::cerr << "Inconsistent vertex count " << nv << " != "
		  << numVertices << std::endl;
    
    int numElements = 0;
    for (unsigned int i = 0; i < entities.size(); i++)
	numElements += entities[i]->getNumMeshElements();

    std::cerr << "Elements: " << numElements << std::endl;
    
    nlist.New(numVertices);
    std::vector<int> vtxId(numVertices);
    int numMax = 0;
    for (i = 0; i < numVertices; i++)
	nlist[i].New(3);

    std::cerr << "Node list allocated." << std::endl;
    
    for (i = k = 0; i < entities.size(); i++)
	for (j = 0; j < entities[i]->mesh_vertices.size(); j++) {
	    nlist[k][0] = entities[i]->mesh_vertices[j]->x();
    	    nlist[k][1] = entities[i]->mesh_vertices[j]->y();
	    nlist[k][2] = entities[i]->mesh_vertices[j]->z();
	    int num = entities[i]->mesh_vertices[j]->getNum();
	    vtxId[k] = num;
	    if (num > numMax) numMax = num;
	    k++;
	}

    std::cerr << "Vertices initialised." << std::endl;
    
    double minz = nlist[0][2];
    double maxz = nlist[0][2];
    for (i = 1; i < numVertices; i++) {
	if (nlist[i][2] > maxz)
	    maxz = nlist[i][2];
	else if (nlist[i][2] < minz)
	    minz = nlist[i][2];
    }
    bool is2D = (maxz-minz) < 1e-6;
    if (is2D) { // remove z-components from node coordinates
	for (i = 0; i < numVertices; i++) {
	    double x = nlist[i][0];
	    double y = nlist[i][1];
	    nlist[i].New(2);
	    nlist[i][0] = x;
	    nlist[i][1] = y;
	}
    }
    
    // Inverse vertex id locator
    std::vector<int> r_vtxId(numMax+1);
    for (i = 0; i < numVertices; i++)
	r_vtxId[vtxId[i]] = i;

    int nregion = 0, idx = 0;
    int ntets = 0, ntri = 0;
    for (GModel::riter it = model->firstRegion(); it != model->lastRegion();
	 ++it) {
	ntets += (*it)->tetrahedra.size();
    }
    for (GModel::fiter it = model->firstFace(); it != model->lastFace(); ++it) {
	ntri  += (*it)->triangles.size();
    }

    if (ntets) { // set up tets for volume mesh
	
	elist.New(ntets);
	for (GModel::riter it = model->firstRegion(); it != model->lastRegion();
	     ++it) {
	    for (i = 0; i < (*it)->tetrahedra.size(); i++) {
		elist[idx] = new Tetrahedron4;
		std::vector<int> vtx;
		(*it)->tetrahedra[i]->getVerticesIdForMSH(vtx);
		if (vtx.size() < 4)
		    xERROR("Inconsistent vertex list");
		for (j = 0; j < 4; j++) {
		    elist[idx]->Node[j] = r_vtxId[vtx[j]];
		}
		elist[idx]->SetRegion(nregion);
		idx++;
	    }
	    nregion++;
	}
	std::cerr << "Tetrahedra initialised." << std::endl;
	
    } else if (ntri) { // setup up triangles for surface mesh

	elist.New(ntri);
	for (GModel::fiter it = model->firstFace(); it != model->lastFace();
	     ++it) {
	    for (i = 0; i < (*it)->triangles.size(); i++) {
		if (is2D) elist[idx] = new Triangle3;
		else      elist[idx] = new Triangle3D3;
		std::vector<int> vtx;
		(*it)->triangles[i]->getVerticesIdForMSH(vtx);
		for (j = 0; j < 3; j++)
		    elist[idx]->Node[j] = r_vtxId[vtx[j]];
		idx++;
	    }
	}
	std::cerr << "Triangles initialised." << std::endl;

    } else {
	xERROR("No elements found in mesh definition.");
    }
    
    mesh->Setup();
    
    std::cerr << "Mesh is set up." << std::endl;
    
    return mesh;
}
