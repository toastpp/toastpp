#define STOASTLIB_IMPLEMENTATION
#include "stoastlib.h"

#include "projector.h"

void Projector::getBoundaryNodes()    
{
    // Store boundary list
    //FEMMesh.MarkBoundary();
    bndNodes = new int[nBndNodes];
    for (int i=0, j=0; i<FEMMesh->nlen(); ++i)
    {
	if (FEMMesh->nlist[i].isBnd())
	{
	    bndNodes[j] = i;
	    ++j;
	}
    }
}

Projector::Projector(Camera * cam, QMMesh * mesh)
    :	camera(cam), FEMMesh(mesh), 
	nBndNodes(mesh->nbnd()), nNodes(mesh->nlen()), 
        w(cam->getImageWidth()), h(cam->getImageHeight()), nImagePts(w*h),
	toImage(nImagePts, nNodes), toField(nNodes, nImagePts)//, toData(nNodes, nNodes)
{
    getBoundaryNodes();
}

void Projector::init(Camera * cam, QMMesh * mesh)
{
    camera=cam; FEMMesh=mesh;
    nBndNodes = mesh->nbnd(); nNodes = mesh->nlen(); 
    w = cam->getImageWidth(); h = cam->getImageHeight(); nImagePts = w*h;
    toImage.New(nImagePts, nNodes); toField.New(nNodes, nImagePts);// toData.New(nNodes, nNodes);

    getBoundaryNodes();    
}
