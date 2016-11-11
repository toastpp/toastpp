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



RayProjector::RayProjector (Camera *cam, QMMesh *mesh)
    : Projector (cam, mesh)
{
    constructMatrix();
}

void RayProjector::constructMatrix()
{
    double ps = camera->getPixelSize();
    int    iw = camera->getImageWidth();
    int    ih = camera->getImageHeight();
    double w = iw*ps;
    double h = ih*ps;

    RVector p0 = camera->getPos();
    RVector up = camera->getUp();
    RVector vd = camera->getViewDir();
    RVector rt = cross (up, vd);

    PinholeCamera *ph_camera = dynamic_cast<PinholeCamera*>(camera);
    if (ph_camera) {

	// this is for an (upright) "perspective" camera
	RVector ph = p0;
	p0 = ph + vd*ph_camera->getFocalLength();

	// this is for an actual (inverted) pinhole camera
	//RVector ph = p0 + vd*ph_camera->getFocalLength();

	int i, j, k, r, el;
	for (j = r = 0; j < ih; j++) {
	    for (i = 0; i < iw; i++) {
		RVector p = p0 + rt*((i+0.5)*ps-w*0.5) +
		    up*((j+0.5)*ps-h*0.5);
		Point s = FEMMesh->BndIntersect(p, ph, &el);
		if (el >= 0) {
		    RVector fun = FEMMesh->elist[el]->GlobalShapeF(FEMMesh->nlist, s);
		    RVector row(FEMMesh->nlen());
		    for (k = 0; k < fun.Dim(); k++)
			row[FEMMesh->elist[el]->Node[k]] = fun[k];
		    toImage.SetRow(r, row);
		}
		r++;
	    }
	}
	return;
    }

    OrthoCamera *ot_camera = dynamic_cast<OrthoCamera*>(camera);
    if (ot_camera) {
	int i, j, k, r, el;
	for (j = r = 0; j < ih; j++) {
	    for (i = 0; i < iw; i++) {
		RVector p = p0 + rt*((i+0.5)*ps-w*0.5) +
		    up*((j+0.5)*ps-h*0.5);
		Point s = FEMMesh->BndIntersect(p, p+vd, &el);
		if (el >= 0) {
		    RVector fun = FEMMesh->elist[el]->GlobalShapeF(FEMMesh->nlist, s);
		    RVector row(FEMMesh->nlen());
		    for (k = 0; k < fun.Dim(); k++)
			row[FEMMesh->elist[el]->Node[k]] = fun[k];
		    toImage.SetRow(r, row);
		}
		r++;
	    }
	}
	return;
    }
    
    xERROR("RayProjector: invalid camera type.");
}
