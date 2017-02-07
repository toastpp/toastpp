#ifndef PROJECTOR_H
#define PROJECTOR_H


#include "stoastlib.h"
#include "slu_zdefs.h"
//#include "supermatrix.h"
#include "camera.h"

//#if defined(_WIN32)
//#include <windows.h>
//#endif
//#include "GL/gl.h"
//#include "GL/osmesa.h"

class STOASTLIB Projector  
{
    public:
	Projector() {}
	Projector(Camera * cam, QMMesh * mesh);
	virtual void init(Camera * cam, QMMesh * mesh);

	// Masked
	void projectFieldToImage(const RVector & field, RVector & image, RCompRowMatrix & mask)
	{
	    //RVector fld(field.Dim(), 1.0);
	    image = mask * toImage * field;
	}
	void projectImageToField(const RVector & image, RVector & field, RCompRowMatrix & maskT)
	{
	    //RVector img(image.Dim(), 1.0);
	    field = toField * (maskT * image);
	}

	// Unmasked
	void projectFieldToImage(const RVector & field, RVector & image)
	{
	    //RVector fld(field.Dim(), 1.0);
	    image = toImage * field;
	}
	void projectImageToField(const RVector & image, RVector & field)
	{
	    //RVector img(image.Dim(), 1.0);
	    field = toField * image;
	}
	void mapToData(const RVector & field, RVector & data)
	{
	    //data = toData * field;
	}
	void getBoundaryNodes();
	int getImageWidth() { return w; }
	int getImageHeight() { return h; }
	IVector getImageDim() { IVector i(2); i[0]=w; i[1]=h; return i; }
	int getImageSize() {return nImagePts;}

    protected:
	QMMesh * FEMMesh;
	Camera * camera;
	int * bndNodes;
	int w, h, nImagePts, nBndNodes, nNodes;
	RCompRowMatrix toImage, toField; //, toData;
	virtual void constructMatrix()=0;
};

class STOASTLIB RayProjector : public Projector
{
 public:
    RayProjector() {}
    RayProjector(Camera * cam, QMMesh * mesh);

 protected:
    int intersectBoundary(const RVector & ray);
    void constructMatrix();
};


#endif
