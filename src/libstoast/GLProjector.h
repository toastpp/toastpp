#ifndef GLPROJECTOR_H
#define GLPROJECTOR_H

#ifdef MESA_SUPPORT

#include "projector.h"
#include "GL/osmesa.h"

class STOASTLIB GLProjector : public Projector
{
    public:
	GLProjector() {}
	GLProjector(Camera * cam, QMMesh * mesh, int nBndElems_, int * bndellist_, int * bndsdlist_);
	GLProjector(Camera * cam, QMMesh * mesh);
	~GLProjector();
	void init(Camera * cam, QMMesh * mesh);
	void setupMesaGL();
	void constructMatrix();
	void constructMatrix_SF();
	RVector getImageBuffer();
	void renderImage();
	void renderImageElem();
    protected:
	OSMesaContext mesaGLContext;
	GLubyte * imageBuffer;
	int * bndellist, * bndsdlist;
	int nBndElems;
};

#endif // MESA_SUPPORT
#endif //!GLPROJECTOR_H
