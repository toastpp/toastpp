#define STOASTLIB_IMPLEMENTATION

#ifdef MESA_SUPPORT

#if defined(_WIN32)
#include <windows.h>
#include <WinGDI.h>
#endif
#include "stoastlib.h"
#include "GLProjector.h"
#ifdef __APPLE__
#include <OpenGl/OpenGl.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <OpenGL/CGLTypes.h>
#include <OpenGL/CGLCurrent.h>
#else
#include "GL/gl.h"
#include "GL/glu.h"
#endif
#include <iostream>
#include "qmmesh.h"

using namespace std;

GLProjector::GLProjector(Camera * cam, QMMesh * mesh, int nBndElems_, int * bndellist_, int * bndsdlist_)
    : Projector(cam, mesh), bndellist(bndellist_), bndsdlist(bndsdlist_), nBndElems(nBndElems_)
{
    setupMesaGL();
    renderImageElem();
    constructMatrix_SF();
    return;

}

GLProjector::GLProjector(Camera * cam, QMMesh * mesh)
 : Projector(cam, mesh)
{
    nBndElems = FEMMesh->BoundaryList (&bndellist, &bndsdlist);
    setupMesaGL();
    renderImageElem();
    constructMatrix_SF();
}

void GLProjector::init(Camera * cam, QMMesh * mesh)
{
    Projector::init(cam, mesh);
    nBndElems = FEMMesh->BoundaryList (&bndellist, &bndsdlist);
    setupMesaGL();
    renderImageElem();
    constructMatrix_SF();
}

void GLProjector::setupMesaGL()
{
#ifdef __APPLE__
    // TODO
#else
    // Set up MesaGL context here
#if OSMESA_MAJOR_VERSION * 100 + OSMESA_MINOR_VERSION >= 305
    /* specify Z, stencil, accum sizes */
    OSMesaContext mesaGLContext = OSMesaCreateContextExt( GL_RGBA, 16, 0, 0, NULL );
#else
    OSMesaContext mesaGLContext = OSMesaCreateContext( GL_RGBA, NULL );
#endif
    if (!mesaGLContext) {
	cout<<"OSMesaCreateContext failed!\n"<<endl;
    }   

    // Allocate buffer
    imageBuffer = new GLubyte[w * h * 4];
    /* Bind the buffer to the context and make it current */
    if (!OSMesaMakeCurrent( mesaGLContext, imageBuffer, GL_UNSIGNED_BYTE, w, h )) {
	cout<<"OSMesaMakeCurrent failed!\n"<<endl;
    }
    
#endif

}

GLProjector::~GLProjector()
{
    delete[] imageBuffer;
    delete[] bndellist;
    delete[] bndsdlist;
#ifdef __APPLE__
    // TODO
#else
    OSMesaDestroyContext(mesaGLContext);
#endif
}

RVector GLProjector::getImageBuffer()
{
    RVector result(w*h);
    for (int i=0, j=0; i<(w*h); ++i, j+=4)
    {
	int idx = imageBuffer[j] 
	    + imageBuffer[j+1]*256
	    + imageBuffer[j+2]*256*256 -1;
	result[i] = (double)(idx);
    }
    return result;
}

void GLProjector::renderImage()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    
    // Perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    GLdouble zFar = l2norm(camera->getPos())*2.0;
    double fovy = camera->getFoVy();
    if (!fovy) { // assume ortho projection
	glOrtho(-50.0, 50.0, -50.0, 50.0
		/*-18.3084, 18.3084, -18.3084, 18.3084*/, 10.0, zFar);
	//glScalef(2.5, 2.5, 2.5);
    } else {
	double aspect = camera->getAspect();
	gluPerspective( fovy, aspect, 5.0, zFar );
    }

    // Viewing angle, viewport
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0,0,w,h);

    RVector pos = camera->getPos();

    RVector cnt = pos + camera->getViewDir();
    RVector y = camera->getUp();
    gluLookAt( pos[0], pos[1], pos[2],
	    cnt[0], cnt[1], cnt[2],
	    y[0], y[1], y[2]);

    // Loop over boundary sides, and render as quads with nearest-neighbour node idx as colour
   //int nBndElems = FEMMesh->BoundaryList (&bndellist, &bndsdlist);
   glBegin(GL_QUADS);
   for (int e=0; e<nBndElems; ++e)
   {
	// Get element
	int elemIdx = bndellist[e];
	Element * elem = FEMMesh->elist[elemIdx];
	// Get side info
	int side = bndsdlist[e];
	int nSideNodes = elem->nSideNode(side);
	// Centre of element in global (ie. not local to element) coords
	Point elemCnt = elem->Global(FEMMesh->nlist, elem->SideCentre(side)); 

	// Loop over nodes and render a quad for each
	for (int n=0; n<nSideNodes; ++n)
	{
	    int nodeIdx = elem->Node[elem->SideNode(side, n)];  // index in mesh->s list

	    // Generate colour from node index (+1 so that can distinguish bg as 0)
	    unsigned char b1 = (nodeIdx +1) % 256;
	    unsigned char b2 = (nodeIdx + 1)/256 % 256;
	    unsigned char b3 = (nodeIdx + 1)/(256*256) % 256;
	    glColor4ub(b1, b2, b3, 255);

	    // Get mid points of adjacent edges
	    int nextNodeIdx = elem->Node[elem->SideNode(side, (n+1)%nSideNodes)];
	    int prevNodeIdx;
	    if (n==0)
	    {
		prevNodeIdx = elem->Node[elem->SideNode(side, nSideNodes-1)];
	    }else{
		prevNodeIdx = elem->Node[elem->SideNode(side, n-1)];
	    }
	    Node & node = FEMMesh->nlist[nodeIdx];
	    Node & nextNode = FEMMesh->nlist[nextNodeIdx];
	    Node & prevNode = FEMMesh->nlist[prevNodeIdx];
	    Point midNext = 0.5*(node + nextNode);	
	    Point midPrev = 0.5*(node + prevNode);

	    // Render quad
	    glVertex3f(node[0], node[1], node[2]);
	    glVertex3f(midNext[0], midNext[1], midNext[2]);
	    glVertex3f(elemCnt[0], elemCnt[1], elemCnt[2]);
	    glVertex3f(midPrev[0], midPrev[1], midPrev[2]);
	}
    }
    glEnd();
    glFlush();

}
void GLProjector::renderImageElem()
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    
    // Perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    GLdouble zFar = l2norm(camera->getPos())*2.0;
    double fovy = camera->getFoVy();
    double iw = camera->getImageWidth() * camera->getPixelSize();
    double ih = camera->getImageHeight() * camera->getPixelSize();
    if (!fovy) { // assume ortho projection
	glOrtho(-iw*0.5, iw*0.5, -ih*0.5, ih*0.5, 5.0, zFar);
	//glOrtho(-18.3084, 18.3084, -18.3084, 18.3084, 10.0, zFar);
    } else {
	double aspect = camera->getAspect();
	gluPerspective( fovy, aspect, 5.0, zFar );
    }
    // Viewing angle, viewport
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0,0,w,h);

    RVector pos = camera->getPos();

    RVector cnt = pos + camera->getViewDir();
    RVector y = camera->getUp();
    gluLookAt( pos[0], pos[1], pos[2],
	    cnt[0], cnt[1], cnt[2],
	    y[0], y[1], y[2]);

    // Loop over boundary sides, and render as quads with nearest-neighbour node idx as colour
   //int nBndElems = FEMMesh->BoundaryList (&bndellist, &bndsdlist);
   for (int e=0; e<nBndElems; ++e)
   {
        glBegin(GL_POLYGON);

	// Get element
	int elemIdx = bndellist[e];
	Element * elem = FEMMesh->elist[elemIdx];
	// Get side info
	int side = bndsdlist[e];
	int nSideNodes = elem->nSideNode(side);
	// Centre of element in global (ie. not local to element) coords
	Point elemCnt = elem->Global(FEMMesh->nlist, elem->SideCentre(side)); 

	// Generate colour from element index (+1 so that can distinguish bg as 0)
	unsigned char b1 = (e + 1) % 256;
	unsigned char b2 = (e + 1)/256 % 256;
	unsigned char b3 = (e + 1)/(256*256) % 256;
	glColor4ub(b1, b2, b3, 255);
	//cout << "Drawing colour "<< int(b1) << ", " << int(b2) << ", " << int(b3) << " for boundary elem "<< e << endl;

	// Loop over nodes and render a quad for each
	for (int n=0; n<nSideNodes; ++n)
	{
	    int nodeIdx = elem->Node[elem->SideNode(side, n)];  // index in mesh->s list
	    Node & node = FEMMesh->nlist[nodeIdx];
	    // Render quad
	    glVertex3f(node[0], node[1], node[2]);
	}

	glEnd(); //end polygon
    }
    glFlush();

}

void GLProjector::constructMatrix()
{
    RVector rowvec(nNodes);
    for (int r=0, j=0; r<(w*h); ++r, j+=4)
    {
	int nodeIdx = imageBuffer[j] 
	    + imageBuffer[j+1]*256
	    + imageBuffer[j+2]*256*256 - 1;
	if (nodeIdx>=0)
	{
	    rowvec[nodeIdx] = 1.0;
	    toImage.SetRow(r, rowvec);
	    rowvec[nodeIdx] = 0.0;
	}
    }

    // Scale each column by 1/l1norm(column) so that P^T*P=I
    toField = transp(toImage);
    RVector scale(nNodes);
    for (int c=0; c<nBndNodes; ++c)
    {
	 double s = l1norm(toField.Row(bndNodes[c]));
	 if (s>0.0) scale[bndNodes[c]] = 1.0 / s;
    }
    toField.RowScale(scale);
}

void GLProjector::constructMatrix_SF()
{
    RDenseMatrix jacMat; // not used, just for argument
    RVector rowvec(nNodes), frowvec(nImagePts); //, dataRowVec(nNodes);
    int cnt=0, tot=0;
    for (int r=0, j=0; r<(w*h); ++r, j+=4)
    {
	int ix = r%w, iy = r/w;
/*	if ( (ix>60) && (ix<80) && (iy>40) && (iy<80) )
	{*/
	    int e = imageBuffer[j] 
		+ imageBuffer[j+1]*256
		+ imageBuffer[j+2]*256*256 - 1;
	    //cout << "Reading colour " << int(imageBuffer[j]) << ", " << int(imageBuffer[j+1]) << ", " << int(imageBuffer[j+2]) << " for boundary elem "<< e << endl;
	    if (e>=0)
	    {
		int elemIdx = bndellist[e];
		int sideIdx = bndsdlist[e];
		RVector rayDir, rayPos;
		camera->getRayVector(ix, iy, rayPos, rayDir);	// ray through pixel
		Element * elem = FEMMesh->elist[elemIdx];
		int nSideNodes = elem->nSideNode(sideIdx);
		/*cout << "Getting direction cosine for side "<<sideIdx
		    << " of element "<<elemIdx
		    << ", with boundary elem idx "<<e
		    << " / "<<nBndElems<<endl;*/
		RVector normal = elem->DirectionCosine(sideIdx, jacMat);
		RVector elemPos = FEMMesh->nlist[elem->Node[elem->SideNode(sideIdx, 0)]];
		Point intersect = rayPos + rayDir * dot(elemPos-rayPos, normal) / dot(rayDir, normal);
		/*cout << "intersect error = "<<dot(intersect-elemPos, normal)<<endl;
		  cout << "intersect dist to elem[0] "<<l2norm(intersect-elemPos)<<endl;
		  cout << "Side size "<< elem->SideSize(sideIdx, FEMMesh->nlist) <<endl;*/
		RVector shapeFun = elem->GlobalShapeF(FEMMesh->nlist, intersect);
		shapeFun /= l1norm(shapeFun);
		tot+=nSideNodes;
		if (vmin(shapeFun)<-1e-2)
		{
//		    cout <<"Shape fun is "<<shapeFun<<endl;
		    for (int i=0; i<nSideNodes; ++i)
			if (shapeFun[i]<-1e-2) cnt++; 
		}
		//shapeFun.Clip(0.0, 1.0);
		for (int n=0; n<nSideNodes; ++n)
		{
		    int elemNodeIdx = elem->SideNode(sideIdx, n);
		    int meshNodeIdx = elem->Node[elemNodeIdx];
		   
            // To field matrix maps node to 1 pixel only
		    /*double npx, npy;
		    camera->getPixelCoords(node, npx, npy);
		    double pixel = npx + npy*w;
		    frowvec[(int)(pixel)] = 1.0;
		    toField.SetRow(meshNodeIdx, frowvec);
		    frowvec[(int)(pixel)] = 0.0;*/

		    // Shape function based mapping to image
		    rowvec[meshNodeIdx] = shapeFun[elemNodeIdx];
		}
		toImage.SetRow(r, rowvec);
		rowvec.Clear();
	    }
/*	}*/
    }

    // Scale each column by 1/l1norm(column) so that P^T*P=I
    toField = transp(toImage);
    RVector scale(nNodes);
    double m=0;
    for (int c=0; c<nBndNodes; ++c)
    {
	int nidx = bndNodes[c];
	RVector row = toField.Row(nidx);
	double s = l1norm(row);
        //cout <<s<<endl;
	if (m<s) m=s;
    }



    for (int c=0; c<nBndNodes; ++c)
    {
	
	if (m>0.0)
	{
 	//cout <<"s = "<<s<<endl;
	    scale[bndNodes[c]] = 1.0/m;
	}
    }

    toField.RowScale(scale);
}

#endif // MESA_SUPPORT