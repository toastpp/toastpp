// some goo

#ifndef __SURFACENET_H
#define __SURFACENET_H

// Symbol import/export direction
#ifdef SHARMLIB_IMPLEMENTATION
#define SHARMLIB DLLEXPORT
#else
#define SHARMLIB DLLIMPORT
#endif



class Point;

class SHARMLIB SurfaceNet {
	

public:
	
	SurfaceNet();        //        Constructor
	~SurfaceNet();       //        Destructor


	void Vox2surf(IDenseMatrix &voxels, IDenseMatrix &Faces, RDenseMatrix &Vertex);

	//Point Vertex(int i){return VertexList[i];}//The list of vertices
	//IVector Faces(int i) {return AntiClocFaces[i];}//The faces Anticlockwise
	//int NoVertices (void) {return  NeWNoVRT;}
	//int NoFaces (void) {return LengoFaces;}
	//void ReadVoxels(char *fname);
	
 // friend istream& operator>> (istream& i, SurfaceNet& sn);
 // friend ostream& operator<< (ostream& o, SurfaceNet& sn);
 // void WriteVoxels(char *ofname);
 // void WriteFaces(char *ofname);  
	
  
private:
	Point *VertexList;

	IVector *AntiClocFaces;
	int NeWNoVRT;
	int LengoFaces;

};

bool inpointlist(Point *list,int length,int i,int j, int k);
//int  drop(int *list,int length,int chop,int * newlist);


#endif
