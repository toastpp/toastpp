#ifndef __SURFACENET_H
#define __SURFACENET_H

// Symbol import/export direction
#ifdef HARMONICLIB_IMPLEMENTATION
#define HARMONICLIB DLLEXPORT
#else
#define HARMONICLIB DLLIMPORT
#endif



class Point;

class HARMONICLIB SurfaceNet {
 public:
  SurfaceNet();        //        Constructor
  ~SurfaceNet();      //         Destructor

  Point Vertex(int i){return VertexList[i];}//The list of vertices
 
  IVector Faces(int i) {return AntiClocFaces[i];}//The faces Anticlockwise
 
  int NoVertices (void) {return  NeWNoVRT;}
  int NoFaces (void) {return LengoFaces;}
  friend std::istream &operator>> (std::istream& i, SurfaceNet& sn);
  friend std::ostream &operator<< (std::ostream& o, SurfaceNet& sn);
  void ReadVoxels(char *fname);
  void WriteVoxels(char *ofname);
  void WriteFaces(char *ofname);  
 
  
 
 private:
  Point *VertexList;

  IVector *AntiClocFaces;
  int NeWNoVRT;
  int LengoFaces;

};

bool inpointlist(Point *list,int length,int i,int j, int k);
//int  drop(int *list,int length,int chop,int * newlist);


#endif
