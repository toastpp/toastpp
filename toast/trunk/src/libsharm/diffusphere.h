#ifndef __DIFFUSPHERE_H
#define __DIFFUSPHERE_H

// Symbol import/export direction
#ifdef SHARMLIB_IMPLEMENTATION
#define SHARMLIB DLLEXPORT
#else
#define SHARMLIB DLLIMPORT
#endif

class SHARMLIB diffusphere {
 
	public:

	diffusphere();  //  Constructor
	~diffusphere();  // Destructor

  void  DoDif (IDenseMatrix &Faces, RDenseMatrix &Vertices, RVector &Thita, RVector &Fi);
  //IVector DirNeighbours(int i) {return DirNBR[i];}//The list of Direct Neighbours anticlockwise
  
 // IVector DirNeighboursClock(int i) {return DirNBRClock[i];}// The list of direct neighbours clockwise
 // void read_net(char *name);
 // friend istream& operator>> (istream& i, DiffuSphere& df);
 // friend ostream& operator<< (ostream& o, DiffuSphere& df);

 // double GiveThita(int i) {return Thita[i];}
//  double GiveFi(int i) {return Fi[i];}
 // int NoDirNeighbours(int i){return NoDNBR[i];}
  // int NoVRTin(){return NoVRTforOutput;}
  //void DoDifussion (char *fname);

  
 private:

	int NoVRT , NoFac;
    IDenseMatrix node;//the mesh net. 


  IVector  *DirNBRClock;
  
  //SurfaceNet surfnet;
  
  RVector  Thita, Fi, zerv;
  //RVector  Fi;
  IVector  *DirNBR;
  int *NoDNBR;


};
#endif
