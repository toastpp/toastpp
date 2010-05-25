#ifndef __DIFFUSPHERE_H
#define __DIFFUSPHERE_H



class DiffuSphere {
 public:

  DiffuSphere();  //  Constructor
  ~DiffuSphere();  // Destructor

  IVector DirNeighbours(int i) {return DirNBR[i];}//The list of Direct Neighbours anticlockwise
  
  IVector DirNeighboursClock(int i) {return DirNBRClock[i];}// The list of direct neighbours clockwise
  void read_net(char *name);
  friend std::istream& operator>> (std::istream& i, DiffuSphere& df);
  friend std::ostream& operator<< (std::ostream& o, DiffuSphere& df);

  double GiveThita(int i) {return Thita[i];}
  double GiveFi(int i) {return Fi[i];}
  int NoDirNeighbours(int i){return NoDNBR[i];}
  // int NoVRTin(){return NoVRTforOutput;}
  void DoDifussion (char *fname);

  
 private:
  IVector  *DirNBRClock;
  
  SurfaceNet surfnet;
  
  RVector  Thita;
  RVector  Fi;
  IVector  *DirNBR;
  int *NoDNBR;


};
#endif
