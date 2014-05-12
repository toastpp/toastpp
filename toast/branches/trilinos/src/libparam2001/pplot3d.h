#ifndef __PPLOT3D_H
#define __PPLOT3D_H



class Pplot3d
{

 public:

  Pplot3d();   //Constructor
  ~Pplot3d();  //Destructor
  void CalcSh (char *name, char *name1, char *name2, int deg1);// 
  //  friend ostream& operator<< (ostream& o,Pplot3d& PP);
  //  friend istream& operator>> (istream& i,Pplot3d& PP);

 private:
  //some useful parts of other classes
  
  SurfaceNet surfnet;//Member of the SurfaceNet class
  DiffuSphere diffsphr;//Member of the  DiffuSphere class
  OptimiZZZ optimiz;// Member of the OptimiZZZ class


  RVector *OTF;//the optimised Thita-Phi
  RVector *Vertice;//the initial vertice possitions
  CVector *Coef;//The complex coefficients vector;
  CVector *CVertice;//the initial vertice possitions
  int NoVRT;//the number of vertices


  int maxdegree, len ;
  IVector mm, ll, NsH ;
  

};

#endif
