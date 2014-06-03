#ifndef __OPTIMIZZZ_H
#define __OPTIMIZZZ_H

// Symbol import/export direction
#ifdef HARMONICLIB_IMPLEMENTATION
#define HARMONICLIB DLLEXPORT
#else
#define HARMONICLIB DLLIMPORT
#endif


class HARMONICLIB OptimiZZZ {
 public:
  
  OptimiZZZ();  //  Constructor
  ~OptimiZZZ();  // Destructor

 
  void optim ( char *name, char *name2);// The main optimisation function
  friend std::ostream& operator<< (std::ostream& o,OptimiZZZ& TF);
  friend std::istream& operator>> (std::istream& i,OptimiZZZ& TF);
  double FinalThita(int i) {return ParaVec[i][0];}
  double FinalPhi(int i) {return ParaVec[i][1];}


 private:
  
  // some useful parts of the other classes

  SurfaceNet surfnet;//Member of the SurfaceNet class
  DiffuSphere diffsphr;//Member of the  DiffuSphere class


  //The functions of the optimiZZZ class

  int NofIneq (RVector *Vec1,int Len);// Returns the number of active inequalities
  double CVecConstruction (RVector *Vec1,int Len);// creates the constraint vector
  void CVecJacobian (void);// creates the Jacobian of the constraint vector
  void CVecDiJacobian(void);//The discrete Jacobian calculation
  double  Fobj_calc(RVector *Vec1,int Len);// calculates the objective function at the current point
  void Vec_normaZ(RVector *Vec1,int Len);// normalises the nodes to be on the sphere surface
  void CalcDeltaF(void);//Calculate DeltaF : Gradient of objective function
  void F_min(void);// minimises the F objective function using Lagrange multipliers
  void Con_Min(void);//Minimises the C constraint vector using Newton Raphson
  void Springs(void);//Uses springs 
  // Some global variables and vectors


  int NoVRT;//the number of vertices
  int NoFac;// The number of faces

  RVector *Vec;// the current possition of the nodes
  RVector *ParaVec;//The parameters for Thita and Fi

  int LoC;//The length of the constraints vector
  RVector ConVec;//The constraints vector
  RCompRowMatrix  JConVec;//The jacobian of the constraint vector


  double Fobj_value;// The current value of the objective function
  RVector DeltaF;//The Gradient of the objective function F
  RVector Lamda;// The vector of Lagrange multipliers

 
  RVector GFstore1, GFstore2, Pstore2, Pstore1;//Used for the conjugate gradient
  RVector P;// The direction of search
  double gamma;//the gamma factor for CG 


  int Iter;// the current iteration
  int NoIter;//The number of total iterations
  RVector MonitorCost, MonitorP,MonitorF;// vectors for monitoring purpose
  
  double steplarge;//the max step
  double His;// The threshold for the inequalities, from 0 to 1 with 1 the strongest.
  int ToloIneq;//Number of  inequalities to tolerate
  double linetol;// The line tolerance for the line search
  int fail;//indicator for failure to procede
};





#endif
