#ifdef __BORLANDC__
#include <strstrea.h>
#else
#include <sstream>
#endif
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include "harmoniclib.h"
#include "surfacenet.h" 
#include "diffusphere.h"
#include "optimizzz.h"
#include "usefulfan.h"
#include "legendre.h"
#include "pplot3d.h"


//        Constructor Definition

Pplot3d::Pplot3d()
{
  OTF=NULL;
  Vertice=NULL;
  Coef=NULL;
  CVertice=NULL;
}

//         Destructor Definition

Pplot3d::~Pplot3d()
{
  if (OTF) delete []OTF ;
  if (Vertice) delete []Vertice;
  if (Coef) delete []Coef;
  if (CVertice) delete []CVertice;  
}



// *****************************************************************************
// *       Few Additions to the RCompRowmatrix and RCoordMatrix                *
// *****************************************************************************
//*                                                                            *

//** Use of RCompRowMatrix from Toast to create the RCompMatrix ATA  from the RCompmatrix A
//** The increase in speed is significant 

class myCCompRowMatrix: public CCompRowMatrix
{
public:
  myCCompRowMatrix(const CCoordMatrix &cm): CCompRowMatrix (cm) {}
  myCCompRowMatrix(const CCompRowMatrix  &cm): CCompRowMatrix (cm) {}
  CCompRowMatrix AttA ();
  CVector DiagATA ();
};


CVector  myCCompRowMatrix::DiagATA ()
{
  // cerr << "start AttA" << endl;
  myCCompRowMatrix At(*this);
  At.Transpone();
  CVector Diag(At.nRows());

  for(int g=0;g<At.nRows();g++)Diag[g]=0.0;
  
  for(int i=0;i<At.nRows();i++)
    {
      for(int j=0;j<At.nCols();j++)
	{
	  //  cerr<<i<<"  " <<At(i,j)<<endl;
	  if(At.Get(i,j)!=0)  Diag[i]+=At(i,j)*At(i,j); 
	}
    }
 
  return Diag;
}




CCompRowMatrix myCCompRowMatrix::AttA ()
{
  // cerr << "start AttA" << endl;

  myCCompRowMatrix At(*this);

  At.Transpone();

  int i, j,k;//, n = nCols();
  complex at_ik, at_jk;
  // cerr<<"A"<<endl;
  //CCoordMatrix ata(At.nRows(), At.nRows());
  CDenseMatrix ata2(At.nRows(), At.nRows());
  for (i = 0; i < At.nRows(); i++)
    {
    for (j = 0; j < At.nRows(); j++)
      {
      for(k=0;k<At.nCols();k++)
	{
	  if((at_ik = At.Get(i,k))!=0 && (at_jk = At.Get(j,k))!=0)
	    {
	      ata2(i,j)+=at_ik*at_jk;
	    }
	}
      }
    }

  // cerr<<"B"<<endl;
  CCompRowMatrix ata(ata2);

  return ata;
}
 

//************************************************************ 
//**************   Pplot3d functions   *********************** 
//************************************************************ 




void  Pplot3d::CalcSh(char *name, char *name2, char *name3, int deg1)

{
  
  //************************************************************ 
  //     Loads the class object of the class surfacenet
  //************************************************************
    
  ifstream fin(name);
  if (!fin)
    cout << "Couldn;t open fin " <<endl;
  fin >>  surfnet; 
  //************************************************************ 
  //     Loads the class object of the class DiffuSphere
  //************************************************************
    
  ifstream f2in(name2);
  if (!f2in)
    cout << "Couldn;t open f2in " <<endl;
  f2in >> diffsphr; 
  //************************************************************ 
  //     Loads the class object of the class OptimiZZZ
  //************************************************************
    
  ifstream f3in(name3);
  if (!f3in)
    cout << "Couldn;t open f3in " <<endl;
  f3in >> optimiz; 



  //************************************************************ 
  //*****************   Initialisation  ************************ 
  //************************************************************ 



    NoVRT=surfnet.NoVertices();// The number of vertices 


    //Creates the OTF where Thita and Fi are stored 

    OTF =new RVector[NoVRT];
    for(int i=0;i<NoVRT;i++) OTF[i].New(2);
    

    //  Thita and Fi from optimiz

    for(int i=0;i<NoVRT;i++)
      {
	OTF[i][0]= (optimiz.FinalThita (i));
	OTF[i][1]= (optimiz.FinalPhi (i));
      }

    //Creates the Vertice  where the initial coordinates are stored from surfnet

    Vertice =new RVector[NoVRT];
    for(int i=0;i<NoVRT;i++) Vertice[i]=surfnet.Vertex(i);


    //Creates the CVertice  where the initial coordinates are stored as complex

    CVertice =new CVector[NoVRT]; 
    for(int i=0;i<NoVRT;i++) 
      {
	CVertice[i].New(3);
	SetReal( CVertice[i],surfnet.Vertex(i));
      }

 



    //The degree  of the reconstruction with spherical harmonics

    maxdegree = deg1;

    NsH.New(maxdegree);


    // and the length len of the Spherical Harmonic coefficients

    len=1;
    for (int j=1; j<=maxdegree; j++)
      {
	len+=2*j+1;
	NsH[j-1]=len;
      }  

    cout<<len<<endl;
 
    mm.New(len);
    ll.New(len);
 

    int w = 0;
    
    for(int i=0;i<=maxdegree;i++)
      {
	for(int p=-i;p<=i;p++)
	  { 
	    mm[w]=p;
	    ll[w]=i;
	    w++;
	  }
      }


    cout<<mm<<endl;
    cout<<ll<<endl;
    cout <<NsH<<endl;

    //Creates the Coef to store the coefficients
    Coef =new CVector[len]; 
    for(int i=0;i<len;i++) Coef[i].New(3);
    

    // The matrix with the spherical harmonics calculated on the thia and fi
    CDenseMatrix B1(NoVRT,len);
  
    for ( int i=0;i<NoVRT;i++)
      {
	for ( int j =0;j<len;j++)
	  {
	    B1(i,j)=SphericalHarmonic(ll[j], mm[j], OTF[i][0], OTF[i][1]); 
	  }
      }

    //Transform to CompRow Matrix
    CCompRowMatrix B2(B1);
    myCCompRowMatrix B(B2);

    // calculate AtA = BTB
    CCompRowMatrix AtA = B.AttA();
 
    cerr<<"matrix constructed !!"<<endl;

    double tol=1.0e-16;

    for (int  r=0;r<3;r++)//for x/y/z
      {
	CVector CVerX(NoVRT);// The Vertice coordinates for x/y/z
	CVector CoefX(len);// The Coefficients for x/y/z

	for(int i=0;i<NoVRT;i++)
	  {
	    CVerX[i] = CVertice[i][r];
	  }

	CVector BTCVerX(len);
	B.ATx(CVerX, BTCVerX);// Calculate BT.Vertice

	BiCG( AtA, BTCVerX, CoefX, tol);

	cerr<<"BiCG DONE !!!!"<<endl;
	cout<<CoefX<<endl;


	//Copying the values to Coef

     	for(int i=0;i<len;i++)
	  {
	    Coef[i][r] =  CoefX[i];
	  }
      }


}//Main function ends here !!


