#define HARMONICLIB_IMPLEMENTATION
//============================================================//
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mathlib.h"
#include "felib.h"
#include "harmoniclib.h"
#include "surfacenet.h" 
#include "diffusphere.h"
#include "optimizzz.h"
#include "usefulfan.h"

//        Constructor Definition

OptimiZZZ::OptimiZZZ()
{
 
  Vec=NULL;
  ParaVec=NULL;
}


//         Destructor Definition

OptimiZZZ::~OptimiZZZ() 
{
  if (Vec) delete []Vec ;
  if (ParaVec) delete []ParaVec ;
}





 std::ostream& operator<< ( std::ostream& o, OptimiZZZ& TF)// Exporting the results
{
  int Nove= TF.surfnet.NoVertices();   
  o << Nove<< std::endl;// Writing the number of Vertices  
  for(int j=0;j<Nove;j++)
    {
      o<< TF.ParaVec[j][0]<< std::endl;// Writing thita
    }
  o<< std::endl;
  for(int j=0;j<Nove;j++)
    {
      o<< TF.ParaVec[j][1]<< std::endl;// Writing Fi
    }
  return o;
}

 std::istream& operator>> ( std::istream& i, OptimiZZZ& TF)// Exporting the results
{
  int Nove= TF.surfnet.NoVertices();   
  i >> Nove; // reading the number of Vertices
  
  TF.ParaVec=new RVector[Nove];
  RVector tmp(2);
  for(int j=0;j<Nove;j++)
    {
      TF.ParaVec[j]=tmp;
      i>> TF.ParaVec[j][0];// Writing thita
    }
  for(int j=0;j<Nove;j++)
    {
      i>> TF.ParaVec[j][1];// Writing Fi
    }


  return i;
}




//******************************************************************************
//**************   OptimiZZZ::optim(char *name, char *name2)    ****************
//******************************************************************************


void OptimiZZZ::optim(char *name, char *name2)
{
  
//************************************************************
//      Loads the class object of the class SurfaceNet
//************************************************************

  std:: ifstream fin(name);
  if (!fin)
     std::cout << "Couldn;t open fin " << std::endl;
  fin >> surfnet; 

//************************************************************ 
//     Loads the class object of the class DiffuSphere
//************************************************************
    
    std:: ifstream f2in(name2);
    if (!f2in)
      std:: cout << "Couldn;t open f2in " << std::endl;
    f2in >> diffsphr; 


//************************************************************
//* Initiallisation of the nessecary vectors and matrices;   *
//************************************************************
 

  
 NoVRT=surfnet.NoVertices();
 NoFac=surfnet.NoFaces();
  
 Vec =new RVector[NoVRT]; //Creates the Vec where all the unit vectors are stored
 for(int i=0;i<NoVRT;i++) Vec[i].New(3);
  
 ParaVec =new RVector[NoVRT]; //Creates the Vec where Thita and Fi are stored
 for(int i=0;i<NoVRT;i++) ParaVec[i].New(2);
  

 //  from Spherical Coordinates to XYZ for Vec

 for(int i=0;i<NoVRT;i++)
    {
      Vec[i][0]= cos( diffsphr.GiveFi (i)) * sin(diffsphr.GiveThita (i));
      Vec[i][1]= sin( diffsphr.GiveFi (i)) * sin(diffsphr.GiveThita (i));
      Vec[i][2]= cos( diffsphr.GiveThita (i));
    }
  
 
 GFstore1.New(3*NoVRT);
 GFstore2.New(3*NoVRT);

 Pstore1.New(3*NoVRT);
 Pstore2.New(3*NoVRT);

 P.New(3*NoVRT);
 Lamda.New(NoFac);

 ToloIneq=0;//tolerance to inequalities set to zero
 linetol=0.00000001;

 
 steplarge=0;// Step set to Zero 
 gamma=0.0;// set to zero
 His=-0.1;// The tolerance for sin of angles
 Iter=0;// set to zero

 fail=0;
 
 for(int d=0;d<(3*NoVRT);d++)// Initialising to zero
   { 
     GFstore1[d]=0.0;
     GFstore2[d]=0.0;
   }
 for(int i=0;i<NoFac;i++) Lamda[i]=0.0;


 
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 NoIter=150;// The number of the Iterations to be performed

 //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 MonitorP.New(NoIter);
 MonitorCost.New(NoIter+1);
 MonitorF.New(NoIter+1);
 MonitorCost[0]= CVecConstruction((RVector *) 0,0); 
 MonitorF[0]=Fobj_calc((RVector *) 0,0);
 
  std::cout<<"No of inequalities = "<<NofIneq((RVector *) 0,0)<< std::endl;

//************************************************************
//*               The iteration steps                        *
//************************************************************


  for (Iter=0;Iter<NoIter;Iter++)
    {
      Con_Min();
      F_min();
     Springs();

       std::cout<<"--------------------------------------------------------------------------------------"<< std::endl;
       std::cout<<"! Iter # "<<Iter<<" F_Cost= "<<Fobj_calc((RVector *) 0,0)<<" C_Cost= "
	  <<CVecConstruction((RVector *) 0,0)<<"Noineq = "<<NofIneq((RVector *) 0,0)<< std::endl;
       std::cout<<"--------------------------------------------------------------------------------------"<< std::endl;
       std::cout<< std::endl;

      Vec_normaZ((RVector *) 0,0);

      MonitorCost[Iter+1]=CVecConstruction((RVector *) 0,0);

      MonitorF[Iter+1]=Fobj_calc((RVector *) 0,0); 
    }
  


//************************************************************
//******            The Results Display                 ****** 	 
//************************************************************


  std::cout<<"C: ";
 for(int k=0;k<NoIter+1;k++) std::cout<<MonitorCost[k]<<" , ";
  std::cout<< std::endl<< std::endl<<"F: ";
 for(int k=0;k<NoIter+1;k++) std::cout<<MonitorF[k]<<" , ";
  std::cout<< std::endl<< std::endl<<"P: ";
 for(int k=0;k<NoIter;k++) std::cout<<MonitorP[k]<<" , ";
  std::cout<< std::endl<< std::endl<<"convec: ";
 for(int k=0;k<NoFac;k++) std::cout<<ConVec[k]<<" , ";
 
 for(int i = 0; i < NoVRT; i++)
   {
     double r=NormV(Vec[i]);
     double fi=atan3(Vec[i][1],Vec[i][0]);
     double thita = acos(Vec[i][2]/r);
     ParaVec[i][0] = thita;
     ParaVec[i][1]= fi;
   }

  std::cout<< std::endl<< std::endl;
  std::cout<< " Thita = {";
 for(int i = 0; i < NoVRT; i++) std::cout<<ParaVec[i][0]<<" , ";
  std::cout<<" } "<< std::endl<< std::endl;
  std::cout<< " Fi = {";
 for(int i = 0; i < NoVRT; i++) std::cout<<ParaVec[i][1]<<" , ";
  std::cout<<" } "<< std::endl;
}//end of Function 


//*                                                                            *
//*                                                                            *
//******************************************************************************
//*                    The other functions of the class OptimiZZZ              *
//******************************************************************************
//*                                                                            *



//*                                                                            *
//******************************************************************************
//* NoFIneq : returns the number of active inequalities according to His       * 
//******************************************************************************



int OptimiZZZ::NofIneq (RVector *Vec1,int Len)
{
  int ineq=0;
 
  if(Len==0)
    {
      for (int j=0;j<NoFac;j++)
	{
	  for (int ve=0;ve<4;ve++)
	    {
	      int Tv=surfnet.Faces(j)[ve];
	      int Pve=ve-1;
	      int Nve=ve+1;
	      if (ve==3){Pve=2;Nve=0;}
	      if (ve==0){Pve=3;Nve=1;}
	      int Pv =surfnet.Faces(j)[Pve];
	      int Nv =surfnet.Faces(j)[Nve];
	      double angle = LDAB(Vec[Pv],Vec[Tv],Vec[Nv]);
	      double sinangle=sin(angle);
	      if(sinangle < His){ineq++;}
	    }
	}
    }
  else
    {
      for (int j=0;j<NoFac;j++)
	{
	  for (int ve=0;ve<4;ve++)
	    {
	      int Tv=surfnet.Faces(j)[ve];
	      int Pve=ve-1;
	      int Nve=ve+1;
	      if (ve==3){Pve=2;Nve=0;}
	      if (ve==0){Pve=3;Nve=1;}
	      int Pv =surfnet.Faces(j)[Pve];
	      int Nv =surfnet.Faces(j)[Nve];
	      double angle = LDAB(Vec1[Pv],Vec1[Tv],Vec1[Nv]);
	      double sinangle=sin(angle);
	      if(sinangle < His){ineq++;}
	    }
	}
    } 
  return ineq;
}

//*                                                                            *
//******************************************************************************
//* CVecConstruction : Constructs the vector of the constraints : ConVec       *  
//******************************************************************************
//*                                                                            *

double OptimiZZZ::CVecConstruction (RVector *Vec1,int Len)
{
  if(Len==0)
    {      
      int NoIneq=NofIneq((RVector*)0,0);
      double Final_Area= (4 * Pi)/(NoFac);
     
      LoC= NoFac+NoIneq;
      ConVec.New(LoC);
      
      //First the Area constraints
      
      for(int k=0;k<NoFac;k++)
	{
	  double Current_Area=0;
	  int NoDA=  surfnet.Faces(k)[0];
	  int NoDB=  surfnet.Faces(k)[1];
	  int NoDC=  surfnet.Faces(k)[2];
	  int NoDD=  surfnet.Faces(k)[3];
	  Current_Area= Area(Vec[NoDA],Vec[NoDB],Vec[NoDC],Vec[NoDD] );
	  ConVec[k]=Final_Area-Current_Area;
	}  
  
      // The Euclidean constraints to keep the nodes on the surface of the sphere
      //Not neccesary since will be enforced
  
      //The angle constraints
      
      if(NoIneq>0) 
	{
	  int IFI= NoFac;
	  for (int j=0;j<NoFac;j++)
	    {
	      for (int ve=0;ve<4;ve++)
		{
		  int Tv=surfnet.Faces(j)[ve];
		  int Pve=ve-1;
		  int Nve=ve+1;
		  if (ve==3){Pve=2;Nve=0;}
		  if (ve==0){Pve=3;Nve=1;}
		  int Pv =surfnet.Faces(j)[Pve];
		  int Nv =surfnet.Faces(j)[Nve];
		  double angle = LDAB(Vec[Pv],Vec[Tv],Vec[Nv]);
		  double sinangle=sin(angle);
		  // cout<<angle<<"  " <<sinangle<< endl;
		  if(sinangle < His){ ConVec[IFI]=sinangle;IFI++;}
		}
	    }
	}
    }
  else// For any other vector
    {
    
      int NoIneq=NofIneq(Vec1,NoVRT);
      double Final_Area= (4 * Pi)/(NoFac);
      LoC= NoFac+NoIneq;
      ConVec.New((LoC));
      
      //First the Area constraints
      
      for(int k=0;k<NoFac;k++)
	{
	  double Current_Area=0;
	  int NoDA=  surfnet.Faces(k)[0];
	  int NoDB=  surfnet.Faces(k)[1];
	  int NoDC=  surfnet.Faces(k)[2];
	  int NoDD=  surfnet.Faces(k)[3];
	  Current_Area= Area(Vec1[NoDA],Vec1[NoDB],Vec1[NoDC],Vec1[NoDD] );
	  ConVec[k]=Final_Area-Current_Area;
	}  
  
      // The Euclidean constraints to keep the nodes on the surface of the sphere
      //Not neccesary since will be enforced


  
      //The angle constraints
      
      if(NoIneq>0) 
	{
	  int IFI= NoFac;
	  for (int j=0;j<NoFac;j++)
	    {
	      for (int ve=0;ve<4;ve++)
		{
		  int Tv=surfnet.Faces(j)[ve];
		  int Pve=ve-1;
		  int Nve=ve+1;
		  if (ve==3){Pve=2;Nve=0;}
		  if (ve==0){Pve=3;Nve=1;}
		  int Pv =surfnet.Faces(j)[Pve];
		  int Nv =surfnet.Faces(j)[Nve];
		  double angle = LDAB(Vec1[Pv],Vec1[Tv],Vec1[Nv]);
		  double sinangle=sin(angle);
		  //  cout<<angle<<"  " <<sinangle<< endl;
		  if(sinangle < His){ ConVec[IFI]=sinangle;IFI++;}
		}
	    }
	}
    

    }

  double cost=0;
  for(int k=0;k<LoC;k++)
    {
      cost +=ConVec[k]*ConVec[k];
      
    }

  return cost;

}//End here 



//*                                                                            *
//******************************************************************************
//*             The jacobian of the constraint  : JConVec using DC/Dx          *
//******************************************************************************
//*                                                                            *

void OptimiZZZ::CVecDiJacobian (void)
{
 
  int NoIneq = NofIneq((RVector*)0,0);
  LoC = NoFac+NoIneq;
  int LoCu=LoC;
  
   JConVec.New(3*NoVRT,LoC);
  
  RCoordMatrix ConJConVec(3*NoVRT,LoC);
  
 

  RVector *InVec = 0;
  InVec = new RVector[NoVRT];
  for(int j=0; j < NoVRT; j++) InVec[j] = RVector(3);
    
  CVecConstruction((RVector *) 0,0);
     
  RVector Cv1(LoC),Cv2(LoC);//create the Cv1 to map the constraint vector at the beginning

 
  for(int s=0;s<LoC;s++) Cv1[s] = ConVec[s];//Cv1 gets the values from the constraint vector ConVec

     
  for(int j=0;j<NoVRT;j++)
    {
      InVec[j][0] = Vec[j][0];
      InVec[j][1] = Vec[j][1];
      InVec[j][2] = Vec[j][2];
    }
  double m_delta;
  
  for(int i = 0 ; i < NoVRT; i++)// For the first NoVRT  columns of JConVec
    {
      if(InVec[i][0]!=0.0) m_delta= (InVec[i][0] / 10.0);
      else m_delta= 0.05;

	  InVec[i][0] += m_delta;  
	  
	  CVecConstruction(InVec,NoVRT);
	
	  RVector Cv3(LoCu);

	  int LoCc=ConVec.Dim();

	  int  LoC_correct = LoCu;

	  if(LoCc<LoCu)
	    {
	      LoC_correct = LoCc;

	      for(int f = LoC_correct; f < LoCu; f++)
		{
		  Cv3[f] = ( - Cv1[f]) / m_delta ;
		  if(Cv3[f])   ConJConVec(i,f)=Cv3[f];
		}
	    }
	  
	  for(int f = 0; f < LoC_correct; f++)
	    {
	      Cv3[f] = (ConVec[f] - Cv1[f]) / m_delta ;
	      if(Cv3[f])   ConJConVec(i,f)=Cv3[f];
	    }

	  InVec[i][0] = Vec[i][0];
    }
      
  
  
  for(int i = 0 ; i < NoVRT; i++)// For the second  NoVRT  columns of JConVec
    {
      if(InVec[i][1]!=0.0) m_delta= (InVec[i][1] / 10.0);
	
      else  m_delta= 0.05;

      InVec[i][1] += m_delta;  
      
      CVecConstruction(InVec,NoVRT);
             
      RVector Cv3(LoCu);

      int LoCc=ConVec.Dim();
      
      int  LoC_correct = LoCu;

      if(LoCc<LoCu)
	{
	  LoC_correct = LoCc;
	  
	  for(int f = LoC_correct; f < LoCu; f++)
	    {
	      Cv3[f] = ( - Cv1[f]) / m_delta ;
	      if(Cv3[f])   ConJConVec(i+NoVRT,f)=Cv3[f];
	    }
	}
                
      for(int f = 0; f < LoC_correct; f++)
	{
	  Cv3[f] = (ConVec[f] - Cv1[f]) / m_delta ;
	  if(Cv3[f])   ConJConVec(i+NoVRT,f)=Cv3[f];
	}
      
      
      InVec[i][1] = Vec[i][1];
    }
  
  
  
  
  for(int i = 0 ; i < NoVRT; i++)// For the third  NoVRT  columns of JConVec
    {
      if(InVec[i][2]!=0.0)  m_delta= (InVec[i][2] / 10.0);
      else  m_delta= 0.05;

      
      InVec[i][2] += m_delta;  
      
      CVecConstruction(InVec,NoVRT);
      
	  
      RVector Cv3(LoCu);
	  
      int LoCc=ConVec.Dim();
      
      int  LoC_correct = LoCu;

      if(LoCc<LoCu)
	{
	  LoC_correct = LoCc;
	  
	  for(int f = LoC_correct; f < LoCu; f++)
	    {
	      Cv3[f] = ( - Cv1[f]) / m_delta ;
	      if(Cv3[f])   ConJConVec(i+2*NoVRT,f)=Cv3[f];
	    }
	}
                
      for(int f = 0; f < LoC_correct; f++)
	{
	  Cv3[f] = (ConVec[f] - Cv1[f]) / m_delta ;
	  if(Cv3[f])   ConJConVec(i+2*NoVRT,f)=Cv3[f];
	}
      
	  InVec[i][2] = Vec[i][2];
    }
  
  JConVec=ConJConVec;
  JConVec.Transpone();
  
  delete []InVec;
  
  
}
//*                                                                            *
//******************************************************************************
//*             The jacobian of the constraint vector : JConVec : A            *
//******************************************************************************
//*                                                                            
void OptimiZZZ::CVecJacobian (void)
{
  
   
  int NoIneq=NofIneq((RVector*)0,0);
  double Final_Area= (4 * Pi)/(NoFac);
  LoC= NoFac+NoIneq;
  RCoordMatrix ConJConVec(LoC,3*NoVRT);
  JConVec.New(LoC,3*NoVRT);
  double ForX,ForY,ForZ,term1,term2,term3,NORMtx,tgtx;
  
  //First the Area constraints, Normalised !
      
  for(int k=0;k<NoFac;k++)
    {
      for(int w=0;w<4;w++)
	{
	  
	  int NoDA=  surfnet.Faces(k)[(0+w)%4];
	  int NoDB=  surfnet.Faces(k)[(1+w)%4];
	  int NoDC=  surfnet.Faces(k)[(2+w)%4];
	  int NoDD=  surfnet.Faces(k)[(3+w)%4];
	  
	  
	  ForX=DAreaDX(Vec[NoDA],Vec[NoDB],Vec[NoDC],Vec[NoDD]);
	  ForY=DAreaDY(Vec[NoDA],Vec[NoDB],Vec[NoDC],Vec[NoDD]);
	  ForZ=DAreaDZ(Vec[NoDA],Vec[NoDB],Vec[NoDC],Vec[NoDD]);
	  term1 =  ForX * Vec[NoDA][0];
	  term2 =  ForY * Vec[NoDA][1];
	  term3 =  ForZ * Vec[NoDA][2];
	  NORMtx = Vec[NoDA][0]* Vec[NoDA][0] +  Vec[NoDA][1]* Vec[NoDA][1] +  Vec[NoDA][2]* Vec[NoDA][2] ;
	  tgtx = (term1 + term2 + term3) / NORMtx ;
	  

	  ConJConVec(k,NoDA) = ForX ;// Not normalized
	  ConJConVec(k,NoDA + NoVRT) = ForY;;
	  ConJConVec(k,NoDA + 2 * NoVRT) = ForZ;
	  /* 
	  ConJConVec(k,NoDA) = ForX - (tgtx * Vec[NoDA][0]);//Normalised
	  ConJConVec(k,NoDA + NoVRT) = ForY- (tgtx * Vec[NoDA][1]);
	  ConJConVec(k,NoDA + 2 * NoVRT) = ForZ - (tgtx * Vec[NoDA][2]);
	  */
	}   
    }  
 
  // The Euclidean constraints to keep the nodes on the surface of the sphere
  //Not neccesary since will be enforced
  
 

  //The angle constraints, Normalised
 
  if(NoIneq>0) 
    {
      int IFI= NoFac;
      for (int j=0;j<NoFac;j++)
	{
	  for (int ve=0;ve<4;ve++)
	    {
	      int Tv=surfnet.Faces(j)[ve];
	      int Pve=ve-1;
	      int Nve=ve+1;
	      if (ve==3){Pve=2;Nve=0;}
	      if (ve==0){Pve=3;Nve=1;}
	      int Pv =surfnet.Faces(j)[Pve];
	      int Nv =surfnet.Faces(j)[Nve];
	      double angle = LDAB(Vec[Pv],Vec[Tv],Vec[Nv]);
	      double sinangle=sin(angle);
	     
	      if(sinangle < His)
		{
		      
		  ForX=DAx(Vec[Pv],Vec[Tv],Vec[Nv]);//For Pv Vector
		  ForY=DAy(Vec[Pv],Vec[Tv],Vec[Nv]);
		  ForZ=DAz(Vec[Pv],Vec[Tv],Vec[Nv]);

		  term1 =  ForX * Vec[Pv][0];
		  term2 =  ForY * Vec[Pv][1];
		  term3 =  ForZ * Vec[Pv][2];
		 
		  NORMtx = Vec[Pv][0]* Vec[Pv][0] +  Vec[Pv][1]* Vec[Pv][1] +  Vec[Pv][2]* Vec[Pv][2] ;
		  tgtx = (term1 + term2 + term3) / NORMtx ;
		  
		  ConJConVec(IFI,Pv) = ForX;// not Normalised
		  ConJConVec(IFI,Pv+NoVRT) = ForY;
		  ConJConVec(IFI,Pv+2*NoVRT) = ForZ;
		  
		  /*
		  ConJConVec(IFI,Pv) = ForX- (tgtx * Vec[Pv][0]); //Normalised
		  ConJConVec(IFI,Pv+NoVRT) = ForY-  (tgtx * Vec[Pv][1]);
		  ConJConVec(IFI,Pv+2*NoVRT) = ForZ- (tgtx * Vec[Pv][2]);
		  */

		  ForX=DBx(Vec[Pv],Vec[Tv],Vec[Nv]);//For Tv Vector
		  ForY=DBy(Vec[Pv],Vec[Tv],Vec[Nv]);
		  ForZ=DBz(Vec[Pv],Vec[Tv],Vec[Nv]);
		  
		  term1 =  ForX * Vec[Tv][0];
		  term2 =  ForY * Vec[Tv][1];
		  term3 =  ForZ * Vec[Tv][2];
		      
		  NORMtx = Vec[Tv][0]* Vec[Tv][0] +  Vec[Tv][1]* Vec[Tv][1] +  Vec[Tv][2]* Vec[Tv][2] ;
		  tgtx = (term1 + term2 + term3) / NORMtx ;
	
		  ConJConVec(IFI,Tv) = ForX - (tgtx * Vec[Tv][0]);
		  ConJConVec(IFI,Tv+NoVRT) = ForY - (tgtx * Vec[Tv][1]);
		  ConJConVec(IFI,Tv+2*NoVRT) = ForZ - (tgtx * Vec[Tv][2]);

		  ForX=DCx(Vec[Pv],Vec[Tv],Vec[Nv]);//For Nv Vector
		  ForY=DCy(Vec[Pv],Vec[Tv],Vec[Nv]);
		  ForZ=DCz(Vec[Pv],Vec[Tv],Vec[Nv]);

		  term1 =  ForX * Vec[Nv][0];
		  term2 =  ForY * Vec[Nv][1];
		  term3 =  ForZ * Vec[Nv][2];
		 
		  NORMtx = Vec[Nv][0]* Vec[Nv][0] +  Vec[Nv][1]* Vec[Nv][1] +  Vec[Nv][2]* Vec[Nv][2] ;
		  tgtx = (term1 + term2 + term3) / NORMtx ;
	
		  /*
		  ConJConVec(IFI,Nv) = ForX - (tgtx * Vec[Nv][0]);// Normalised
		  ConJConVec(IFI,Nv+NoVRT) = ForY - (tgtx * Vec[Nv][1]);
		  ConJConVec(IFI,Nv+2*NoVRT) = ForZ - (tgtx * Vec[Nv][2]);
		  */
		  ConJConVec(IFI,Nv) = ForX;// not normalised
		  ConJConVec(IFI,Nv+NoVRT) = ForY; 
		  ConJConVec(IFI,Nv+2*NoVRT) = ForZ;
		  IFI++;
		}
	    }
	}
    }
  JConVec=ConJConVec;


}//End here 


//*                                                                            *
//******************************************************************************
//*                            Objective Function : f(x)                       * 
//******************************************************************************
//*                                                                            *

double OptimiZZZ::Fobj_calc(RVector *Vec1,int Len)
{
  
  if (Len==0)
    {
    Fobj_value=0;
    double CoSij=0;

    for(int l=0;l<NoVRT;l++)
      {
	int stelos= diffsphr.NoDirNeighbours(l);
	
	for(int s=0;s<stelos;s++)
	  {
	    int Dnl=diffsphr.DirNeighboursClock(l)[s];
	    
	    CoSij = (Vec[l] & Vec[Dnl]);
	    Fobj_value +=(-0.5)*CoSij;
	    // cerr<<"l="<<l<<" d= "<<Dnl<<"Cos="<<CoSij<<"F="<<Fobj_value<<endl;
	  }
      }
    }
  else
    {
      Fobj_value=0;
      double CoSij=0;
      for(int l=0;l<NoVRT;l++)
	{
	  int stelos= diffsphr.NoDirNeighbours(l);
	  
	  for(int s=0;s<stelos;s++)
	    {
	      int Dnl=diffsphr.DirNeighboursClock(l)[s];
	      CoSij = (Vec1[l] & Vec1[Dnl]);
	      Fobj_value +=(-0.5)*CoSij;
	      
	    }
	}
    }
  //  cerr<<" F = "<<Fobj_value<<endl;
  return Fobj_value;
}//end here

//*                                                                            *
//******************************************************************************
//*               Sents the nodes back on the sphere surface                   *
//******************************************************************************
//*                                                                            *

void OptimiZZZ::Vec_normaZ(RVector *Vec1,int Len)
{ 

  if(Len==0)
    {      
      
      for(int i=0;i<NoVRT;i++)
	{
	  double metro= NormV(Vec[i]);
	  Vec[i][0]=Vec[i][0]/metro;
	  Vec[i][1]=Vec[i][1]/metro;
	  Vec[i][2]=Vec[i][2]/metro;
	}
    }
  else
    {
     for(int i=0;i<NoVRT;i++)
	{
	  double metro= NormV(Vec1[i]);
	  Vec1[i][0]=Vec1[i][0]/metro;
	  Vec1[i][1]=Vec1[i][1]/metro;
	  Vec1[i][2]=Vec1[i][2]/metro;
	 
	}
    }
}//end



//*                                                                            *
//******************************************************************************
//*     The Gradient of the  Objective function  DeltaF                        *
//******************************************************************************
//*                                                                            *

 void OptimiZZZ::CalcDeltaF(void)
{
// The gradient of f(x) DeltaF

  DeltaF.New(3 *NoVRT);
  
  //Calculating the gradient of F
  
  for(int l=0;l<NoVRT;l++) //for all nodes
  {
    int stelos= diffsphr.NoDirNeighbours(l);

    for(int s=0;s<stelos;s++)//for all the neighbours
      {
	int Dnl=diffsphr.DirNeighboursClock(l)[s];
	DeltaF[l]+= (-0.5)*Vec[Dnl][0];
	DeltaF[l+NoVRT]+= (-0.5)*Vec[Dnl][1]; 
	DeltaF[l+ 2*NoVRT]+= (-0.5)*Vec[Dnl][2];
      }
  }

  //Projecting DeltaF orthogonal to Vec

  double metro=0.0;
  for(int i=0;i<NoVRT;i++)
  {
    metro = ( DeltaF[i] * Vec[i][0] + DeltaF[i+NoVRT] * Vec[i][1] + DeltaF[i+2*NoVRT] * Vec[i][2] ) 
      / ( NormV(Vec[i]) * NormV(Vec[i]) );
    DeltaF[i] -= metro * Vec[i][0];
    DeltaF[i+ NoVRT] -= metro * Vec[i][1];	
    DeltaF[i+ 2 * NoVRT] -= metro * Vec[i][2]; 
  }
  
}


//*                                                                            *
// *****************************************************************************
// *       Few Additions to the RCompRowmatrix and RCoordMatrix                *
// *****************************************************************************
//*                                                                            *

//** Use of RCompRowMatrix from Toast to create the RCompMatrix ATA  from the RCompmatrix A
//** The increase in speed is significant 

class myRCompRowMatrix: public RCompRowMatrix
{
public:
  myRCompRowMatrix(const RCoordMatrix &cm): RCompRowMatrix (cm) {}
  myRCompRowMatrix(const RCompRowMatrix  &cm): RCompRowMatrix (cm) {}
  // myRCompRowMatrix(int rows, int cols): RCompRowMatrix (rows, cols) {}
  RCoordMatrix AttA ();
  RVector DiagATA ();
};


RVector  myRCompRowMatrix::DiagATA ()
{
  //cerr << "start AttA" << endl;
  myRCompRowMatrix At(*this);
  At.Transpone();
  RVector Diag(At.nRows());

 

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




RCoordMatrix myRCompRowMatrix::AttA ()
{
  //  cerr << "start AttA" << endl;
  myRCompRowMatrix At(*this);
  At.Transpone();

  int i, j,k;//, n = nCols();
 
  RCoordMatrix ata(At.nRows(), At.nRows());
  for (i = 0; i < At.nRows(); i++)
    {
    for (j = 0; j < At.nRows(); j++)
      {
      for(k=0;k<At.nCols();k++)
	{
	  if(At.Get(i,k)!=0 && At.Get(j,k)!=0)
	    {
	      ata(i,j)+=At(i,k)*At(j,k);
	    }
	}
      }
    }
  return ata;
}



//*                                                                            *
//******************************************************************************
//*     minimisation of the Objective function L(x,lamda)=F(x)-lamda.C(x)      *
//******************************************************************************
//*                                                                            *

void OptimiZZZ::F_min(void)
{
  
  double F1= Fobj_calc((RVector *) 0,0); // the value of the objective function at the beginning

  
  CalcDeltaF();// the Gradient of the objective Function

   
  CVecJacobian ();// The jacobian of the constraint vector

  RVector ATDeltaF(LoC);

 //Calculate AT.DeltaF

  for(int i=0;i<LoC;i++)
  {
    for (int f=0;f<3*NoVRT;f++)
      {
   	 if(JConVec.Get(i,f)!=0)  ATDeltaF[i]+= DeltaF[f]*JConVec(i,f);
      }
  }


  JConVec.Transpone();

  myRCompRowMatrix myA (JConVec);// Create myRCompRowMatrix myA for the calculation of the ATA

  
  RVector diagAtA = myA.DiagATA();// using the diagonal of the ATA only
   
  RVector m_lamda(LoC),dx(3*NoVRT);//dv and dx
  for(int g=0;g<LoC;g++) m_lamda[g]=0.0;// set to 0.0


  double dcxTemp;
  for(int h=0;h<LoC;h++)
    {
      if(diagAtA[h]!=0)
	{
	  dcxTemp=(-1.0)* ATDeltaF[h] / diagAtA[h];//-????
	  m_lamda[h]=dcxTemp;
	}
      else
	{
	   m_lamda[h]=0;
	}
    }

 JConVec.Transpone();


  //Calculate  A.m_lamda 

  RVector Alamda(3*NoVRT);

  for (int f=0;f<3*NoVRT;f++)
  {
    for(int i=0;i<LoC;i++)
    {
   
      if(JConVec.Get(i,f)!=0)  Alamda[f]+= m_lamda[i]*JConVec(i,f);
    }
  }

  //Caclulate Gz=DeltaF-A.m_lamda

  for(int t=0;t<3*NoVRT;t++)
  {
    GFstore2[t]=DeltaF[t]-Alamda[t];
    
  }






 
  

  //Refreshing gamma : applies only if iteration counter >1

  if (Iter>1)
  { 
    double dotGF2=0.0;
    double dotGF1=0.0;
    for(int k=0;k<(3*NoVRT);k++)
      {
	dotGF1 += ( ( GFstore2[k] - GFstore1[k] ) * GFstore2[k]);
	dotGF2 += ( GFstore2[k] * GFstore2[k] );
      }
    gamma=dotGF1/dotGF2;
  }
  
  // The Conjugate Gradient scheme
  
  double res;
  for(int kk=0;kk<(3*NoVRT);kk++)
  {
    res= gamma * Pstore1[kk];
    Pstore2[kk] = GFstore2[kk]+ res;
  }
  
  //Create the search direction P
  
  double P_cost=0.0;  
  
  for(int kkk=0;kkk<(3*NoVRT);kkk++)
  {
    P[kkk]=Pstore2[kkk]; 
    P_cost+= P[kkk]* P[kkk];
  }
  MonitorP[Iter]= P_cost;
  
  
  //Storing for the next step of conjugate gradient

  for(int klk=0;klk<(3*NoVRT);klk++)
  {
    Pstore1[klk]=Pstore2[klk];
    GFstore1[klk]=GFstore2[klk];
  }
   
  //The penalty factor
  
  double penalty_limit=0.8;
  double Penalty;
  double Min_Cost=0.0;
  if(Iter>1)
  {
    Min_Cost=MonitorCost[0];
    for(int g=0;g<Iter;g++)
      {
	if(MonitorCost[g]<Min_Cost)  Min_Cost=MonitorCost[g];
      }
  }
  double CurrentCost= CVecConstruction((RVector *) 0,0);
  double Increase;
  Increase= CurrentCost-Min_Cost;
  if (Increase<0.001){ Increase=CurrentCost;}
  
  Penalty =(1.5 *( CurrentCost / Increase))+ penalty_limit;

  //cout<<"Penalty factor "<<Penalty<<endl;
  












  // Line Search

  double step;
  step=0.5 + steplarge;
  double a=0.0;
  int GoBack=1;
 
 
  RVector *CentralVec=0;
  CentralVec = new RVector[NoVRT];
  RVector *LeftVec=0;
  LeftVec = new RVector[NoVRT];
  RVector *RightVec=0;
  RightVec = new RVector[NoVRT];
 
  for(int j=0;j<NoVRT;j++)
  {
    CentralVec[j] =RVector(3);
    LeftVec[j] =RVector(3);
    RightVec[j] =RVector(3);
  }
 
  while(GoBack==1)
  {

     std::cout<< "(((((((((  a=  "<<a<<"   step = "<<step<<" steplarge= "<<steplarge<< std::endl;


    for(int j=0;j<NoVRT;j++)
    {
      CentralVec[j][0]=Vec[j][0] + a * P[j];
      LeftVec[j][0]=Vec[j][0]+(a - step)* P[j];
      RightVec[j][0]=Vec[j][0] + (a + step)* P[j];
      
      CentralVec[j][1]=Vec[j][1]+ a * P[j+NoVRT];
      LeftVec[j][1]=Vec[j][1]+(a-step)* P[j+NoVRT];
      RightVec[j][1]=Vec[j][1]+(a+step)* P[j+NoVRT];
      
      CentralVec[j][2]=Vec[j][2]+ a * P[j+2*NoVRT];
      LeftVec[j][2]=Vec[j][2]+(a-step)* P[j+2*NoVRT];
      RightVec[j][2]=Vec[j][2]+(a+step)* P[j+2*NoVRT];
    }
 
    // Normalisation of vectors to get on the sphere surface
    
    
    for(int i=0;i<NoVRT;i++)
    {
      double metro= NormV(CentralVec[i]);
      CentralVec[i][0]=CentralVec[i][0]/metro;
      CentralVec[i][1]=CentralVec[i][1]/metro;
      CentralVec[i][2]=CentralVec[i][2]/metro;
    }
    for(int i=0;i<NoVRT;i++)
    {
      double metro= NormV(RightVec[i]);
      RightVec[i][0]=RightVec[i][0]/metro;
      RightVec[i][1]=RightVec[i][1]/metro;
      RightVec[i][2]=RightVec[i][2]/metro;
    }
    for(int i=0;i<NoVRT;i++)
    {
      double metro= NormV(LeftVec[i]);
      LeftVec[i][0]=LeftVec[i][0]/metro;
      LeftVec[i][1]=LeftVec[i][1]/metro;
      LeftVec[i][2]=LeftVec[i][2]/metro;
    }

   
    
    //Evaluate the Lagrangian for the L R and C
    
    double CostCentral=CVecConstruction (CentralVec,NoVRT);
    double LaCcentral=0.0;
    for(int f=0;f<NoFac;f++)
      {
	LaCcentral += dabs(m_lamda[f]* ConVec[f]);
      }
    double Fcentral =Fobj_calc(CentralVec,NoVRT);
    double Lcentral =  -Fcentral + LaCcentral + Penalty * CostCentral;
    
    double CostRight=CVecConstruction (RightVec,NoVRT);
    double LaCright=0.0;
    for(int f=0;f<NoFac;f++)
       {
	 LaCright += dabs(m_lamda[f]* ConVec[f]);
       }
    double Fright =  Fobj_calc(RightVec,NoVRT);
    double Lright =  -Fright + LaCright + Penalty * CostRight;
     
    double CostLeft=CVecConstruction (LeftVec,NoVRT);
    double LaCleft=0.0;
    for(int f=0;f<NoFac;f++)
      {
	LaCleft+=dabs(m_lamda[f]* ConVec[f]);
      }
    double Fleft =  Fobj_calc(LeftVec,NoVRT);
    double Lleft=   -Fleft + LaCleft+ Penalty *CostLeft;
      
   
    int Next1Right=0;
    int Next1Left=0;
    
    //checking the active inequalities
    if((NofIneq(CentralVec,NoVRT)+ ToloIneq)> NofIneq(RightVec,NoVRT)) Next1Right=1;
    if((NofIneq(CentralVec,NoVRT)+ ToloIneq)> NofIneq(LeftVec,NoVRT)) Next1Left=1;
    
    // printing for monitoring
     std::cout<<"+------------------------------------------------------------------------------"<< std::endl;
     std::cout<<"|  Lan  <-: "<< Lleft<<"     :-|-: "<<Lcentral<<"      :-> "<<Lright<< std::endl;
     std::cout<<"|  LaC  <-: "<< LaCleft<<"   :-|-: "<<LaCcentral<<"    :-> "<<LaCright<< std::endl;
     std::cout<<"|   C   <-: "<< CostLeft<<"  :-|-: "<<CostCentral<<"     :-> "<<CostRight<< std::endl;
     std::cout<<"|   F   <-: "<< Fleft<<"  ```:-|-: "<<Fcentral<<"   `  :-> "<<Fright<< std::endl;

 
    
    
    //Find the minimun Lagrangian
    
    int Next2;
    if(Lcentral<Lright && Lcentral <Lleft){Next2=0;}
    if(Lcentral>Lright && Lright <Lleft){Next2=1;}
    if(Lleft < Lright && Lcentral >Lleft){Next2=-1;}

    //Find the minimun objective f
    /* 
   int Next2;
    if(Fcentral<Fright && Fcentral <Fleft){Next2=0;}
    if(Fcentral>Fright && Fright <Fleft){Next2=1;}
    if(Fleft < Fright && Fcentral >Fleft){Next2=-1;}
 */
    
    // making the decision and moving a
   
    
    if(Next1Right==1 && Next2==1){a = a + step;  std::cout<<"---- Right ----"<< std::endl; ToloIneq=0;}
    if(Next1Left==1 && Next2==-1){a = a - step;  std::cout<<"---- Left ----"<< std::endl; ToloIneq=0;}
    
    //Step halved unless the line tollarance is reached
    
     step = step * 0.5;
    
     if(step < linetol)
    {
      fail=1;
      GoBack = 0;
      for ( int g = 0; g < NoVRT; g++ )
	{
	  Vec[g][0] =CentralVec[g][0];// Vec[g][0] + a * P[g];
	  Vec[g][1] =CentralVec[g][1];// Vec[g][1] + a * P[g];
	  Vec[g][2] =CentralVec[g][2];// Vec[g][2] + a * P[g];
	}
    }
  
  }//End of while loop
  
  
  //Normalise Vec
  
  Vec_normaZ((RVector*)0,0);
 
 
  double F2= Fobj_calc((RVector*)0,0);
 
  if(dabs(F1-F2)<0.000001)
  {
    steplarge +=1.6;
    ToloIneq=1;
     std::cout<<" Failed to improve 8-( "<< std::endl;
    fail=1;
  
  } 	
  else
  {
    steplarge=0;

  } 
  
  delete []CentralVec;
  delete []LeftVec;
  delete []RightVec;
  
}//End of F_min function



//*                                                                            *
//******************************************************************************
//*                   Springs instead of Objective function                    *
//******************************************************************************
//*                                                                            *

void OptimiZZZ::Springs(void)
{

  double normf;
  double eps=1.0;
  RVector sumZ(NoVRT),sumY(NoVRT),sumX(NoVRT);

  for ( int g = 0; g < NoVRT; g++ )
    {

      for ( int j = 0; j < NoVRT; j++ )
	{
	if(j!=g)
	  {
	    normf= NormS((Vec[g][0]-Vec[j][0]),(Vec[g][1]-Vec[j][1]),(Vec[g][2]-Vec[j][2]));
	    sumX[g]+=( Vec[g][0]-Vec[j][0])/normf;
	    sumY[g]+= (Vec[g][1]-Vec[j][1])/normf;
	    sumZ[g]+=( Vec[g][2]-Vec[j][2])/normf;
	  }
	}
      sumX[g]= eps * sumX[g]+ Vec[g][0];
      sumY[g]= eps * sumY[g]+ Vec[g][1];
      sumZ[g]= eps * sumZ[g]+ Vec[g][2];

    }


  for ( int gg = 0; gg < NoVRT; gg++ )
    {
      Vec[gg][0]= sumX[gg];
      Vec[gg][1]= sumY[gg];
      Vec[gg][2]= sumZ[gg];
    }


 Vec_normaZ((RVector*)0,0);
 


}//end of springs

//*                                                                            *
//******************************************************************************
//*                     Minimisation of the constraint vector                  *
//******************************************************************************
//*                                                                            *




void OptimiZZZ::Con_Min(void)
{
  
  double Con_tol=0.1;

  CVecJacobian();// Create the jacobian of the constraint vector ConVec
 
  int NiQ1=NofIneq ((RVector*)0,0);// number of inequalities at the beginning 

  double CCost1= CVecConstruction((RVector*)0,0);//The cost for the constraints at the beginning
 
  RVector Cv1(LoC);//create the Cv1 to map the constraint vector at the beginning
  
  int  LoC_local=LoC;
  

  Cv1=ConVec;//Cv1 gets the values from the constraint vector ConVec
    
 
  JConVec.Transpone();

  myRCompRowMatrix myA (JConVec);// Create myRCompRowMatrix myA for the calculation of the ATA

  // RCompRowMatrix AtA = myA.AttA();// calculate ATA ->AtA
  RVector diagAtA = myA.DiagATA();// using the diagonal of the ATA only
   
  RVector dv(LoC_local),dx(3*NoVRT);//dv and dx
  for(int g=0;g<LoC_local;g++) dv[g]=0.0;// set to 0.0

//Calculate dx step using the conjugate gradient solver from Toast
//  double tolcg=0.008;//tolerance for the solver
// CG(AtA, dv, Cv2, tolcg);//The solver

  double dcxTemp;
  for(int h=0;h<LoC_local;h++)
    {
      if(diagAtA[h]!=0)
	{
	  dcxTemp=(-1.0)* Cv1[h] / diagAtA[h];//-????
	  dv[h]=dcxTemp;
	}
      else
	{
	  dv[h]=0;
	}
    }



 JConVec.Ax(dv,dx);// calculate dx

 JConVec.Transpone();




 //********************************************
 //search for the optimal along th dx direction
 //******************************************** 
    
    RVector *XVec=0;// XVec will store the positions of the nodes after each step
    XVec = new RVector[NoVRT];
    for(int j=0;j<NoVRT;j++)	XVec[j] =RVector(3);
   

    int Goback=1;// For the while loop 
    int iN=0;// counter 
  
    while(Goback==1)
      {
	iN++;//+1 every time it halves the step
	
	// Apply the step dx 
	for(int f=0;f<NoVRT;f++)
	  {
	    XVec[f][0] = Vec[f][0] + dx[f];
	    XVec[f][1] = Vec[f][1] + dx[f+ NoVRT];
	    XVec[f][2] = Vec[f][2] + dx[f+ 2 * NoVRT];
	  }
		
	//Projection on the Sphere surface 
	for(int i=0;i<NoVRT;i++)
	  {
	    double metro= NormV(XVec[i]);
	    XVec[i][0]=XVec[i][0]/metro;
	    XVec[i][1]=XVec[i][1]/metro;
	    XVec[i][2]=XVec[i][2]/metro;
	  }


	double CCost2= CVecConstruction (XVec,NoVRT);	// The cost after the step
	int NiQ2=NofIneq(XVec,NoVRT);// number of active inequalities after the step
	RVector	Cv2=ConVec; // The constraint vector after the step

	//Monitoring

	 std::cout <<" | "<< iN <<" |  C= " <<CCost2<<" |  inEq = "<<NiQ2<<" |  F = "<<Fobj_calc(XVec,NoVRT)<< std::endl;
	
	int GoingBad = 0;

	for(int t=0;t<NoFac;t++)
	  {
	    if( (dabs(Cv2[t]) - dabs(Cv1[t])) > Con_tol) GoingBad = 1;//dabs->absolute value
	  }
	

	//  Checks the step to make sure it reduces the constraints and does not increases the inequalities
       	int Cineq_tol=0;
	if(fail==1 && Cineq_tol < 3  )Cineq_tol += 1;

	if(NiQ2>NiQ1+Cineq_tol || CCost2>=CCost1 || GoingBad==1)

	  {
	    Goback=1;//Go back 

	    double dxdx=0.0;

	    for(int r=0;r<3*NoVRT;r++) // Halve the step !!
	      {
		dx[r] -= dx[r]/3; 
		dxdx += dx[r]*dx[r];
	      }
	    //    cout<<"dx"<<dxdx<<endl;
  	    if(dxdx<linetol)//Checks that the step is not smaller than is useful
	      {
		Goback = 0;
		 std::cout<<" Step too small  "<< std::endl;
		fail=1;
	      }

	  }

	else   Goback=0;//finished	   


	  }//end while loop 


	for(int f=0;f<NoVRT;f++)//store the new position of the nodes
	  {
	    Vec[f][0]=XVec[f][0];
	    Vec[f][1]=XVec[f][1];
	    Vec[f][2]=XVec[f][2];
	  } 
   
    

    delete []XVec;
    
}// End of C_Min

