#define SHARMLIB_IMPLEMENTATION

//#include <stdlib.h>
//#include <string.h>

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "mathlib.h"
#include "felib.h"
#include "sharmlib.h"
#include "diffusphere.h"
using namespace std;


diffusphere::diffusphere()
{
}
diffusphere::~diffusphere() 
{
}
 

void diffusphere::DoDif (IDenseMatrix &Faces, RDenseMatrix &Vertices, RVector &Thita, RVector &Fi)
{
  
  NoVRT = Vertices.nRows();               // Number of Elements in the  region
  NoFac = Faces.nRows();                  // Number of Nodes in the  region
 

  RCoordMatrix A (NoVRT,NoVRT);

  
 
  // Setting up the poles

  A(0,0)=1;//South Pole
  A(NoVRT-1,NoVRT-1)=1;//North Pole

  cerr<<'01'<<endl;

 
  // From the faces find the direct neighbours


  for(int v=1;v< NoVRT-1;v++)
    {
		
      for(int j=0;j<NoFac;j++)
	{
	
	  if(Faces(j,0) == v)
	    {
	      A(v,Faces(j,1))=-1;
	      A(v,Faces(j,3))=-1;
	    }

	  if(Faces(j,1) == v)
	    {
	      A(v,Faces(j,0))=-1;
	      A(v,Faces(j,2))=-1;
	    }

	  if(Faces(j,2) == v)
	    {
	      A(v,Faces(j,1))=-1;
	      A(v,Faces(j,3))=-1;
	    }

	  if(Faces(j,3) == v)
	    {
	      A(v,Faces(j,2))=-1;
	      A(v,Faces(j,0))=-1;
	    }
	}
    }

  cerr<<'02'<<endl;

  RCompRowMatrix A1(A);

  // The number of direct neighbours in the diagonal

  for(int  nn=1 ;nn < NoVRT ; nn++)
    {
      double  diag=0;
      RVector The_nn_Row=A1.Row(nn); 
      for(int g=0;g<NoVRT ; g++)
	{
	  // diag -= A(nn,g);
	  diag -=The_nn_Row[g];
	}
      

      A(nn,nn)=diag;
    }

  cerr<<'03'<<endl;

  //Printing the A matrix for debuging
  /* for(int n=0;n<NoVRT;n++)
     {
      cerr<<endl;
      for(int g=0;g<NoVRT;g++)
	{
	  cerr<<"   " <<A(n,g);
	}
    }
    cerr<<endl;
*/
   

   //Creating the reduced A Matrix



  RCompRowMatrix  ReduceA(A);
 
  //ReduceA.RemoveRow(0);
  //RVector 
	  zerv.New(ReduceA.nRows());
  ReduceA.SetRow (0, zerv);

   cerr<<'04'<<endl;
  //ReduceA.RemoveRow(NoVRT-2);
 // RVector zerv(ReduceA.nRows());
	  zerv.New(ReduceA.nRows());
  ReduceA.SetRow (NoVRT-2, zerv);
 
  ReduceA. Transpone();
    cerr<<'05'<<endl;
  //ReduceA.RemoveRow(0);
  //RVector zerv(ReduceA.nRows());
	  zerv.New(ReduceA.nRows());
  ReduceA.SetRow (0, zerv);

  //ReduceA.RemoveRow(NoVRT-2);
 // RVector zerv(ReduceA.nRows());
	  zerv.New(ReduceA.nRows());
  ReduceA.SetRow (NoVRT-2, zerv);
 
  ReduceA. Transpone();
  
  cerr<<'06'<<endl;



 

   //Printing the AReduced matrix for debuging
 /* 
    for(int n=0;n<NoVRT-2;n++)
     {
       cerr<<endl;
       for(int g=0;g<NoVRT-2;g++)
	 {
	   cerr<<"   " <<ReduceA(n,g);
	 }
     }*/
   
 
   //Creating the Thita vector and the B constand vector

   RVector RedThita(NoVRT-2),BConVeCTR(NoVRT-2);
   Thita.New (NoVRT);

   //filling the b constand vector

   for(int d=1;d<NoVRT-1;d++)
     {
       if (A(d,0)== -1)BConVeCTR[d-1]=0;
       if (A(d,NoVRT-1)== -1)BConVeCTR[d-1]=M_PI;
     } 

   cerr<<'07'<<endl;

   double tol=1.0e-16;
//   cout<<"  Calculating Thita  "<<endl;
   CG( ReduceA, BConVeCTR, RedThita,tol);
 //  cout<<"   Thita Done  "<<endl;
   // delete B;

    cerr<<'08'<<endl;

 // Create  the vector Thita

  Thita[0]=0;
  Thita[NoVRT-1]=M_PI;

   for(int k=1;k<NoVRT-1;k++)
     {
       Thita[k]=RedThita[k-1];
     }

 
    cerr<<'09'<<endl;

 
   //##########################################################################
   //          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
   //                                The Longitude
   //          ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   //##########################################################################

    RCoordMatrix  RedA=ReduceA.MakeCoordMatrix();

   //modify reduced A for the longitude calculation

   for(int n=0;n<NoVRT-2;n++)
     {
       double diagon=0;
       RedA(n,n)=0;
  
      RVector The_n_Row=ReduceA.Row(n);
      The_n_Row[n]=0;
     
      for(int g=0;g<NoVRT-2 ; g++)
	{
	  diagon -=The_n_Row[g];
	}
     
     RedA(n,n)=diagon;
     }
   RedA(0,0)+=2;


  



    cerr<<'10'<<endl;






//Printing the AReduced matrix for debuging
/* 
   
   for(int n=0;n<NoVRT-2;n++)
     {
     cerr<<endl;
     for(int g=0;g<NoVRT-2;g++)
	 {
	   cerr <<RedA(n,g)<<"   ";
	 }
     }
*/






   //Set Up the Reduced Fi Vector

   RVector RedFi(NoVRT-2);
   Fi.New(NoVRT);
   //double RedBlong[NoVRT-2];
	 RVector RedBlong(NoVRT-2);

  
     cerr<<'11'<<endl;


  //  Find the direct neighbours anticlockwise from the faces 
   
   DirNBR = new IVector[NoVRT];
   for(int i=0;i<NoVRT;i++) DirNBR[i].New(6);


   DirNBRClock = new IVector[NoVRT];
   for(int i=0;i<NoVRT;i++) DirNBRClock[i].New(6);
   NoDNBR= new int[NoVRT];
   int l,CurFace[6][4],ThisNGBR, NextNGBR,s;
  
   //double  Th1[NoVRT];
   RVector Th1(NoVRT);
   int NorthPole=0;
   int SouthPole=NoVRT-1;
   int PrevNode,t, neighbors_ID;
   int here,NextPos,PrevPos,kk;
   double maxTh;


   
  cerr<<'12'<<endl;
  
   //initialise DirNBR an DirNBRClock

   for(int h=0;h<NoVRT-2;h++)
     {
       RedBlong[h]=0;
     }

   for(int g=0;g<NoVRT;g++)
     {
       for(int k=0;k<6;k++)
	 {
	   DirNBRClock[g][k]=-1;
	   DirNBR[g][k]=-1;
	   
	 }
     }

    cerr<<'13'<<endl;

   for(int v=0;v<NoVRT;v++)// For all the vertices
     {
       //cerr<<v<<endl;
       l=0;//initialisation
      
       for(int in=0;in<6;in++){// up to six faces for each vertex
	 for(int on=0;on<4;on++){//four vertices  define each of the six faces
	   CurFace[in][on]=-1;//initialising CurFace to accept the temporary faces 
	 }
       }


       for(int j=0;j<NoFac;j++)//For all the faces 
	 {
	   for(int g=0;g<4;g++) //For the four vertices  on each face
	     { 
	       if(Faces(j,g)==v)// If vertex v participates in face j ...
		 {
		   // cerr<<" o vetrex :"<<v<<" einai sto "<<j <<endl;
		   CurFace[l][0]= Faces(j,0);// Temporary table of patricipating faces
		   CurFace[l][1]= Faces(j,1);
		   CurFace[l][2]= Faces(j,2);
		   CurFace[l][3]= Faces(j,3);
		   l++;
		   
		 } 
	     }
	 }

  cerr<<'14'<<endl;
 

       for(int n=0;n<4;n++) //The First Two Neighbours
	 {
 

	   if(CurFace[0][n]==v)
	     {
	       ThisNGBR=CurFace[0][(n+1)%4];
	       NextNGBR=CurFace[0][(n+3)%4]; 
	       CurFace[0][0]=-1;
	       CurFace[0][1]=-1;
	       CurFace[0][2]=-1;
	       CurFace[0][3]=-1;
	      
	       } 
	  
	 }
	   DirNBR[v][0]= NextNGBR;
	   s=1;


  cerr<<'15'<<endl;

	   while(ThisNGBR != NextNGBR)//loops until it returns to the first one
	     {
	       for (int i=0;i<6;i++)
		 {
		   for(int jf=0;jf<4;jf++)
		     { 
		     
		      
		       if(CurFace[i][jf] == NextNGBR)
			 {
			 
			  

			   kk=(jf+2)%4;
			  
			   NextNGBR = CurFace[i][kk];
		
			   DirNBR[v][s++]=NextNGBR;
			  
			   NoDNBR[v]=s-1;
			  
			   CurFace[i][0]=-1;
			   CurFace[i][1]=-1;
			   CurFace[i][2]=-1;
			   CurFace[i][3]=-1; 
			   //  cerr<<" mesa "<<kk<<" Next Ngbr " << NextNGBR<<endl;
			 }
		     }
		 }
	     }
 	 
     }
 
  
    cerr<<'16'<<endl;

 //Direct Neighbours clockwise
   
   int rh;

   for(int g=0;g<NoVRT;g++)
     {
      
       for(int h=0;h<NoDNBR[g]+1;h++)
	 
	 {
  	   rh=NoDNBR[g]-h;  
	   DirNBRClock[g][h]= DirNBR[g][rh];
	  
	 }
     }

  
   cerr<<'17'<<endl;
 
   //Printing the matrix 
   /*
  
   for(int n=0;n<NoVRT;n++)
     {
       cerr<<endl;
       for(int g=0;g<6;g++)
	 {
	   cerr <<DirNBRClock[n][g] <<"   ";
	 }
     }
  
   */










       // Setup Bvector for longitude
       
       NorthPole=0;

       SouthPole=NoVRT-1;

       for(int i=0;i<NoVRT;i++)
	 {
	   Th1[i]=Thita[i];
	   
	   
	 }

     ;
     
       //Initialisation

       PrevNode= NorthPole;
       here=1;
       maxTh=0;
       PrevPos=0;
      
       while(here != SouthPole)
	 {
	   
	  
	   for(int w=0;w<NoDNBR[here]+1;w++)// for every node in the neighbourhood of here
	     {
	       neighbors_ID=DirNBR[here][w];
	       //      cerr<<"neighbors_ID "<<neighbors_ID<<"w= "<<w<<endl;

	       if(Th1[neighbors_ID]>maxTh)
		 {
		   // cout <<maxTh<<" is less than "<< Th1[neighbors_ID]<<" w= "<<w<<" here is "<<here <<endl;
		   maxTh=Th1[neighbors_ID];
		   NextPos=w;
		   //	   cerr<<"NextPos "<<NextPos<<endl;
		 }
	       
	       if(neighbors_ID == PrevNode)
		 {
	
		   PrevPos=w;
		   //	   cerr<<"PrevPos "<<PrevPos<<endl;
		 }
	       
	     }
	 

	   // All the direct Neighbours clockwise from PrevPos to NextPos

	   if(PrevPos < NextPos)

	     {
	       //    cerr<<"1     :  PrevPos  "<<PrevPos<<"  NextPos  " <<NextPos<<endl;
	       for(int ii=PrevPos+1;ii<NextPos;ii++)
		 {
		   // cerr <<" sto prwto ii=  "<<ii<<endl; 
		   //ii=ii%4;
		   //cerr <<" sto prwto ii2 =  "<<ii<<endl; 

		   RedBlong[DirNBR[here][ii]-1 ] += 2* M_PI;// add 2 Pi to b[neighbour]
		   RedBlong[here-1] -= 2 * M_PI ;//Subtrack2Pi from b[here]
		 }
	       
	     }  
	  

	   //cerr<<"2out     :  PrevPos  "<<PrevPos<<"  NextPos  " <<NextPos<<endl;

	   
	   if(PrevPos > NextPos)

	     {
	       // cerr<<"2    :  PrevPos  "<<PrevPos<<"  NextPos  " <<NextPos<<"  n odn "<<NoDNBR[here]<<endl;

	       for(int iii=PrevPos+1;iii<(NoDNBR[here]+NextPos+1);iii++)
		 {

		   //  cerr <<" sto deutero  iii=  "<<iii<<endl; 

		   t=iii% (NoDNBR[here]+1);
		   // cerr <<" sto deutero  t=  "<<t<<" <DirNBRClock[here][t]=   "<<DirNBR[here][t]<<"here   "<<here<<endl; 
		   RedBlong[(DirNBR[here][t]-1)] += 2 * M_PI;// add 2 Pi to b[neighbour]
		   RedBlong[here-1] -= 2 *  M_PI;//Subtrack2Pi from b[here]
		 }
	     }
	  
	  
	   

	   PrevNode=here;
	 

	   here=DirNBR[here][NextPos];
	  
	 }
     
      
     
       //The vector  b for longitude
       RVector VectorBLong(NoVRT-2);
    
        for (int s=0;s<NoVRT-2;s++)
	  {
	    VectorBLong[s]=RedBlong[s];
	    
	  }
         cerr<<'19'<<endl;


	//The solution with congugate gradiant
	//cout<<"   Calculating Fi  "<<endl;
	CG(RedA,VectorBLong, RedFi,tol);
	//cout<<"    Fi done  "<<endl;
	Fi[0]=0;
	Fi[NoVRT-1]=0;
	for (int s=0;s<NoVRT-2;s++)
	  {
	    Fi[s+1]=RedFi[s];
	  }

  cerr<<'20'<<endl;
}//end of Function 
