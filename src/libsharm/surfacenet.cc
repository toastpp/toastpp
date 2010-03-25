
//============================================================//

#define SHARMLIB_IMPLEMENTATION

#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>
#include "felib.h"
#include "sharmlib.h"
#include "surfacenet.h" 
using namespace toast;


SurfaceNet::SurfaceNet()
{
  VertexList =0;
  AntiClocFaces=0;
  
  NeWNoVRT=0;
  LengoFaces=0;
}


SurfaceNet::~SurfaceNet() 
{
  if (VertexList) delete []VertexList ;
  if (AntiClocFaces) delete [] AntiClocFaces ;
  //if (voxels) delete [] voxels;
}


//============================================================//

/*
istream& operator>>  (istream& i, SurfaceNet& sn)//  ?
{
  
  i>> sn.NeWNoVRT;
  
  sn.VertexList = new Point[sn.NeWNoVRT];
  for(int j=0;j<sn.NeWNoVRT;j++)
    { 
      
      Point x(3);
      i>>x;
      
      sn.VertexList[j]=x;
    }
  
  i>> sn.LengoFaces;
  
  sn.AntiClocFaces = new IVector[sn.LengoFaces];
  for(int j=0;j<sn.LengoFaces;j++)
    {
      IVector x(4);
      i>>x;
      sn.AntiClocFaces[j]=x;
    }
  
  return i;
}*/
//============================================================//

/* ostream& operator<< (ostream& o, SurfaceNet& sn)//   ?
{
  o<< sn.NeWNoVRT<<endl;;
  for(int j=0;j<sn.NeWNoVRT;j++)o<< sn.VertexList[j]<<endl;
  o<< sn.LengoFaces<<endl;
  for(int j=0;j<sn.LengoFaces;j++)o<<sn.AntiClocFaces[j]<<endl;
  return o;
}
*/
//============================================================//
/*
void SurfaceNet::WriteVoxels(char *ofname)
{
  ofstream is(ofname);
  
  int NoVo= NeWNoVRT;
  
  is<<NoVo<<endl;
  
  for(int i = 0;i< NoVo ;i++)
    { 
      is<< VertexList[i] <<endl;
    }
  
}
*/
//============================================================//
/*
void SurfaceNet::WriteFaces(char *ofname)
{
  ofstream outFac(ofname);
  
  int NoFa= LengoFaces;
  
  outFac<<NoFa<<endl;
  
  for(int i = 0;i< NoFa ;i++)
    { 
      outFac<< AntiClocFaces[i] <<endl;
      
    }
}
*/

//============================================================//

//void SurfaceNet::ReadVoxels(char *fname)
//{
//  ifstream is(fname);
 // int NoVoxels;// number of voxels
 // is >> NoVoxels;
 // Point voxels[NoVoxels];
  
 // for (int i=0;i<NoVoxels;i++)// reads the voxel data from a file 
 //   {
 //     is >> voxels[i];
//    } 
 

void SurfaceNet::Vox2surf(IDenseMatrix &voxels, IDenseMatrix &Faces, RDenseMatrix &Vertex);
{
  // .................................................................... 
  // finds the borders of the voxel space
  // ....................................................................
  
  
  int maxX =  (int)voxels[0][0];
  int maxY =  (int)voxels[0][1];
  int maxZ =  (int)voxels[0][2];
  
  for (int i=1;i<NoVoxels;i++)
    {
      if(voxels[i][0]>maxX)maxX =  (int) voxels[i][0];
      if(voxels[i][1]>maxY)maxY =  (int) voxels[i][1];
      if(voxels[i][2]>maxZ)maxZ =  (int) voxels[i][2];
    } 
  
  
  maxX+=2;
  maxY+=2;
  maxZ+=2;
  
  
  cout<<NoVoxels<<"  "<<maxX<<" "<<maxY<<" "<<maxZ<<endl;
  
  int bitmap[maxX][maxY][maxZ];
  int bitmap2[maxX][maxY][maxZ];
  
  
  for(int i=0;i<maxX;i++)
    {
      for(int j=0;j<maxY;j++)
	{
	  for(int k=0;k<maxZ;k++)
	    {
	      if(inpointlist(voxels,NoVoxels,i,j,k))
		{
		  bitmap[i][j][k] =1;
		}
	      else
		{
		  bitmap[i][j][k] =0;
		}
	      bitmap2[i][j][k] =0;
	    }
	}
    }
  
  
  /*Printing bitmap for debuging  */
  /*
    for(int k=0;k<maxZ;k++)
    {
    cout<<endl<<"z= "<<k; 
    for(int i=0;i<maxY;i++)
    {
    cout<<endl;
    for(int j=0;j<maxX;j++)
    {
    cout<<bitmap[j][i][k];
    }
    } 
    } 
  */
  
  // ....................................................................
  // fills in the rest of the voxel vertices and creates bitmap2
  // ....................................................................
  
  for(int i=0;i<maxX;i++)
    {
      for(int j=0;j<maxY;j++)
	{
	  for(int k=0;k<maxZ;k++)
	    {
	      if( bitmap[i][j][k] ==1)
		{
		  // cout<<"jas  "<<i<<j<<k<<endl;
		  bitmap2[i][j][k]=1;
		  bitmap2[i+1][j][k]=1;
		  bitmap2[i+1][j+1][k]=1;
		  bitmap2[i][j+1][k]=1;
		  bitmap2[i][j][k+1]=1;
		  bitmap2[i+1][j][k+1]=1;
		  bitmap2[i][j+1][k+1]=1;
		  bitmap2[i+1][j+1][k+1]=1;
		}
	    }
	}
    }
  
  
  
  
  
  
  // ....................................................................
  //Printing bitmap2 for debuging 
  // ....................................................................
  
  /*
    for(int k= 0; k < maxZ; k++)
    { 
    cout<<endl<<"z= "<<k; 
    for(int j=0 ;j < maxY ; j++)
    {
    cout<<endl;
    for(int i=0 ;i < maxX ; i++)
    {
    cout<<bitmap2[i][j][k];
    }
    } 
    }
  */
  
  // ....................................................................
  // creates the list of the existing vertices
  // ....................................................................
  
  
  Point verts1[NoVoxels*8];
  for(int rr=0;rr<NoVoxels*8;rr++)
    { 
      Point np(3);//Sets verts to 3d point
      verts1[rr]=np;
      verts1[rr][0] = -1;
      verts1[rr][1] = -1;
      verts1[rr][2] = -1;
    }
  int kk=0;
  
  for(int k=0;k<maxZ;k++)
    {
      for(int j=0;j<maxY;j++)
	{
	  for(int i=0;i<maxX;i++)
	    {//cerr<<"k= "<<k<<"j= "<<j<<"i="<<i<<endl;
	      if( bitmap2[i][j][k] ==1)
		{	//cerr<<" kk = "<<kk<<endl;
		  verts1[kk][0] = i;
		  verts1[kk][1] = j;
		  verts1[kk][2] = k;
		  kk++;
		  
		}
	    }
	}
    }
  Point verts[kk];
  
  //cerr<<" kkfinal = "<<kk<<endl;
  
  for(int i = 0;i<kk;i++)// 
    {
      
      verts[i]=verts1[i];
      // cerr<<" -- "<<verts [i][0]<<" "<<verts [i][1]<<" "<<verts [i][2]<<endl;
    }
  
  
  
  // ....................................................................
  // reconstruction of the voxel shapes
  // ....................................................................
  
  
  int NoVRT=kk;
  
  int ver1[NoVoxels];
  int ver2[NoVoxels];
  int ver3[NoVoxels];
  int ver4[NoVoxels];
  
  for(int jj=0;jj<NoVoxels;jj++)
    {
      ver1[jj]=0;
      ver2[jj]=0;
      ver3[jj]=0; 
      ver4[jj]=0;
    }
  //cerr<<"Novrt=== "<<NoVRT<<endl;
  
  // ***  For every voxel we retrieve the one side parallel to yOz
  
  for(int k = 0;k <NoVoxels;k++)
    { //cerr<<voxels[k][0]<<"  "<<voxels[k][1]<<"  "<<voxels[k][2]<<endl;
      for( int j=0;j < NoVRT ;j++)
	{
	  
	  
	  if( verts[j][0] == voxels[k][0] && verts[j][1] == voxels[k][1] &&
	      verts[j][2] ==  voxels[k][2])
	    {  ver1[k]= j;// the vertice has same coordinates with the voxel
	    // cerr<<" ver1 = "<<ver1[k]<<endl;
	    }
	  
	  if( verts[j][0] == voxels[k][0] && verts[j][1] == voxels[k][1] &&
	      verts[j][2] ==  ((voxels[k][2])+1)) 
	    {  ver2[k]= j; // Z vertice = Z voxel+1
	    //  cerr<<" ver2 = "<<ver2[k]<<endl;
	    }
	  
	  if( verts[j][0] == voxels[k][0] && verts[j][1] == voxels[k][1]+1 &&
	      verts[j][2] ==  voxels[k][2]+1) 
	    {  ver3[k]= j;// y vertice = y voxel+1..  Z vertice = Z voxel+1
	    // cerr<<" ver3 = "<<ver3[k]<<endl;
	    }
	  
	  if( verts[j][0] == voxels[k][0] && verts[j][1] == (voxels[k][1])+1 &&
	      verts[j][2] ==  voxels[k][2]) 
	    { ver4[k]= j;// y vertice = y voxel+1
	    // cerr<<" ver4 = "<<ver4[k]<<endl;
	    }
	  
	}
    }
  
  // just printing 
  /*
    for(int g=0;g<NoVoxels;g++)
    {
    cerr<<ver1[g]<<"   "<<ver2[g]<<"   "<<ver3[g]<<"   "<<ver4[g]<<endl;
    } */
  
  // ....................................................................
  //(*Stores all the faces of the shape Inner and outer*)
  // ....................................................................
  
  
  int InitialFaces[6*NoVoxels][4];
  
  
  for ( int f=0;f<NoVoxels;f++)
    {
      InitialFaces [f][0]=ver1[f];//face1
      InitialFaces [f][1]=ver2[f];
      InitialFaces [f][2]=ver3[f];
      InitialFaces [f][3]=ver4[f];
      
      InitialFaces [f+NoVoxels][0]=ver1[f];//face2
      InitialFaces [f+NoVoxels][1]=ver1[f]+1;
      InitialFaces [f+NoVoxels][2]=ver2[f]+1;
      InitialFaces [f+NoVoxels][3]=ver2[f];
      
      InitialFaces [f+ 2*NoVoxels][0]=ver1[f]+1;//face3
      InitialFaces [f+ 2*NoVoxels][1]=ver4[f]+1;
      InitialFaces [f+ 2*NoVoxels][2]=ver3[f]+1;
      InitialFaces [f+ 2*NoVoxels][3]=ver2[f]+1;
      
      InitialFaces[f + 3*NoVoxels][0]=ver4[f];//face4
      InitialFaces[f + 3*NoVoxels][1]=ver3[f];
      InitialFaces[f + 3*NoVoxels][2]=ver3[f]+1;
      InitialFaces[f + 3*NoVoxels][3]=ver4[f]+1;
      
      InitialFaces [f+ 4*NoVoxels][0]=ver1[f];//face5
      InitialFaces [f+ 4*NoVoxels][1]=ver4[f];
      InitialFaces [f+ 4*NoVoxels][2]=ver4[f]+1;
      InitialFaces [f+ 4*NoVoxels][3]=ver1[f]+1;
      
      InitialFaces[f+ 5*NoVoxels][0]=ver2[f];//face6
      InitialFaces[f+ 5*NoVoxels][1]=ver2[f]+1;
      InitialFaces[f+ 5*NoVoxels][2]=ver3[f]+1;
      InitialFaces[f+ 5*NoVoxels][3]=ver3[f];
      
    }
  
  
  //Printing for debugging
  /*
    for ( int f=0;f<6*NoVoxels;f++)
    {
    cout<<InitialFaces[f][0]<<"  "<<InitialFaces[f][1]<<"  "<<InitialFaces[f][2]<<"  "<<InitialFaces[f][3]<<endl;
    }
  */
  
  
  // ....................................................................
  //Copy initialFaces to Surfaces
  // ....................................................................
  
  int Surfaces[6*NoVoxels][4];
  
  
  
  for(int rr=0;rr<6*NoVoxels;rr++)
    
    
    {
      
      
      for(int ff=0;ff<4;ff++)
	
	{  
	  
	  
	  Surfaces[rr][ff]=InitialFaces[rr][ff];
	  // cerr<< rr<<" "<<6*NoVoxels<<" "<<ff<<" "<<4 <<endl;
	}
    }
  //cerr<<"endloop"<<endl;
  /*
    ofstream the("desauto.dat");
    
    int NoVo= 6*NoVoxels;
    
    the<<NoVo<<endl;
    
       for(int i = 0;i< NoVo ;i++)
       { 
       the<<"["<< InitialFaces[i][0]<<" "<< InitialFaces[i][1]<<" "<< InitialFaces[i][2]<<" "<< InitialFaces[i][3]<<" ] " <<endl;
       }
  */
  
  // ....................................................................
  //  Puts zeros in the place of the interior faces
  // ....................................................................
  
  
  for(int re=0;re<6*NoVoxels;re++)//puts zeros in the place of the interior faces
    {
      for(int fe=0;fe<6*NoVoxels;fe++)
	{
	  if(InitialFaces[re][0]==InitialFaces[fe][0])
	    {
	      if(InitialFaces[re][1]==InitialFaces[fe][3])
		{
		  if(InitialFaces[re][3]==InitialFaces[fe][1])
		    {
		      for(int in=0;in<4;in++)
			{
			  Surfaces[fe][in]=-1;
			  Surfaces[re][in]=-1;
			}
		    }
		}
	    }
	}
    }
  

  
  
  
  // ....................................................................
  //     Finds the length of the non zero faces
  // ....................................................................
  
  
  //int LengoFaces=0;

  for(int gy=0;gy<NoVoxels*6;gy++)
    {
      if(Surfaces[gy][0]!=-1 || Surfaces[gy][1]!=-1 || Surfaces[gy][2]!=-1 || Surfaces[gy][3]!=-1)
	{
	  LengoFaces++;
	}
    }
  
  
  //int AntiClocFaces[LengoFaces][4];
  AntiClocFaces = new IVector[LengoFaces];
  
  for(int i=0;i<LengoFaces;i++)AntiClocFaces[i].New(4);
  
  
  LengoFaces=0;
  
  
  
  // ....................................................................
  //     Makes a new array Surfaces with the surface faces
  // ....................................................................
  
  
  
  for(int ga=0;ga<NoVoxels*6;ga++)
    {
      if(Surfaces[ga][0]!=-1 ||  Surfaces[ga][1]!=-1 || Surfaces[ga][2]!=-1 || Surfaces[ga][3]!=-1)
	{
	  AntiClocFaces[LengoFaces][0]=Surfaces[ga][0];
	  AntiClocFaces[LengoFaces][1]=Surfaces[ga][1];
	  AntiClocFaces[LengoFaces][2]=Surfaces[ga][2];
	  AntiClocFaces[LengoFaces][3]=Surfaces[ga][3];
	  LengoFaces++;
	}
    }
  
  
  
  
  
  
  // ....................................................................
  //      Find the inside vertices
  // ....................................................................
  
  
  int TobeRemov[NoVRT];
  int map0[NoVRT];
  int li=0;
  int la=0;
  
  for(int ji=0;ji<NoVRT;ji++)
    {
      if(inpointlist(voxels,NoVoxels,((int)verts1[ji][0]),((int)verts1[ji][1]),((int)verts1[ji][2]))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]-1),((int)verts1[ji][1]),((int)verts1[ji][2]))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]-1),((int)verts1[ji][1]-1),((int)verts1[ji][2]))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]-1),((int)verts1[ji][1]-1),((int)verts1[ji][2]-1))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]),((int)verts1[ji][1]),((int)verts1[ji][2]-1))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]-1),((int)verts1[ji][1]),((int)verts1[ji][2]-1))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]),((int)verts1[ji][1]-1),((int)verts1[ji][2]))&&
	 inpointlist(voxels,NoVoxels,((int)verts1[ji][0]),((int)verts1[ji][1]-1),((int)verts1[ji][2]-1)))
	{ 
	  TobeRemov[li]=ji;
	  //	 cerr<<"TobeRemov "<<li<<" - "  <<TobeRemov[li]<<endl;
	  li++;
	}
      else
	{
	  map0[la]=ji;
	  
	  //	 cerr<<"map  "<<map0[la]<<endl;
	  la++;
	}
      
    }
  
  
  
  // ....................................................................
  //Remove inside vertices from the list
  // ....................................................................
  
  int hh=0;
  int jf=0;
  
  if (li>0) 
    {
      NeWNoVRT=NoVRT-li;
      VertexList= new Point[NeWNoVRT];
      
      // ....................................................................
      //       Updates the vertex list
      // ....................................................................
      
      for (int jj=0;jj<NeWNoVRT;jj++)
	{
	  VertexList[jj]=verts[map0[jj]];
	  //	 cout<<VertexList[jj]<<endl;
	} 
      
      //  hh=0;
      // jf=0;
      

      
      // ....................................................................
      //   Updates the names of the vertices in the faces lists
      // ....................................................................
      
      
      for (int gr=0;gr<LengoFaces;gr++)
	{
	  for (int gk=0;gk<4;gk++)
	    {
	      for(int fg=0;fg<NeWNoVRT;fg++)
		{
		 if(AntiClocFaces[gr][gk]==map0[fg])
		   {
		     AntiClocFaces[gr][gk]=fg;
		   }
		}
	    }
	}
    }
  else
    {
      NeWNoVRT=NoVRT;
      VertexList= new Point[NeWNoVRT];
      // ....................................................................
      for (int jj=0;jj<NeWNoVRT;jj++)
	{
	 VertexList[jj]=verts[jj];
	} 
      
    }
  
  
  // ....................................................................
  //     Printing The final surface lego faces
  // ....................................................................   
     
  // for(int nnn=0;nnn<LengoFaces;nnn++)
  //   { 
  //     cout<< AntiClocFaces[nnn][0]<<" - "<<AntiClocFaces[nnn][1]<<" - "<<AntiClocFaces[nnn][2]<<" - "<<AntiClocFaces[nnn][3]<<endl;
  //   }  /**/
  
  
  
  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>> EDW <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

}//End of function 


// boolean for the existance or not of a point in the list.


bool inpointlist(Point *list,int length,int i,int j, int k)
{
  bool inlist = false;
  Point p(3);
  p[0]=i;
   p[1]=j;
   p[2]=k;
   for(int i = 0;(i<length)&&!inlist;i++)
     {
       if((list[i][0]==((int)p[0])) &&(list[i][1]==((int)p[1]))&&(list[i][2]==((int)p[2]))) inlist=true;
     }
   return inlist;
} 

//before calling do int *newlist; 
//    pass newlist as fourth variable
// when finished with new list don't forget delete [] newlist;
/*
  int  drop(int *list,int length,int chop, int *newlist)//drops part of a list ,
  {
  int newlength=length-chop;
  //int newlist[newlength];
  newlist = new int[newlength];
  if (chop>0)//if chop >0 drops the first chop pieces of the list
  {
  for (int jj=0;jj<newlength;jj++)
  {
  newlist[jj]=list[jj+chop];
  }
  }
  if (chop<=0)//if chop >0 drops the last chop pieces of the list
  {
  for (int jf=0;jf<newlength;jf++)
  {
  newlist[jf]=list[jf];
  }
    }
    return newlength;//newlist;
    }*/

