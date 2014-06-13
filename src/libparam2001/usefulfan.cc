#define HARMONICLIB_IMPLEMENTATION
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
 

using namespace std;
//******************************************************************************
//*                          useful functions                                  *
//******************************************************************************


//Super atan function for angles from 0 to 2 Pi.

double atan3(double hmit, double syn)
{
double gwnia;
int mhden=0;
if(syn==0 && hmit==1)
{ 
  gwnia =Pi/2;
  mhden=1;
}
if(syn==0 && hmit==-1)
{
  gwnia =3 * Pi / 2 ;
  mhden=1;
}
if (mhden==0)
{
  double efa =atan2(hmit, syn);
  if(efa <0)
    {
      gwnia= efa + 2 * Pi;
    }
  else
    {
      gwnia = efa;
    }
}
return gwnia;
}


//Cross Product
RVector cross_product(RVector edge1,RVector edge2)
{
  RVector Snormal(3);
   Snormal[0] = edge1[1] * edge2[2] - edge1[2] * edge2[1];
   Snormal[1] = edge1[2] * edge2[0] - edge1[0] * edge2[2];
   Snormal[2] = edge1[0] * edge2[1] - edge1[1] * edge2[0];

   return Snormal;
}


//Side of Quadrilateral Sij----RVectors
 double  SijV(RVector  x1,RVector x2)
{
  double dot=x1[0]*x2[0]+x1[1]*x2[1]+x1[2]*x2[2];
  double s=acos(dot);
  return s;
}



//Side of Quadrilateral Sij--- coordinates
 double  SijS(double x1, double y1,double z1,double x2,double y2,double z2)
{
  double dot=x1*x2+y1*y2+z1*z2;
  double s=acos(dot);
  return s;
}


//Magnitude of vector x
double NormV(RVector x)
{
  double n=sqrt( pow(x[0],2.0)+pow(x[1],2.0)+pow(x[2],2.0) );
  return n;
}


//Magnitude of vector {a,b,c}
double NormS(double a, double b, double c)
{
  double n=sqrt( pow(a,2.0)+pow(b,2.0)+pow(c,2.0) );
  return n;
}


//Sine of the side using Cross product
double  Sin_Side(RVector edge1,RVector edge2)
{
  RVector cropro(3);
  cropro = cross_product( edge1, edge2);
  double N1=NormV(edge1);
  double N2=NormV(edge2);
  double NC=NormV(cropro);
  double sin_res= NC/(N1*N2);
  return sin_res;
}


//The CosLDAB

double CosLDAB(RVector VecD, RVector VecA, RVector VecB)
{
  double res= (VecB & VecD)- (VecA & VecD)*(VecA & VecB);
  return res; 
}

// The SinLDAB 
double SinLDAB(RVector VecD, RVector VecA, RVector VecB)
{
  double res= VecA & cross_product(VecB, VecD) ;
  return res; 
}

//The angle LDAB 

 double  LDAB (RVector VecD, RVector VecA, RVector VecB)
{
  double resa= atan3(SinLDAB( VecD, VecA, VecB), CosLDAB( VecD, VecA, VecB ));
    return resa;
}

//Find the area of the quadrilateral ABCD

double Area( RVector VecA, RVector VecB, RVector  VecC, RVector VecD)
{
  double area ;
  area= ( LDAB (VecB, VecC, VecD)+ LDAB (VecC, VecD, VecA) +LDAB (VecD, VecA, VecB)+ LDAB (VecA, VecB, VecC) )- 2* Pi;
  return area;
}


double DAreaDX( RVector VecA, RVector VecB, RVector  VecC, RVector VecD)
{

  double res=(-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) + 
           VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*(VecC[0] - VecB[0]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))\
       + (-(VecB[2]*VecC[1]) + VecB[1]*VecC[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] - 
         (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))/
    (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) + 
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) + 
      pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] - 
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0)) + 
   ((-(VecB[2]*VecD[1]) + VecB[1]*VecD[2])*(VecB[0]*VecD[0] + VecB[1]*VecD[1] + VecB[2]*VecD[2] - 
         (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecA[0]*VecD[0] + VecA[1]*VecD[1] + VecA[2]*VecD[2])) - 
      (-((VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*VecD[0]) - 
         VecB[0]*(VecA[0]*VecD[0] + VecA[1]*VecD[1] + VecA[2]*VecD[2]))*
       (VecA[2]*(-(VecB[1]*VecD[0]) + VecB[0]*VecD[1]) + VecA[1]*(VecB[2]*VecD[0] - VecB[0]*VecD[2]) + 
         VecA[0]*(-(VecB[2]*VecD[1]) + VecB[1]*VecD[2])))/
    (pow(VecB[0]*VecD[0] + VecB[1]*VecD[1] + VecB[2]*VecD[2] - 
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecA[0]*VecD[0] + VecA[1]*VecD[1] + VecA[2]*VecD[2]),2.0) + 
      pow(VecA[2]*(-(VecB[1]*VecD[0]) + VecB[0]*VecD[1]) + VecA[1]*(VecB[2]*VecD[0] - VecB[0]*VecD[2]) + 
        VecA[0]*(-(VecB[2]*VecD[1]) + VecB[1]*VecD[2]),2.0)) + 
   (-(((-(VecA[2]*VecC[1]) + VecA[1]*VecC[2])*VecD[0] + (VecA[2]*VecC[0] - VecA[0]*VecC[2])*VecD[1] + 
           (-(VecA[1]*VecC[0]) + VecA[0]*VecC[1])*VecD[2])*
         (VecC[0] - VecD[0]*(VecC[0]*VecD[0] + VecC[1]*VecD[1] + VecC[2]*VecD[2]))) + 
      (-(VecC[2]*VecD[1]) + VecC[1]*VecD[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] - 
         (VecA[0]*VecD[0] + VecA[1]*VecD[1] + VecA[2]*VecD[2])*(VecC[0]*VecD[0] + VecC[1]*VecD[1] + VecC[2]*VecD[2])))/
    (pow((-(VecA[2]*VecC[1]) + VecA[1]*VecC[2])*VecD[0] + (VecA[2]*VecC[0] - VecA[0]*VecC[2])*VecD[1] + 
        (-(VecA[1]*VecC[0]) + VecA[0]*VecC[1])*VecD[2],2.0) + 
      pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] - 
        (VecA[0]*VecD[0] + VecA[1]*VecD[1] + VecA[2]*VecD[2])*(VecC[0]*VecD[0] + VecC[1]*VecD[1] + VecC[2]*VecD[2]),2.0));
 return res;
}

double DAreaDY( RVector a, RVector b, RVector  c, RVector d)
{

  double res=(-((b[2]*(a[1]*c[0] - a[0]*c[1]) + b[1]*(-(a[2]*c[0]) + a[0]*c[2]) + b[0]*(a[2]*c[1] - a[1]*c[2]))*
         (c[1] - b[1]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]))) + 
      (b[2]*c[0] - b[0]*c[2])*(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - 
         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])))/
    (pow(b[2]*(a[1]*c[0] - a[0]*c[1]) + b[1]*(-(a[2]*c[0]) + a[0]*c[2]) + b[0]*(a[2]*c[1] - a[1]*c[2]),2.0) + 
      pow(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]),2.0))\
    + ((b[2]*d[0] - b[0]*d[2])*(b[0]*d[0] + b[1]*d[1] + b[2]*d[2] - 
         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(a[0]*d[0] + a[1]*d[1] + a[2]*d[2])) - 
      (-((a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*d[1]) - b[1]*(a[0]*d[0] + a[1]*d[1] + a[2]*d[2]))*
       (a[2]*(-(b[1]*d[0]) + b[0]*d[1]) + a[1]*(b[2]*d[0] - b[0]*d[2]) + a[0]*(-(b[2]*d[1]) + b[1]*d[2])))/
    (pow(b[0]*d[0] + b[1]*d[1] + b[2]*d[2] - (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(a[0]*d[0] + a[1]*d[1] + a[2]*d[2]),2.0) + 
      pow(a[2]*(-(b[1]*d[0]) + b[0]*d[1]) + a[1]*(b[2]*d[0] - b[0]*d[2]) + a[0]*(-(b[2]*d[1]) + b[1]*d[2]),2.0)) + 
   (-(((-(a[2]*c[1]) + a[1]*c[2])*d[0] + (a[2]*c[0] - a[0]*c[2])*d[1] + (-(a[1]*c[0]) + a[0]*c[1])*d[2])*
         (c[1] - d[1]*(c[0]*d[0] + c[1]*d[1] + c[2]*d[2]))) + 
      (c[2]*d[0] - c[0]*d[2])*(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - 
         (a[0]*d[0] + a[1]*d[1] + a[2]*d[2])*(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])))/
    (pow((-(a[2]*c[1]) + a[1]*c[2])*d[0] + (a[2]*c[0] - a[0]*c[2])*d[1] + (-(a[1]*c[0]) + a[0]*c[1])*d[2],2.0) + 
      pow(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - (a[0]*d[0] + a[1]*d[1] + a[2]*d[2])*(c[0]*d[0] + c[1]*d[1] + c[2]*d[2]),2.0));

 return res;
}


double DAreaDZ( RVector a, RVector b, RVector  c, RVector d)
{

  double res=(-((b[2]*(a[1]*c[0] - a[0]*c[1]) + b[1]*(-(a[2]*c[0]) + a[0]*c[2]) + b[0]*(a[2]*c[1] - a[1]*c[2]))*
         (c[2] - b[2]*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]))) + 
      (-(b[1]*c[0]) + b[0]*c[1])*(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - 
         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2])))/
    (pow(b[2]*(a[1]*c[0] - a[0]*c[1]) + b[1]*(-(a[2]*c[0]) + a[0]*c[2]) + b[0]*(a[2]*c[1] - a[1]*c[2]),2.0) + 
      pow(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(b[0]*c[0] + b[1]*c[1] + b[2]*c[2]),2.0))\
    + ((-(b[1]*d[0]) + b[0]*d[1])*(b[0]*d[0] + b[1]*d[1] + b[2]*d[2] - 
         (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(a[0]*d[0] + a[1]*d[1] + a[2]*d[2])) - 
      (-((a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*d[2]) - b[2]*(a[0]*d[0] + a[1]*d[1] + a[2]*d[2]))*
       (a[2]*(-(b[1]*d[0]) + b[0]*d[1]) + a[1]*(b[2]*d[0] - b[0]*d[2]) + a[0]*(-(b[2]*d[1]) + b[1]*d[2])))/
    (pow(b[0]*d[0] + b[1]*d[1] + b[2]*d[2] - (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])*(a[0]*d[0] + a[1]*d[1] + a[2]*d[2]),2.0) + 
      pow(a[2]*(-(b[1]*d[0]) + b[0]*d[1]) + a[1]*(b[2]*d[0] - b[0]*d[2]) + a[0]*(-(b[2]*d[1]) + b[1]*d[2]),2.0)) + 
   (-(((-(a[2]*c[1]) + a[1]*c[2])*d[0] + (a[2]*c[0] - a[0]*c[2])*d[1] + (-(a[1]*c[0]) + a[0]*c[1])*d[2])*
         (c[2] - d[2]*(c[0]*d[0] + c[1]*d[1] + c[2]*d[2]))) + 
      (-(c[1]*d[0]) + c[0]*d[1])*(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - 
         (a[0]*d[0] + a[1]*d[1] + a[2]*d[2])*(c[0]*d[0] + c[1]*d[1] + c[2]*d[2])))/
    (pow((-(a[2]*c[1]) + a[1]*c[2])*d[0] + (a[2]*c[0] - a[0]*c[2])*d[1] + (-(a[1]*c[0]) + a[0]*c[1])*d[2],2.0) + 
     pow(a[0]*c[0] + a[1]*c[1] + a[2]*c[2] - (a[0]*d[0] + a[1]*d[1] + a[2]*d[2])*(c[0]*d[0] + c[1]*d[1] + c[2]*d[2]),2.0));
      
 return res;
}



double DAx(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3(
       VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) , VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))*
     (-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*(VecC[0] - VecB[0]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))
         + (-(VecB[2]*VecC[1]) + VecB[1]*VecC[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }


double DAy(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3( VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) , VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
      ))*
     (-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*(VecC[1] - VecB[1]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))
         + (VecB[2]*VecC[0] - VecB[0]*VecC[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

double DAz(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3(  VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) , VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))*
     (-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*(VecC[2] - VecB[2]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])))
         + (-(VecB[1]*VecC[0]) + VecB[0]*VecC[1])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

//done
double DBx(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3(  VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) , VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
     ))*
     (-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*
          (-((VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*VecC[0]) -
            VecA[0]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))) +
       (VecA[2]*VecC[1] - VecA[1]*VecC[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

//done
double DBy(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3(     VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) , VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
   ))*
     (-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*
          (-((VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*VecC[1]) -
            VecA[1]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))) +
       (-(VecA[2]*VecC[0]) + VecA[0]*VecC[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

//done
double DBz(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3(   VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2])  ,  VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
    ))*
     (-((VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))*
          (-((VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*VecC[2]) -
            VecA[2]*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))) +
       (VecA[1]*VecC[0] - VecA[0]*VecC[1])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

//done
double DCx(RVector VecA,RVector VecB,RVector VecC)
{



double res= (cos(atan3(  VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) ,  VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
      ))*
     (-((VecA[0] - VecB[0]*(VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2]))*
          (VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))) +
       (-(VecA[2]*VecB[1]) + VecA[1]*VecB[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

//done

double DCy(RVector VecA,RVector VecB,RVector VecC)
{



double res= (cos(atan3(   VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) ,  VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
     ))*
     (-((VecA[1] - VecB[1]*(VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2]))*
          (VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))) +
       (VecA[2]*VecB[0] - VecA[0]*VecB[2])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

//done

double DCz(RVector VecA,RVector VecB,RVector VecC)
{


double res= (cos(atan3(  VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
        VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]) , VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
        (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2])
     ))*
     (-((VecA[2] - VecB[2]*(VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2]))*
          (VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
            VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]))) +
       (-(VecA[1]*VecB[0]) + VecA[0]*VecB[1])*(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
          (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]))))/
   (pow(VecB[2]*(VecA[1]*VecC[0] - VecA[0]*VecC[1]) + VecB[1]*(-(VecA[2]*VecC[0]) + VecA[0]*VecC[2]) +
       VecB[0]*(VecA[2]*VecC[1] - VecA[1]*VecC[2]),2.0) +
     pow(VecA[0]*VecC[0] + VecA[1]*VecC[1] + VecA[2]*VecC[2] -
       (VecA[0]*VecB[0] + VecA[1]*VecB[1] + VecA[2]*VecB[2])*(VecB[0]*VecC[0] + VecB[1]*VecC[1] + VecB[2]*VecC[2]),2.0));
       return res;
       }

double dabs(double in)
{
  if (in <0)in=-in;
  return in;
}
