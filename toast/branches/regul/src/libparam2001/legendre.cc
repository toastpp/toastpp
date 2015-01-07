#define HARMONICLIB_IMPLEMENTATION
#include "mathlib.h"
#include "felib.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>

//#include "harmoniclib.h"
//#include "surfacenet.h" 
//#include "diffusphere.h"
//#include "optimizzz.h"
//#include "usefulfan.h"
//#include "complex.h"


using namespace std;
//const double M_PI=3.141592653589793;
int factorial(int i)
{
  if (i==2) return 2;
  if (i==1) return 1;
  if (i==0) return 1;
  int result=i*factorial(i-1);
  return result;
}


double plgndr(int l, int m, double x) //Computes the associated Legendre polynomial 
                                    //Pm,l(x). Here m and l are integers satisfying
                                    // 0 < m < l, while x lies in the range -1 < x < 1.
{
  double fact,pll,pmm,pmmp1,somx2; 
  int i,ll;

  if (m < 0 || m > l || fabs(x) > 1.0)
    {
      cout<<" Bad arguments in routine plgndr ";
      return 0;
    }
 
  pmm= 1.0;					//Compute Pm,m

  if ( m > 0 )
    {
      somx2 = sqrt ( ( 1.0 - x ) * ( 1.0 + x ) );
      fact = 1.0;
      for ( i = 1 ; i <= m ; i++)
	{
	  pmm *= - fact * somx2 ;
	  fact += 2.0 ;
	}
    }
  
  if ( l == m )
    return pmm;
  
  else                                            //Compute Pm,m+1
    {   
      pmmp1 = x * ( 2 * m + l ) * pmm;
      
      if (l == ( m + 1 ) ) return pmmp1;
      
      else 
	{		                          //Compute Pm,l  l>m+1
	  for( ll= (m+2); ll<= l; ll++) 
	    {
	      pll = ( x * ( 2 * ll - 1 ) * pmmp1 - ( ll + m - 1 ) * pmm ) / ( ll - m );  
	      pmm = pmmp1;
	      pmmp1 = pll;
	    }
	  return pll;
	}
    }
}



std::complex<double> SphericalHarmonic(int l,int m, double thi, double fi)// calculates the spherical harmonic
                                                         //of degree l order m at the (thi,fi) point
{ 

  double d1, d2, p, r,i,re,im;
  double ll= ((double) l );
  double mm = ((double) m );
 
 if(m >= 0)
    {
    
  
    

      d1 = sqrt ( ( 2 * ll + 1 ) / ( 4 * Pi ) );
     
      d2 = sqrt ( ((double)factorial ( l - m )) / ((double)factorial ( l + m )) );
     
      p = plgndr (l, m, cos( thi ) );
      
      i = sin ( mm * fi);
     
      r = cos ( mm * fi);
     
      re = d1 * d2 * p * r;
      im = d1 * d2 * p * i;
      //cout <<re<<"  "<<im<<" i";
    }

  if(m < 0)
    {
      double mon;
     
      m = (-1) * m;
      mm = (-1.0) * mm;

      d1 = sqrt ( ( 2 * ll + 1 ) / ( 4 * Pi ) );
     
      d2 = sqrt ( ((double)factorial ( l - m )) / ((double)factorial ( l + m )) );
      
      p = plgndr (l, m, cos( thi ) );
      
      i = sin ( mm * fi);
     
      r = cos ( mm * fi);
      mon = pow((-1.0),2.0);
      re = d1 * mon * d2 * p * r;
      im = (-1)*d1 * mon * d2 * p * i;
      //  cout <<re<<"  "<<im<<" i";
    }  
  std::complex<double> z(re,im);
  //  cout << " complex z= " << z << endl;
  return z;   
}





