
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <iostream.h>


double factorial(double i)
{
  if (i==2) return 2.0;
  if (i==1) return 1.0;
  if (i==0) return 1.0;
  double result=i*factorial(i-1);
  double result2=result*1.0;
  return result2;
}


double plgndr(int l, int m, double x) //Computes the associated Legendre polynomial 
                                    //Pm,l(x). Here m and l are integers satisfying
                                    // 0 < m < l, while x lies in the range -1 < x < 1.
{
  long double fact,pll,pmm,pmmp1,somx2; 
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
    {
      return pmm;
    }
  else                                            //Compute Pm,m+1
    {   
      pmmp1 = x * ( 2 * m + 1 ) * pmm;
      //  cout<<"===="<<pmmp1<<endl;
      
      if (l == ( m + 1 ) ) 
	{
	  return pmmp1;
	}
      else 
	{		                          //Compute Pm,l  l>m+1
	  for( ll= (m+2); ll<= l; ll++) 
	    {

	      pll = ( x * ( 2 * ll - 1 ) * pmmp1 - ( ll + m - 1 ) * pmm ) / ( ll - m );  
	      //  cout<<pll<<"----"<<pmmp1<<endl;
	      pmm = pmmp1;
	      pmmp1 = pll;
	    }
	  return pll;
	}
    }
}



double SphericalHarmonic(int l,int m, double thi, double fi)// calculates the spherical harmonic
                                                         //of degree l order m at the (thi,fi) poi99nt
{
  cout.precision(100);
  long double d1, d2, p, r,i,re,im, ll , mm;
  ll=l*1.0;
  mm=m*1.0;
  d1 = sqrt ( ( 2 * ll + 1 ) / ( 4 * M_PI ) );

  d2=1.0;
  for (int rr=(l-m+1);rr<(l+m+1);rr++)
    {
     long  double rrd=(1.0/rr);
     d2=(rrd)*d2;
      //cout<<"--> "<<rr<<"d2="<<d2<<endl;
    }
  d2 = sqrt(d2);
 
  // d2 = sqrt ( factorial ( ll - mm ) / factorial ( ll + mm ) );
 
  p = plgndr (l, m, cos ( thi ) );
   cout<<"p="<<p<<endl;
  i = sin ( mm * fi);
  // cout<<"r="<<r<<endl;
  r = cos ( mm * fi);
  // cout<<"i="<<i<<endl;
  re = d1 * d2 * p * r;
  im = d1 * d2 * p * i;
  cout <<re<<" + i "<<im<<endl;
  return 0;

}







int main(void)
{
  /*  for(int i=0;i<10;i++)
    {
      for(int j=0;j<=i;j++)
	{
	  cout<<"l= "<<i<<" m= "<<j<<" "<<SphericalHarmonic(i,j,9,6)<<endl;
	}
	}*/
   // cout<<plgndr(2,0,0.5)<<endl;
  SphericalHarmonic(9,6,1,1);
  return 0;
}
