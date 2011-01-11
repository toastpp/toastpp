/***************************************************************************
 * test_genPN_integrals.cc         Simon Arridge           26.07.10        *
 *                                                                         *
 ***************************************************************************/

#ifdef __BORLANDC__
#include <strstrea.h>
#include <conio.h>
#include <process.h>

typedef unsigned pid_t;
#else
#include <sstream>
#include <unistd.h>

#endif
#include <time.h>
#include <fstream.h>
#include <iomanip.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mathlib.h>
#include <felib.h>
#include "toast.h"

#include "rte3D_math.h"
#include "sphints.h"
#define Alm(l, m) (sqrt(complex((((l)+(m))*((l)-(m)))/((double)(2*(l)+1)*(2*(l)-1)), 0)))
#define Blm(l, m) (sqrt(complex((((l)+(m))*((l)+(m)-1))/((double)(2*(l)+1)*(2*(l)-1)), 0)))



#define MIN(A,B) ( (A) < (B) ? (A) : (B))

QMMesh qmmesh;
NodeList &nlist=qmmesh.nlist;
ElementList &elist=qmmesh.elist;



// error handler for FE library routines *************************************

void LocalErrorhandler (char *msg)
{
    cerr << "\nread_toastbemmesh (PID " << getpid() << ")\n" << msg << endl << flush;
    cerr << "Aborted.\n";
    //    logfile << msg << endl << "Aborted." << endl;
    exit (1);
}

void sincosY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c)
{

	if(m==0)
	{
		a[0] = Blm(l, 0)/complex(2.0, 0); b[0] = complex(-1, 0)*Blm(l+1, 1)/complex(2.0, 0); 
		c[0] = complex(-1, 0)*Blm(l, 0)/complex(2.0, 0); d[0] = Blm(l+1, 1)/complex(2.0, 0);
		
		a[1] = 0; b[1] = 0; 
		c[1] = 0; d[1] = 0;
		
		a1c(0, 0) = l-1; a1c(0, 1) = 1; b1c(0, 0) = l+1; b1c(0, 1) = 1;
		c1c(0, 0) = l-1; c1c(0, 1) = -1; d1c(0, 0) = l+1; d1c(0, 1) = -1; 
		
		a1c(1, 0) = -1; a1c(1, 1) = -1; b1c(1, 0) = -1; b1c(1, 1) = -1;
		c1c(1, 0) = -1; c1c(1, 1) = -1; d1c(1, 0) = -1; d1c(1, 1) = -1; 				
	}
	else if(m>0)
	{
	    	a[0] = Blm(l, -m)/(complex(2.0*sqrt(2), 0)); b[0] =  complex(-1, 0)*Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0)); 
	    	c[0] = complex(-1, 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); d[0] = Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0));
	    	a[1] = complex(sign(m), 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); b[1] = complex(-1*sign(m), 0)*Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0));
		c[1] =  complex(-1*sign(m), 0)*Blm(l, -m)/(complex(2.0*sqrt(2), 0)); d[1] = complex(sign(m), 0)*Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0));

	   	a1c(0, 0) = l-1; a1c(0, 1) = m+1; b1c(0, 0) = l+1; b1c(0, 1) = m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = m-1; d1c(0, 0) = l+1; d1c(0, 1) = m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = -m+1; b1c(1, 0) = l+1; b1c(1, 1) = -m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = -m-1; d1c(1, 0) = l+1; d1c(1, 1) = -m-1; 
 
	}
	else
	{
	   	a[0] =  Blm(l, m)/(complex(0, 2.0*sqrt(2))); b[0] =  complex(-1, 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2))); 
	    	c[0] =   complex(-1, 0)*Blm(l, -m)/(complex(0, 2.0*sqrt(2))); d[0] = Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2)));
	    	a[1] = complex(sign(m+1), 0)*Blm(l, -m)/(complex(0, 2.0*sqrt(2))); b[1] = complex(-1*sign(m+1), 0)*Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2)));
		c[1] =  complex(-1*sign(m+1), 0)*Blm(l, m)/(complex(0, 2.0*sqrt(2))); d[1] =  complex(sign(m+1), 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2)));


	   	a1c(0, 0) = l-1; a1c(0, 1) = -m+1; b1c(0, 0) = l+1; b1c(0, 1) = -m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = -m-1; d1c(0, 0) = l+1; d1c(0, 1) = -m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = m+1; b1c(1, 0) = l+1; b1c(1, 1) = m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = m-1; d1c(1, 0) = l+1; d1c(1, 1) = m-1; 
	}

}


void sinsinY(const int l, const int m, CVector& a, CVector& b, CVector& c, CVector& d, IDenseMatrix& a1c, IDenseMatrix& b1c, IDenseMatrix& c1c, IDenseMatrix& d1c)
{

	if(m==0)
	{
		a[0] = Blm(l, 0)/complex(0, 2.0); b[0] = complex(-1, 0)*Blm(l+1, 1)/complex(0, 2.0); 
		c[0] = Blm(l, 0)/complex(0, 2.0); d[0] = complex(-1, 0)*Blm(l+1, 1)/complex(0, 2.0);
		
		a[1] = 0; b[1] = 0; 
		c[1] = 0; d[1] = 0;
		
		a1c(0, 0) = l-1; a1c(0, 1) = 1; b1c(0, 0) = l+1; b1c(0, 1) = 1;
		c1c(0, 0) = l-1; c1c(0, 1) = -1; d1c(0, 0) = l+1; d1c(0, 1) = -1; 
		
		a1c(1, 0) = -1; a1c(1, 1) = -1; b1c(1, 0) = -1; b1c(1, 1) = -1;
		c1c(1, 0) = -1; c1c(1, 1) = -1; d1c(1, 0) = -1; d1c(1, 1) = -1; 				
	}
	else if(m>0)
	{
	    	a[0] = Blm(l, -m)/(complex(0, 2.0*sqrt(2))); b[0] =  complex(-1, 0)*Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2))); 
	    	c[0] = Blm(l, m)/(complex(0, 2.0*sqrt(2))); d[0] = complex(-1, 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2)));
	    	a[1] = complex(sign(m), 0)*Blm(l, m)/(complex(0, 2.0*sqrt(2))); b[1] = complex(-1*sign(m), 0)*Blm(l+1, -m+1)/(complex(0, 2.0*sqrt(2)));
		c[1] =  complex(sign(m), 0)*Blm(l, -m)/(complex(0, 2.0*sqrt(2))); d[1] = complex(-1*sign(m), 0)*Blm(l+1, m+1)/(complex(0, 2.0*sqrt(2)));

	   	a1c(0, 0) = l-1; a1c(0, 1) = m+1; b1c(0, 0) = l+1; b1c(0, 1) = m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = m-1; d1c(0, 0) = l+1; d1c(0, 1) = m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = -m+1; b1c(1, 0) = l+1; b1c(1, 1) = -m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = -m-1; d1c(1, 0) = l+1; d1c(1, 1) = -m-1; 
 
	}
	else
	{
	   	a[0] =  complex(-1, 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); b[0] =  Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0)); 
	    	c[0] =   complex(-1, 0)*Blm(l, -m)/(complex(2.0*sqrt(2), 0)); d[0] = Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0));
	    	a[1] = complex(sign(m), 0)*Blm(l, -m)/(complex(2.0*sqrt(2), 0)); b[1] = complex(-1*sign(m), 0)*Blm(l+1, m+1)/(complex(2.0*sqrt(2), 0));
		c[1] =  complex(sign(m), 0)*Blm(l, m)/(complex(2.0*sqrt(2), 0)); d[1] =  complex(-1*sign(m), 0)*Blm(l+1, -m+1)/(complex(2.0*sqrt(2), 0));


	   	a1c(0, 0) = l-1; a1c(0, 1) = -m+1; b1c(0, 0) = l+1; b1c(0, 1) = -m+1;
	   	c1c(0, 0) = l-1; c1c(0, 1) = -m-1; d1c(0, 0) = l+1; d1c(0, 1) = -m-1;

	   	a1c(1, 0) = l-1; a1c(1, 1) = m+1; b1c(1, 0) = l+1; b1c(1, 1) = m+1;
	   	c1c(1, 0) = l-1; c1c(1, 1) = m-1; d1c(1, 0) = l+1; d1c(1, 1) = m-1; 
	}
}

void cosY(const int l, const int m, CVector& e, CVector& f, IDenseMatrix& e1c, IDenseMatrix& f1c)
{
	if(m==0)
	{
		e[0] = Alm(l, 0); f[0] = Alm(l+1, 0);
		e[1] = 0; f[1] = 0;
		
		e1c(0, 0) = l-1; e1c(0, 1) = 0;
		f1c(0, 0) = l+1; f1c(0, 1) = 0;
		
		e1c(1, 0) = -1; e1c(1, 1) = -1;
		f1c(1, 0) = -1; f1c(1, 1) = -1;

	}	
	else if(m>0)
	{
		e[0] = Alm(l, m)/complex(sqrt(2), 0); f[0] = Alm(l+1, m)/complex(sqrt(2), 0);
		e[1] = complex(sign(m), 0)*Alm(l, -m)/complex(sqrt(2), 0); f[1] = complex(sign(m), 0)*Alm(l+1, -m)/complex(sqrt(2), 0);
	
		e1c(0, 0) = l-1; e1c(0, 1) = m;
		f1c(0, 0) = l+1; f1c(0, 1) = m;

		e1c(1, 0) = l-1; e1c(1, 1) = -m;
		f1c(1, 0) = l+1; f1c(1, 1) = -m;
	}
	else
	{
		e[0] = Alm(l, -m)/complex(0, sqrt(2)); f[0] = Alm(l+1, -m)/complex(0, sqrt(2));
		e[1] = complex(-1*sign(m), 0)*Alm(l, m)/complex(0, sqrt(2)); f[1] = complex(-1*sign(m), 0)*Alm(l+1, m)/complex(0, sqrt(2));
	
		e1c(0, 0) = l-1; e1c(0, 1) = -m;
		f1c(0, 0) = l+1; f1c(0, 1) = -m;

		e1c(1, 0) = l-1; e1c(1, 1) = m;
		f1c(1, 0) = l+1; f1c(1, 1) = m;

	}

}

void sphY(const int l, const int m, CVector& p, IDenseMatrix& p1c)
{
	if(m==0)
	{
		p[0] = 1.0; p[1] = 0;
		
		p1c(0, 0) = l; p1c(0, 1) = m;
		p1c(1, 0) = -1; p1c(1, 1) = -3;
		
	}	
	else if(m>0)
	{
		p[0] = complex(1.0, 0)/complex(sqrt(2), 0); p[1] = complex(sign(m), 0)/complex(sqrt(2), 0);
	
		p1c(0, 0) = l; p1c(0, 1) = m;
		p1c(1, 0) = l; p1c(1, 1) = -m;

	}
	else
	{
		p[0] = complex(1.0, 0)/complex(0, sqrt(2)); p[1] = complex(-1*sign(m), 0)/complex(0, sqrt(2));
	
		p1c(0, 0) = l; p1c(0, 1) = -m;
		p1c(1, 0) = l; p1c(1, 1) = m;

	}
}

double Integrate2(CVector &a, CVector &b, IDenseMatrix &a1c, IDenseMatrix &b1c)
{
   complex temp;
   temp =  a[0]*b[0]*kronD(a1c.Row(0), b1c.Row(0)) + a[0]*b[1]*kronD(a1c.Row(0), b1c.Row(1)) + a[1]*b[0]*kronD(a1c.Row(1), b1c.Row(0)) + a[1]*b[1]*kronD(a1c.Row(1), b1c.Row(1));
   return(temp.re);
}

void genmat_angint_3D(const int sphOrder, const int angN, RCompRowMatrix& Aint, RCompRowMatrix& Aintsc, RCompRowMatrix& Aintss, RCompRowMatrix& Aintc, RCompRowMatrix& Aintscsc, RCompRowMatrix& Aintscss, RCompRowMatrix& Aintscc, RCompRowMatrix& Aintssss, RCompRowMatrix& Aintssc, RCompRowMatrix& Aintcc)
{
	RDenseMatrix dnsAint(angN, angN), dnsAintsc(angN, angN), dnsAintss(angN, angN), dnsAintc(angN, angN);
	RDenseMatrix dnsAintscsc(angN, angN), dnsAintscss(angN, angN), dnsAintcc(angN, angN);
	RDenseMatrix dnsAintscc(angN, angN), dnsAintssss(angN, angN), dnsAintssc(angN, angN);

	CVector a1(2), b1(2), c1(2), d1(2), e1(2), f1(2), p1(2), a2(2), b2(2), c2(2), d2(2), e2(2), f2(2), p2(2);
        IDenseMatrix a1c(2, 2), b1c(2, 2), c1c(2, 2), d1c(2, 2), e1c(2, 2), f1c(2, 2), p1c(2, 2);
	IDenseMatrix a2c(2, 2), b2c(2, 2), c2c(2, 2), d2c(2, 2), e2c(2, 2), f2c(2, 2), p2c(2, 2);
	 

	int is, js, indl1, indl2;

	for(int l1 = 0; l1 <= sphOrder; l1++){
		indl1 = getPos(l1, -1*l1);
		cout<<"angint: order: "<<l1<<endl;
		for(int m1 = -1*l1; m1 <= l1; m1++){
		        			
			is = indl1 + l1 + m1;    	
			for(int l2 = 0; l2 <= sphOrder; l2++){
				indl2 = getPos(l2, -1*l2);
				for(int m2 = -1*l2; m2 <= l2; m2++){
					
					js = indl2 + l2 + m2;
				
					sphY(l1, m1, p1, p1c);
					sphY(l2, m2, p2, p2c);
					dnsAint(is, js) = Integrate2(p1, p2, p1c, p2c);

					sincosY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintsc(is, js) = Integrate2(p1, a2, p1c, a2c) + Integrate2(p1, b2, p1c, b2c);
					dnsAintsc(is, js) += Integrate2(p1, c2, p1c, c2c) + Integrate2(p1, d2, p1c, d2c);

					sinsinY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintss(is, js) = Integrate2(p1, a2, p1c, a2c) + Integrate2(p1, b2, p1c, b2c);
					dnsAintss(is, js) += Integrate2(p1, c2, p1c, c2c) + Integrate2(p1, d2, p1c, d2c);

					cosY(l2, m2, e2, f2, e2c, f2c);
					dnsAintc(is, js) = Integrate2(p1, e2, p1c, e2c) + Integrate2(p1, f2, p1c, f2c);

					sincosY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
					sincosY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintscsc(is, js) = Integrate2(a1, a2, a1c, a2c) + Integrate2(a1, b2, a1c, b2c);
					dnsAintscsc(is, js) += Integrate2(a1, c2, a1c, c2c) + Integrate2(a1, d2, a1c, d2c);
					dnsAintscsc(is, js) += Integrate2(b1, a2, b1c, a2c) + Integrate2(b1, b2, b1c, b2c);
					dnsAintscsc(is, js) += Integrate2(b1, c2, b1c, c2c) + Integrate2(b1, d2, b1c, d2c);
					dnsAintscsc(is, js) += Integrate2(c1, a2, c1c, a2c) + Integrate2(c1, b2, c1c, b2c);
					dnsAintscsc(is, js) += Integrate2(c1, c2, c1c, c2c) + Integrate2(c1, d2, c1c, d2c);
					dnsAintscsc(is, js) += Integrate2(d1, a2, d1c, a2c) + Integrate2(d1, b2, d1c, b2c);
					dnsAintscsc(is, js) += Integrate2(d1, c2, d1c, c2c) + Integrate2(d1, d2, d1c, d2c);

					
					sinsinY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintscss(is, js) = Integrate2(a1, a2, a1c, a2c) + Integrate2(a1, b2, a1c, b2c);
					dnsAintscss(is, js) += Integrate2(a1, c2, a1c, c2c) + Integrate2(a1, d2, a1c, d2c);
					dnsAintscss(is, js) += Integrate2(b1, a2, b1c, a2c) + Integrate2(b1, b2, b1c, b2c);
					dnsAintscss(is, js) += Integrate2(b1, c2, b1c, c2c) + Integrate2(b1, d2, b1c, d2c);
					dnsAintscss(is, js) += Integrate2(c1, a2, c1c, a2c) + Integrate2(c1, b2, c1c, b2c);
					dnsAintscss(is, js) += Integrate2(c1, c2, c1c, c2c) + Integrate2(c1, d2, c1c, d2c);
					dnsAintscss(is, js) += Integrate2(d1, a2, d1c, a2c) + Integrate2(d1, b2, d1c, b2c);
					dnsAintscss(is, js) += Integrate2(d1, c2, d1c, c2c) + Integrate2(d1, d2, d1c, d2c);

					cosY(l2, m2, e2, f2, e2c, f2c);
					dnsAintscc(is, js) = Integrate2(a1, e2, a1c, e2c) + Integrate2(a1, f2, a1c, f2c);
					dnsAintscc(is, js) += Integrate2(b1, e2, b1c, e2c) + Integrate2(b1, f2, b1c, f2c);
					dnsAintscc(is, js) += Integrate2(c1, e2, c1c, e2c) + Integrate2(c1, f2, c1c, f2c);
					dnsAintscc(is, js) += Integrate2(d1, e2, d1c, e2c) + Integrate2(d1, f2, d1c, f2c);

					cosY(l1, m1, e1, f1, e1c, f1c);
					dnsAintcc(is, js) = Integrate2(e1, e2, e1c, e2c) + Integrate2(e1, f2, e1c, f2c);
					dnsAintcc(is, js) += Integrate2(f1, e2, f1c, e2c) + Integrate2(f1, f2, f1c, f2c);


					sinsinY(l1, m1, a1, b1, c1, d1, a1c, b1c, c1c, d1c);
					dnsAintssc(is, js) = Integrate2(a1, e2, a1c, e2c) + Integrate2(a1, f2, a1c, f2c);
					dnsAintssc(is, js) += Integrate2(b1, e2, b1c, e2c) + Integrate2(b1, f2, b1c, f2c);
					dnsAintssc(is, js) += Integrate2(c1, e2, c1c, e2c) + Integrate2(c1, f2, c1c, f2c);
					dnsAintssc(is, js) += Integrate2(d1, e2, d1c, e2c) + Integrate2(d1, f2, d1c, f2c);

					sinsinY(l2, m2, a2, b2, c2, d2, a2c, b2c, c2c, d2c);
					dnsAintssss(is, js) = Integrate2(a1, a2, a1c, a2c) + Integrate2(a1, b2, a1c, b2c);
					dnsAintssss(is, js) += Integrate2(a1, c2, a1c, c2c) + Integrate2(a1, d2, a1c, d2c);
					dnsAintssss(is, js) += Integrate2(b1, a2, b1c, a2c) + Integrate2(b1, b2, b1c, b2c);
					dnsAintssss(is, js) += Integrate2(b1, c2, b1c, c2c) + Integrate2(b1, d2, b1c, d2c);
					dnsAintssss(is, js) += Integrate2(c1, a2, c1c, a2c) + Integrate2(c1, b2, c1c, b2c);
					dnsAintssss(is, js) += Integrate2(c1, c2, c1c, c2c) + Integrate2(c1, d2, c1c, d2c);
					dnsAintssss(is, js) += Integrate2(d1, a2, d1c, a2c) + Integrate2(d1, b2, d1c, b2c);
					dnsAintssss(is, js) += Integrate2(d1, c2, d1c, c2c) + Integrate2(d1, d2, d1c, d2c);
					}
			}
		}
	}
	Aint = shrink(dnsAint);Aintsc = shrink(dnsAintsc); Aintss = shrink(dnsAintss);Aintc = shrink(dnsAintc);
	Aintscsc = shrink(dnsAintscsc);Aintscss = shrink(dnsAintscss);Aintscc = shrink(dnsAintscc); Aintssc = shrink(dnsAintssc);
	Aintssss = shrink(dnsAintssss);Aintcc = shrink(dnsAintcc);
}


// main routine **************************************************************

int main (int argc, char *argv[])
{
    char cbuf[200];
    int i,j,el,is,js;
    const int sporder    = 2;
    const int nso = sporder+1;
    int nsp = nso*nso;

    CCompRowMatrix RY;

    RCompRowMatrix Aint1, Aintsc1, Aintss1, Aintc1, Anvec1;
    RCompRowMatrix Aintscsc1,  Aintscss1, Aintscc1,  Aintssss1,  Aintssc1,  Aintcc1;
    RCompRowMatrix Anvec_sc1, Anvec_ss1,  Anvec_c1;

    
    /* explicit RRz matrix */
    RDenseMatrix expRRx(nsp,nsp);
    expRRx(0,1)= expRRx(1,0) = 1/sqrt(3);
    expRRx(1,4) = 1/sqrt(5); expRRx(1,6) = -1/sqrt(15);
    expRRx(2,5) = 1/sqrt(5);
    expRRx(3,8) = -1/sqrt(5);
    expRRx(4,1) = 1/sqrt(5);
    expRRx(5,2) = 1/sqrt(5);
    expRRx(6,1) = -1/sqrt(15);
    expRRx(8,4) = -1/sqrt(5);
   

    //cout << "------------------Calling genmat_angint_3D_PN------------\n";
    genmat_angint_3D_PN(Aint1, Aintsc1, Aintss1, Aintc1, Anvec1, RY, sporder);
    //cout << "------------------returned from genmat_angint_3D_PN OK---\n";
    //cerr << "Aint and Anvec are the same ?\n" << (Aint1 - Anvec1) << endl;
    //    cout << "zintegral \n"<< Aintc << endl;

    /*RDenseMatrix RRz = Sparse2Dense(Aintc1);

    cout << "RRz\n" << RRz << endl;
    //   cout << "xintegal\n" << Aintsc << endl;
    RDenseMatrix RRx = Sparse2Dense(Aintsc1);   

    cout << "RRx\n" << RRx << endl;
    cout << "Hand built RRx\n" << expRRx << endl;

    RDenseMatrix RRy = Sparse2Dense(Aintss);   

    cout << "RRy\n" << RRy << endl;*/
 
 
    //cout << "--------------Calling genmat_angint_sdm_3D_PN-------------\n";
    genmat_angint_sdm_3D_PN(Aintscsc1,Aintscss1, Aintscc1, Aintssss1, Aintssc1, Aintcc1,  Anvec_sc1, Anvec_ss1,  Anvec_c1, RY, sporder);
    //cout << "------------returned from genmat_angint_sdm_3D_PN OK------\n";

    RCompRowMatrix Aint2, Aintsc2, Aintss2, Aintc2;
    RCompRowMatrix Aintscsc2,  Aintscss2, Aintscc2,  Aintssss2,  Aintssc2,  Aintcc2;


    genmat_angint_3D(sporder, (sporder+1)*(sporder+1), Aint2, Aintsc2, Aintss2, Aintc2, Aintscsc2, Aintscss2, Aintscc2, Aintssss2, Aintssc2, Aintcc2);

    for(int i=0; i<9; i++){
	for(int j=0; j<9; j++)
		cout<<i<<" "<<j<<" "<<Aintsc1.Get(i, j)<<"  "<<Aintsc2.Get(i, j)<<endl;;
	}
   // cout<<"Aint comparison: "<<(Aint1-Aint2)<<endl;
    /*cout<<"Aintsc1 comparison: "<<Sparse2Dense(Aintsc1)<<endl;
    cout<<"Aintsc2 comparison: "<<Sparse2Dense(Aintsc2)<<endl;
    cout<<"Aintss1 comparison: "<<Sparse2Dense(Aintss1)<<endl;
    cout<<"Aintss2 comparison: "<<Sparse2Dense(Aintss2)<<endl;*/

   

    //cout<<"Aintss1 comparison: "<<Aintss1<<endl;
    //cout<<"Aintss2 comparison: "<<Aintss2<<endl;

    //cout<<"Aintc comparison: "<<(Aintc1-Aintc2)<<endl;


	
   /* RDenseMatrix RRxx = Sparse2Dense(Aintscsc);

    cout << "RRxx\n" << RRxx << endl;

    RDenseMatrix RRyy = Sparse2Dense(Aintssss);

    cout << "RRyy\n" << RRyy << endl;

    RDenseMatrix RRzz = Sparse2Dense(Aintcc);

    cout << "RRzz\n" << RRzz << endl;

    RDenseMatrix RRxy = Sparse2Dense(Aintscss);

    cout << "RRxy\n" << RRxy << endl;

    RDenseMatrix RRxz = Sparse2Dense(Aintscc);

    cout << "RRxz\n" << RRxz << endl;

    RDenseMatrix RRyz = Sparse2Dense(Aintssc);

    cout << "RRyz\n" << RRyz << endl;*/

    /*--------------- do local basis matrix, and rotate -------------*/

    // set up angular "mesh"

   /* Mesh S2Mesh;

    ifstream ifm(argv[1]);
    ifm >> S2Mesh;
    ifm.close();
    cout << "Angular  " << S2Mesh.elen() << " elements, " << S2Mesh.nlen()
         << " nodes\n";
    S2Mesh.Setup();
    const int& SN =  S2Mesh.nlen();       // dimensions are size of nodes.

    cout<<"Calling genmat_angint_3D ..."<<endl;

    genmat_angint_3D(Aint, Aintsc, Aintss, Aintc, Anvec, S2Mesh);
    cerr << "Aint and Anvec are the same ?\n" << (Aint - Anvec) << endl;

    cout<<"Calling genmat_angint_sdm_3D ..."<<endl;
 
    genmat_angint_sdm_3D(Aintscsc,  Aintscss,   Aintscc,  Aintssss,  Aintssc,  Aintcc,  Anvec_sc, Anvec_ss,  Anvec_c, S2Mesh);

    cerr << "Aintc and Anvec_c are the same ?\n" << (Aintc - Anvec_c) << endl;
    cerr << "Aintsc and Anvec_sc are the same ?\n"<<(Aintsc - Anvec_sc)<< endl;
    cerr << "Aintss and Anvec_ss are the same ?\n"<<(Aintss - Anvec_ss)<< endl;

    cout<<"Calling genmat_apu ..."<<endl;*/
    //    genmat_apu(S2Mesh, Anvec, Anvec_sc, Anvec_ss, Anvec_c, apu1, apu1sc, apu1ss, apu1c, sigma, sktyp);
    /* create tables of spherical Harmonic samples */
    /*RDenseMatrix YS(SN,nsp);
    RDenseMatrix YST(nsp,SN);
    for(i = 0; i < SN; i++) {
        YST(0,i) = YS(i,0) = irtpi/2;
        double x = S2Mesh.nlist[i][0];
        double y = -S2Mesh.nlist[i][1]; // seems necessary to change sign..
        double z = S2Mesh.nlist[i][2];
        double rad = sqrt( SQR(x) +SQR(y) + SQR(z)  );
	//        cout << rad << endl;
        x = x/rad; y = y/rad; z = z/rad;
	YST(1,i) = YS(i,1)  = (sqrt(3)*irtpi/2) * x;
	YST(2,i) = YS(i,2)  = (sqrt(3)*irtpi/2) * z;
	YST(3,i) = YS(i,3)  = (sqrt(3)*irtpi/2) * y;
	YST(4,i) = YS(i,4)  = (sqrt(15)*irtpi/4) * (SQR(x) - SQR(y));
	YST(5,i) = YS(i,5)  = (sqrt(15)*irtpi/2) * x*z;
	YST(6,i) = YS(i,6)  = (sqrt(5)*irtpi/4) * (3*SQR(z) - 1);
	YST(7,i) = YS(i,7)  = (sqrt(15)*irtpi/2) * y*z;
	YST(8,i) = YS(i,8)  = -(sqrt(15)*irtpi/2) *x*y;
    }
    cerr << "spherical harmonic samples \n" << YS << endl;
    //   cout <<  Aint;
    RDenseMatrix YY = YST*YS;
    YY = Chop(YY);
    cout << "spherical harmonic inner product:\n" <<   YY << endl;

    RDenseMatrix AAD = Sparse2Dense(Aint);
    RDenseMatrix YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aint :\n" <<   YAY << endl;

    AAD = Sparse2Dense(Aintsc);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintsc :\n" <<   YAY << endl;
    cout << "direct integral RRx\n" << Chop(RRx) << endl;
    cout << "Hand built RRx\n" << expRRx << endl;

    //   cout << "direct integral YYx\n" << Chop(YYx) << endl;

    AAD = Sparse2Dense(Aintss);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintss :\n" <<   YAY << endl;
    cout << "direct integral RRy\n" << Chop(RRy) << endl;
    //    cout << "direct integral YYy\n" << Chop(YYy) << endl;

    AAD = Sparse2Dense(Aintc);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintc :\n" <<   YAY << endl;
    cout << "direct integral RRz\n" << Chop(RRz) << endl;
    //    cout << "direct integral YYz\n" << Chop(YYz) << endl;


    AAD = Sparse2Dense(Aintscsc);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintscsc :\n" << YAY << endl;
    cout << "direct integral RRxx\n" << Chop(RRxx) << endl;
 
    AAD = Sparse2Dense(Aintssss);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintssss :\n" << YAY << endl;
    cout << "direct integral RRyy\n" << Chop(RRyy) << endl;
 
    AAD = Sparse2Dense(Aintcc);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintcc :\n" <<   YAY << endl;
    cout << "direct integral RRzz\n" << Chop(RRzz) << endl;
 

    AAD = Sparse2Dense(Aintscss);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintscss :\n" << YAY << endl;
    cout << "direct integral RRxy\n" << Chop(RRxy) << endl;
 
    AAD = Sparse2Dense(Aintscc);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintscc :\n" << YAY << endl;
    cout << "direct integral RRxz\n" << Chop(RRxz) << endl;
 
    AAD = Sparse2Dense(Aintssc);
    YAY = YST*AAD*YS;
    YAY = Chop(YAY);
    cout << "spherical harmonic product integrals Aintssc :\n" << YAY << endl;
    cout << "direct integral RRyz\n" << Chop(RRyz) << endl;
 */
}




















