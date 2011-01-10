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

#include "rte3D.h"
#include "sphints.h"



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

// main routine **************************************************************

int main (int argc, char *argv[])
{
    char cbuf[200];
    int i,j,el,is,js;
    const int sporder    = 2;
    const int nso = sporder+1;
    int nsp = nso*nso;

    CCompRowMatrix RY;

    RCompRowMatrix Aint, Aintsc, Aintss, Aintc, Anvec;
    RCompRowMatrix Aintscsc,  Aintscss, Aintscc,  Aintssss,  Aintssc,  Aintcc;
    RCompRowMatrix Anvec_sc, Anvec_ss,  Anvec_c;
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
   

    cout << "------------------Calling genmat_angint_3D_PN------------\n";
    genmat_angint_3D_PN(Aint, Aintsc, Aintss, Aintc, Anvec, RY, sporder);
    cout << "------------------returned from genmat_angint_3D_PN OK---\n";
    cerr << "Aint and Anvec are the same ?\n" << (Aint - Anvec) << endl;
    //    cout << "zintegral \n"<< Aintc << endl;

    RDenseMatrix RRz = Sparse2Dense(Aintc);

    cout << "RRz\n" << RRz << endl;
    //   cout << "xintegal\n" << Aintsc << endl;
    RDenseMatrix RRx = Sparse2Dense(Aintsc);   

    cout << "RRx\n" << RRx << endl;
    cout << "Hand built RRx\n" << expRRx << endl;

    RDenseMatrix RRy = Sparse2Dense(Aintss);   

    cout << "RRy\n" << RRy << endl;
 
 
    cout << "--------------Calling genmat_angint_sdm_3D_PN-------------\n";
    genmat_angint_sdm_3D_PN(Aintscsc,Aintscss, Aintscc, Aintssss, Aintssc, Aintcc,  Anvec_sc, Anvec_ss,  Anvec_c, RY, sporder);
    cout << "------------returned from genmat_angint_sdm_3D_PN OK------\n";


    RDenseMatrix RRxx = Sparse2Dense(Aintscsc);

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

    cout << "RRyz\n" << RRyz << endl;

    /*--------------- do local basis matrix, and rotate -------------*/

    // set up angular "mesh"

    Mesh S2Mesh;

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

    cout<<"Calling genmat_apu ..."<<endl;
    //    genmat_apu(S2Mesh, Anvec, Anvec_sc, Anvec_ss, Anvec_c, apu1, apu1sc, apu1ss, apu1c, sigma, sktyp);

    /* create tables of spherical Harmonic samples */
    RDenseMatrix YS(SN,nsp);
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
 
}




















